#include "OsdActor.hpp"
#include "CephCommon.hpp"
#include "MonActor.hpp"
#include <sstream>
#include <variant>

XBT_LOG_NEW_DEFAULT_CATEGORY(s4u_ceph_sim_osd,
                             "Messages specific for OsdActors");

std::string Osd::primary_pgs_to_string() const {
  std::ostringstream ss;

  if (!my_primary_pgs.empty())
    ss << "primary for pgs ";
  std::string sep = "";
  for (const auto &pg : my_primary_pgs) {
    ss << sep << pg->get_id();
    sep = ", ";
  }

  if (!needs_backfill_pgs.empty())
    ss << " needs backfill";
  for (const auto &pg : needs_backfill_pgs) {
    ss << sep << pg->get_id();
  }

  return ss.str();
}

void Osd::on_pgmap_change() {
  auto my_pgs = pgmap->primary_osd_get_pgs(id);
  my_primary_pgs.clear();
  needs_backfill_pgs.clear();

  for (const auto &pg : my_pgs) {
    my_primary_pgs.insert(pg);
    if (pg->needs_backfill())
      needs_backfill_pgs.insert(pg);
  }

  std::string pg_str = primary_pgs_to_string();
  if (pg_str.empty()) {
    XBT_INFO("No primary pgs");
  } else {
    XBT_INFO("%s", pg_str.c_str());
  }
}

void Osd::maybe_reserve_backfill() {
  if (backfilling_pg)
    return;
  backfilling_pg =
      needs_backfill_pgs.empty() ? nullptr : *needs_backfill_pgs.begin();
  if (backfilling_pg == nullptr)
    return;
  XBT_INFO("backfilling pg %u", backfilling_pg->get_id());
}

void Osd::maybe_schedule_object_backfill() {
  if (used_recovery_threads >= max_recovery_threads || !backfilling_pg ||
      !backfilling_pg->schedule_recovery())
    return;

  int op_id = last_op_id++;
  auto a = disk->read_async(backfilling_pg->get_object_size());
  activities.push(a);
  OpContext *oc = new OpContext{
      .id = op_id,
      .type = OpType::BACKFILL,
      .pgid = backfilling_pg->get_id(),
      .sender = id,
      .size = backfilling_pg->get_object_size(),
      .state = OpState::OP_WAITING_DISK,
  };
  op_context_map[a] = oc;
  op_contexts[op_id] = oc;

  used_recovery_threads++;
}

// todo: make it do stuff, currently instant ack
void Osd::on_osd_op_message(int sender, const OsdOpMsg &osd_op_msg) {
  Op *op = osd_op_msg.op;

  XBT_DEBUG("received op msg of type %d from osd.%u for pg %u op_id %u",
            static_cast<int>(op->type), sender, op->pgid, op->id);

  switch (op->type) {
  case OpType::REPLICA_WRITE: {
    sg4::Mailbox *target_osd_mb = pgmap->get_osd_mailbox(sender);
    Message *ack_msg = make_message<OsdOpAckMsg>(op->id);
    target_osd_mb->put_async(ack_msg, 0).detach(); // note detached()
    XBT_DEBUG("received replica write op %u from osd.%u. Acking immediately",
              op->id, sender);
    break;
  }
  default:
    xbt_die("osd.%u received op message with unknown type %d", id,
            static_cast<int>(op->type));
  }
}

void Osd::on_osd_op_ack_message(int sender, const OsdOpAckMsg &msg) {
  OpContext *context = op_contexts[msg.op_id];
  XBT_DEBUG("received ack for op %u", msg.op_id);
  xbt_assert(context, "op context is null");

  switch (context->type) {
  case OpType::BACKFILL:
    advance_backfill_op(context, sender);
    break;
  default:
    xbt_die("osd.%u received ack message with unknown type %d", id,
            static_cast<int>(context->type));
  }
}

void Osd::process_message(Message *msg) {
  std::visit(overloaded{[&](const PGMapNotification &) {
                          XBT_INFO("Received PGMapNotification");
                          on_pgmap_change();
                        },
                        [&](const KillMsg &) {
                          XBT_INFO("Received KillMsg, exiting");
                          simgrid::s4u::this_actor::exit();
                        },
                        [&](const OsdOpMsg &) {
                          on_osd_op_message(sender_str_to_int(msg->sender),
                                            std::get<OsdOpMsg>(msg->payload));
                        },
                        [&](const OsdOpAckMsg &) {
                          on_osd_op_ack_message(
                              sender_str_to_int(msg->sender),
                              std::get<OsdOpAckMsg>(msg->payload));
                        },
                        [&](const auto &unknown_payload) {
                          xbt_die("OSD %s received unexpected message type: %s",
                                  my_host->get_name().c_str(),
                                  typeid(unknown_payload).name());
                        }},
             msg->payload);

  delete msg;
}

void Osd::send_op(Op *op) {
  // TODO: store the mailboxes somewhere always looking up by name seems like a
  // bad idea.
  sg4::Mailbox *target_osd_mb = pgmap->get_osd_mailbox(op->recipient);

  Message *msg = make_message<OsdOpMsg>(op);
  target_osd_mb->put_async(msg, op->size).detach();
}

void Osd::advance_backfill_op(OpContext *context, int peer_osd_id) {
  // peer osd id is the sender of the ack message, used only for OP_WAITING_PEER
  switch (context->state) {
  case OpState::OP_WAITING_DISK:
    // send to peers
    for (auto peer_osd_shard : backfilling_pg->get_up().members) {

      if (peer_osd_shard->get_osd_id() == id || peer_osd_shard->is_acting())
        continue;

      Op *op = new Op{
          .type = OpType::REPLICA_WRITE,
          .id = context->id,
          .recipient = peer_osd_shard->get_osd_id(),
          .pgid = backfilling_pg->get_id(),
          .size = backfilling_pg->get_object_size(),
      };
      send_op(op); // network send is not recorded in activities. I don't care
                   // when network send is done.
      XBT_DEBUG("sending backfill op %u to osd.%u", context->id, op->recipient);
      context->state = OpState::OP_WAITING_PEER;
      context->pending_peers.insert(op->recipient);
    }
    break;
  case OpState::OP_WAITING_PEER:
    context->pending_peers.erase(peer_osd_id);
    if (context->pending_peers.empty()) {
      backfilling_pg->on_object_recovered();
      op_contexts.erase(context->id);
      delete context;
      used_recovery_threads--;
      XBT_DEBUG("backfill op %u completed", context->id);

      if (!backfilling_pg->needs_backfill()) {
        XBT_INFO("Backfill complete for pg %u", backfilling_pg->get_id());
        needs_backfill_pgs.erase(backfilling_pg);
        backfilling_pg = nullptr;
      }
    }
    break;
  default:
    xbt_die("osd.%u advancing backfill op in unknown state", id);
  }
}

void Osd::process_finished_activity(sg4::ActivityPtr activity) {

  if (op_context_map.find(activity) == op_context_map.end()) {
    xbt_die("osd.%u finished unknown activity", id);
  }

  OpContext *context = op_context_map[activity];

  switch (context->type) {
  case OpType::BACKFILL:
    op_context_map.erase(activity);
    advance_backfill_op(context, id);
    break;
  default:
    xbt_die("osd.%u finished activity with unknown type", id);
  }
}

void Osd::make_progress() {
  maybe_reserve_backfill();
  maybe_schedule_object_backfill();
}

void Osd::operator()() {
  XBT_INFO("Started");

  CephActor::main_loop();
  xbt_die("should never reach here");
}

Osd::Osd(PGMap *pgmap, int osd_id, std::string disk_name)
    : CephActor(osd_id, pgmap) {

  xbt_assert(osd_id >= 0, "osd_id must be non-negative");
  disk = my_host->get_disk_by_name(disk_name);
  on_pgmap_change();
}