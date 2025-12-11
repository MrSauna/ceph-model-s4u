#include "OsdActor.hpp"
#include "CephCommon.hpp"
#include "MonActor.hpp"
#include <sstream>
#include <variant>

constexpr size_t OBJECT_SIZE(4 * 1024 * 1024); // 4MB

XBT_LOG_NEW_DEFAULT_CATEGORY(s4u_ceph_sim_osd,
                             "Messages specific for OsdActors");

std::string Osd::primary_pgs_to_string() const {
  std::ostringstream ss;

  if (!my_primary_pgs.empty())
    ss << "osd." << osd_id << ": primary for pgs ";
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
  auto my_pgs = pgmap->primary_osd_get_pgs(osd_id);
  my_primary_pgs.clear();
  needs_backfill_pgs.clear();

  for (const auto &pg : my_pgs) {
    my_primary_pgs.insert(pg);
    if (pg->needs_backfill())
      needs_backfill_pgs.insert(pg);
  }

  XBT_INFO("%s", primary_pgs_to_string().c_str());
}

void Osd::kill_self() {
  XBT_INFO("%s is killing itself", my_host->get_name().c_str());
  simgrid::s4u::this_actor::exit();
}

void Osd::maybe_reserve_backfill() {
  if (backfilling_pg)
    return;
  backfilling_pg =
      needs_backfill_pgs.empty() ? nullptr : *needs_backfill_pgs.begin();
  if (backfilling_pg == nullptr)
    return;
  XBT_INFO("osd.%u is backfilling pg %u", osd_id, backfilling_pg->get_id());
}

void Osd::maybe_schedule_object_backfill() {
  if (used_recovery_threads >= max_recovery_threads || !backfilling_pg)
    return;

  auto a = disk->read_async(OBJECT_SIZE);
  activities.push(a);
  op_context_map[a] = OpContext{
      .id = 0, // TODO
      .type = OpType::BACKFILL,
      .pgid = backfilling_pg->get_id(),
      .sender = osd_id,
      .size = OBJECT_SIZE,
      .state = OpState::OP_WAITING_DISK,
  };
  used_recovery_threads++;
}

void Osd::on_osd_op_message(const OsdOpMsg &osd_op_msg) {
  Op *op = osd_op_msg.op;

  XBT_INFO("osd.%u received op message of type %d from osd.%u for pg %u",
           osd_id, static_cast<int>(op->type), op->sender, op->pgid);

  switch (op->type) {
  case OpType::REPLICA_WRITE: {
    sg4::Mailbox *target_osd_mb =
        simgrid::s4u::Mailbox::by_name("osd." + std::to_string(op->recipient));
    Message *ack_msg = make_message<OsdOpAckMsg>(op->id);
    target_osd_mb->put_async(ack_msg, 0).detach(); // note detached()
    break;
  }
  default:
    xbt_die("osd.%u received op message with unknown type %d", osd_id,
            static_cast<int>(op->type));
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
                          on_osd_op_message(std::get<OsdOpMsg>(msg->payload));
                        },
                        [&](const auto &unknown_payload) {
                          xbt_die("OSD %s received unexpected message type: %s",
                                  my_host->get_name().c_str(),
                                  typeid(unknown_payload).name());
                        }},
             msg->payload);

  delete msg;
}

// todo: should this wait for awknowledgment?
// Do I wait for acks or do I wait for network activity to finish?
void Osd::send_op(Op *op) {
  // TODO: store the mailboxes somewhere always looking up by name seems like a
  // bad idea.
  sg4::Mailbox *target_osd_mb =
      simgrid::s4u::Mailbox::by_name("osd." + std::to_string(op->recipient));

  Message *msg = make_message<OsdOpMsg>(op);
  target_osd_mb->put_async(msg, op->size).detach();
}

void Osd::advance_backfill_op(OpContext &context) {
  switch (context.state) {
  case OpState::OP_WAITING_DISK:
    context.state = OpState::OP_WAITING_PEER;
    // send to peers
    for (auto peer_osd_shard : backfilling_pg->get_up().members) {

      if (peer_osd_shard->get_osd_id() == osd_id || peer_osd_shard->is_acting())
        continue;

      Op *op = new Op{
          .type = OpType::REPLICA_WRITE,
          .id = 0, // TODO
          .sender = osd_id,
          .recipient = peer_osd_shard->get_osd_id(),
          .pgid = backfilling_pg->get_id(),
          .size = OBJECT_SIZE,
      };
      send_op(op); // network send is not recorded in activities. I don't care
                   // when network send is done.
    }
  case OpState::OP_WAITING_PEER:
    // process_ack_message
    context.state = OpState::OP_COMPLETED;
    xbt_assert(context.pending_peers.empty(),
               "osd.%u backfill op still has pending peers", osd_id);
    used_recovery_threads--;
    break;
  default:
    xbt_die("osd.%u advancing backfill op in unknown state", osd_id);
  }
  needs_backfill_pgs.erase(backfilling_pg);
  backfilling_pg = nullptr;
}

void Osd::process_finished_activity(sg4::ActivityPtr activity) {

  if (op_context_map.find(activity) == op_context_map.end()) {
    xbt_die("osd.%u finished unknown activity", osd_id);
  }

  OpContext context = op_context_map[activity];

  switch (context.type) {
  case OpType::BACKFILL:
    advance_backfill_op(context);
    break;
  default:
    xbt_die("osd.%u finished activity with unknown type", osd_id);
  }

  // fixme: activity type based processing
  // get rid of the hardcoded backfill send
}

void Osd::main_loop() {

  // start off listener
  Message *message;
  sg4::CommPtr listener = mb->get_async(&message);
  activities.push(listener);

  while (true) {

    maybe_reserve_backfill();
    maybe_schedule_object_backfill(); // dies because pgmap can be null. If
                                      // another osd is faster to start it can
                                      // send a backfill operation before I have
                                      // received the pgmap.

    // wait any activity
    sg4::ActivityPtr finished = activities.wait_any();

    // new message arrived
    if (finished == listener) {
      XBT_INFO("received message");
      process_message(message);

      // restart listener
      listener = mb->get_async(&message);
      activities.push(listener);
      continue;
    }

    // something else finished
    else {
      process_finished_activity(finished);
    }

    simgrid::s4u::this_actor::sleep_for(0.1); // todo: careful with this
  }
}

void Osd::operator()() {
  XBT_INFO("I was created %s", my_host->get_cname());

  main_loop();
  xbt_assert(false, "should never reach here");
}

Osd::Osd(PGMap *pgmap, int osd_id) : pgmap(pgmap), osd_id(osd_id) {

  my_host = simgrid::s4u::this_actor::get_host();
  mb = simgrid::s4u::Mailbox::by_name(my_host->get_name());
  mon_mb = simgrid::s4u::Mailbox::by_name("mon");

  disk = my_host->get_disks().at(0);
  on_pgmap_change();
}