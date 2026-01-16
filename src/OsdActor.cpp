#include "OsdActor.hpp"
#include "CephCommon.hpp"
#include "MonActor.hpp"
#include "simgrid/s4u/Engine.hpp"
#include "simgrid/s4u/Mailbox.hpp"
#include "xbt/log.h"
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
  std::set<PG *> my_pgs;
  try {
    my_pgs = pgmap->primary_osd_get_pgs(id);
  } catch (const std::out_of_range &) {
    XBT_WARN("OSD %d not in pgmap (likely extra OSD from topology)", id);
    return;
  }
  my_primary_pgs.clear();
  needs_backfill_pgs.clear();

  for (const auto &pg : my_pgs) {
    my_primary_pgs.insert(pg);
    if (pg->needs_backfill()) {
      needs_backfill_pgs.insert(pg);
    }
  }

  std::string pg_str = primary_pgs_to_string();
  if (pg_str.empty()) {
    XBT_DEBUG("No primary pgs");
  } else {
    XBT_DEBUG("%s", pg_str.c_str());
  }

  // todo check if pending backfills can be advanced
}

void Osd::maybe_reserve_backfill() {

  if (pending_retry) {
    return;
  }
  if (backfill_reservation_local &&
      backfill_reservation_remote_pending.empty()) {
    return;
  }

  backfilling_pg =
      needs_backfill_pgs.empty() ? nullptr : *needs_backfill_pgs.begin();

  if (backfilling_pg == nullptr)
    return;

  xbt_assert(!(backfill_reservation_local && !pending_retry &&
               backfill_reservation_remote.size() == 0 &&
               backfill_reservation_remote_pending.size() == 0));
  backfill_reservation_local = true;

  // start greedy remote backfill slot reservation algorithm
  if (backfill_reservation_remote.empty() &&
      backfill_reservation_remote_pending.empty()) {

    std::vector<int> targets = backfilling_pg->get_backfill_targets();
    backfill_reservation_remote_pending =
        std::set<int>(targets.begin(), targets.end());

    for (int pending_osd_id : backfill_reservation_remote_pending) {
      BackfillReservationOp *op = new BackfillReservationOp{
          id, pending_osd_id, backfilling_pg->get_id(),
          BackfillReservationOpType::REQUEST_SLAVE};

      Message *msg = make_message<BackfillReservationMsg>(op);
      sg4::Mailbox *peer_mb = pgmap->get_osd_mailbox(pending_osd_id);
      peer_mb->put_async(msg, 0).detach();
    }
    xbt_assert(!(backfill_reservation_local && !pending_retry &&
                 backfill_reservation_remote.empty() &&
                 backfill_reservation_remote_pending.empty()));
  }
}

void Osd::maybe_schedule_object_backfill() {
  // We check if we *can* schedule (threads, reservation, pg status),
  // but instead of running immediately, we push to dmclock.

  while (used_recovery_threads < max_recovery_threads &&
         backfill_reservation_local && !backfill_reservation_remote.empty() &&
         backfill_reservation_remote_pending.empty() &&
         backfilling_pg->maybe_schedule_recovery()) {

    int op_id = last_op_id++;
    // Create OpContext immediately
    // We allocate on heap to avoid large stack frames, but we dereference to
    // move into queue. Ideally we would construct in place or on stack if small
    // enough. OpContext contains a set, so it's not tiny but manageable. Let's
    // create on stack for simplicity and move.

    OpContext oc = {
        .local_id = op_id,
        .client_op_id = op_id, // strictly internal
        .type = OpType::BACKFILL,
        .pgid = backfilling_pg->get_id(),
        .sender = id,
        .size = backfilling_pg->get_object_size(),
        .state = OpState::OP_QUEUED, // New initial state
    };
    // Do NOT add to op_contexts yet. The callback will do it with the stable
    // address.

    dmc::ReqParams params(1, 1);
    queue->add_request_time(
        std::move(oc), CLIENT_ID_BACKFILL, params,
        sg4::Engine::get_clock()); // todo cost is default 1u
    used_recovery_threads++;
  }
}

void Osd::on_osd_op_message(int sender, const OsdOpMsg &osd_op_msg) {
  Op *op = osd_op_msg.op;

  XBT_DEBUG("received op msg of type %d from osd.%u for pg %u op_id %u",
            static_cast<int>(op->type), sender, op->pgid, op->id);

  switch (op->type) {

  case OpType::CLIENT_WRITE: {
    int op_id = last_op_id++;
    PG *pg = pgmap->get_pg(op->pgid);
    auto acting = pg->get_acting_ids();
    std::set<int> peers(acting.begin(), acting.end());
    OpContext oc = {
        .local_id = op_id,
        .client_op_id = op->id,
        .type = OpType::CLIENT_WRITE,
        .pgid = op->pgid,
        .sender = sender,
        .size = op->size,
        .state = OpState::OP_QUEUED,
        .pending_peers = peers,
    };
    // op_contexts[op_id] = oc; // Don't track yet

    dmc::ReqParams params(1, 1);
    queue->add_request_time(std::move(oc), CLIENT_ID_USER, params,
                            sg4::Engine::get_clock());
    break;
  }

  case OpType::CLIENT_READ: {
    int op_id = last_op_id++;
    OpContext oc = {
        .local_id = op_id,
        .client_op_id = op->id,
        .type = OpType::CLIENT_READ,
        .pgid = op->pgid,
        .sender = sender,
        .size = op->size,
        .state = OpState::OP_QUEUED,
    };
    // op_contexts[op_id] = oc; // Don't track yet

    dmc::ReqParams params(1, 1);
    queue->add_request_time(std::move(oc), CLIENT_ID_USER, params,
                            sg4::Engine::get_clock());
    break;
  }

  case OpType::REPLICA_WRITE: {
    // create opcontext for replica write
    int op_id = last_op_id++;
    OpContext *oc = new OpContext{
        .local_id = op_id,
        .client_op_id = op->id,
        .type = OpType::REPLICA_WRITE,
        .pgid = op->pgid,
        .sender = sender,
        .size = op->size,
        .state = OpState::OP_WAITING_DISK,
    };
    op_contexts[op_id] = oc;
    auto a = disk->write_async(op->size);
    activities.push(a);
    op_context_map[a] = oc;
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

  // this is the context type when we receive an ack from a peer (this is
  // primary's op type, peer's is replica_write)
  case OpType::CLIENT_WRITE:
    context->pending_peers.erase(sender);
    if (context->pending_peers.empty()) {
      op_contexts.erase(msg.op_id);
      // ack client
      sg4::Mailbox *client_mb =
          sg4::Mailbox::by_name("client." + std::to_string(-context->sender));
      Message *ack_msg = make_message<OsdOpAckMsg>(context->client_op_id);
      client_mb->put_async(ack_msg, 0).detach();
      delete context;
    }
    break;

  default:
    xbt_die("osd.%u received ack message with unknown type %d", id,
            static_cast<int>(context->type));
  }
}

void Osd::on_backfill_reservation_message(int sender,
                                          const BackfillReservationMsg &msg) {
  switch (msg.op->type) {

  case BackfillReservationOpType::REQUEST_SLAVE: {
    if (!backfill_reservation_local && backfill_reservation_remote.empty()) {
      backfill_reservation_remote.insert(sender);
      // send accept
      BackfillReservationOp *accept_op = new BackfillReservationOp{
          .type = BackfillReservationOpType::ACCEPT,
      };
      Message *accept_msg = make_message<BackfillReservationMsg>(accept_op);
      pgmap->get_osd_mailbox(sender)->put_async(accept_msg, 0).detach();

    } else {
      // send reject
      BackfillReservationOp *reject_op = new BackfillReservationOp{
          .type = BackfillReservationOpType::REJECT,
      };
      Message *reject_msg = make_message<BackfillReservationMsg>(reject_op);
      pgmap->get_osd_mailbox(sender)->put_async(reject_msg, 0).detach();
    }
    break;
  }

  case BackfillReservationOpType::RELEASE_SLAVE:

    if (!backfill_reservation_local &&
        backfill_reservation_remote.find(sender) !=
            backfill_reservation_remote.end()) {
      backfill_reservation_remote.erase(sender);
    }
    break;

  // For primary from slave
  case BackfillReservationOpType::ACCEPT:
    if (!pending_retry) {
      backfill_reservation_remote.insert(sender);
      backfill_reservation_remote_pending.erase(sender);
      if (backfill_reservation_remote_pending.empty()) {
        XBT_INFO("Starting backfill for pg %i", backfilling_pg->get_id());
        backfilling_pg->set_state(PGState::BACKFILL);
        Message *msg = make_message<PGNotification>(backfilling_pg->get_id());
        mon_mb->put_async(msg, 0).detach();
      }
    }
    break;

  case BackfillReservationOpType::REJECT: {
    // triggers retry to self, if no retry is scheduled
    if (pending_retry)
      break;

    backfill_reservation_remote.clear();
    backfill_reservation_remote_pending.clear();
    backfill_reservation_local = false;

    // immediately release all slaves
    for (int peer_osd_id : backfilling_pg->get_backfill_targets()) {
      if (peer_osd_id == sender) {
        continue;
      }
      BackfillReservationOp *release_op = new BackfillReservationOp{
          .primary_osd_id = id,
          .target_osd_id = peer_osd_id,
          .pg_id = backfilling_pg->get_id(),
          .type = BackfillReservationOpType::RELEASE_SLAVE,
      };
      Message *release_msg = make_message<BackfillReservationMsg>(release_op);
      pgmap->get_osd_mailbox(peer_osd_id)->put_async(release_msg, 0).detach();
    }

    // retry to self
    pending_retry = true;
    sg4::Mailbox *return_mb = mb;
    std::string sender_str = "osd." + std::to_string(id);
    double random_delay = random.uniform_real(0L, 1L);
    my_host->add_actor(
        "backfill_reservation", [return_mb, sender_str, random_delay]() {
          sg4::this_actor::sleep_for(random_delay);
          // create retry op
          BackfillReservationOp *retry_op = new BackfillReservationOp{
              .type = BackfillReservationOpType::RETRY,
          };
          Message *retry_msg = make_message<BackfillReservationMsg>(retry_op);
          retry_msg->sender = sender_str;
          return_mb->put(retry_msg, 0);
        });

    break;
  }

  // for async backoff workaround
  case BackfillReservationOpType::RETRY:
    backfill_reservation_local = false;
    pending_retry = false;
    XBT_DEBUG("retrying backfill for pg %i", backfilling_pg->get_id());
    maybe_reserve_backfill();
    break;

  default:
    xbt_die("to be implemented");
  }
  xbt_assert(!(backfill_reservation_local && !pending_retry &&
               backfill_reservation_remote.size() == 0 &&
               backfill_reservation_remote_pending.size() == 0));
  delete msg.op;
}

void Osd::process_message(Message *msg) {
  xbt_assert(
      !(!backfill_reservation_local && backfill_reservation_remote.size() > 1),
      "illegal state, receiver should only have 1 peer");
  std::visit(
      overloaded{[&](const PGMapNotification &) {
                   XBT_DEBUG("Received PGMapNotification");
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
                   on_osd_op_ack_message(sender_str_to_int(msg->sender),
                                         std::get<OsdOpAckMsg>(msg->payload));
                 },
                 [&](const BackfillReservationMsg &) {
                   int peer = sender_str_to_int(msg->sender);
                   on_backfill_reservation_message(
                       peer, std::get<BackfillReservationMsg>(msg->payload));
                 },
                 [&](const auto &unknown_payload) {
                   xbt_die("OSD %s received unexpected message type: %s",
                           my_host->get_name().c_str(),
                           typeid(unknown_payload).name());
                 }},
      msg->payload);

  xbt_assert(
      !(!backfill_reservation_local && backfill_reservation_remote.size() > 1),
      "illegal state, receiver should only have 1 peer");
  delete msg;
}

void Osd::send_op(Op *op) {
  sg4::Mailbox *target_osd_mb = pgmap->get_osd_mailbox(op->recipient);

  Message *msg = make_message<OsdOpMsg>(op);
  target_osd_mb->put_async(msg, op->size).detach();
}

void Osd::advance_backfill_op(OpContext *context, int peer_osd_id) {
  // peer osd id is the sender of the ack message, used only for OP_WAITING_PEER
  xbt_assert(context->type == OpType::BACKFILL);
  switch (context->state) {

  case OpState::OP_QUEUED: {
    // read from disk, start disk activity
    context->state = OpState::OP_WAITING_DISK;
    auto a = disk->read_async(backfilling_pg->get_object_size());

    activities.push(a);
    op_context_map[a] = context;
    break;
  }

  case OpState::OP_WAITING_DISK:
    // send to peers
    for (auto peer_osd_id : backfilling_pg->get_backfill_targets()) {

      Op *op = new Op{
          .type = OpType::REPLICA_WRITE,
          .id = context->local_id,
          .pgid = backfilling_pg->get_id(),
          .size = backfilling_pg->get_object_size(),
          .recipient = peer_osd_id,
      };
      send_op(op); // network send is not recorded in activities. I don't care
                   // when network send is done.
      XBT_DEBUG("sending backfill op %u to osd.%u", context->local_id,
                op->recipient);
      context->state = OpState::OP_WAITING_PEER;
      context->pending_peers.insert(op->recipient);
    }
    break;

  case OpState::OP_WAITING_PEER:
    context->pending_peers.erase(peer_osd_id);
    if (context->pending_peers.empty()) {
      backfilling_pg
          ->on_object_recovered(); // if it finished this updates acting = up
      op_contexts.erase(context->local_id);
      delete context;
      used_recovery_threads--;
      XBT_DEBUG("backfill op %u completed", context->local_id);

      if (!backfilling_pg->needs_backfill()) {
        XBT_INFO("Backfill complete for pg %u", backfilling_pg->get_id());
        backfilling_pg->set_state(PGState::ACTIVE_CLEAN);

        // notify monitor of backfill completion
        Message *msg = make_message<PGNotification>(backfilling_pg->get_id());
        mon_mb->put_async(msg, 0).detach();

        // free remote reservations, comment is wrong???
        for (auto peer_osd_id : backfilling_pg->get_backfill_targets()) {
          BackfillReservationOp *reservation_op = new BackfillReservationOp{
              .primary_osd_id = id,
              .target_osd_id = peer_osd_id,
              .pg_id = backfilling_pg->get_id(),
              .type = BackfillReservationOpType::RELEASE_SLAVE,
          };
          Message *msg = make_message<BackfillReservationMsg>(reservation_op);
          sg4::Mailbox *peer_mb = pgmap->get_osd_mailbox(peer_osd_id);
          peer_mb->put_async(msg, 0).detach();
        }

        // free local reservations
        needs_backfill_pgs.erase(backfilling_pg);
        backfilling_pg = nullptr;
        backfill_reservation_remote.clear();
        backfill_reservation_local = false;
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

  case OpType::REPLICA_WRITE: {
    // disk write finished; send ack to peer
    sg4::Mailbox *peer_mb =
        sg4::Mailbox::by_name("osd." + std::to_string(context->sender));
    Message *peer_msg = make_message<OsdOpAckMsg>(context->client_op_id);
    peer_mb->put_async(peer_msg, 0).detach();

    op_contexts.erase(context->local_id);
    op_context_map.erase(activity);
    delete context;
    break;
  }

  case OpType::CLIENT_READ: {
    op_context_map.erase(activity);
    Message *ack_msg = make_message<OsdOpAckMsg>(context->client_op_id);
    sg4::Mailbox *target_mb =
        sg4::Mailbox::by_name("client." + std::to_string(-context->sender));
    target_mb->put_async(ack_msg, context->size).detach();
    break;
  }

  default:
    xbt_die("osd.%u finished activity with unknown type", id);
  }
}

void Osd::opcontext_dispatch(OpContext *context) {
  switch (context->type) {

  case OpType::BACKFILL:
    advance_backfill_op(context, id);
    break;

  case OpType::CLIENT_READ: {
    // start disk read activity, update op context state
    auto a = disk->read_async(context->size);
    activities.push(a);
    op_context_map[a] = context;
    break;
  }

  case OpType::CLIENT_WRITE:
    // send replica write ops to peers (including self), update op context state
    for (int peer_osd_id : pgmap->get_pg(context->pgid)->get_acting_ids()) {
      Op *op = new Op{
          .type = OpType::REPLICA_WRITE,
          .id = context->local_id,
          .pgid = context->pgid,
          .size = context->size,
          .recipient = peer_osd_id,
      };
      send_op(op);
    }
    break;

  default:
    xbt_die("Dispatched opcontext of unknown type %d", context->type);
  }
}

void Osd::make_progress() {
  maybe_reserve_backfill();
  maybe_schedule_object_backfill();
  auto result = queue->pull_request(sg4::Engine::get_clock());
  if (result.is_retn()) {
    auto &retn = result.get_retn();
    OpContext *oc = retn.request.release();
    op_contexts[oc->local_id] = oc;
    opcontext_dispatch(oc);
  }
}

void Osd::operator()() {
  XBT_INFO("Started");

  CephActor::main_loop();
  xbt_die("should never reach here");
}

Osd::Osd(PGMap *pgmap, int osd_id, std::string disk_name, double iops,
         SchedulerProfile profile)
    : CephActor(osd_id, pgmap) {

  xbt_assert(osd_id >= 0, "osd_id must be non-negative");
  disk = my_host->get_disk_by_name(disk_name);
  on_pgmap_change();
  init_scheduler(iops, profile);
}

void Osd::init_scheduler(double iops, SchedulerProfile profile) {
  // following ceph implementation where bytes are used as the unit
  // the cost is calcuated as: base_cost + io_size, where base_cost =
  // bandwidth/random_iops
  double avg_bandwidth =
      (disk->get_read_bandwidth() + disk->get_write_bandwidth()) / 2;
  base_cost = avg_bandwidth / iops;
  double limit = base_cost + avg_bandwidth;

  // Client Configurations (reservation, weight, limit)
  user_info = std::make_unique<dmc::ClientInfo>(0, 0, 0);
  backfill_info = std::make_unique<dmc::ClientInfo>(0, 0, 0);
  // profiles copied from
  // https://docs.ceph.com/en/latest/rados/configuration/mclock-config-ref/
  switch (profile) {
  case SchedulerProfile::BALANCED:
    user_info->update(limit * 0.5, 1, limit);
    backfill_info->update(0, 1, limit * 0.9);
    break;
  case SchedulerProfile::HIGH_CLIENT_OPS:
    user_info->update(limit * 0.6, 2, limit);
    backfill_info->update(0, 1, limit * 0.7);
    break;
  case SchedulerProfile::HIGH_RECOVERY_OPS:
    user_info->update(limit * 0.3, 1, limit);
    backfill_info->update(0, 1, limit);
    break;
  default:
    xbt_die("Unknown scheduler profile %d", profile);
  }

  auto client_info_f = [this](ClientId c) -> const dmc::ClientInfo * {
    if (c == CLIENT_ID_USER)
      return user_info.get();
    if (c == CLIENT_ID_BACKFILL)
      return backfill_info.get();
    return nullptr;
  };

  queue = std::make_unique<Queue>(
      client_info_f, std::chrono::seconds(100), std::chrono::seconds(200),
      std::chrono::seconds(50),
      dmc::AtLimit::Wait); // only second time arg means anything for the
                           // simulation (idle age)
}
