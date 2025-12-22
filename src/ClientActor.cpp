#include "ClientActor.hpp"

XBT_LOG_NEW_DEFAULT_CATEGORY(s4u_ceph_sim_client,
                             "Messages specific for ClientActors");

namespace sg4 = simgrid::s4u;

Client::Client(PGMap *pgmap, int client_id) : CephActor(client_id, pgmap) {
  xbt_assert(client_id < 0, "client_id must be negative");
}

void Client::make_progress() {

  // Create Op
  int op_id = last_op_id++;
  int pg_id = last_op_pg++ % pgmap->size();
  int target_osd_id = pgmap->get_pg(pg_id)->get_primary();

  Op *op = new Op{
      .type = OpType::CLIENT_WRITE,
      .id = op_id,
      .recipient = target_osd_id,
      .pgid = pg_id,
      .size = (size_t)communication_cost,
  };

  // Track Context
  OpContext *oc = new OpContext{
      .id = op_id,
      .type = OpType::CLIENT_WRITE,
      .pgid = pg_id,
      .sender = id,
      .size = op->size,
      .state = OpState::OP_WAITING_PEER,
  };
  op_contexts[op_id] = oc;

  // Send
  // Send
  sg4::Mailbox *target_osd_mb = pgmap->get_osd_mailbox(target_osd_id);
  Message *msg = make_message<OsdOpMsg>(op);
  // Overwrite sender with client name/ID if needed, but make_message uses
  // actor name

  target_osd_mb->put_async(msg, communication_cost).detach();

  XBT_INFO("Sent op %d to %s", op_id, target_osd_mb->get_cname());
}

void Client::process_message(Message *msg) {
  std::visit(overloaded{[&](const OsdOpAckMsg &ack) {
                          if (op_contexts.count(ack.op_id)) {
                            XBT_DEBUG("Received Ack for op %d", ack.op_id);
                            delete op_contexts[ack.op_id];
                            op_contexts.erase(ack.op_id);
                            tasks_acked++;
                          } else {
                            XBT_WARN("Received Ack for unknown op %d",
                                     ack.op_id);
                          }
                        },
                        [&](const auto &) {
                          XBT_WARN("Client received unexpected message");
                        }},
             msg->payload);
  delete msg;
}

void Client::process_finished_activity(sg4::ActivityPtr activity) {
  // Only listener activities are tracked in main loop if we don't push others.
  // Client put_async are detached.
  // implementation needed for abstract base class
}

void Client::operator()() { CephActor::main_loop(); }
