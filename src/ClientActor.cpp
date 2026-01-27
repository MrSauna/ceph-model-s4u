#include "ClientActor.hpp"
#include "CephCommon.hpp"

XBT_LOG_NEW_DEFAULT_CATEGORY(s4u_ceph_sim_client,
                             "Messages specific for ClientActors");

namespace sg4 = simgrid::s4u;

std::ofstream Client::metrics_stream;
std::mutex Client::metrics_mutex;

void Client::set_metrics_output(const std::string &filename) {
  metrics_stream.open(filename);
  if (metrics_stream.is_open()) {
    metrics_stream << "time,client_id,op_type,op_size,duration\n";
  } else {
    XBT_ERROR("Failed to open client metrics file: %s", filename.c_str());
  }
}

Client::Client(PGMap *pgmap, int client_id, int read_queue, int write_queue)
    : CephActor(client_id, pgmap), max_concurrent_reads(read_queue),
      max_concurrent_writes(write_queue) {
  xbt_assert(client_id < 0, "client_id must be negative");
  random.set_seed(client_id);
}

void Client::gen_op(OpType type) {

  if (type == OpType::CLIENT_READ) {
    in_flight_reads++;
  } else if (type == OpType::CLIENT_WRITE) {
    in_flight_writes++;
  }
  // Create Op
  int op_id = last_op_id++;
  int pg_id = last_op_pg++ % pgmap->size();
  int target_osd_id = pgmap->get_pg(pg_id)->get_primary();

  Op *op = new Op{
      .type = type,
      .id = op_id,
      .pgid = pg_id,
      .size = pgmap->get_object_size(),
      .recipient = target_osd_id,
  };

  // Track Context
  OpContext *oc = new OpContext{
      .local_id = op_id,
      .client_op_id = op_id,
      .type = type,
      .pgid = pg_id,
      .sender = id,
      .size = op->size,
      .state = OpState::OP_WAITING_PEER,
      .start_time = sg4::Engine::get_clock(),
  };
  op_contexts[op_id] = oc;

  // Send
  sg4::Mailbox *target_osd_mb = pgmap->get_osd_mailbox(target_osd_id);
  Message *msg = make_message<OsdOpMsg>(op);

  target_osd_mb
      ->put_async(msg,
                  type == OpType::CLIENT_READ ? 0 : pgmap->get_object_size())
      .detach();

  // XBT_INFO("Sent op %d to %s", op_id, target_osd_mb->get_cname());
}

std::optional<double> Client::make_progress() {
  if ((in_flight_reads >= max_concurrent_reads &&
       in_flight_writes >= max_concurrent_writes) ||
      shutting_down)
    return std::nullopt;

  // Fill both pipelines
  while (in_flight_reads < max_concurrent_reads) {
    gen_op(OpType::CLIENT_READ);
  }
  while (in_flight_writes < max_concurrent_writes) {
    gen_op(OpType::CLIENT_WRITE);
  }

  return std::nullopt;
}

void Client::on_osd_op_ack_message(int sender, const OsdOpAckMsg &msg) {

  OpContext *context = op_contexts[msg.op_id];
  XBT_DEBUG("received ack for op %u", msg.op_id);
  xbt_assert(context, "op context is null");

  switch (context->type) {
  case OpType::CLIENT_WRITE: {
    double dt = sg4::Engine::get_clock() - context->start_time;
    XBT_DEBUG("wr op %d took %f seconds", msg.op_id, dt);
    {
      std::lock_guard<std::mutex> lock(metrics_mutex);
      if (metrics_stream.is_open()) {
        metrics_stream << context->start_time << "," << id << ","
                       << "write" << "," << context->size << "," << dt << "\n";
      }
    }
    break;
  }
  case OpType::CLIENT_READ: {
    double dt = sg4::Engine::get_clock() - context->start_time;
    XBT_DEBUG("rd op %d took %f seconds", msg.op_id, dt);
    {
      std::lock_guard<std::mutex> lock(metrics_mutex);
      if (metrics_stream.is_open()) {
        metrics_stream << context->start_time << "," << id << ","
                       << "read" << "," << context->size << "," << dt << "\n";
      }
    }
    break;
  }
  default:
    xbt_die("osd.%u received ack message with unknown type %d", id,
            static_cast<int>(context->type));
  }

  if (context->type == OpType::CLIENT_READ) {
    in_flight_reads--;
  } else if (context->type == OpType::CLIENT_WRITE) {
    in_flight_writes--;
  }
  op_contexts.erase(msg.op_id);
  delete context;

  xbt_assert(in_flight_reads + in_flight_writes == op_contexts.size());

  if (shutting_down) {
    XBT_INFO("Shutting down, %d read ops and %d write ops remaining",
             in_flight_reads, in_flight_writes);
    // print op contexts with field names
    for (auto const &[key, val] : op_contexts) {
      XBT_INFO("Op %d: local_id=%d, client_op_id=%d, type=%d, pgid=%d, "
               "sender=%d, size=%zu",
               key, val->local_id, val->client_op_id, val->type, val->pgid,
               val->sender, val->size);
    }
  }
  if (shutting_down && (in_flight_reads + in_flight_writes) == 0) {
    auto mon_mb = simgrid::s4u::Mailbox::by_name("mon");
    auto msg = make_message<KillAckMsg>();
    mon_mb->put(msg, 0);
    XBT_INFO("Shutting down now");
    CephActor::kill_self();
  }
}

void Client::process_message(Message *msg) {
  std::visit(
      overloaded{
          [&](const OsdOpAckMsg &ack) {
            XBT_DEBUG("Received Ack for op %d", ack.op_id);
            on_osd_op_ack_message(sender_str_to_int(msg->sender), ack);
          },
          [&](const KillMsg &kill) {
            shutting_down = true;
            if ((in_flight_reads + in_flight_writes) == 0) {
              auto mon_mb = simgrid::s4u::Mailbox::by_name("mon");
              auto msg = make_message<KillAckMsg>();
              mon_mb->put(msg, 0);
              XBT_INFO("Shutting down now");
              CephActor::kill_self();
            } else {
              XBT_INFO("Received KillMsg, draining %d ops (%d read, %d write)",
                       in_flight_reads + in_flight_writes, in_flight_reads,
                       in_flight_writes);
            }
          },
          [&](const auto &) { xbt_die("Client received unexpected message"); }},
      msg->payload);
  delete msg;
}

void Client::process_finished_activity(sg4::ActivityPtr activity) {
  xbt_die("Clients should not have activities");
}

void Client::operator()() { CephActor::main_loop(); }
