#include "ClientActor.hpp"
#include "CephCommon.hpp"

XBT_LOG_NEW_DEFAULT_CATEGORY(s4u_ceph_sim_client,
                             "Messages specific for ClientActors");

namespace sg4 = simgrid::s4u;

// Static member definitions
std::ofstream Client::metrics_stream;
std::mutex Client::metrics_mutex;
std::map<int, ThroughputBucket> Client::throughput_buckets;
digestible::tdigest<float, uint32_t>
    Client::latency_digest(100); // 100 compression factor
bool Client::aggregate_mode = false;
static std::string aggregate_output_dir;

void Client::set_metrics_output(const std::string &filename) {
  metrics_stream.open(filename);
  if (metrics_stream.is_open()) {
    metrics_stream << "time,client_id,op_type,op_size,duration\n";
  } else {
    XBT_ERROR("Failed to open client metrics file: %s", filename.c_str());
  }
}

Client::Client(PGMap *pgmap, int client_id, int read_queue, int write_queue,
               int op_size)
    : CephActor(client_id, pgmap), max_concurrent_reads(read_queue),
      max_concurrent_writes(write_queue), op_size(op_size) {
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
  int pg_id = random.uniform_int(0, pgmap->size() - 1);
  int target_osd_id = pgmap->get_pg(pg_id)->get_primary();

  Op *op = new Op{
      .type = type,
      .id = op_id,
      .pgid = pg_id,
      .size = static_cast<size_t>(op_size),
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

  activities.push(
      target_osd_mb->put_async(msg, type == OpType::CLIENT_READ ? 0 : op_size));

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

  double dt = sg4::Engine::get_clock() - context->start_time;
  int time_bucket = static_cast<int>(context->start_time);

  // Record metrics based on mode
  {
    std::lock_guard<std::mutex> lock(metrics_mutex);

    // Aggregated mode: update T-Digest and throughput buckets
    if (aggregate_mode) {
      latency_digest.insert(static_cast<float>(dt));

      ThroughputBucket &bucket = throughput_buckets[time_bucket];
      if (context->type == OpType::CLIENT_READ) {
        bucket.read_ops++;
        bucket.read_bytes += context->size;
      } else if (context->type == OpType::CLIENT_WRITE) {
        bucket.write_ops++;
        bucket.write_bytes += context->size;
      }
    }

    // Legacy per-op CSV (if stream is open)
    if (metrics_stream.is_open()) {
      const char *op_name =
          (context->type == OpType::CLIENT_READ) ? "read" : "write";
      metrics_stream << context->start_time << "," << id << "," << op_name
                     << "," << context->size << "," << dt << "\n";
    }
  }

  if (context->type == OpType::CLIENT_READ) {
    in_flight_reads--;
  } else if (context->type == OpType::CLIENT_WRITE) {
    in_flight_writes--;
  }
  op_contexts.erase(msg.op_id);
  delete context;

  xbt_assert(in_flight_reads + in_flight_writes == op_contexts.size());

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
              XBT_INFO("op_context.size() = %zu", op_contexts.size());
              auto some_op = op_contexts.begin()->second;
              XBT_INFO("Received KillMsg, draining %d ops (%d read, %d write)",
                       in_flight_reads + in_flight_writes, in_flight_reads,
                       in_flight_writes);
            }
          },
          [&](const auto &) { xbt_die("Client received unexpected message"); }},
      msg->payload);
  delete msg;
}

void Client::on_finished_activity(sg4::ActivityPtr activity) {
  // Async sends complete here - no special handling needed
}

void Client::operator()() { CephActor::main_loop(); }

void Client::set_aggregate_output(const std::string &output_dir) {
  aggregate_mode = true;
  aggregate_output_dir = output_dir;
  XBT_INFO("Client metrics aggregation enabled, output dir: %s",
           output_dir.c_str());
}

void Client::write_aggregated_metrics() {
  if (!aggregate_mode) {
    return;
  }

  // Merge any buffered data in T-Digest
  latency_digest.merge();

  // Write throughput aggregates: client_metrics_agg.csv
  std::string throughput_path =
      aggregate_output_dir + "/client_metrics_agg.csv";
  std::ofstream tp_out(throughput_path);
  if (tp_out.is_open()) {
    tp_out << "time,read_ops,read_bytes,write_ops,write_bytes\n";
    for (const auto &[time, bucket] : throughput_buckets) {
      tp_out << time << "," << bucket.read_ops << "," << bucket.read_bytes
             << "," << bucket.write_ops << "," << bucket.write_bytes << "\n";
    }
    tp_out.close();
    XBT_INFO("Wrote throughput aggregates to %s (%zu buckets)",
             throughput_path.c_str(), throughput_buckets.size());
  }

  // Write T-Digest centroids: client_latency_digest.csv
  std::string digest_path = aggregate_output_dir + "/client_latency_digest.csv";
  std::ofstream digest_out(digest_path);
  auto centroids = latency_digest.get(); // Returns vector<pair<mean, weight>>

  if (digest_out.is_open()) {
    digest_out << "mean,weight\n";
    for (const auto &[mean, weight] : centroids) {
      digest_out << mean << "," << weight << "\n";
    }
    digest_out.close();
    XBT_INFO("Wrote latency T-Digest to %s (%zu centroids)",
             digest_path.c_str(), centroids.size());
  }

  // Compute mean from centroids (weighted average)
  double weighted_sum = 0.0;
  size_t total_weight = 0;
  for (const auto &[mean, weight] : centroids) {
    weighted_sum += mean * weight;
    total_weight += weight;
  }
  double avg = (total_weight > 0) ? weighted_sum / total_weight : 0.0;

  // Also write computed percentiles for convenience (summary file)
  std::string summary_path =
      aggregate_output_dir + "/client_latency_summary.csv";
  std::ofstream summary_out(summary_path);
  if (summary_out.is_open()) {
    summary_out << "metric,value\n";
    summary_out << "count," << latency_digest.size() << "\n";
    summary_out << "min," << latency_digest.min() << "\n";
    summary_out << "max," << latency_digest.max() << "\n";
    summary_out << "avg," << avg << "\n";
    summary_out << "p50," << latency_digest.quantile(50) << "\n";
    summary_out << "p95," << latency_digest.quantile(95) << "\n";
    summary_out << "p99," << latency_digest.quantile(99) << "\n";
    summary_out << "p99.5," << latency_digest.quantile(99.5) << "\n";
    summary_out.close();
    XBT_INFO("Wrote latency summary to %s", summary_path.c_str());
  }
}
