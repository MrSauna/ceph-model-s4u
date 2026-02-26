#include "ClientActor.hpp"
#include "CephCommon.hpp"

XBT_LOG_NEW_DEFAULT_CATEGORY(s4u_ceph_sim_client,
                             "Messages specific for ClientActors");

namespace sg4 = simgrid::s4u;

// Static member definitions
std::ofstream Client::metrics_stream;
std::mutex Client::metrics_mutex;
std::map<int, ThroughputBucket> Client::throughput_buckets;
digestible::tdigest<float, uint32_t> Client::read_latency_digest(100);
digestible::tdigest<float, uint32_t> Client::write_latency_digest(100);
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
               int op_size, double read_bandwidth, double write_bandwidth)
    : CephActor(client_id, pgmap), max_concurrent_reads(read_queue),
      max_concurrent_writes(write_queue), op_size(op_size),
      read_bandwidth(read_bandwidth), write_bandwidth(write_bandwidth) {
  xbt_assert(client_id < 0, "client_id must be negative");
  random.set_seed(client_id);
  read_tokens =
      read_bandwidth > 0 ? 2.0 * op_size : 0; // max 2 ops burst capacity
  write_tokens =
      write_bandwidth > 0 ? 2.0 * op_size : 0; // max 2 ops burst capacity
  last_token_time = get_mock_epoch();
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

  // double jitter = random.uniform_real(0, 0.01);
  // sg4::this_actor::sleep_for(jitter);

  // XBT_INFO("Sent op %d to %s", op_id, target_osd_mb->get_cname());
}

std::optional<double> Client::make_progress() {
  if (shutting_down)
    return std::nullopt;

  double now = get_mock_epoch();
  double dt = now - last_token_time;
  last_token_time = now;

  if (dt > 0) {
    if (read_bandwidth > 0) {
      read_tokens += dt * read_bandwidth;
      read_tokens = std::min(read_tokens, 2.0 * op_size); // max 2 ops buffer
    }
    if (write_bandwidth > 0) {
      write_tokens += dt * write_bandwidth;
      write_tokens = std::min(write_tokens, 2.0 * op_size); // max 2 ops buffer
    }
  }

  // Fill both pipelines
  while (in_flight_reads < max_concurrent_reads) {
    if (read_bandwidth > 0) {
      if (read_tokens >= op_size) {
        read_tokens -= op_size;
        gen_op(OpType::CLIENT_READ);
      } else {
        break; // not enough tokens
      }
    } else {
      gen_op(OpType::CLIENT_READ);
    }
  }

  while (in_flight_writes < max_concurrent_writes) {
    if (write_bandwidth > 0) {
      if (write_tokens >= op_size) {
        write_tokens -= op_size;
        gen_op(OpType::CLIENT_WRITE);
      } else {
        break; // not enough tokens
      }
    } else {
      gen_op(OpType::CLIENT_WRITE);
    }
  }

  if (in_flight_reads >= max_concurrent_reads &&
      in_flight_writes >= max_concurrent_writes) {
    return std::nullopt;
  }

  double next_wakeup = -1.0;

  if (in_flight_reads < max_concurrent_reads && read_bandwidth > 0) {
    double needed = op_size - read_tokens;
    double time_needed = needed / read_bandwidth;
    next_wakeup = now + time_needed;
  }

  if (in_flight_writes < max_concurrent_writes && write_bandwidth > 0) {
    double needed = op_size - write_tokens;
    double time_needed = needed / write_bandwidth;
    if (next_wakeup < 0 || now + time_needed < next_wakeup) {
      next_wakeup = now + time_needed;
    }
  }

  if (next_wakeup > 0) {
    // main loop expects epoch time, it's used as timeout
    // Prevent 0s sleep / hot loop from FP inaccuracies. clamp to minimum 1
    // microsecond.
    return std::max(next_wakeup, now + 1e-6);
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
      if (context->type == OpType::CLIENT_READ) {
        read_latency_digest.insert(static_cast<float>(dt));
      } else if (context->type == OpType::CLIENT_WRITE) {
        write_latency_digest.insert(static_cast<float>(dt));
      }

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
  read_latency_digest.merge();
  write_latency_digest.merge();

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

  if (digest_out.is_open()) {
    digest_out << "op_type,mean,weight\n";
    auto read_centroids = read_latency_digest.get();
    for (const auto &[mean, weight] : read_centroids) {
      digest_out << "read," << mean << "," << weight << "\n";
    }
    auto write_centroids = write_latency_digest.get();
    for (const auto &[mean, weight] : write_centroids) {
      digest_out << "write," << mean << "," << weight << "\n";
    }
    digest_out.close();
    XBT_INFO("Wrote latency T-Digest to %s (%zu + %zu centroids)",
             digest_path.c_str(), read_centroids.size(),
             write_centroids.size());
  }

  // Also write computed percentiles for convenience (summary file)
  std::string summary_path =
      aggregate_output_dir + "/client_latency_summary.csv";
  std::ofstream summary_out(summary_path);
  if (summary_out.is_open()) {
    summary_out << "op_type,metric,value\n";

    auto write_summary = [&](const std::string &type,
                             digestible::tdigest<float, uint32_t> &digest) {
      if (digest.size() == 0)
        return;

      auto centroids = digest.get();
      double weighted_sum = 0.0;
      size_t total_weight = 0;
      for (const auto &[mean, weight] : centroids) {
        weighted_sum += mean * weight;
        total_weight += weight;
      }
      double avg = (total_weight > 0) ? weighted_sum / total_weight : 0.0;

      summary_out << type << ",count," << digest.size() << "\n";
      summary_out << type << ",min," << digest.min() << "\n";
      summary_out << type << ",max," << digest.max() << "\n";
      summary_out << type << ",avg," << avg << "\n";
      summary_out << type << ",p50," << digest.quantile(50) << "\n";
      summary_out << type << ",p95," << digest.quantile(95) << "\n";
      summary_out << type << ",p99," << digest.quantile(99) << "\n";
      summary_out << type << ",p99.5," << digest.quantile(99.5) << "\n";
    };

    write_summary("read", read_latency_digest);
    write_summary("write", write_latency_digest);

    summary_out.close();
    XBT_INFO("Wrote latency summary to %s", summary_path.c_str());
  }
}
