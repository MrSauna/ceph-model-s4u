#include "CephActor.hpp"
#include <fstream>
#include <map>
#include <mutex>
#include <xbt/random.hpp>

// T-Digest for streaming percentile calculation
#include "digestible/digestible.h"

// Per-second throughput aggregation bucket
struct ThroughputBucket {
  int64_t read_ops = 0;
  int64_t read_bytes = 0;
  int64_t write_ops = 0;
  int64_t write_bytes = 0;
};

class Client : public CephActor {
  int max_concurrent_reads;
  int max_concurrent_writes;
  int op_size;
  double read_bandwidth;
  double write_bandwidth;
  double read_tokens = 0;
  double write_tokens = 0;
  double last_token_time = 0;
  int in_flight_reads = 0;
  int in_flight_writes = 0;
  bool shutting_down = false;

  long tasks_sent = 0;
  long tasks_acked = 0;

  simgrid::xbt::random::XbtRandom random;

  void process_message(Message *msg) override;
  void gen_op(OpType type);
  void on_osd_op_ack_message(int sender, const OsdOpAckMsg &msg);
  void on_finished_activity(sg4::ActivityPtr activity) override;
  std::optional<double> make_progress() override;

  int last_op_pg = 0;

  // Legacy per-op output (optional)
  static std::ofstream metrics_stream;
  static std::mutex metrics_mutex;

  // Aggregated metrics
  static std::map<int, ThroughputBucket> throughput_buckets;
  static digestible::tdigest<float, uint32_t> read_latency_digest;
  static digestible::tdigest<float, uint32_t> write_latency_digest;
  static bool aggregate_mode;

public:
  explicit Client(PGMap *pgmap, int client_id, int read_queue, int write_queue,
                  int op_size, double read_bandwidth = 0,
                  double write_bandwidth = 0);
  void operator()();

  // Legacy per-op CSV output
  static void set_metrics_output(const std::string &filename);

  // Aggregated output (memory-efficient)
  static void set_aggregate_output(const std::string &output_dir);
  static void write_aggregated_metrics();
};
