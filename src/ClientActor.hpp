#include "CephActor.hpp"
#include <fstream>
#include <mutex>
#include <xbt/random.hpp>

class Client : public CephActor {
  int max_concurrent_reads;
  int max_concurrent_writes;
  int in_flight_reads = 0;
  int in_flight_writes = 0;
  bool shutting_down = false;

  long tasks_sent = 0;
  long tasks_acked = 0;

  simgrid::xbt::random::XbtRandom random;

  void process_message(Message *msg) override;
  void gen_op(OpType type);
  void on_osd_op_ack_message(int sender, const OsdOpAckMsg &msg);
  void process_finished_activity(sg4::ActivityPtr activity) override;
  std::optional<double> make_progress() override;

  int last_op_pg = 0;

  static std::ofstream metrics_stream;
  static std::mutex metrics_mutex;

public:
  explicit Client(PGMap *pgmap, int client_id, int read_queue, int write_queue);
  void operator()();

  static void set_metrics_output(const std::string &filename);
};
