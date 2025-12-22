#include "CephActor.hpp"

class Client : public CephActor {
  int max_concurrent_ops = 3;
  int in_flight_ops = 0;

  long tasks_sent = 0;
  long tasks_acked = 0;

  void process_message(Message *msg) override;
  void gen_op(OpType type);
  void on_osd_op_ack_message(int sender, const OsdOpAckMsg &msg);
  void process_finished_activity(sg4::ActivityPtr activity) override;
  void make_progress() override;

  int last_op_pg = 0;

public:
  explicit Client(PGMap *pgmap, int client_id);
  void operator()();
};
