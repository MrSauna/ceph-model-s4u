#include "CephActor.hpp"

class Client : public CephActor {
  long tasks_count;
  double compute_cost;
  double communication_cost;
  std::vector<simgrid::s4u::Mailbox *> workers; // These are OSD mailboxes?

  long tasks_sent = 0;
  long tasks_acked = 0;

  void process_message(Message *msg) override;
  void process_finished_activity(sg4::ActivityPtr activity) override;
  void make_progress() override;

  int last_op_pg = 0;

public:
  explicit Client(PGMap *pgmap, int client_id);
  void operator()();
};
