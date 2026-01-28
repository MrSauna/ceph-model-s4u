#include "CephActor.hpp"
#include "CephCommon.hpp"
#include <fstream>
#include <mutex>

#include <simgrid/s4u.hpp>

class Mon : public CephActor {
  std::map<int, simgrid::s4u::Mailbox *> osd_mailboxes;
  std::vector<std::string> client_names;

  void process_message(Message *msg) override;
  void on_finished_activity(sg4::ActivityPtr activity) override;
  bool is_cluster_balanced();
  void kill_all_osds();
  void kill_self() override;
  static PGMap create_initial_pgmap(std::vector<std::string> args);
  void on_pgmap_change(int pg_id);
  void dispatch_kill_cluster();

public:
  explicit Mon(PGMap *pgmap, std::vector<std::string> client_names,
               long start_up_delay, long shut_down_delay);
  void operator()();
  void on_subscribe_pgmap_change(const std::string &sender,
                                 const SubscribeToPGMapChangeMsg &payload);

  static void set_metrics_output(const std::string &filename);

private:
  static std::ofstream metrics_stream;
  static std::mutex metrics_mutex;
  long start_up_delay;
  long shut_down_delay;
};
