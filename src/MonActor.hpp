#include "CephCommon.hpp"

#include <simgrid/s4u.hpp>

class Mon {
  simgrid::s4u::Mailbox *mailbox;
  PGMap *pgmap;
  std::map<int, simgrid::s4u::Mailbox *> osd_mailboxes;
  std::vector<std::string> client_names;

  void process_message(Message *msg);
  bool is_cluster_balanced();
  void kill_all_osds();
  void kill_self();
  static PGMap create_initial_pgmap(std::vector<std::string> args);
  void main_loop();
  void on_pgmap_change(int pg_id);

public:
  explicit Mon(PGMap *pgmap, std::vector<std::string> client_names);
  void operator()();
  void on_subscribe_pgmap_change(const std::string &sender,
                                 const SubscribeToPGMapChangeMsg &payload);
};
