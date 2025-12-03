#include "CephCommon.hpp"

#include <simgrid/s4u.hpp>

class Mon
{ 
    simgrid::s4u::Mailbox* mailbox;
    PGMap pg_map;
    std::vector<simgrid::s4u::Mailbox*> osds;

    void process_message(Message* msg);
    bool is_cluster_balanced();
    void kill_all_osds();
    static PGMap create_initial_pg_map(std::vector<std::string> args);

public:
    explicit Mon(std::vector<std::string> args);
    void operator()();
    void on_subscribe_pgmap_change(const std::string& sender, const SubscribeToPGMapChangeMsg& payload);
};
