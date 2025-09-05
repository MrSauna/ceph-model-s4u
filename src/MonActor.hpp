#include "CephCommon.hpp"

#include <simgrid/s4u.hpp>

class Mon
{ 
    simgrid::s4u::Mailbox* mailbox;
    std::vector<PGMap> pg_history;
    unsigned int pool_id;
    std::vector<simgrid::s4u::Mailbox*> osds;

    void process_message(Message* msg);
    bool is_cluster_balanced();
    void kill_all_osds();

public:
    explicit Mon(std::vector<std::string> args);
    void operator()();
};
