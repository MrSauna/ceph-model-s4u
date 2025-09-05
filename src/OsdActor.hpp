#include "CephCommon.hpp"

#include <simgrid/s4u.hpp>

class Osd {
    simgrid::s4u::Host* my_host;
    simgrid::s4u::Mailbox* mailbox;
    simgrid::s4u::Mailbox* mon_mb;
    PGMap* pg_map;
public:
    explicit Osd(std::vector<std::string> args);
    void operator()();
};