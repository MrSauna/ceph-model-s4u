#include <simgrid/s4u.hpp>

class Client
{
    long tasks_count;
    double compute_cost;
    double communication_cost;
    std::vector<simgrid::s4u::Mailbox*> workers;
public:
    explicit Client(std::vector<std::string> args);
    void operator()();
};
