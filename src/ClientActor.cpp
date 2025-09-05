#include "ClientActor.hpp"

XBT_LOG_NEW_DEFAULT_CATEGORY(s4u_ceph_sim_client, "Messages specific for ClientActors");

Client::Client(std::vector<std::string> args) {
  xbt_assert(args.size() > 4, "The master function expects at least 3 arguments");

  tasks_count        = std::stol(args[1]);
  compute_cost       = std::stod(args[2]);
  communication_cost = std::stod(args[3]);

  for (unsigned int i = 4; i < args.size(); i++)
    workers.push_back(simgrid::s4u::Mailbox::by_name(args[i]));

  XBT_INFO("Got %zu workers and %ld tasks to process", workers.size(), tasks_count);
}

void Client::operator()()
{
  for (int i = 0; i < tasks_count; i++) { /* For each task to be executed: */
    /* - Select a worker in a round-robin way */
    simgrid::s4u::Mailbox* mailbox = workers[i % workers.size()];

    /* - Send the computation cost to that worker */
    XBT_INFO("Sending task %d of %ld to mailbox '%s'", i, tasks_count, mailbox->get_cname());
    mailbox->put(new double(compute_cost), communication_cost);
  }

  XBT_INFO("All tasks have been dispatched. Request all workers to stop.");
  for (unsigned int i = 0; i < workers.size(); i++) {
    /* The workers stop when receiving a negative compute_cost */
    simgrid::s4u::Mailbox* mailbox = workers[i % workers.size()];

    mailbox->put(new double(-1.0), 0);
  }
}
