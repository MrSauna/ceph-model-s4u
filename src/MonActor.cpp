#include "MonActor.hpp"

#include <filesystem>

namespace fs = std::filesystem;

XBT_LOG_NEW_DEFAULT_CATEGORY(s4u_ceph_sim_mon, "Messages specific for MonActors");

Mon::Mon(std::vector<std::string> args)
{
  xbt_assert(args.size() == 3, "The master function expects 3 argument");
  xbt_assert(fs::is_regular_file(args[1]), "args[1] not a file path");

  // Initial pg map
  pool_id = std::stoul(args[2]);
  PGMap initial_pg_map(std::string(args[1]), pool_id);
  pg_history.push_back(initial_pg_map);
  XBT_INFO("PGMap initialized. Size is %lu", initial_pg_map.size());

  mailbox = simgrid::s4u::Mailbox::by_name("mon");

  for (unsigned int i = 0; i < 3; i++)
    osds.push_back(simgrid::s4u::Mailbox::by_name("osd."+std::to_string(i)));

  XBT_INFO("Got %zu osds", osds.size());

}

void Mon::process_message(Message* msg) {

  std::visit([&](auto&& payload) {
      // 'msg' has the correct, concrete type here! No casting needed.
      using T = std::decay_t<decltype(payload)>;

      // osd should have this, a monitor should never get kill signal
      if constexpr (std::is_same_v<T, KillMsg>) {
          simgrid::s4u::this_actor::exit();
      } else if constexpr (std::is_same_v<T, SubscribeToPGMapChangeMsg>) {
          // 'payload' is a SubscribeToPGMapChangeMsg object
          XBT_INFO("%s subscribed to pg map changes", msg->sender.c_str());
          PGMap* pg_map_ptr = &pg_history.back();
          auto response_msg = make_message<PGMapNotification>(pg_map_ptr);
          auto return_mb = simgrid::s4u::Mailbox::by_name(msg->sender);
          return_mb->put(response_msg, 0);
      } else {
        xbt_die("Monitor received unexpected message");
      }
  }, msg->payload);

  delete msg;
}

// TODO implement actual logic
bool Mon::is_cluster_balanced() {

  if (simgrid::s4u::Engine::get_clock() >= 100)
    return true;
  return false;
}

void Mon::kill_all_osds() {
  for (auto osd_mb : osds) {
    Message* msg = make_message<KillMsg>();
    osd_mb->put(msg, 0);
  }
}

void Mon::operator()()
{
  XBT_INFO("Monitor was created");
  // main loop?
  while (true) {
    // subscriptions and migration updates
    if (!mailbox->empty()) {
      Message* msg = mailbox->get<Message>();
      process_message(msg);
    }
    // stop condition for simulation
    if(is_cluster_balanced()) {
      break;
    }
    simgrid::s4u::this_actor::sleep_for(0.1);
  }
  XBT_INFO("Sending kill to all osds");
  kill_all_osds();
  XBT_INFO("ALL OSDs exited");
  XBT_INFO("Simulation over");

}
