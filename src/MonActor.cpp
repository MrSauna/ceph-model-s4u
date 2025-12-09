#include "MonActor.hpp"
#include <filesystem>
#include <typeinfo>

namespace fs = std::filesystem;

XBT_LOG_NEW_DEFAULT_CATEGORY(s4u_ceph_sim_mon,
                             "Messages specific for MonActors");

PGMap Mon::create_initial_pg_map(std::vector<std::string> args) {
  xbt_assert(args.size() == 4, "The master function expects 4 argument");
  xbt_assert(fs::is_regular_file(args[2]), "args[2] not a file path");
  xbt_assert(fs::is_regular_file(args[3]), "args[3] not a file path");

  // Initial pg map
  unsigned int pool_id = std::stoul(args[1]);
  PGMap initial_pgmap(std::string(args[2]), pool_id); // acting
  PGMap up_pgmap(std::string(args[3]), pool_id);      // up

  // todo combine maps
  for (size_t i = 0; i < initial_pgmap.size(); i++) {
    std::vector<unsigned int> up_members = up_pgmap.get_pg(i)->get_up_ids();
    initial_pgmap.get_pg(i)->set_up(up_members);
  }
  initial_pgmap.rebuild_primary_osd_to_pg_index();
  return initial_pgmap;
}

Mon::Mon(std::vector<std::string> args) : pg_map(create_initial_pg_map(args)) {

  XBT_INFO("Initial PGMap\n%s", pg_map.to_string().c_str());
  XBT_INFO("Initial map's reverse index\n%s",
           pg_map.primary_osds_to_pgs_string().c_str());

  mailbox = simgrid::s4u::Mailbox::by_name("mon");

  for (unsigned int i = 0; i < 3; i++)
    osds.push_back(simgrid::s4u::Mailbox::by_name("osd." + std::to_string(i)));

  XBT_INFO("Got %zu osds", osds.size());
}

void Mon::process_message(Message *msg) {

  std::visit(
      overloaded{[&](const KillMsg &) {
                   xbt_die("Monitor received unexpected KillMsg");
                 },
                 [&](const SubscribeToPGMapChangeMsg &) {
                   this->on_subscribe_pgmap_change(
                       msg->sender,
                       std::get<SubscribeToPGMapChangeMsg>(msg->payload));
                 },
                 [&](const auto &unknown_payload) {
                   xbt_die("Monitor received unexpected message type: %s",
                           typeid(unknown_payload).name());
                 }},
      msg->payload);

  delete msg;
}

void Mon::on_subscribe_pgmap_change(const std::string &sender,
                                    const SubscribeToPGMapChangeMsg &payload) {
  XBT_INFO("%s subscribed to pg map changes", sender.c_str());
  PGMap *pg_map_ptr = &pg_map;
  auto response_msg = make_message<PGMapNotification>(pg_map_ptr);
  auto return_mb = simgrid::s4u::Mailbox::by_name(sender);
  return_mb->put(response_msg, 0);
}

// TODO implement actual logic
bool Mon::is_cluster_balanced() {

  if (simgrid::s4u::Engine::get_clock() >= 100)
    return true;
  return false;
}

void Mon::kill_all_osds() {
  for (auto osd_mb : osds) {
    Message *msg = make_message<KillMsg>();
    osd_mb->put(msg, 0);
  }
}

void Mon::operator()() {
  XBT_INFO("Monitor was created");
  // main loop?
  while (true) {
    // subscriptions and migration updates
    if (!mailbox->empty()) {
      Message *msg = mailbox->get<Message>();
      process_message(msg);
    }
    // stop condition for simulation
    if (is_cluster_balanced()) {
      break;
    }
    simgrid::s4u::this_actor::sleep_for(0.1);
  }
  XBT_INFO("Sending kill to all osds");
  kill_all_osds();
  XBT_INFO("ALL OSDs exited");
  XBT_INFO("Simulation over");
}
