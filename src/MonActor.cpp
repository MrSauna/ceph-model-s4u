#include "MonActor.hpp"
#include "CephCommon.hpp"
#include <typeinfo>

XBT_LOG_NEW_DEFAULT_CATEGORY(s4u_ceph_sim_mon,
                             "Messages specific for MonActors");

Mon::Mon(PGMap *pgmap, std::vector<std::string> client_names)
    : pgmap(pgmap), client_names(client_names) {

  XBT_INFO("Initial PGMap\n%s", pgmap->to_string().c_str());
  XBT_INFO("Initial map's reverse index\n%s",
           pgmap->primary_osds_to_pgs_string().c_str());

  mailbox = simgrid::s4u::Mailbox::by_name("mon");

  // dynamically discover osds
  std::vector<int> osds = pgmap->get_osds();
  for (int i = 0; i < osds.size(); i++) {
    osd_mailboxes[i] =
        simgrid::s4u::Mailbox::by_name("osd." + std::to_string(osds[i]));
  }

  // XBT_INFO("Got %zu osds", osds.size());
  sg4::Engine *engine = sg4::Engine::get_instance();
  engine->get_filtered_actors([](const sg4::ActorPtr &actor) {
    return actor->get_name().find("osd.") != std::string::npos;
  });
}

void Mon::main_loop() {
  XBT_INFO("Monitor entered main loop");
  while (true) {
    // subscriptions and migration updates
    if (!mailbox->empty()) {
      Message *msg = mailbox->get<Message>();
      process_message(msg);
    }
    // stop condition for simulation
    if (!pgmap->needs_backfill()) {
      XBT_INFO("Cluster does not need backfill, exiting main loop");
      break;
    }
    simgrid::s4u::this_actor::sleep_for(0.1);
  }
  XBT_INFO("Monitor exiting main loop");
}

void Mon::on_pgmap_change(int pg_id) {
  XBT_INFO("PG %d changed", pg_id);
  PG *pg = pgmap->get_pg(pg_id);
  xbt_assert(!pg->needs_backfill(),
             "currently only pgs that have finished backfilling should change");

  // collect everybody that should be notified
  auto up = pg->get_up_ids();
  auto acting = pg->get_acting_ids();
  // create a set of them
  std::set<int> notify(up.begin(), up.end());
  notify.insert(acting.begin(), acting.end());

  // set acting = up
  pg->update_acting(pg->get_up_ids()); // runs prune_shards

  // update the index stuff
  pgmap->update_primary_osd_to_pg_index();

  // notify
  for (int osd_id : notify) {
    auto mb = osd_mailboxes[osd_id];
    Message *msg = make_message<PGMapNotification>(pgmap);
    mb->put(msg, 0);
  }
}

void Mon::process_message(Message *msg) {

  std::visit(
      overloaded{[&](const KillMsg &) {
                   xbt_die("Monitor received unexpected KillMsg");
                 },
                 [&](const SubscribeToPGMapChangeMsg &) {
                   on_subscribe_pgmap_change(
                       msg->sender,
                       std::get<SubscribeToPGMapChangeMsg>(msg->payload));
                 },
                 [&](const KillAckMsg &) {
                   XBT_WARN("Monitor received unexpected KillAckMsg");
                 },
                 [&](const PGNotification &payload) {
                   on_pgmap_change(payload.pg_id);
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
  auto response_msg = make_message<PGMapNotification>(pgmap);
  auto return_mb = simgrid::s4u::Mailbox::by_name(sender);
  return_mb->put(response_msg, 0);
}

void Mon::kill_all_osds() {
  XBT_INFO("Killing all clients");
  sg4::ActivitySet kill_activities;

  // 1. Kill Clients
  for (auto client_name : client_names) {
    auto mb = simgrid::s4u::Mailbox::by_name(client_name);
    Message *msg = make_message<KillMsg>();
    auto a = mb->put_async(msg, 0);
    kill_activities.push(a);
  }
  kill_activities.wait_all(); // Wait for KillMsg to be put into mailboxes

  // 2. Wait for clients to finish draining
  int alive_clients = client_names.size();
  while (alive_clients > 0) {
    Message *msg = mailbox->get<Message>();
    std::visit(
        overloaded{[&](const KillAckMsg &) {
                     alive_clients--;
                     XBT_INFO("Client %s finished", msg->sender.c_str());
                   },
                   [&](const auto &) {
                     XBT_WARN(
                         "Monitor received unexpected message while waiting "
                         "for clients to die: %s",
                         msg->sender.c_str());
                   }},
        msg->payload);
    delete msg;
  }
  XBT_INFO("All clients dead");

  // 3. Kill OSDs
  XBT_INFO("Sending kill to all osds");
  for (auto mb_tuple : osd_mailboxes) {
    auto mb = mb_tuple.second;
    Message *msg = make_message<KillMsg>();
    auto a = mb->put_async(msg, 0);
    kill_activities.push(a);
  }
  kill_activities.wait_all();
}

void Mon::kill_self() {
  XBT_INFO("Monitor is killing itself");
  simgrid::s4u::this_actor::exit();
}

void Mon::operator()() {
  XBT_INFO("Monitor was created");
  main_loop();
  kill_all_osds();
  kill_self();
}
