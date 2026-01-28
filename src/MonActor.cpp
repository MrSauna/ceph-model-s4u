#include "MonActor.hpp"
#include "CephCommon.hpp"
#include <typeinfo>

XBT_LOG_NEW_DEFAULT_CATEGORY(s4u_ceph_sim_mon,
                             "Messages specific for MonActors");

std::ofstream Mon::metrics_stream;
std::mutex Mon::metrics_mutex;

void Mon::set_metrics_output(const std::string &filename) {
  metrics_stream.open(filename);
  if (metrics_stream.is_open()) {
    metrics_stream << "time,active_clean,backfill,backfill_wait\n";
  } else {
    XBT_ERROR("Failed to open mon metrics file: %s", filename.c_str());
  }
}

Mon::Mon(PGMap *pgmap, std::vector<std::string> client_names,
         long start_up_delay, long shut_down_delay)
    : CephActor(-1, pgmap), client_names(client_names),
      start_up_delay(start_up_delay), shut_down_delay(shut_down_delay) {

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

void Mon::on_pgmap_change(int pg_id) {
  PG *pg = pgmap->get_pg(pg_id);
  XBT_INFO("PG %d state changed to %d", pg_id, pg->get_state());

  {
    std::lock_guard<std::mutex> lock(metrics_mutex);
    if (metrics_stream.is_open()) {
      int active_clean = 0;
      int backfill = 0;
      int backfill_wait = 0;
      for (size_t i = 0; i < pgmap->size(); i++) {
        switch (pgmap->get_pg(i)->get_state()) {
        case PGState::ACTIVE_CLEAN:
          active_clean++;
          break;
        case PGState::BACKFILL:
          backfill++;
          break;
        case PGState::BACKFILL_WAIT:
          backfill_wait++;
          break;
        }
      }
      metrics_stream << sg4::Engine::get_clock() << "," << active_clean << ","
                     << backfill << "," << backfill_wait << "\n";
    }
  }

  if (pg->get_state() == PGState::ACTIVE_CLEAN) {
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

  // Check if backfill is complete
  if (!pgmap->needs_backfill()) {
    XBT_INFO("Cluster is balanced");
    dispatch_kill_cluster();
  }
}

void Mon::dispatch_kill_cluster() {
  if (shut_down_delay > 0) {
    sg4::this_actor::sleep_for(shut_down_delay);
    XBT_INFO("Waited %ld seconds before shutting down cluster",
             shut_down_delay);
  }
  XBT_INFO("Shutting down cluster");
  kill_all_osds();
  kill_self();
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

void Mon::on_finished_activity(sg4::ActivityPtr activity) {
  // Mon currently doesn't track background activities in the main loop
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
    Message *msg = mb->get<Message>();
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
  // CephActor::kill_self(); // This calls exit(), but we want to log first as
  // above
  simgrid::s4u::this_actor::exit();
}

void Mon::operator()() {
  XBT_INFO("Monitor was created");

  // start up delay
  sg4::this_actor::sleep_for(start_up_delay);

  // Send PGMapNotification to all OSDs
  for (auto mb_tuple : osd_mailboxes) {
    auto mb = mb_tuple.second;
    Message *msg = make_message<PGMapNotification>(pgmap);
    mb->put(msg, 0);
  }
  if (start_up_delay > 0) {
    XBT_INFO("Monitor waited %ld seconds before sending PGMapNotification",
             start_up_delay);
  }
  // Calls the base class main_loop which is event driven
  CephActor::main_loop();
}
