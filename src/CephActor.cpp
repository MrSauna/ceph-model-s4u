#include "CephActor.hpp"

XBT_LOG_NEW_DEFAULT_CATEGORY(s4u_ceph_actor,
                             "Messages specific for CephActors");

CephActor::CephActor(int id, PGMap *pgmap) : id(id), pgmap(pgmap) {
  my_host = simgrid::s4u::this_actor::get_host();
  my_actor = sg4::Actor::self();
  mb = simgrid::s4u::Mailbox::by_name(my_actor->get_name());
  mon_mb = simgrid::s4u::Mailbox::by_name("mon");
}

void CephActor::kill_self() {
  XBT_INFO("%s is killing itself", my_host->get_name().c_str());

  simgrid::s4u::this_actor::exit();
}

void CephActor::main_loop() {
  // start off listener
  Message *message;
  listener = mb->get_async(&message);
  activities.push(listener);

  while (true) {
    std::optional<double> next_event = make_progress();
    double timeout = 1000000.0; // Default large timeout (infinity)

    if (next_event.has_value()) {
      timeout = next_event.value() - sg4::Engine::get_clock();
      if (timeout < 0)
        timeout = 0;
    }

    // wait any activity
    sg4::ActivityPtr finished = nullptr;
    try {
      if (next_event.has_value()) {
        finished = activities.wait_any_for(timeout);
      } else {
        finished = activities.wait_any();
      }
    } catch (const simgrid::TimeoutException &) {
      // Timeout occurred, just loop back to make_progress
      continue;
    }

    // new message arrived
    if (finished == listener) {
      process_message(message);

      // restart listener
      listener = mb->get_async(&message);
      activities.push(listener);
      continue;
    }

    // something else finished
    else if (finished != nullptr) {
      on_finished_activity(finished);
    }

    // timeout
    else {
      continue;
    }
  }
}
