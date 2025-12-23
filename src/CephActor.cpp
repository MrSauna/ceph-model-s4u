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
    make_progress();

    // wait any activity
    sg4::ActivityPtr finished = activities.wait_any();

    // new message arrived
    if (finished == listener) {
      process_message(message);

      // restart listener
      listener = mb->get_async(&message);
      activities.push(listener);
      continue;
    }

    // something else finished
    else {
      process_finished_activity(finished);
    }
  }
}
