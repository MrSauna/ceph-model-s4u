#include "OsdActor.hpp"
#include "MonActor.hpp"

XBT_LOG_NEW_DEFAULT_CATEGORY(s4u_ceph_sim_osd, "Messages specific for OsdActors");

Osd::Osd(std::vector<std::string> args)
{
  xbt_assert(args.size() == 1, "The worker expects no argument");

  my_host = simgrid::s4u::this_actor::get_host();
  mailbox = simgrid::s4u::Mailbox::by_name(my_host->get_name());
  mon_mb  = simgrid::s4u::Mailbox::by_name("mon");
}

void Osd::operator()()
{
  XBT_INFO("I was created %s", my_host->get_cname());
  // Subscribe to pgmap
  auto* init_msg = make_message<SubscribeToPGMapChangeMsg>();
  mon_mb->put(init_msg, 0);

  // Wait for pg map
  auto msg = mailbox->get<Message>();
  PGMapNotification notification = std::get<PGMapNotification>(msg->payload);
  pg_map = notification.pg_map;

  // start stuff
  XBT_INFO("hello my pid is=%ld", simgrid::s4u::this_actor::get_pid());
  XBT_INFO("received pg size=%lu", pg_map->size());

  // receive anything
  msg  = mailbox->get<Message>();
  XBT_INFO("got a message from %s", msg->sender.c_str());
  delete msg;

  XBT_INFO("Exiting now.");
}
