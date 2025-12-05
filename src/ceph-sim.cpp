#include "CLI11.hpp"
#include "ClientActor.hpp"
#include "MonActor.hpp"
#include "OsdActor.hpp"
#include <simgrid/s4u.hpp>

XBT_LOG_NEW_DEFAULT_CATEGORY(s4u_ceph_sim, "Messages specific for ceph-sim");

int main(int argc, char *argv[]) {
  simgrid::s4u::Engine e(&argc, argv);
  xbt_assert(argc > 2, "Usage: %s platform_file deployment_file\n", argv[0]);

  /* Register the functions representing the actors */
  e.register_actor<Mon>("mon");
  e.register_actor<Osd>("osd");
  e.register_actor<Client>("client");

  /* Load the platform description and then deploy the application */
  e.load_platform(argv[1]);
  e.load_deployment(argv[2]);

  /* Run the simulation */
  e.run();

  XBT_INFO("Simulation is over");

  return 0;
}