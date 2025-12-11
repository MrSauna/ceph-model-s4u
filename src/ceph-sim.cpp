#include "CLI11.hpp"
#include "ClientActor.hpp"
#include "MonActor.hpp"
#include "OsdActor.hpp"
#include <simgrid/s4u.hpp>

#define READ_BANDWIDTH (120 * 1024 * 1024)
#define WRITE_BANDWIDTH (80 * 1024 * 1024)
XBT_LOG_NEW_DEFAULT_CATEGORY(s4u_ceph_sim, "Messages specific for ceph-sim");

int main(int argc, char *argv[]) {
  simgrid::s4u::Engine e(&argc, argv);
  xbt_assert(argc > 2, "Usage: %s platform_file deployment_file\n", argv[0]);

  /* Register the functions representing the actors */
  // e.register_actor<Mon>("mon");
  // e.register_actor<Osd>("osd");
  // e.register_actor<Client>("client");

  /* Load the platform description and then deploy the application */
  e.load_platform(argv[1]);

  // initialize pg map
  int pool_id = 1;                                         // hardcoded for now
  PGMap *pgmap = new PGMap(pool_id, std::string{argv[3]}); // acting
  PGMap up_pgmap(pool_id, std::string{argv[4]});           // up
  for (size_t i = 0; i < pgmap->size(); i++) {
    std::vector<int> up_members = up_pgmap.get_pg(i)->get_up_ids();
    pgmap->get_pg(i)->init_up(up_members);
  }
  pgmap->init_primary_osd_to_pg_index();
  XBT_INFO("%s", pgmap->to_string().c_str());

  // do deployment
  // e.load_deployment(argv[2]);
  std::vector<sg4::Host *> hosts = e.get_all_hosts();
  for (auto host : hosts) {
    std::string host_name = host->get_name();
    if (host_name.find("osd.") != std::string::npos) {
      size_t pos = host_name.find(".");
      int osd_id = std::stoi(host_name.substr(pos + 1));
      e.add_actor(host_name, host, [pgmap, osd_id]() {
        Osd osd(pgmap, osd_id);
        osd();
      });
    } else if (host_name == "mon") {
    } else {
      XBT_INFO("Unknown host %s in platform file", host_name.c_str());
    }
  }

  // temp add disks
  for (auto host : e.get_all_hosts()) {
    host->add_disk("disk0", READ_BANDWIDTH, WRITE_BANDWIDTH);
  }

  /* Run the simulation */
  e.run();

  XBT_INFO("Simulation is over");

  return 0;
}