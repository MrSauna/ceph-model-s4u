#include "CephActor.hpp"

#include <memory>
#include <simgrid/s4u.hpp>
#include <unordered_map>

// dmclock includes
#include "dmclock_server.h"
#include "dmclock_util.h"

namespace sg4 = simgrid::s4u;
namespace dmc = crimson::dmclock;

class Osd : public CephActor {
  sg4::Disk *disk;

  int max_recovery_threads = 3; // real hdd default
  int used_recovery_threads = 0;

  std::set<PG *> my_primary_pgs;
  std::set<PG *> needs_backfill_pgs;
  bool backfill_reservation_local = false;
  bool backfill_reservation_remote = false;
  std::set<int> backfill_reservation_remote_pending;
  PG *backfilling_pg = nullptr;
  std::unordered_map<int, sg4::Mailbox *> peer_osd_mailboxes;

  // Scheduler types
  using ClientId = int;
  using Queue = dmc::PushPriorityQueue<ClientId, OpContext>;

  std::unique_ptr<Queue> queue;

  // QoS Clients
  // 1 = User Generic (Reservation 50%, Weight 1, Unlimited)
  // 2 = Backfill (Reservation 0%, Weight 1, Limit 90%)
  static constexpr ClientId CLIENT_ID_USER = 1;
  static constexpr ClientId CLIENT_ID_BACKFILL = 2;

  // Helper to init QoS
  void init_scheduler();

  void send_op(Op *op);
  void process_message(Message *msg) override;
  void maybe_reserve_backfill();
  void maybe_schedule_object_backfill();
  void process_finished_activity(sg4::ActivityPtr activity) override;
  void make_progress() override;

  void on_pgmap_change();
  void on_osd_op_message(int sender, const OsdOpMsg &osd_op_msg);
  void on_osd_op_ack_message(int sender, const OsdOpAckMsg &msg);
  void advance_backfill_op(OpContext *context, int peer_osd_id);
  void advance_write_op(int op_id);
  void advance_read_op(int op_id);

public:
  explicit Osd(PGMap *pgmap, int osd_id, std::string disk_name);
  void operator()();

  std::string primary_pgs_to_string() const;
};