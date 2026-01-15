#pragma once
#include "CephActor.hpp"

#include <memory>
#include <simgrid/s4u.hpp>
#include <unordered_map>

#include "dmclock_server.h"

namespace sg4 = simgrid::s4u;
namespace dmc = crimson::dmclock;

enum struct SchedulerProfile {
  BALANCED,
  HIGH_CLIENT_OPS,
  HIGH_RECOVERY_OPS,
};

struct NoOpJob {
  // Matches constructor signature of RunEvery from dmclock
  template <typename... Args> NoOpJob(Args &&...args) {}

  // Method called by PriorityQueueBase to update timer
  void try_update(std::chrono::milliseconds) {}
};

#include <xbt/random.hpp>

class Osd : public CephActor {
  sg4::Disk *disk;
  simgrid::xbt::random::XbtRandom random;

  int max_recovery_threads = 3; // real hdd default
  int used_recovery_threads = 0;

  double base_cost;

  std::set<PG *> my_primary_pgs;
  std::set<PG *> needs_backfill_pgs;
  bool backfill_reservation_local = false;
  std::set<int> backfill_reservation_remote;
  std::set<int> backfill_reservation_remote_pending;
  bool pending_retry = false;
  PG *backfilling_pg = nullptr;
  std::unordered_map<int, sg4::Mailbox *> peer_osd_mailboxes;

  // Scheduler types
  using ClientId = int;
  using Queue =
      dmc::PullPriorityQueue<ClientId, OpContext, false, false, 2, NoOpJob>;

  std::unique_ptr<dmc::ClientInfo> user_info;
  std::unique_ptr<dmc::ClientInfo> backfill_info;

  std::unique_ptr<Queue> queue;

  // QoS Clients
  static constexpr ClientId CLIENT_ID_USER = 1;
  static constexpr ClientId CLIENT_ID_BACKFILL = 2;

  // Helper to init QoS
  void init_scheduler(double iops, SchedulerProfile profile);

  void send_op(Op *op);
  void process_message(Message *msg) override;
  void maybe_reserve_backfill();
  void maybe_schedule_object_backfill();
  void process_finished_activity(sg4::ActivityPtr activity) override;
  void make_progress() override;

  void on_pgmap_change();
  void on_osd_op_message(int sender, const OsdOpMsg &osd_op_msg);
  void on_osd_op_ack_message(int sender, const OsdOpAckMsg &msg);
  void on_backfill_reservation_message(int sender,
                                       const BackfillReservationMsg &msg);
  void opcontext_dispatch(OpContext *context);
  void advance_backfill_op(OpContext *context, int peer_osd_id);

public:
  explicit Osd(PGMap *pgmap, int osd_id, std::string disk_name, double iops,
               SchedulerProfile profile);
  void operator()();

  std::string primary_pgs_to_string() const;
};