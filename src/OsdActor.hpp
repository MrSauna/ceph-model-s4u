#include "CephActor.hpp"

#include <simgrid/s4u.hpp>
#include <unordered_map>
namespace sg4 = simgrid::s4u;

class Osd : public CephActor {
  sg4::Disk *disk;

  int max_recovery_threads = 3; // real hdd default
  int used_recovery_threads = 0;

  std::set<PG *> my_primary_pgs;
  std::set<PG *> needs_backfill_pgs;
  PG *backfilling_pg = nullptr;
  std::unordered_map<int, sg4::Mailbox *> peer_osd_mailboxes;

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