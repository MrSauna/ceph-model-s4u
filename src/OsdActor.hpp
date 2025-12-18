#include "CephCommon.hpp"

#include <simgrid/s4u.hpp>
#include <unordered_map>
namespace sg4 = simgrid::s4u;

// activities started by osd that have to be tracked
// enum class ActivityType {
//   READ_BACKFILL_OBJECT, // read from disk before sending to new peers
//   WRITE_OBJECT,         // immediate write
//   READ_OBJECT,          // immediate read
//   SEND_OBJECT,          // immediate write order
//   CLIENT_READ,          // goes throush scheduler
//   CLIENT_WRITE          // goes through scheduler
// };
// => can I do with OpContext only?

class Osd {
  simgrid::s4u::Host *my_host;
  sg4::Actor *my_actor;
  sg4::Disk *disk;
  simgrid::s4u::Mailbox *mb;
  simgrid::s4u::Mailbox *mon_mb;
  PGMap *pgmap;
  int osd_id;

  int max_recovery_threads = 3; // real hdd default
  int used_recovery_threads = 0;

  // I might want to wrap this in a class later
  sg4::ActivitySet activities; // all activities, need to keep track of what
                               // they are seperately
  // std::unordered_map<simgrid::s4u::ActivityPtr, ActivityType>
  // activity_type_map; op context tracking
  int last_op_id = 0;
  std::unordered_map<int,
                     OpContext *>
      op_contexts; // if I need to do this ownership gets annying
                   // => do I need to use pointers again?
  std::unordered_map<sg4::ActivityPtr, OpContext *> op_context_map;

  std::set<PG *> my_primary_pgs;
  std::set<PG *> needs_backfill_pgs;
  PG *backfilling_pg = nullptr;
  std::unordered_map<int, sg4::Mailbox *> peer_osd_mailboxes;

  void send_op(Op *op);
  void process_message(Message *msg);
  void kill_self();
  void maybe_reserve_backfill();
  void maybe_schedule_object_backfill();
  void process_finished_activity(sg4::ActivityPtr activity);
  void main_loop();
  void on_pgmap_change();
  void on_osd_op_message(int sender, const OsdOpMsg &osd_op_msg);
  void on_osd_op_ack_message(int sender, const OsdOpAckMsg &msg);
  void advance_backfill_op(OpContext *context, int peer_osd_id);

public:
  explicit Osd(PGMap *pgmap, int osd_id, std::string disk_name);
  void operator()();

  std::string primary_pgs_to_string() const;
};