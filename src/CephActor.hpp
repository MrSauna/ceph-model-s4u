#pragma once

#include "CephCommon.hpp"
#include <simgrid/s4u.hpp>
#include <unordered_map>

namespace sg4 = simgrid::s4u;

class CephActor {
protected:
  simgrid::s4u::Host *my_host;
  sg4::Actor *my_actor;
  simgrid::s4u::Mailbox *mb;
  simgrid::s4u::Mailbox *mon_mb;
  sg4::ActivitySet activities;
  sg4::CommPtr listener;

  std::unordered_map<sg4::ActivityPtr, OpContext *>
      op_context_map; // activity to op context map
  std::unordered_map<int, OpContext *> op_contexts; // op id to op context map
  int last_op_id = 0;

  int id;
  PGMap *pgmap;

public:
  CephActor(int id, PGMap *pgmap);
  virtual ~CephActor() = default;
  virtual void main_loop();
  virtual void process_message(Message *msg) = 0;
  virtual void process_finished_activity(sg4::ActivityPtr activity) = 0;
  virtual void kill_self();
  virtual void make_progress() {}
};
