#include <fstream>

#include <llvm/Support/raw_ostream.h>
#include <llvm/Support/raw_os_ostream.h>
#include <llvm/IR/Instruction.h>
#include <llvm/IR/Instructions.h>

#include "VCEvent.h"
#include "VCIID.h"
#include "VCAnnotation.h"
#include "VCGraphVclock.h"

llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const VCEvent& ev)
{
  out << ev.cpid << "-" << ev.instruction_order;
  if (ev.instruction)
    out << *ev.instruction;
  else
    out << " <null>";

  return out;
}

void VCEvent::dump() const {
  llvm::errs() << *this << "\n";
}

llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const VCIID& iid)
{
  out << iid.cpid << "-" << iid.instruction_order;
  if (iid.instruction)
    out << *iid.instruction;
  else
    out << " <null>";

  return out;
}

void VCIID::dump() const {
	llvm::errs() << *this << "\n";
}

llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const std::pair<int,CPid> annotvalue)
{
  out << "v" << annotvalue.first << "-" << annotvalue.second;

  return out;
}

llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const VCAnnotation& annot) {
  out << "Annotation+ {\n";
  for (auto& pr : annot) {
    out << "("  << pr.first << ",\n"
                << " "  << pr.second << ")\n";
  }
  out << "}\n";

  return out;
}

void VCAnnotation::dump() const {
  llvm::errs() << *this;
}

void VCGraphVclock::to_dot(const char *edge_params) const
{
	/*
  llvm::errs() << "digraph {\n";
  for (auto& it : *this) {
    llvm::errs() << "NODE" << &it << " [label=\"";
    if (!it.getEvent())
      llvm::errs() << "init";
    else
      llvm::errs() << *it.getEvent();

    llvm::errs() << "\"]\n";

  }

  llvm::errs() << "\n";
  for (auto& it : *this) {
    for (unsigned idx = 0; idx < it.successors.size(); ++idx) {
      unsigned eid = it.successors[idx];
      // no edge to this thread
      if (eid == INT_MAX)
        continue;
      auto succ = processes[idx][eid];
      llvm::errs() << "NODE" << &it << " -> NODE" << succ;
      if (edge_params) {
        llvm::errs() << "[" << edge_params << "]\n";
      } else {
        llvm::errs() << "\n";
      }
    }
  }

  llvm::errs() << "}\n";
	*/
}

void VCGraphVclock::dump() const
{
	/*
  llvm::errs() << " --- annotation:\n";
  annotation.dump();
  llvm::errs() << " --- HB:\n";
  happens_before.dump();
  llvm::errs() << " --- processes:\n";
  unsigned idx = 0;
  llvm::errs() << "[";
  for (int x : initial_node->successors) {
    if (x == INT_MAX)
      llvm::errs() << "x, ";
    else
      llvm::errs() << x << ", ";
  }
  llvm::errs() << "] ";
  llvm::errs() << " -- init\n";

  for (auto& proc : processes) {
    llvm::errs() << "process " << idx++ << "\n";
    unsigned id = 0;
    for (const Node *nd : proc) {
      llvm::errs() << id++ << " [";
      for (int x : nd->successors) {
        if (x == INT_MAX)
          llvm::errs() << "x, ";
        else
          llvm::errs() << x << ", ";
      }
      llvm::errs() << "] ";
      llvm::errs() << *nd->getEvent() << "\n";
   }
  }
  llvm::errs() << " -----------\n";
	*/
}
