#include <fstream>

#include <llvm/Support/raw_ostream.h>
#include <llvm/Support/raw_os_ostream.h>
#include <llvm/IR/Instruction.h>
#include <llvm/IR/Instructions.h>

#include "VCEvent.h"
#include "VCIID.h"
#include "VCAnnotation.h"
#include "VCBasis.h"
#include "VCHappensBeforeGraph.h"

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

void VCBasis::dump() const
{
  llvm::errs() << "Basis {\n";
  for (auto& process : processes) {
    llvm::errs() << process[0]->cpid << ":\n";
    for (const VCEvent *ev : process) {
      llvm::errs() << "  " << *ev << "\n";
    }
  }
}

void VCHappensBeforeGraph::to_dot(const char *edge_params) const
{
  llvm::errs() << "digraph {\n";
  for (auto& it : nodes) {
    llvm::errs() << "NODE" << it.second << " [label=\"";
    if (!it.second->event)
      llvm::errs() << "init";
    else
      llvm::errs() << *it.second->event;

    llvm::errs() << "\"]\n";

  }

  llvm::errs() << "\n";
  for (auto& it : nodes) {
    // succ is a pointer
    for (auto succ : it.second->successors) {
      llvm::errs() << "NODE" << it.second << " -> NODE" << succ;
      if (edge_params) {
        llvm::errs() << "[" << edge_params << "]\n";
      } else {
        llvm::errs() << "\n";
      }
    }
  }

  llvm::errs() << "}\n";
}

void VCHappensBeforeGraph::dump() const
{
  llvm::errs() << " --- annotation:\n";
  annotation.dump();
  llvm::errs() << " --- basis:\n";
  basis.dump();
}
