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

  return out;
}

void VCIID::dump() const {
	llvm::errs() << *this << "\n";
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

  llvm::errs() << "digraph {\n";
	for (auto it = nodes_begin(); it != nodes_end(); ++it) {
    const Node *nd = *it;
		llvm::errs() << "NODE"
								 << nd->getProcessID() * 100000 + nd->getEventID()
								 << " [label=\"";
		if (!nd->getEvent())
      llvm::errs() << "init";
    else
      llvm::errs() << *nd->getEvent();
		llvm::errs() << "\"]\n";
	}

  llvm::errs() << "\n";
  for (unsigned ti=0; ti<succ_original.size(); ++ti)
		for (unsigned tj=0; tj<succ_original[ti].size(); ++tj) {
      int curval = -1;
			for (unsigned tiev=0; tiev<succ_original[ti][tj].size(); ++tiev) {
				
				if (succ_original[ti][tj][tiev] == INT_MAX)
					break;
				
				if (succ_original[ti][tj][tiev] > curval) {
          if (curval != -1) {
						llvm::errs() << "NODE"
												 << ti * 100000 + tiev - 1
												 << " -> NODE"
												 << tj * 100000 + curval;
						if (edge_params) {
							llvm::errs() << "[" << edge_params << "]\n";
						} else {
							llvm::errs() << "\n";
						}
					}
					curval = succ_original[ti][tj][tiev];
				}
				
			}
		}

  llvm::errs() << "}\n";

}
