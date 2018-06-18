#include <fstream>

#include <llvm/Support/raw_ostream.h>
#include <llvm/Support/raw_os_ostream.h>
#include <llvm/IR/Instruction.h>
#include <llvm/IR/Instructions.h>

#include "VCEvent.h"
#include "VCIID.h"
#include "VCAnnotation.h"
#include "VCGraphVclock.h"

void removeSubstrings(std::string& s, std::string&& p) { 
  std::string::size_type n = p.length();
  for (std::string::size_type i = s.find(p);
      i != std::string::npos;
      i = s.find(p))
      s.erase(i, n);
}

llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const VCEvent& ev)
{
  //out << "ins-" << ev.instruction_order << "  ";
  
  switch(ev.kind) {
   case VCEvent::Kind::DUMMY : {
		 out << "<>";
	 } break;
			
   case VCEvent::Kind::LOAD  : {
		 out << "read";
		 assert(ev.instruction);
		 std::string str;
		 llvm::raw_string_ostream rso(str);
		 ev.instruction->print(rso);
		 std::string::size_type i1 = str.find("@");
		 if (i1 != std::string::npos) {
			 std::string::size_type i2 = str.find(",", i1);
			 out << " " << str.substr(i1, i2 - i1);
		 } else out << str;
	 } break;
      
   case VCEvent::Kind::STORE  : {
		 out << "write " << ev.value;
		 assert(ev.instruction);
		 std::string str;
		 llvm::raw_string_ostream rso(str);
		 ev.instruction->print(rso);
		 std::string::size_type i1 = str.find("@");
		 if (i1 != std::string::npos) {
			 std::string::size_type i2 = str.find(",", i1);
			 out << " " << str.substr(i1, i2 - i1);
			 //out << " v-" << ev.value;
		 } else out << str;
	 } break;
      
   case VCEvent::Kind::SPAWN  : {
		 out << "spawn " << ev.childs_cpid;
     assert(ev.instruction);
	 } break;
      
   case VCEvent::Kind::JOIN  : {
     out << "join " << ev.childs_cpid;
     assert(ev.instruction);		 
   } break;

   case VCEvent::Kind::M_INIT  : {
     out << "mutexinit";
     assert(ev.instruction);
     out << *ev.instruction;
   } break;
   
   case VCEvent::Kind::M_DESTROY  : {
     out << "mutexdestroy";
     assert(ev.instruction);
     out << *ev.instruction;
   } break;
      
   case VCEvent::Kind::M_LOCK  : {
     out << "lock";
     assert(ev.instruction);
     out << *ev.instruction;
   } break;
   
   case VCEvent::Kind::M_UNLOCK  : {
      out << "unlock";
      assert(ev.instruction);
      out << *ev.instruction;
   } break;
   default :
      out << "unknown";
	}
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
                << " "  << pr.second.first;
		char loc = (pr.second.second == VCAnnotation::Loc::LOCAL)?'L':
			((pr.second.second == VCAnnotation::Loc::REMOTE)?'R':'A');
		out << loc << ")\n";
  }
  out << "}\n";

  return out;
}

void VCAnnotation::dump() const {
  llvm::errs() << *this;
}

void VCGraphVclock::dump_po(const PartialOrder& po) const
{
  ThreadPairsVclocks& succ = *(po.first);
	ThreadPairsVclocks& pred = *(po.second);

	llvm::errs() << "\n";
  for (unsigned ti=0; ti<succ.size(); ++ti)
		for (unsigned tj=0; tj<succ[ti].size(); ++tj) {
			llvm::errs() << "succ " << ti << "->" << tj << " size: " << succ[ti][tj].size() << "   content:";
			for (unsigned tiev=0; tiev<succ[ti][tj].size(); ++tiev)
				llvm::errs() << " " << succ[ti][tj][tiev];
			llvm::errs() << "\n";
		}
  for (unsigned ti=0; ti<pred.size(); ++ti)
		for (unsigned tj=0; tj<pred[ti].size(); ++tj) {
			llvm::errs() << "pred " << ti << "->" << tj << " size: " << pred[ti][tj].size() << "   content:";
			for (unsigned tiev=0; tiev<pred[ti][tj].size(); ++tiev)
				llvm::errs() << " " << pred[ti][tj][tiev];
			llvm::errs() << "\n";
		}
	llvm::errs() << "\n";
}

void VCGraphVclock::to_dot(const PartialOrder& po, const char *edge_params) const
{
  ThreadPairsVclocks& succ = *(po.first);

  llvm::errs() << "\ndigraph {\n";
  
  for (unsigned tid = 0; tid < processes.size(); ++tid) {
		assert((int) tid == processes[tid][0]->getEvent()->iid.get_pid() / 2);
		llvm::errs() << "subgraph cluster_" << tid << "{\n";
		llvm::errs() << "label = \"Process" << tid
								 << " " << processes[tid][0]->getEvent()->cpid;
		if (tid == this->starRoot())
			llvm::errs() << " STAR-ROOT";
		llvm::errs() << "\"\n";
		for (unsigned evid = 0; evid < processes[tid].size(); ++evid) {
			const Node *nd = processes[tid][evid];
			llvm::errs() << "NODE"
									 << tid * 100000 + evid
									 << " [label=\"";
			if (!nd->getEvent())
				llvm::errs() << "init";
			else
				llvm::errs() << *nd->getEvent();
			llvm::errs() << "\"]\n";
		}
		llvm::errs() << "}\n";
	}
  llvm::errs() << "\n";
  
  for (unsigned ti=0; ti<succ.size(); ++ti)
		for (unsigned tiev = 0; tiev < processes[ti].size() - 1; ++tiev) {
			llvm::errs() << "NODE"
			             << ti * 100000 + tiev
			             << " -> NODE"
			             << ti * 100000 + tiev + 1 << "[style=bold]\n";
		}
  
  for (unsigned ti=0; ti<succ.size(); ++ti)
		for (unsigned tj=0; tj<succ[ti].size(); ++tj) {
      int curval = INT_MAX;
			for (int tiev = succ[ti][tj].size() - 1; tiev >= 0; --tiev) {
				
				if (succ[ti][tj][tiev] < curval) {
					curval = succ[ti][tj][tiev];
						llvm::errs() << "NODE"
												 << ti * 100000 + tiev
												 << " -> NODE"
												 << tj * 100000 + curval;
					if (edge_params) {
						llvm::errs() << "[" << edge_params << "]\n";
					} else {
						llvm::errs() << "\n";
					}
				}
				
			}
		}

  llvm::errs() << "}\n\n";

}
