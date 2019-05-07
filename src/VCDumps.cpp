#include <fstream>

#include <llvm/Support/raw_ostream.h>
#include <llvm/Support/raw_os_ostream.h>
#include <llvm/IR/Instruction.h>
#include <llvm/IR/Instructions.h>

#include "VCEvent.h"
#include "ZAnnotation.h"
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
  switch(ev.kind) {
   case VCEvent::Kind::DUMMY :
     out << ev.id << "_<size: " << ev.size << ">";
     break;
   case VCEvent::Kind::LOAD :
     out << ev.id << "_read " << ev.ml.addr.to_string() << " <- " << ev.value;
     break;
   case VCEvent::Kind::STORE :
     out << ev.id << "_write " << ev.value << " -> " << ev.ml.addr.to_string();
     break;
   case VCEvent::Kind::SPAWN :
     out << ev.id << "_spawn " << ev.childs_cpid;
     break;
   case VCEvent::Kind::JOIN :
     out << ev.id << "_join " << ev.childs_cpid << "_<size: " << ev.size << ">";
     break;
   case VCEvent::Kind::M_INIT :
     out << ev.id << "_mutexinit " << ev.ml.addr.to_string();
     break;
   case VCEvent::Kind::M_DESTROY :
     out << ev.id << "_mutexdestroy " << ev.ml.addr.to_string();
     break;
   case VCEvent::Kind::M_LOCK :
     out << ev.id << "_lock " << ev.ml.addr.to_string() << "_<size: " << ev.size << ">";
     break;
   case VCEvent::Kind::M_UNLOCK :
     out << ev.id << "_unlock " << ev.ml.addr.to_string();
     break;
   default :
     out << "_unknown";
  }
  return out;
}

void VCEvent::dump() const {
  llvm::errs() << *this << "\n";
}

llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const Node& nd)
{
  out << "[" << nd.getProcessID() << "][" << nd.getEventID() << "]";
  if (nd.getEvent())
    out << "__" << *(nd.getEvent());
  return out;
}

void Node::dump() const {
  llvm::errs() << *this << "\n";
}

llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZAnnotation::Loc& loc) {
  char c = (loc == ZAnnotation::Loc::LOCAL)?'L':
    ((loc == ZAnnotation::Loc::REMOTE)?'R':'A');
  out << c;
  return out;
}

llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZAnnotation::Ann& ann) {
  out << ann.value << "-" << ann.loc << "_";
  if (ann.loc != ZAnnotation::Loc::LOCAL) {
    out << "remotegood:";
    for (auto& iid : ann.goodRemote)
      out << "[" << iid.first << "][" << iid.second << "] ";
  }
  if (ann.loc != ZAnnotation::Loc::REMOTE
      && ann.goodLocal) {
    if (ann.goodLocal->first == INT_MAX)
      out << "_localgood:INIT ";
    else
      out << "_localgood:[" << ann.goodLocal->first
          << "][" << ann.goodLocal->second << "] ";
  }
  return out;
}

void ZAnnotation::Ann::dump() const {
  llvm::errs() << *this << "\n";
}

llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZAnnotation& annot) {
  out << "Annotation+ {\n";
  for (auto& pr : annot) {
    out << "( ["  << pr.first.first << "][" << pr.first.second
    << "] observes:: ";
    out << pr.second;
    out << ")\n";
  }
  out << "}\n";
  return out;
}

void ZAnnotation::dump() const {
  llvm::errs() << *this;
}

void VCGraphVclock::dump_po(const PartialOrder& po) const {
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

void VCGraphVclock::to_dot(const PartialOrder& po, const char *edge_params) const {
  ThreadPairsVclocks& succ = *(po.first);

  llvm::errs() << "\ndigraph {\n";

  for (unsigned tid = 0; tid < processes.size(); ++tid) {
    assert((int) tid == processes[tid][0]->getEvent()->iid.get_pid() / 2);
    llvm::errs() << "subgraph cluster_" << tid << "{\n";
    llvm::errs() << "label = \"Pr" << tid
                 << " " << processes[tid][0]->getEvent()->cpid;
    if (tid == this->starRoot())
      llvm::errs() << " ROOT";
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

  assert(succ.size() == processes.size());

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
