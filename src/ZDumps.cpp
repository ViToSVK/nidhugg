#include <fstream>

#include <llvm/Support/raw_ostream.h>
#include <llvm/Support/raw_os_ostream.h>
#include <llvm/IR/Instruction.h>
#include <llvm/IR/Instructions.h>

#include "ZEvent.h"
#include "ZPartialOrder.h"
#include "ZAnnotation.h"
#include "ZGraph.h"

#include <sstream>
#include <set>


void removeSubstrings(std::string& s, std::string&& p) {
  std::string::size_type n = p.length();
  for (std::string::size_type i = s.find(p);
       i != std::string::npos;
       i = s.find(p))
    s.erase(i, n);
}


llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZEvent& ev)
{
  out << "[" << ev.threadID() << "," << ev.auxID() << "," << ev.eventID() << "]__";
  out << ev.trace_id << "__";
  switch(ev.kind) {
   case ZEvent::Kind::DUMMY :
     "<size: " << ev.size << ">";
     break;
   case ZEvent::Kind::READ :
     out << "read <- " << ev.ml.addr.to_string() << "_<observes: " << ev.observed_trace_id << ">";
     break;
   case ZEvent::Kind::WRITEB :
     out << "writeBuf -> " << ev.ml.addr.to_string();
     break;
   case ZEvent::Kind::WRITEM :
     assert(ev.writeOther);
     out << "writeMem -> " << ev.ml.addr.to_string() << "_<Buf [" << ev.writeOther->threadID()
         << "," << ev.writeOther->auxID() << "," << ev.writeOther->eventID() << "]>";
     break;
   case ZEvent::Kind::SPAWN :
     out << "spawn " << ev.childs_cpid;
     break;
   case ZEvent::Kind::JOIN :
     out << "join " << ev.childs_cpid << "_<size: " << ev.size << ">";
     break;
   case ZEvent::Kind::M_INIT :
     out << "mutexinit " << ev.ml.addr.to_string();
     break;
   case ZEvent::Kind::M_DESTROY :
     out << "mutexdestroy " << ev.ml.addr.to_string();
     break;
   case ZEvent::Kind::M_LOCK :
     out << "lock " << ev.ml.addr.to_string() << "_<size: " << ev.size << ">";
     break;
   case ZEvent::Kind::M_UNLOCK :
     out << "unlock " << ev.ml.addr.to_string();
     break;
   default :
     out << "unknown";
  }
  return out;
}
void ZEvent::dump() const {
  llvm::errs() << *this << "\n";
}


void ZPartialOrder::dump() const {
  std::stringstream res;

  res << "\ndigraph {\n";

  int max_thread = -1;
  int max_aux = -2;
  for (const auto& taux_line : basis.thread_aux_to_line_id) {
    if (taux.first.first > max_thread)
      max_thread = taux.first.first;
    if (taux.first.second > max_aux)
      max_aux = taux.first.second;
  }

  assert(max_thread >= 0);
  assert(max_aux >= -1);
  std::vector<std::set<int>> th_aux;
  for (unsigned i = 0; i < max_thread+1; ++i)
    th_aux.push_back(std::set<int>());

  for (const auto& taux_line : basis.thread_aux_to_line_id)
    th_aux.at(taux.first.first).emplace(taux.first.second);

  // NODES
  for (unsigned tid = 0; tid < th_aux.size(); ++tid) {
    res << "subgraph cluster_" << tid << "{\n";
    res << "label = \"Th" << tid
        << " " << basis[tid, -1][0]->cpid;
    if (basis.isRoot(tid))
      res << " ROOT";
    res << "\"\n";
    for (const auto& aux : th_aux[tid]) {
      unsigned line =  basis.thread_aux_to_line_id.at(std::pair<unsigned,int>(tid, aux));
      for (unsigned evid = 0; evid < basis[tid, aux].size(); ++evid) {
        const ZEvent *ev = basis[tid, aux][evid];
        res << "NODE" << line * 100000 + evid
            << " [label=\"" << *ev << "\"]\n";
      }
    }

    res << "}\n";
  }
  res << "\n";


  // THREAD ORDER
  for (unsigned tid = 0; tid < th_aux.size(); ++tid) {
    for (const auto& aux : th_aux[tid]) {
      unsigned line =  basis.thread_aux_to_line_id.at(std::pair<unsigned,int>(tid, aux));
      for (unsigned evid = 0; evid < basis[tid, aux].size() - 1; ++evid) {
        res << "NODE" << line * 100000 + evid
            << " -> NODE" << line * 100000 + evid + 1 << "[style=bold]\n";
      }
    }
  }

  // REST ORDER
  for (unsigned li=0; ti<_succ.size(); ++li)
    for (unsigned lj=0; tj<_succ[li].size(); ++lj) {
      int curval = INT_MAX;
      for (int liev = _succ[li][lj].size() - 1; liev >= 0; --liev) {
        if (_succ[li][lj][liev] < curval) {
          curval = _succ[li][lj][liev];
            res << "NODE" << li * 100000 + liev
                << " -> NODE" << lj * 100000 + curval << "\n";
        }
      }
    }

  res << "}\n\n";

  llvm::errs() << res.to_string() << "\n";
}


/*
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
*/
