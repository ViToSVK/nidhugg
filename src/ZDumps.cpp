#include <fstream>

#include <llvm/Support/raw_ostream.h>
#include <llvm/Support/raw_os_ostream.h>
#include <llvm/IR/Instruction.h>
#include <llvm/IR/Instructions.h>

#include <sstream>
#include <set>

#include "ZBuilderTSO.h"
#include "ZPartialOrder.h"
#include "ZAnnotation.h"
#include "ZGraph.h"


void removeSubstrings(std::string& s, std::string&& p) {
  std::string::size_type n = p.length();
  for (std::string::size_type i = s.find(p);
       i != std::string::npos;
       i = s.find(p))
    s.erase(i, n);
}


std::string ZEvent::to_string(bool write_cpid = true) const {
  std::stringstream res;

  res << traceID() << "_";
  if (write_cpid)
    res << cpid << "_";
  res << "[" << threadID() << "," << auxID() << "," << eventID() << "]";
  switch(kind) {
   case ZEvent::Kind::DUMMY :
     res << " <s:" << size << ">";
     break;
   case ZEvent::Kind::READ :
     res << " read <- " << ml.addr.to_string() << " <O:" << observed_trace_id << ">";
     break;
   case ZEvent::Kind::WRITEB :
     res << " writeB -> " << ml.addr.to_string();
     break;
   case ZEvent::Kind::WRITEM :
     res << " writeM -> " << ml.addr.to_string()
         << " <B:" << write_other_trace_id << ">";
     break;
   case ZEvent::Kind::SPAWN :
     res << " spawn " << childs_cpid;
     break;
   case ZEvent::Kind::JOIN :
     res << " join " << childs_cpid;
     break;
   case ZEvent::Kind::M_INIT :
     res << " mutexinit " << ml.addr.to_string();
     break;
   case ZEvent::Kind::M_DESTROY :
     res << " mutexdestroy " << ml.addr.to_string();
     break;
   case ZEvent::Kind::M_LOCK :
     res << " lock " << ml.addr.to_string();
     break;
   case ZEvent::Kind::M_UNLOCK :
     res << " unlock " << ml.addr.to_string();
     break;
   default :
     res << " unknown";
  }

  return res.str();
}
llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZEvent& ev)
{
  out << ev.to_string();
  return out;
}
void ZEvent::dump() const {
  llvm::errs() << *this << "\n";
}


void ZBuilderTSO::dump(const std::vector<ZEvent>& trace)
{
  llvm::errs() << "TRACE::: " << trace.size() << " EVENTS\n";
  for (const auto& ev : trace)
    ev.dump();
}


void ZPartialOrder::dump() const {
  std::stringstream res;

  res << "\ndigraph {\n";

  int max_thread = -1;
  int max_aux = -2;
  for (const auto& taux_line : basis.thread_aux_to_line_id) {
    if ((int) taux_line.first.first > max_thread)
      max_thread = taux_line.first.first;
    if (taux_line.first.second > max_aux)
      max_aux = taux_line.first.second;
  }

  assert(max_thread >= 0);
  assert(max_aux >= -1);
  std::vector<std::set<int>> th_aux;
  for (int i = 0; i < max_thread+1; ++i)
    th_aux.push_back(std::set<int>());

  for (const auto& taux_line : basis.thread_aux_to_line_id)
    th_aux.at(taux_line.first.first).emplace(taux_line.first.second);

  // NODES
  for (unsigned tid = 0; tid < th_aux.size(); ++tid) {
    res << "subgraph cluster_" << tid << "{\n";
    res << "label = \"Th" << tid
        << " " << basis(tid, -1)[0]->cpid;
    if (basis.isRoot(tid))
      res << " ROOT";
    res << "\"\n";
    for (const auto& aux : th_aux[tid]) {
      unsigned line =  basis.thread_aux_to_line_id.at(std::pair<unsigned,int>(tid, aux));
      for (unsigned evid = 0; evid < basis(tid, aux).size(); ++evid) {
        const ZEvent *ev = basis(tid, aux)[evid];
        res << "NODE" << line * 100000 + evid
            << " [shape=\"rectangle\", label=\"" << ev->to_string(false) << "\"]\n";
      }
    }

    res << "}\n";
  }
  res << "\n";


  // THREAD ORDER
  for (unsigned tid = 0; tid < th_aux.size(); ++tid) {
    for (const auto& aux : th_aux[tid]) {
      unsigned line =  basis.thread_aux_to_line_id.at(std::pair<unsigned,int>(tid, aux));
      for (unsigned evid = 0; evid < basis(tid, aux).size() - 1; ++evid) {
        res << "NODE" << line * 100000 + evid
            << " -> NODE" << line * 100000 + evid + 1 << "[style=bold]\n";
      }
    }
  }

  // REST ORDER
  for (unsigned li=0; li<_succ.size(); ++li)
    for (unsigned lj=0; lj<_succ[li].size(); ++lj) {
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

  llvm::errs() << res.str();
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
