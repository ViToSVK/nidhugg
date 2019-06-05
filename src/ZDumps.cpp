
#ifndef _Z_DUMPS_
#define _Z_DUMPS_

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

  const auto& th_aux = basis.threads_auxes;

  // NODES
  for (unsigned tid = 0; tid < th_aux.size(); ++tid) {
    res << "subgraph cluster_" << tid << "{\n";
    res << "style=\"bold,rounded\" label = \"Th" << tid
        << " " << basis(tid, -1)[0]->cpid;
    if (basis.isRoot(tid))
      res << " ROOT";
    res << "\"\n";
    for (const auto& aux : th_aux[tid]) {
      unsigned line =  basis.thread_aux_to_line_id.at(std::pair<unsigned,int>(tid, aux));
      res << "subgraph cluster_" << 1001+tid*100+aux << "{\n";
      res << "style=\"invis\" label = \"" << ((aux==-1)?"real":"aux") << "\"\n";
      for (unsigned evid = 0; evid < basis(tid, aux).size(); ++evid) {
        const ZEvent *ev = basis(tid, aux)[evid];
        res << "NODE" << line * 100000 + evid
            << " [shape=\"rectangle\", label=\"" << ev->to_string(false) << "\"]\n";
      }

      res << "}\n";
    }

    res << "}\n";
  }
  res << "\n";

  // THREAD ORDER
  for (unsigned tid = 0; tid < th_aux.size(); ++tid) {
    for (const auto& aux : th_aux[tid]) {
      unsigned line =  basis.thread_aux_to_line_id.at(std::pair<unsigned,int>(tid, aux));
      for (int evid = 0; evid < (int) basis(tid, aux).size() - 1; ++evid) {
        res << "NODE" << line * 100000 + evid
            << " -> NODE" << line * 100000 + evid + 1 << "[style=bold]\n";
      }
    }
  }

  // REST ORDER
  for (unsigned tidI = 0; tidI < th_aux.size(); ++tidI) {
  for (const auto& auxI : th_aux.at(tidI)) {
    unsigned lI = basis.thread_aux_to_line_id.at(std::pair<unsigned,int>(tidI, auxI));
    unsigned realI = (auxI == -1) ? lI
    : basis.thread_aux_to_line_id.at(std::pair<unsigned,int>(tidI, -1));

    for (unsigned tidJ = 0; tidJ < th_aux.size(); ++tidJ) {
    for (const auto& auxJ : th_aux.at(tidJ)) {
      unsigned lJ = basis.thread_aux_to_line_id.at(std::pair<unsigned,int>(tidJ, auxJ));
      unsigned realJ = (auxJ == -1) ? lJ
      : basis.thread_aux_to_line_id.at(std::pair<unsigned,int>(tidJ, -1));

      if (lI != lJ) {
        int curval = INT_MAX;
        for (int liev = _succ.at(lI).at(lJ).size() - 1; liev >= 0; --liev) {
          if (auxI != -1 && realI != lJ) {
            // Do not add edges that transitively follow from
            // lI(aux) -> realI -> lJ
            auto fromev = basis.getEvent(tidI, auxI, liev);
            auto succ_idx = succ(fromev, tidI, -1);
            if (succ_idx.second < (int) _succ.at(realI).at(lJ).size() &&
                _succ.at(realI).at(lJ).at(succ_idx.second) < curval)
              curval = _succ.at(realI).at(lJ).at(succ_idx.second);
          }
          if (auxJ != -1 && lI != realJ) {
            // Do not add edges that transitively follow from
            // lI -> realJ -> lJ(aux)
            auto fromev = basis.getEvent(tidI, auxI, liev);
            auto succ_idx = succ(fromev, tidJ, -1);
            if (succ_idx.second < (int) _succ.at(realJ).at(lJ).size() &&
                _succ.at(realJ).at(lJ).at(succ_idx.second) < curval)
              curval = _succ.at(realJ).at(lJ).at(succ_idx.second);
          }
          if (_succ.at(lI).at(lJ).at(liev) < curval) {
            curval = _succ.at(lI).at(lJ).at(liev);
              res << "NODE" << lI * 100000 + liev
                  << " -> NODE" << lJ * 100000 + curval << "\n";
          }
        }
      }
    }
    }
  }
  }

  res << "}\n\n";

  llvm::errs() << res.str();
}


llvm::raw_ostream& operator<<(llvm::raw_ostream& out, const ZAnnotation& annot) {
  out << "Annotation+ {\n";
  for (auto& an : annot) {
    out << an.first.to_string() << "  observes::  ";
    out << an.second.to_string();
    out << "\n";
  }
  out << "}\n";
  return out;
}
void ZAnnotation::dump() const {
  llvm::errs() << *this;
}


#endif // _Z_DUMPS_
