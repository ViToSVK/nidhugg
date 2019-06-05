/* Copyright (C) 2016-2017 Marek Chalupa
 * Copyright (C) 2017-2019 Viktor Toman
 *
 * This file is part of Nidhugg.
 *
 * Nidhugg is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Nidhugg is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include "ZTrace.h"
#include "ZDumps.cpp"


bool ZTrace::respectsAnnotation() const
{
  for (const ZEvent& evref: trace) {
    const ZEvent *ev = &evref;
    if (isRead(ev) && annotation.defines(ev->threadID(), ev->eventID())) {
      const ZObs& obs = annotation.getObs(ev);
      const ZEvent *obsEv = (obs.thr == INT_MAX)
        ? graph.getBasis().initial() : graph.getBasis().getEvent(obs);
      const ZEvent *obsM = (isInitial(obsEv))
        ? graph.getBasis().initial() : &(trace.at(obsEv->write_other_trace_id));
      assert(obsEv->value == obsM->value);
      assert(isInitial(obsEv) || (isWriteB(obsEv) && isWriteM(obsM) &&
                                  sameMl(obsEv, obsM) && sameMl(ev, obsEv)));
      const ZEvent *realObservation = &(trace.at(ev->observed_trace_id));

      if (realObservation != obsEv && realObservation != obsM) {
        dumpTrace(trace);
        graph.getPo().dump();
        annotation.dump();
        llvm::errs() << "This read         :::  ";
        ev->dump();
        llvm::errs() << "Should observeB   :::  ";
        obsEv->dump();
        llvm::errs() << "Should observeM   :::  ";
        obsM->dump();
        llvm::errs() << "Actually observed :::  ";
        realObservation->dump();
        return false;
      }
    }
  }

  return true;
}
