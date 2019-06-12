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


bool ZTrace::respectsAnnotation() const
{
  for (const ZEvent& evref: trace) {
    const ZEvent *ev = &evref;
    if (isRead(ev) && annotation.defines(ev->threadID(), ev->eventID())) {
      const ZObs& obs = annotation.getObs(ev);
      const ZEvent *obsB = (obs.thr == INT_MAX)
        ? graph.getBasis().initial() : graph.getBasis().getEvent(obs);
      const ZEvent *obsM = (isInitial(obsB))
        ? graph.getBasis().initial() : &(trace.at(obsB->write_other_trace_id));
      assert(obsB->value == obsM->value);
      assert(isInitial(obsB) || (isWriteB(obsB) && isWriteM(obsM) &&
                                 sameMl(obsB, obsM) && sameMl(ev, obsB)));
      const ZEvent *realObservation = (ev->observed_trace_id == -1)
        ? graph.getBasis().initial() : &(trace.at(ev->observed_trace_id));

      if (realObservation != obsB && realObservation != obsM) {
        dumpTrace(trace);
        graph.getPo().dump();
        annotation.dump();
        llvm::errs() << "This read         :::  ";
        ev->dump();
        llvm::errs() << "Should observeB   :::  ";
        obsB->dump();
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
