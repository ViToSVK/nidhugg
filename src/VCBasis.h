#ifndef _VCBASIS_H_
#define _VCBASIS_H_

#include <vector>
#include <unordered_map>

#include "VCEvent.h"
#include "VCIID.h"

//class PositiveAnnotation;

typedef int IPid;

template <typename NodeT>
class Processes {
public:
  typedef std::vector<NodeT> ProcessT;
  typedef std::vector<ProcessT> ProcessesT;

  Processes() = default;
  Processes(Processes&& oth) = default;
  Processes& operator=(Processes&& oth) = default;
  Processes(const Processes& oth) = default;
  Processes& operator=(const Processes& oth) = default;

  size_t size() const { return processes.size(); }
  bool empty() const { return processes.empty(); }

  const ProcessesT& getEvents() const { return processes;}
  const ProcessT& operator[](unsigned idx) const {
    assert(idx < size());
    return processes[idx];
  }

  int cpidToProc(const CPid& cpid) {
    int n = 0;
    for (auto& process: processes) {
      assert(process.size() > 0);
      if (process[0]->cpid == cpid) {
        return n;
      }
      ++n;
    }

    return -1;
  }

  void setTopologyRoot(const CPid& cpid) {
    assert(topology_root_index == -1 && "Already have a topology root");
    topology_root_index = cpidToProc(cpid);
    assert(topology_root_index != -1 && "Do not have a process with such cpid");
  }

  void setTopologyRoot(int idx) {
    assert(topology_root_index == -1 && "Already have a topology root");
    assert(idx >= 0 && (size_t) idx < processes.size());
    topology_root_index = idx;
  }

  bool isTopologyRoot(const CPid& cpid) const {
    assert(topology_root_index != -1);
    return cpid == processes[topology_root_index][0]->cpid;
  }

  bool isTopologyRoot(unsigned idx) const {
      return topology_root_index == (int) idx;
  }

  bool isTopologyRoot(const ProcessT& process) const {
    assert(topology_root_index != -1);
    return &process == &processes[topology_root_index];
  }

  bool hasTopologyRoot() const {
    return topology_root_index != -1;
  }

  ProcessT& getTopologyRoot() {
    assert(topology_root_index != -1);
    return processes[topology_root_index];
  }

  const ProcessT& getTopologyRoot() const {
    assert(topology_root_index != -1);
    return processes[topology_root_index];
  }

  const CPid& getTopologyRootCPid() const {
    assert(topology_root_index != -1);
    return processes[topology_root_index][0]->cpid;
  }

	// typename clarifies that (const_)iterator
	// is a class and not a static member
  typename ProcessesT::iterator begin() { return processes.begin(); }
  typename ProcessesT::iterator end() { return processes.end(); }
  typename ProcessesT::const_iterator begin() const { return processes.begin(); }
  typename ProcessesT::const_iterator end() const { return processes.end(); }

  // iterator to seamlessly iterate over all events in the processes
  // with possibilities to skip whole processes
  class events_iterator {
		const ProcessesT& processes;
		unsigned process_idx = 0;
		unsigned event_idx = 0;

		events_iterator(const ProcessesT& processes, bool end = false)
		: processes(processes), process_idx(end ? processes.size() : 0), event_idx(0) {}

		events_iterator(const ProcessesT& processes, unsigned process, unsigned event)
		: processes(processes), process_idx(process), event_idx(event) {}

		friend class Processes;
		friend class VCBasis;
		
   public:
		events_iterator(const events_iterator& oth)
		: processes(oth.processes), process_idx(oth.process_idx), event_idx(oth.event_idx) {}

		events_iterator& operator=(const events_iterator& oth) {
			processes = oth.processes;
			process_idx = oth.process_idx;
			event_idx = oth.event_idx;
			return *this;
		}

		events_iterator& operator++() {
			assert(event_idx < processes[process_idx].size());
			assert(process_idx < processes.size());

			if (event_idx < processes[process_idx].size() - 1) {
				assert(process_idx < processes.size()
								&& "Moved invalid iterator");
				++event_idx;
			} else {
				event_idx = 0;
				++process_idx;
			}

			return *this;
		}

		events_iterator operator++(int) {
			auto tmp = *this;
			operator++();
			return tmp;
		}

		events_iterator& operator--() {
			assert(event_idx < processes[process_idx].size());
			assert(process_idx < processes.size());

			if (event_idx > 0) {
				assert(process_idx < processes.size()
								&& "Moved invalid iterator");
				--event_idx;
			} else {
				assert(event_idx == 0);

				if (process_idx > 0) {
					--process_idx;
					event_idx = processes[process_idx].size() - 1;
				} else {
					assert(process_idx == 0);
					// we reached the end
					process_idx = processes.size();
				}
			}

			return *this;
		}

		events_iterator operator--(int) {
			auto tmp = *this;
			operator--();
			return tmp;
		}

		unsigned getProcessID() const { return process_idx; }
		unsigned getEventID() const { return event_idx; }
		void shiftProcess() { ++process_idx; event_idx = 0;}

		// void invalidate() { event_idx = 0; process_idx = processes.size(); }
		// returns the processes.end() iterator for the processes from this iterator.
		// This iterator serves as an end in both directions
		events_iterator end() const { return events_iterator(processes, true /* end */); }

		// return an iterator to the first event of the process that the
		// current iterator is in
		events_iterator processStart() const {
			return events_iterator(processes, process_idx, 0);
		}

		// return the iterator to the last event from the process that this
		// iterator currently is in
		events_iterator processEnd() const {
			auto tmp = *this;
			tmp.skipProcess();
			return tmp;
		}

		// return an iterator for the first event of the next process
		// XXX: maybe we could just set the indices instead of doing process end + shift
		events_iterator nextProcess() const {
			auto tmp = processEnd();
			return ++tmp;
		}

		// return an iterator to the first event of the process that the
		// current iterator is in
		events_iterator prevProcess() const {
			auto tmp = processStart();
			return --tmp;
		}

		// calling this results in jumping on the next process
		// upon call to operator++ (it just jumps on the last
		// event of the current process)
		void skipProcess() { event_idx = processes[process_idx].size() - 1; }

    bool isAtProcessEnd() const { return event_idx == processes[process_idx].size() - 1; }
		
		const NodeT operator*() const {
			assert(process_idx < processes.size());
			assert(event_idx < processes[process_idx].size());

			return processes[process_idx][event_idx];
		}

		bool operator==(const events_iterator& oth) const {
			return process_idx == oth.process_idx && event_idx == oth.event_idx;
		}
		bool operator!=(const events_iterator& oth) const {
				return !operator==(oth);
		}
		bool operator<(const events_iterator& oth) const {
			return (process_idx == oth.process_idx) ?
							event_idx < oth.event_idx : process_idx < oth.process_idx;
		}
  };

  events_iterator events_begin() const { return events_iterator(processes); }
  events_iterator events_end() const { return events_iterator(processes, true); }
  events_iterator getIterator(unsigned process = 0, unsigned event = 0) const {
    return events_iterator(processes, process, event);
  }	

protected:
  int topology_root_index = -1;
  ProcessesT processes;
};

class VCBasis : public Processes<const VCEvent *> {
public:
  VCBasis() = default;

  VCBasis(ProcessesT&& procs) {
		processes.swap(procs);
  }

  void swap(ProcessesT& oth) {
		processes.swap(oth);
  }

  // TraceT is a sequence of VCEvents.
  // (Whatever the TraceT is, the iterators must
  // return the VCEvent& on dereference.)
  template <typename TraceT>
  VCBasis(TraceT& trace, const VCEvent *till = nullptr)
  {
    // mapping from already discovered threads to our processes indexes
    std::unordered_map<IPid, unsigned> mapping;
    mapping.reserve(trace.size() / 2);

    for (const VCEvent& event : trace) {
      IPid pid = event.iid.get_pid();
      unsigned idx;

      auto it = mapping.find(pid);
      if (it == mapping.end()) {
        idx = processes.size();
        mapping[pid] = idx;
        processes.emplace_back(std::vector<const VCEvent *>());
      } else
        idx = it->second;

      processes[idx].push_back(&event);

      // we want processes only till this event?
      // (we know that anything that occured after this event
      //  is irrelevant for us)
      if (till == &event)
        break;
    }
  }

	// Get an iterator pointing to the given VCEvent
  events_iterator getIterator(const VCEvent& ev) const {
    unsigned pidx = 0;
    for (auto& proc : processes) {
			if (proc[0]->cpid != ev.cpid) {
				++pidx;
				continue;
			}

			// we found the right process,
			// now find the right event
			unsigned n = 0;
			for (const VCEvent *pev : proc) {
				if (pev == &ev)
					return events_iterator(processes, pidx, n);

				++n;
			}

			break;
    }

    return events_end();
  }

  // Get an iterator pointing to the given VCIID
  events_iterator getIterator(const VCIID& vciid) const {
    unsigned pidx = 0;
    for (auto& proc : processes) {
			if (proc[0]->cpid != vciid.cpid) {
				++pidx;
				continue;
			}

			// we found the right process,
			// now find the right event
			unsigned n = 0;
			for (const VCEvent *pev : proc) {
				if (pev->instruction_order == vciid.instruction_order)
					return events_iterator(processes, pidx, n);

				++n;
			}

			break;
    }

    return events_end();
  }

	// Get an iterator pointing to (the beginning of) the given CPid
  events_iterator getIterator(const CPid& cpid, unsigned event = 0) const {
		unsigned pidx = 0;
		for (auto& proc : processes) {
			if (proc[0]->cpid == cpid)
				return events_iterator(processes, pidx, event);

			++pidx;
		}

		return events_end();
	}

  void dump() const;
};

#endif // _VCBASIS_H_
