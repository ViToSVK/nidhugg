/* Copyright (C) 2016-2017 Marek Chalupa
 * Copyright (C) 2017-2018 Viktor Toman
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

#ifndef _VC_BASIS_H_
#define _VC_BASIS_H_

#include <vector>
#include <unordered_map>

#include "VCIID.h"

class Node {
	// Indices into the basis
	// (INT_MAX, INT_MAX) -> the initial node
	unsigned process_id;
	unsigned event_id;
	// Pointer to the Event
	// nullptr -> the initial node
	const VCEvent *event;

 public:
	Node(unsigned pid, unsigned evid, const VCEvent* ev)
	: process_id(pid), event_id(evid), event(ev) {}
	Node(const Node& oth)
	: process_id(oth.process_id),
		event_id(oth.event_id),
		event(oth.event) {}

	bool operator==(const Node& oth) const {
    return (process_id == oth.process_id &&
						event_id == oth.event_id &&
						event == oth.event);
	}
	
	unsigned getProcessID() const { return process_id; }
	unsigned getEventID() const { return event_id; }
	const VCEvent *getEvent() const { return event; }
	
};

class VCBasis {
 public:
  typedef std::vector<Node *> ProcessT;
  typedef std::vector<ProcessT> ProcessesT;
	
 protected:
  ProcessesT processes;
	std::unordered_map<CPid, unsigned> cpid_to_processid;
	std::unordered_map<const VCEvent *, Node *> event_to_node;
	std::unordered_map<VCIID, Node *> read_vciid_to_node;
	// make use of copy elision: emplace VCIID using copy constructor instead of move
	// std::unordered_map<const llvm::Instruction *, std::vector<Node *>> instr_to_nodeset;
	
 public:
	// Basis is not responsible for any resources
	// All Node* allocation and deletion happens
	// in the graph classes that inherit Basis
  VCBasis() = default;
  VCBasis(VCBasis&& oth) = default;
  VCBasis& operator=(VCBasis&& oth) = default;
  VCBasis(const VCBasis& oth) = default;
  VCBasis& operator=(const VCBasis& oth) = default;

  const ProcessT& operator[](unsigned idx) const {
    assert(idx < size());
    return processes[idx];
  }
	
  bool empty() const { return processes.empty(); }
  size_t size() const { return processes.size(); }
	size_t nodes_size() const {
		size_t res = 0;
		for (auto& pr : processes)
			res += pr.size();
		return res;
	}

  int cpidToProcessID(const CPid& cpid) const {
		auto it = cpid_to_processid.find(cpid);
		if (it == cpid_to_processid.end())
			return -1;
		return it->second;
  }
	
	// typename clarifies that (const_)iterator
	// is a class and not a static member
  typename ProcessesT::iterator begin() { return processes.begin(); }
  typename ProcessesT::iterator end() { return processes.end(); }
  typename ProcessesT::const_iterator begin() const { return processes.begin(); }
  typename ProcessesT::const_iterator end() const { return processes.end(); }

  // Iterator to seamlessly iterate over all events in the processes
  // with possibilities to skip whole processes
  class events_iterator {
		const ProcessesT& processes;
		unsigned process_idx = 0;
		unsigned event_idx = 0;

		events_iterator(const ProcessesT& processes, bool end = false)
		: processes(processes), process_idx(end ? processes.size() : 0), event_idx(0) {}

		events_iterator(const ProcessesT& processes, unsigned process, unsigned event)
		: processes(processes), process_idx(process), event_idx(event) {}

		friend class VCBasis;
		
   public:
		events_iterator(events_iterator&& oth)
			: processes(std::move(oth.processes)), process_idx(oth.process_idx), event_idx(oth.event_idx) {}
		events_iterator& operator=(events_iterator&& oth) = delete;
		
		events_iterator(const events_iterator& oth)
		  : processes(oth.processes), process_idx(oth.process_idx), event_idx(oth.event_idx) {}
		events_iterator& operator=(const events_iterator& oth) = delete;

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

		// First event of the current process
    void jumpToProcessStart() {
			event_idx = 0;
		}
		events_iterator processStart() const {
			auto tmp = *this;
			tmp.jumpToProcessStart();
			return tmp;
		}

		// Last event of the current process
		void jumpToProcessEnd() {
			if (process_idx < processes.size())
			  event_idx = processes[process_idx].size() - 1;
		}
		events_iterator processEnd() const {
			auto tmp = *this;
			tmp.jumpToProcessEnd();
			return tmp;
		}

		// First event of the next process
		void jumpToNextProcess() {
			if (process_idx < processes.size()) {
				++process_idx; event_idx = 0;
			}
		}
		events_iterator nextProcess() const {
			auto tmp = *this;
			tmp.jumpToNextProcess();
			return ++tmp;
		}

    // The end (in both directions)
		void jumpToEnd() {
			process_idx = processes.size(); event_idx = 0;
		}
		events_iterator end() const {
			return events_iterator(processes, true /* end */);
	  }

		// Can not change NodeT through the iterator
		const Node *operator*() const {
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

  events_iterator nodes_begin() const { return events_iterator(processes); }
  events_iterator nodes_end() const { return events_iterator(processes, true); }
  events_iterator nodes_iterator(unsigned process = 0, unsigned event = 0) const {
    return events_iterator(processes, process, event);
  }	

	// Iterator pointing to the given Node
  events_iterator nodes_iterator(const Node *nd) const {
		assert(nd && "Do not have such node");
    return nodes_iterator(nd->getProcessID(), nd->getEventID());
	}

	// Node corresponding to the given VCEvent
	const Node *getNode(const VCEvent& ev) const {
    auto it = event_to_node.find(&ev);
		assert(it != event_to_node.end()
					 && "Given event is not tied to any node");
    return it->second;
  }
	
	// Iterator corresponding to the given VCEvent
  events_iterator nodes_iterator(const VCEvent& ev) const {
    return nodes_iterator(getNode(ev));
  }

	// Node corresponding to the given VCIID
	const Node *getNode(const VCIID& vciid) const {
    auto it = read_vciid_to_node.find(vciid);
		assert(it != event_to_node.end()
					 && "Given VCIID is not tied to any read event");
    return it->second;
  }
	
	// Iterator corresponding to the given VCIID
  events_iterator nodes_iterator(const VCIID& vciid) const {
    return nodes_iterator(getNode(vciid));
  }

	// Get an iterator pointing to (the beginning of) the given CPid
  events_iterator nodes_iterator(const CPid& cpid, unsigned evidx = 0) const {
    int pidx = cpidToProcessID(cpid);
		assert(pidx >= 0 && pidx < processes.size()
					 && "Given cpid is not tied to any process");
		assert(evidx < processes[pidx].size()
					 && "Given evidx is out of range for the corresponding process");
		return nodes_iterator(pidx, evidx);
	}

};

/*
CODE IF WE EVER NEED TO CONSIDER A STAR ARCHITECTURE WITHIN BASIS

  void setTopologyRoot(int idx) {
    assert(topology_root_index == -1 && "Already have a topology root");
    assert(idx >= 0 && (size_t) idx < processes.size());
    topology_root_index = idx;
  }

  void setTopologyRoot(const CPid& cpid) {
    assert(topology_root_index == -1 && "Already have a topology root");
    topology_root_index = cpidToProcessID(cpid);
    assert(topology_root_index != -1 && "Do not have a process with such cpid");
  }

  bool hasTopologyRoot() const {
    return topology_root_index != -1;
  }

  bool isTopologyRoot(unsigned idx) const {
      return topology_root_index == (int) idx;
  }

  bool isTopologyRoot(const ProcessT& process) const {
    assert(topology_root_index != -1);
    return &process == &processes[topology_root_index];
  }

  bool isTopologyRoot(const CPid& cpid) const {
    assert(topology_root_index != -1);
    return cpid == cpidvector[topology_root_index];
  }

  const ProcessT& getTopologyRoot() const {
    assert(topology_root_index != -1);
    return processes[topology_root_index];
  }

  const CPid& getTopologyRootCPid() const {
    assert(topology_root_index != -1);
    return cpidvector[topology_root_index];
  }

 private:
	int topology_root_index = -1;

 */

#endif // _VC_BASIS_H_
