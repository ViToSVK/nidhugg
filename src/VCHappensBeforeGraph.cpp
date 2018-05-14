#include <unordered_map>

#include "VCTBHelpers.h" // isWrite
#include "VCHappensBeforeGraph.h"
#include "SCC.h"

VCHappensBeforeGraph::VCHappensBeforeGraph(VCHappensBeforeGraph&& oth)
: basis(std::move(oth.basis)),
  annotation(std::move(oth.annotation)),
  max_search_index(oth.max_search_index),
  nodes(std::move(oth.nodes)),
  initial_node(oth.initial_node),  
  instr_to_nodeset(std::move(oth.instr_to_nodeset))
{
}

VCHappensBeforeGraph& VCHappensBeforeGraph::operator=(VCHappensBeforeGraph&& oth) {
  if (&oth != this)
    *this = std::move(oth);
  return *this;
}

bool VCHappensBeforeGraph::isAcyclic()
{
  AcyclicTopoOrder<Node> detect(max_search_index,
                                true /* only detect cycles */);
  auto start = begin();
  auto en = end();
  bool ret = detect.compute(start, en);
  max_search_index = detect.getIndex();

  return ret;
}

std::set<VCHappensBeforeGraph::Node *>&
VCHappensBeforeGraph::getNodesReachableFromCompute(Node *n)
{
  auto& ret = reachable_nodes[n];
  assert(ret.first != edge_relation_stage);
  ret.second.clear();
  
  std::vector<Node *> to_process;
  to_process.push_back(n);
  
  while(!to_process.empty()) {
    std::vector<Node *> newly_discovered;
    newly_discovered.reserve(to_process.size());

    for (Node *nd : to_process) {
			if ((ret.second.insert(nd)).second) {
				if (reachable_nodes[nd].first == edge_relation_stage) {
					// nd has up-to-date info, use it
					ret.second.insert(reachable_nodes[nd].second.begin(),
					                  reachable_nodes[nd].second.end());
				} else {
					for (Node *succ : nd->successors)
						newly_discovered.push_back(succ);					
				}
			}
    }

    to_process.swap(newly_discovered);
  }

  // we do not want to be reflexive
  ret.second.erase(n);
  // n now has up-to-date info
  ret.first = edge_relation_stage;
  
  return ret.second;
}

std::set<VCHappensBeforeGraph::Node *>&
VCHappensBeforeGraph::getNodesReachingToCompute(Node *n)
{
  auto& ret = rev_reachable_nodes[n];
  assert(ret.first != edge_relation_stage);
  ret.second.clear();
  
  std::vector<Node *> to_process;
  to_process.push_back(n);
  
  while(!to_process.empty()) {
    std::vector<Node *> newly_discovered;
    newly_discovered.reserve(to_process.size());

    for (Node *nd : to_process) {
			if ((ret.second.insert(nd)).second) {
				if (rev_reachable_nodes[nd].first == edge_relation_stage) {
					// nd has up-to-date info, use it
					ret.second.insert(rev_reachable_nodes[nd].second.begin(),
					                  rev_reachable_nodes[nd].second.end());
				} else {
					for (Node *pred : nd->predecessors)
						newly_discovered.push_back(pred);					
				}
			}
    }

    to_process.swap(newly_discovered);
  }

  // we do not want to be reflexive
  ret.second.erase(n);
  // n now has up-to-date info
  ret.first = edge_relation_stage;
  
  return ret.second;
}

// this code supposes that the graph is acyclic, so it does
// not do a fixpoint (since on acyclic graph the fixpoint
// is reached after one iteration).
void VCHappensBeforeGraph::makeTransitiveClosureSimple()
{
    for (auto& it : nodes) {
        Node *nd = it.second;
        for (Node *s : getNodesReachableFrom(nd))
          nd->addSuccessor(s);
    }
}

void VCHappensBeforeGraph::makeTransitiveClosure()
{
  AcyclicTopoOrder<Node> topo(max_search_index);

  auto start = begin();
  auto en = end();
  bool ret = topo.compute(start, en);
  assert(ret && "The graph is not acyclic");
  ((void)(ret)); // so ret does not appear unused on release

  // XXX: maybe we should do just simple recursive DFS?
  max_search_index = topo.getIndex();
  auto& ord = topo.getOrder();
  assert(ord.size() == size());
  for (unsigned i = 0, e = ord.size(); i < e; ++i) {
    // make a copy of current successors
    auto succs = ord[i]->getSuccessors();
    for (Node *succ : succs)
      for (Node *succ_of_succ : succ->getSuccessors())
          ord[i]->addSuccessor(succ_of_succ);
  }

  assert(isTransitivelyClosed());
}

bool VCHappensBeforeGraph::isTransitivelyClosed()
{
  for (Node& nd1 : *this) {
    for (Node *succ1 : nd1.getSuccessors()) {
        for (Node *nd2 : succ1->getSuccessors()) {
            // we have edge nd1->succ1 && succ1->nd2
            // so check whether we have nd1->nd2
            if (!hasEdge(&nd1, nd2))
              return false;
        }
    }
  }

  return true;
}

// This works since the graph is acyclic
std::vector<VCEvent> VCHappensBeforeGraph::linearize()
{
  std::vector<VCEvent> trace;
  AcyclicTopoOrder<Node> detect(max_search_index);

  auto start = begin();
  auto en = end();
  bool ret = detect.compute(start, en);
  assert(ret && "The graph is not acyclic");
	((void)(ret)); // so ret does not appear unused on release

  max_search_index = detect.getIndex();
  // the order is reversed
  auto& ord = detect.getOrder();
  for (unsigned i = ord.size(); i > 0; --i)
    if (ord[i - 1]->event) {
      trace.push_back(ord[i - 1]->event->blank_copy());
      trace.back().id = trace.size() - 1;
    } else
      // the initial event must be the first event
      assert(i == ord.size());

  return trace;
}

bool VCHappensBeforeGraph::addNecessaryEdges()
{
	// TODO
	
  // makeTransitiveClosure();

  return true;
}

void VCHappensBeforeGraph::computeAnnotationInfo()
{
  // For each annotated read, find the
	// corresponding good and bad writes
  for (auto& ann : annotation) {
    Node *annot_read = getNode(ann.first);
		assert(annot_read && annot_read != initial_event && "Malformed annotation");
		assert(isRead(annot_read->event) && "Event is not a read");

    
    // TODO handle good/bad local write
		
		
		// Find good and bad remote writes
		// Start at the beginning of the corresponding process
    auto it = basis.getIterator(ann.second.second);
		auto event_id = 0;
		bool prEnd = false;
		
		while (!prEnd) {
			// Check if this is the last event of the process
      if (it.isAtProcessEnd())
				prEnd = true;

			if ((*it)->may_conflict && isWrite(**it) &&
					(*it)->ml == annot_read->event->ml) {
				if ((*it)->value == ann.second.first) {
          // Found a good write
				} else {
          // Found a bad write
				}
			}

			++it;
			++event_id;
		}
	  
	}
		
	/*
  // fill in the annotation edges
  for (auto& rw_pair : annotation) {
    Node *read = mapAnnotToNode(rw_pair.first);
    Node *write = mapAnnotToNode(rw_pair.second);

    assert(read);
    if (!write)
      write = initial_node;

    write->addSuccessor(read);
  }
  */
}

void VCHappensBeforeGraph::_constructGraph()
{
  // This should be more efficient than to dynamically adjust the hash tables
  unsigned reserve_size = 0;
  for (auto&b : basis)
    reserve_size += b.size();
  instr_to_nodeset.reserve(reserve_size);
  nodes.reserve(reserve_size);
  reachable_nodes.reserve(reserve_size);
  rev_reachable_nodes.reserve(reserve_size);

  // Create the initial node corresponding to the initial event
  initial_node = new Node(nullptr);
  addNode(initial_node);

  // Set of nodes that call pthread_create and pthread_join
  std::vector<Node *> spawns;
  std::vector<Node *> joins;

	// Create the nodes and add the program structure
	assert(!basis.empty());
  unsigned idx = 0;
  for (auto& base : basis) {
    Node *node;
    Node *last = nullptr;

		// Create for this process
    assert(!base.empty());
    unsigned idx2 = 0;
    for (const VCEvent *event : base) {
      assert(event->iid.get_pid() >= 0);
      node = new Node(event);
      node->process_id = idx;
      node->event_id = idx2;
      addNode(node);
      if (event->instruction)
        instr_to_nodeset[event->instruction].push_back(node);

      if (last)
        last->addSuccessor(node);
      else
        // this is the first node/event in the base,
        // make it a successor of the inital node/event
        initial_node->addSuccessor(node);

      last = node;

      // XXX: what about function pointer calls?
      if (event->instruction) {
        if (is_function_call(event->instruction, "pthread_create"))
          spawns.push_back(node);
        else if (is_function_call(event->instruction, "pthread_join"))
          joins.push_back(node);
      }
      ++idx2;
    }
    ++idx;
  }

	
  // Add edges following from the threads creation
  for (Node *spwn : spawns) {
		bool a = false;
    // find the first event of the thread in the basis
    for (auto& base : basis) {
      assert(!base.empty());
      if (base[0]->cpid == spwn->event->childs_cpid) {
				addSuccessor(spwn, nodes[base[0]]);
				a = true;
        break;
      }
    }
		assert(a && "Invalid pthread_create");
		((void)(a)); // so a does not appear unused on release
  }
  for (Node *jn : joins) {
		bool a = false;		
    // find the first event of the thread in the basis
    for (auto& base : basis) {
      assert(!base.empty());
      if (base[0]->cpid == jn->event->childs_cpid) {
				addSuccessor(nodes[base[base.size() - 1]], jn);
				a = true;
        break;
      }
    }
		assert(a && "Invalid pthread_join");
		((void)(a)); // so a does not appear unused on release		
  }

  // Compute info following from the annotation
	computeAnnotationInfo();
	
}
