#ifndef _VCHAPPENSBEFOREGRAPH_H_
#define _VCHAPPENSBEFOREGRAPH_H_

#include <set>
#include <vector>
#include <unordered_map>

#include "VCEvent.h"
#include "VCAnnotation.h"
#include "VCBasis.h"
#include "VCIID.h"

class VCHappensBeforeGraph {

	typedef VCIID AnnotationKeyT;

	const VCBasis& basis;
  const VCAnnotation& annotation;
  unsigned max_search_index = 0; // for topoorders and sccs
  unsigned edge_relation_stage = 0; // for (re)computing reachable nodes

 public:
  VCHappensBeforeGraph(const VCBasis& b,
	                     const VCAnnotation& annot)
  : basis(b), annotation(annot)
  {
    _constructGraph();
    reachable_nodes.reserve(nodes.size());
    rev_reachable_nodes.reserve(nodes.size());
  }

  VCHappensBeforeGraph(VCHappensBeforeGraph&& oth);
  VCHappensBeforeGraph& operator=(VCHappensBeforeGraph&& oth);

  ~VCHappensBeforeGraph() {
		for (auto& it : nodes)
    delete it.second;
	}



  const VCAnnotation& getAnnotation() const { return annotation; }
  const VCBasis& getBasis() const { return basis; }

  class Node {
    typedef std::set<Node *> EdgesT;
    EdgesT successors;
    EdgesT predecessors;

    // unique id of the node,
    // serves as the index to the vector of nodes
    unsigned id = 0;
    // indices into the basis
    unsigned process_id = 0;
    unsigned event_id = 0;

    // for searching cycles and SCC
    // FIXME: move this out of the node to auxiliary mappings
    unsigned dfsid, lowpt;
    bool is_on_stack;

    bool addSuccessor(Node *nd) {
      assert(nd != this && "We do not want self-loops");
      assert(nd && "Given node is nullptr");
			bool a = (successors.insert(nd)).second;
      bool b = (nd->predecessors.insert(this)).second;
			assert(a == b && "Inconsistent successor and predecessor relations");
			((void)(b)); // so b does not appear unused on release
			return a;
    }

    bool removeSuccessor(Node *nd) {
			assert(nd && "Given node is nullptr");
      bool a = successors.erase(nd);
      bool b = nd->predecessors.erase(this);
      assert(a == b && "Inconsistent successor and predecessor relations");
			((void)(b)); // so b does not appear unused on release
			return a;
    }

   public:
    const VCEvent *event;
    Node(const VCEvent *ev)
      : dfsid(0), lowpt(0), is_on_stack(0), event(ev) {}

    bool operator==(const Node& oth) const { return this->id == oth.id; }
		bool operator!=(const Node& oth) const { return !(*this == oth); }
		
    void set_on_stack(bool a) { is_on_stack = a; }
    void set_dfsid(unsigned id) { dfsid = id; }
    void set_lowpt(unsigned lp) { lowpt = lp; }

    unsigned get_dfsid() const { return dfsid; }
    unsigned get_lowpt() const { return lowpt; }
    bool on_stack() const { return is_on_stack; }

    EdgesT& getSuccessors() { return successors; }
    const EdgesT& getSuccessors() const { return successors; }
    EdgesT& getPredecessors() { return predecessors; }
    const EdgesT& getPredecessors() const { return predecessors; }

    friend class VCHappensBeforeGraph;
  };

 private:
  // Cache for reachable nodes (instead of doing the transitive closure)
  // value.first: the graph relation stage when this was computed
  // value.second: the set of nodes reachable from 'key'
  std::unordered_map<Node *, std::pair<unsigned, std::set<Node *>>> reachable_nodes;
  // Cache for reverse reachable nodes
  // value.first: the graph relation stage when this was computed
  // value.second: the set of nodes reverse reachable from 'key'
  std::unordered_map<Node *, std::pair<unsigned, std::set<Node *>>> rev_reachable_nodes;

 public:

  void addSuccessor(Node *from, Node *to) {
		assert(from && "Given node is nullptr");
		bool a = from->addSuccessor(to);
		if (a) {
			edge_relation_stage++;
			edge_relation_stage %= UINT32_MAX;
		}
	}

  void removeSuccessor(Node *from, Node *to) {
		assert(from && "Given node is nullptr");
		bool a = from->removeSuccessor(to);
		if (a) {
			edge_relation_stage++;
			edge_relation_stage %= UINT32_MAX;
		}
	}

  std::set<Node *>& getNodesReachableFromCompute(Node *n);
  std::set<Node *>& getNodesReachingToCompute(Node *n);

  std::set<Node *>& getNodesReachableFrom(Node *n) {
    auto& ret = reachable_nodes[n];
    if (ret.first == edge_relation_stage)
			return ret.second;

    return getNodesReachableFromCompute(n);
  }

  std::set<Node *>& getNodesReachingTo(Node *n) {
    auto& ret = rev_reachable_nodes[n];
    if (ret.first == edge_relation_stage)
      return ret.second;

    return getNodesReachingToCompute(n);
  }

  std::set<Node *>& getNodesReachableFrom(const VCEvent *e) {
    Node *tmp = getNode(e);
    assert(tmp && "Invalid event");
    return getNodesReachableFrom(tmp);
  }

  const Node *getNode(const VCEvent *e) const {
    auto it = nodes.find(e);
    if (it == nodes.end())
      return nullptr;

    return it->second;
  }

  Node *getNode(const VCEvent *e) {
    auto it = nodes.find(e);
    if (it == nodes.end())
      return nullptr;

    return it->second;
  }

  Node *getNode(const AnnotationKeyT& iid) const {
		if (iid.instruction == nullptr) {
			assert(initial_node);
			return initial_node;
		}	

		// find the read and write nodes for this annotation pair
		auto it = instr_to_nodeset.find(iid.instruction);
		if (it == instr_to_nodeset.end())
			return nullptr;

		for (Node *nd : it->second) {
			if (nd->event->instruction_order == iid.instruction_order &&
					nd->event->cpid == iid.cpid)
				return nd;
		}

		return nullptr;
  }

  unsigned addNode(Node *nd) {
		assert(nodes.find(nd->event) == nodes.end() && "Such node already present");
		nd->id = nodes.size();
		nodes.emplace(nd->event, nd);
		reachable_nodes.emplace(nd,std::make_pair<unsigned, std::set<Node *>>
		                           (UINT32_MAX, {}));
		rev_reachable_nodes.emplace(nd,std::make_pair<unsigned, std::set<Node *>>
		                           (UINT32_MAX, {}));
		return nd->id;
  }

  // check whether ev2 is reachable from ev1
  bool canReach(const VCEvent *ev1, const VCEvent *ev2) {
		Node *ev1_node = getNode(ev1);
		Node *ev2_node = getNode(ev2);
		assert(ev1_node);
		assert(ev2_node);

		auto& reach = getNodesReachableFrom(ev1_node);

		return (reach.find(ev2_node) != reach.end());
	}
	
  bool canReach(const VCEvent& ev1, const VCEvent& ev2) {
    return canReach(&ev1, &ev2);
  }

  bool hasEdge(const VCEvent *e1, const VCEvent *e2) {
    auto n1 = getNode(e1);
    auto n2 = getNode(e2);
    assert(n1 && n2 && "Do not have such node");
    return n1->successors.count(n2) == 1;
  }

  bool hasEdge(Node *n1, Node *n2) {
    assert(n1 && n2 && "Do not have such node");
    return n1->successors.count(n2) == 1;
  }

  size_t size() const { return nodes.size(); }

  bool isAcyclic();
  bool addNecessaryEdges();
  void makeTransitiveClosure();
  void makeTransitiveClosureSimple();
  // for debugging
  bool isTransitivelyClosed();

  // linearize the PO induced by this graph
  std::vector<VCEvent> linearize();

  void to_dot(const char *edge_params=nullptr) const;
  void dump() const;

  class nodes_iterator {
    std::unordered_map<const VCEvent *, Node *>::iterator it;
    nodes_iterator(std::unordered_map<const VCEvent *, Node *>::iterator&& i)
    : it(i) {}

    friend class VCHappensBeforeGraph;
   public:
    nodes_iterator& operator++() { ++it; return *this;}
    nodes_iterator operator++(int) { auto tmp = *this; ++it; return tmp;}
    bool operator==(const nodes_iterator& oth) { return it == oth.it; }
    bool operator!=(const nodes_iterator& oth) { return it != oth.it; }
    Node& operator*() { return *it->second; }
  };

  nodes_iterator begin() { return nodes_iterator(nodes.begin()); }
  nodes_iterator end() { return nodes_iterator(nodes.end()); }

private:
  void _constructGraph();

  void computeAnnotationInfo();
	
  std::unordered_map<const VCEvent *, Node *> nodes;
  
  Node *initial_node = nullptr;

  // create a mapping from instructions to nodes.
  // We could have mapping from (CPid, Instruction) -> Node,
  // but the comparsion of CPid would be too unefficient for our purposes
  // right now...
  std::unordered_map<const llvm::Instruction *, std::vector<Node *>> instr_to_nodeset;

  // Node *mapAnnotToNode(const AnnotationKeyT& iid) const; same as getNode

	std::unordered_map<const VCIID *, const VCEvent *> localGoodWrite;
	std::unordered_map<const VCIID *, const VCEvent *> localBadWrite;

	// Set iterator follows the lexicographic ordering
	std::unordered_map<const VCIID *,
	                   std::set< std::pair<int, const VCEvent *> >> goodWrites;
	std::unordered_map<const VCIID *, std::pair<int, const VCEvent *>> badWrites;
	
};

#endif // _VCHAPPENSBEFOREGRAPH_H_
