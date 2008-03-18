// Projekt, mit dem die Demos getestet werden koennen

#define TEST_BASICS
#define TEST_FIND
#define TEST_ALIGNMENT
#define TEST_MODIFIERS
#define TEST_INDEX
#define TEST_GRAPH


#ifdef TEST_BASICS

	#define main runAllocator
	#include "../../library/demos/allocator.cpp"
	#undef main

	#define main runAlphabet
	#include "../../library/demos/alphabet.cpp"
	#undef main


	#define main runIterator
	#include "../../library/demos/iterator.cpp"
	#undef main


	#define main runRootedIterator
	#include "../../library/demos/rooted_iterator.cpp"
	#undef main

	#define main runString1
	#include "../../library/demos/string_1.cpp"
	#undef main


	#define main runFileFormat
	#include "../../library/demos/file_format.cpp"
	#undef main

#endif

#ifdef TEST_FIND

	#define main runFindExact
	#include "../../library/demos/find_exact.cpp"
	#undef main

	#define main runFindApprox
	#include "../../library/demos/find_approx.cpp"
	#undef main

	#define main runFindWild
	#include "../../library/demos/find_wild.cpp"
	#undef main

#endif

#ifdef TEST_ALIGNMENT

	#define main runAlignment
	#include "../../library/demos/alignment.cpp"
	#undef main

	#define main runAlignmentLocal
	#include "../../library/demos/alignment_local.cpp"
	#undef main

#endif

#ifdef TEST_MODIFIERS

	#define main runModifierModReverse
	#include "../../library/demos/modifier_modreverse.cpp"
	#undef main

	#define main runModifierModView
	#include "../../library/demos/modifier_modview.cpp"
	#undef main

	#define main runModifierNested
	#include "../../library/demos/modifier_nested.cpp"
	#undef main

#endif

#ifdef TEST_INDEX

	#define main runIndexSA
	#include "../../library/demos/index_sufarray.cpp"
	#undef main

	#define main runIndexFind
	#include "../../library/demos/index_find.cpp"
	#undef main

	#define main runIndexMUMs
	#include "../../library/demos/index_mums.cpp"
	#undef main

	#define main runIndexSuperMaxRepeats
	#include "../../library/demos/index_supermaxrepeats.cpp"
	#undef main

	#define main runIndexMaxRepeats
	#include "../../library/demos/index_maxrepeats.cpp"
	#undef main

	#define main runIndexMummy
	#include "../../library/demos/index_mummy.cpp"
	#undef main

	#define main runIndexNodePredicate
	#include "../../library/demos/index_node_predicate.cpp"
	#undef main

#endif

#ifdef TEST_GRAPH

	#define main runGraph0
	#include "../../library/demos/graph_algo_bfs.cpp"
	#undef main

	#define main runGraph1
	#include "../../library/demos/graph_algo_dfs.cpp"
	#undef main

	#define main runGraph2
	#include "../../library/demos/graph_algo_topsort.cpp"
	#undef main

	#define main runGraph3
	#include "../../library/demos/graph_algo_scc.cpp"
	#undef main

	#define main runGraph4
	#include "../../library/demos/graph_algo_tree_prim.cpp"
	#undef main

	#define main runGraph5
	#include "../../library/demos/graph_algo_tree_kruskal.cpp"
	#undef main

	#define main runGraph6
	#include "../../library/demos/graph_algo_path_dag.cpp"
	#undef main

	#define main runGraph7
	#include "../../library/demos/graph_algo_path_bellmanford.cpp"
	#undef main

	#define main runGraph8
	#include "../../library/demos/graph_algo_path_dijkstra.cpp"
	#undef main

	#define main runGraph9
	#include "../../library/demos/graph_algo_path_allpairs.cpp"
	#undef main

	#define main runGraph10
	#include "../../library/demos/graph_algo_path_floydwarshall.cpp"
	#undef main

	#define main runGraph11
	#include "../../library/demos/graph_algo_path_transitive.cpp"
	#undef main

	#define main runGraph12
	#include "../../library/demos/graph_algo_flow_fordfulkerson.cpp"
	#undef main

	#define main runGraph13
	#include "../../library/demos/graph_algo_matching_pathgrowing.cpp"
	#undef main

	#define main runGraph14
	#include "../../library/demos/graph_algo_lis.cpp"
	#undef main

	#define main runGraph15
	#include "../../library/demos/graph_algo_his.cpp"
	#undef main

	#define main runGraph16
	#include "../../library/demos/graph_algo_lcs.cpp"
	#undef main

	#define main runGraph17
	#include "../../library/demos/graph_align_nw.cpp"
	#undef main

	#define main runGraph18
	#include "../../library/demos/graph_align_gotoh.cpp"
	#undef main

	#define main runGraph19
	#include "../../library/demos/graph_align_hirschberg.cpp"
	#undef main

	#define main runGraph25
	#include "../../library/demos/graph_consensus.cpp"
	#undef main

	#define main runGraph20
	#include "../../library/demos/graph_align_sw.cpp"
	#undef main

	#define main runGraph21
	#include "../../library/demos/graph_align_guide_nj.cpp"
	#undef main

	#define main runGraph22
	#include "../../library/demos/graph_align_guide_upgma.cpp"
	#undef main

	#define main runGraph23
	#include "../../library/demos/graph_align_msa.cpp"
	#undef main

	#define main runGraph24
	#include "../../library/demos/seqan_tcoffee.cpp"
	#undef main

	#define main runGraph26
	#include "../../library/demos/graph_hmm.cpp"
	#undef main

#endif


int main(int argc, const char *argv[]) 
{
#ifdef TEST_BASICS
	runAllocator();
	runAlphabet();
	runIterator();
	runRootedIterator();
	runString1();
	runFileFormat();
#endif

#ifdef TEST_FIND
	runFindExact();
	runFindApprox();
	runFindWild();
#endif

#ifdef TEST_ALIGNMENT
	runAlignment();
	runAlignmentLocal();
#endif

#ifdef TEST_MODIFIERS
	runModifierModReverse();
	runModifierModView();
	runModifierNested();
#endif

#ifdef TEST_INDEX
	runIndexSA();
	runIndexFind();
	runIndexMUMs();
	runIndexSuperMaxRepeats();
	runIndexMaxRepeats();
	runIndexMummy(argc, argv);
	runIndexNodePredicate();
#endif

#ifdef TEST_GRAPH
	runGraph0();runGraph1();runGraph2();runGraph3();runGraph4();runGraph5();
	runGraph6();runGraph7();runGraph8();runGraph9();runGraph10();runGraph11();
	runGraph12();runGraph13();runGraph14();runGraph15();runGraph16();runGraph17();
	runGraph18();runGraph19();runGraph20();runGraph21();runGraph22();runGraph23();
	runGraph24(argc, argv);
	runGraph25(argc, argv);
	runGraph26();
#endif

	return 0;
}
