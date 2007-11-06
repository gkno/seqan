// Projekt, mit dem die Demos getestet werden koennen

#define main runAllocator
#include "../../demos/allocator.cpp"
#undef main

#define main runAlphabet
#include "../../demos/alphabet.cpp"
#undef main

/*
//!kompiliert nicht!
#define main runIterator
#include "../../demos/iterator.cpp"
#undef main
*/

#define main runRootedIterator
#include "../../demos/rooted_iterator.cpp"
#undef main

#define main runString
#include "../../demos/string_1.cpp"
#undef main



#define main runIndexSA
#include "../../demos/index_sufarray.cpp"
#undef main

#define main runIndexFind
#include "../../demos/index_find.cpp"
#undef main

#define main runIndexMUMs
#include "../../demos/index_mums.cpp"
#undef main

#define main runIndexSuperMaxRepeats
#include "../../demos/index_supermaxrepeats.cpp"
#undef main

#define main runIndexMaxRepeats
#include "../../demos/index_maxrepeats.cpp"
#undef main



#define main runGraph0
#include "../../demos/graph_algo_bfs.cpp"
#undef main

#define main runGraph1
#include "../../demos/graph_algo_dfs.cpp"
#undef main

#define main runGraph2
#include "../../demos/graph_algo_topsort.cpp"
#undef main

#define main runGraph3
#include "../../demos/graph_algo_scc.cpp"
#undef main

#define main runGraph4
#include "../../demos/graph_algo_tree_prim.cpp"
#undef main

#define main runGraph5
#include "../../demos/graph_algo_tree_kruskal.cpp"
#undef main

#define main runGraph6
#include "../../demos/graph_algo_path_dag.cpp"
#undef main

#define main runGraph7
#include "../../demos/graph_algo_path_bellmanford.cpp"
#undef main

#define main runGraph8
#include "../../demos/graph_algo_path_dijkstra.cpp"
#undef main

#define main runGraph9
#include "../../demos/graph_algo_path_allpairs.cpp"
#undef main

#define main runGraph10
#include "../../demos/graph_algo_path_floydwarshall.cpp"
#undef main

#define main runGraph11
#include "../../demos/graph_algo_path_transitive.cpp"
#undef main

#define main runGraph12
#include "../../demos/graph_algo_flow_fordfulkerson.cpp"
#undef main

#define main runGraph13
#include "../../demos/graph_algo_matching_pathgrowing.cpp"
#undef main

#define main runGraph14
#include "../../demos/graph_algo_lis.cpp"
#undef main

#define main runGraph15
#include "../../demos/graph_algo_his.cpp"
#undef main

#define main runGraph16
#include "../../demos/graph_algo_lcs.cpp"
#undef main

#define main runGraph17
#include "../../demos/graph_align_nw.cpp"
#undef main

#define main runGraph18
#include "../../demos/graph_align_gotoh.cpp"
#undef main

#define main runGraph19
#include "../../demos/graph_align_hirschberg.cpp"
#undef main

#define main runGraph20
#include "../../demos/graph_align_sw.cpp"
#undef main

#define main runGraph21
#include "../../demos/graph_align_guide_nj.cpp"
#undef main

#define main runGraph22
#include "../../demos/graph_align_guide_upgma.cpp"
#undef main

#define main runGraph23
#include "../../demos/graph_align_msa.cpp"
#undef main

#define main runFind
#include "../../demos/find.cpp"
#undef main

#define main runFileFormat
#include "../../demos/file_format.cpp"
#undef main

int main() 
{
	runAllocator();
	runAlphabet();
	//runIterator();
	runRootedIterator();
	runString();

	runIndexSA();
	runIndexFind();
	runIndexMUMs();
	runIndexSuperMaxRepeats();
	runIndexMaxRepeats();

	runGraph0();runGraph1();runGraph2();runGraph3();runGraph4();runGraph5();
	runGraph6();runGraph7();runGraph8();runGraph9();runGraph10();runGraph11();
	runGraph12();runGraph13();runGraph14();runGraph15();runGraph16();runGraph17();
	runGraph18();runGraph19();runGraph20();runGraph21();runGraph22();runGraph23();

	//runFind();
	runFileFormat();
}
