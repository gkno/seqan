#ifndef SEQAN_HEADER_TEST_GRAPH_INTERVAL_TREE_H
#define SEQAN_HEADER_TEST_GRAPH_INTERVAL_TREE_H


namespace SEQAN_NAMESPACE_MAIN
{

	
template<typename TValue, typename TConstructSpec, typename TStoreSpec>
void IntervalTreeTest() {
//____________________________________________________________________________
	
	typedef int TCargo;
	typedef IntervalAndCargo<TValue,TCargo> TInterval;
	
	typedef IntervalTreeNode<TInterval,TStoreSpec> TNode;
	typedef Graph< > TGraph;
	typedef String<TNode> TPropertyMap;
	
	String<TInterval> intervals;
	resize(intervals,10);
	intervals[0].i1 = (TValue)5.0;
	intervals[0].i2 = (TValue)20.0;
	intervals[1].i1 = (TValue)6.0;
	intervals[1].i2 = (TValue)13.0;
	intervals[2].i1 = (TValue)8.0;
	intervals[2].i2 = (TValue)10.0;
	intervals[3].i1 = (TValue)9.0;
	intervals[3].i2 = (TValue)18.0;
	intervals[4].i1 = (TValue)15.0;
	intervals[4].i2 = (TValue)23.0;
	intervals[5].i1 = (TValue)21.0;
	intervals[5].i2 = (TValue)25.0;
	intervals[6].i1 = (TValue)29.0;
	intervals[6].i2 = (TValue)42.0;
	intervals[7].i1 = (TValue)20.0;
	intervals[7].i2 = (TValue)36.0;
	intervals[8].i1 = (TValue)38.0;
	intervals[8].i2 = (TValue)42.0;
	intervals[9].i1 = (TValue)36.0;
	intervals[9].i2 = (TValue)48.0;
	for(unsigned int i = 0; i < length(intervals); ++i) intervals[i].cargo = i;
	
	//construct interval tree graph and its corresponding property map
	TGraph g;
	TPropertyMap pm;
	
	createIntervalTree(g,pm,intervals,(TValue)26.0,TConstructSpec());
	
	TValue query = (TValue)25.0;
	String<TCargo> result;
	findIntervals(g, pm, query, result);
	
	SEQAN_TASSERT(length(result)==1);
	SEQAN_TASSERT(result[0]==7);

	query = (TValue)7.0;
	resize(result,0);
	findIntervals(g, pm, query, result);
	SEQAN_TASSERT(length(result)==2);
	SEQAN_TASSERT((result[0]==0 && result[1]==1) ||(result[0]==1 && result[1]==0));
	
	query = (TValue)37.0;
	resize(result,0);
	findIntervals(g, pm, query, result);
	SEQAN_TASSERT(length(result)==2);
	SEQAN_TASSERT((result[0]==9 && result[1]==6) ||(result[0]==6 && result[1]==9));
	
}

template<typename TValue>
void IntervalTreeRestTest() {
//____________________________________________________________________________

//	std::cout <<"\nIntervalTreeRestTest\n";
	typedef int TCargo;
	typedef IntervalAndCargo<TValue,TCargo> TInterval;

	typedef IntervalTreeNode<TInterval,StoreIntervals> TNode;
	typedef Graph< > TGraph;
	typedef String<TNode> TPropertyMap;


	String<TInterval> intervals;
	resize(intervals,8);
	intervals[0].i1 = (TValue)1.0;
	intervals[0].i2 = (TValue)5.0;
	intervals[1].i1 = (TValue)3.0;
	intervals[1].i2 = (TValue)8.0;
	intervals[2].i1 = (TValue)2.0;
	intervals[2].i2 = (TValue)10.0;
	intervals[3].i1 = (TValue)21.0;
	intervals[3].i2 = (TValue)23.0;
	intervals[4].i1 = (TValue)22.0;
	intervals[4].i2 = (TValue)29.0;
	intervals[5].i1 = (TValue)23.0;
	intervals[5].i2 = (TValue)25.0;
	intervals[6].i1 = (TValue)25.0;
	intervals[6].i2 = (TValue)31.0;
	intervals[7].i1 = (TValue)26.0;
	intervals[7].i2 = (TValue)35.0;
	for(unsigned int i = 0; i < length(intervals); ++i) intervals[i].cargo = i;
	
	TGraph g;
	TPropertyMap pm;
	
	createIntervalTree(g, pm, intervals);


	SEQAN_TASSERT(length(pm)==(TValue)6.0);
	SEQAN_TASSERT(pm[0].center==(TValue)18.0);
	SEQAN_TASSERT(pm[1].center==(TValue)5.5);
	SEQAN_TASSERT(pm[2].center==(TValue)3.0);
	SEQAN_TASSERT(pm[3].center==(TValue)28.0);
	SEQAN_TASSERT(pm[4].center==(TValue)23.0);
	SEQAN_TASSERT(pm[5].center==(TValue)22.0);

	SEQAN_TASSERT(length(pm[0].list1)==0);
	SEQAN_TASSERT(length(pm[1].list2)==2);
	SEQAN_TASSERT(length(pm[2].list1)==1);
	SEQAN_TASSERT(length(pm[3].list2)==3);
	SEQAN_TASSERT(length(pm[4].list1)==1);
	SEQAN_TASSERT(length(pm[5].list2)==1);

	SEQAN_TASSERT(leftBoundary(pm[1].list1[0])==(TValue)2.0);
	SEQAN_TASSERT(leftBoundary(pm[1].list1[1])==(TValue)3.0);
	SEQAN_TASSERT(rightBoundary(pm[1].list2[0])==(TValue)10.0);
	SEQAN_TASSERT(rightBoundary(pm[1].list2[1])==(TValue)8.0);

	SEQAN_TASSERT(leftBoundary(pm[2].list1[0])==(TValue)1.0);
	SEQAN_TASSERT(rightBoundary(pm[2].list2[0])==(TValue)5.0);

	SEQAN_TASSERT(leftBoundary(pm[3].list1[0])==(TValue)22.0);
	SEQAN_TASSERT(leftBoundary(pm[3].list1[1])==(TValue)25.0);
	SEQAN_TASSERT(leftBoundary(pm[3].list1[2])==(TValue)26.0);
	SEQAN_TASSERT(rightBoundary(pm[3].list2[0])==(TValue)35.0);
	SEQAN_TASSERT(rightBoundary(pm[3].list2[1])==(TValue)31.0);
	SEQAN_TASSERT(rightBoundary(pm[3].list2[2])==(TValue)29.0);

	SEQAN_TASSERT(leftBoundary(pm[4].list1[0])==(TValue)23.0);
	SEQAN_TASSERT(rightBoundary(pm[4].list2[0])==(TValue)25.0);

	SEQAN_TASSERT(leftBoundary(pm[5].list1[0])==(TValue)21.0);
	SEQAN_TASSERT(rightBoundary(pm[5].list2[0])==(TValue)23.0);

	//suche nach allen intervallen den query punkt enthalten
	TValue query = (TValue)13.0;
	String<TCargo> result;
	findIntervals(g, pm, query, result);

	SEQAN_TASSERT(length(result)==0);


}


template<typename TValue>
void
testEasyIntervalTree()
{
	String<TValue> begins;
	String<TValue> ends;
	resize(begins,5);
	resize(ends,5);
	begins[0] = 2;
	begins[1] = 5;
	begins[2] = 3;
	begins[3] = 7;
	begins[4] = 1;
	ends[0] = 7;
	ends[1] = 8;
	ends[2] = 5;
	ends[3] = 10;
	ends[4] = 5;

	IntervalTree<TValue> itree(begin(begins),begin(ends),length(begins));

	TValue query = 7;
	//Default cargo int
	String<int> result;
	findIntervals(itree, query, result);
	
	SEQAN_TASSERT(length(result)==2);
	SEQAN_TASSERT(result[0]==1);
	SEQAN_TASSERT(result[1]==3);
	

	typedef IntervalAndCargo<TValue,int> TInterval;
	String<TInterval> intervals;
	resize(intervals,length(begins));
	typename Iterator<String<TValue> >::Type interval_begins = begin(begins);
	typename Iterator<String<TValue> >::Type interval_ends = begin(ends);
	size_t i = 0;
	while(i<length(begins))
	{
		intervals[i].i1 = value(interval_begins);
		++interval_begins;
		intervals[i].i2 = value(interval_ends);
		++interval_ends;
		intervals[i].cargo = i;
		++i;
	}

	IntervalTree<TValue,int> itree2(intervals);
	findIntervals(itree2, query, result);
	SEQAN_TASSERT(length(result)==2);
	SEQAN_TASSERT(result[0]==1);
	SEQAN_TASSERT(result[1]==3);
	//SEQAN_TASSERT(itree == itree2);
	

	//std::cout << "\n\n\n";

	//cout << itree.g;
	//for(int i = 0; i < length(itree.pm); ++i)
	//	std::cout << (itree.pm[i]).center << " ";

	//std::cout << "\n\n\n";

	addInterval(itree,16,20,5);
	addInterval(itree,8,18);
	SEQAN_TASSERT(itree.interval_counter==7);
	

	//addInterval(itree,-10,15,7);
	//addInterval(itree,20,100,8);
	//addInterval(itree,17,18,9);
	//addInterval(itree,0,25,10);
	//addInterval(itree,77,88,11);
	//addInterval(itree,-1,10,12);

	//cout << itree.g;
	//for(int i = 0; i < length(itree.pm); ++i)
	//{
	//	std::cout << "\n"<< (itree.pm[i]).center << "\n";
	//	for(int j = 0; j < length((itree.pm[i]).list1); ++j)
	//		std::cout << leftBoundary((itree.pm[i]).list1[j])<<" ";
	//}
	//std::cout << "\n\n\n";

	//bool chec = removeInterval(itree,-10,15,7);
	//if(chec) std::cout << "gefunden und gelöscht\n";
	//else  std::cout << "interval gibt es nicht\n";
	//chec = removeInterval(itree,20,100,8);
	//if(chec) std::cout << "gefunden und gelöscht\n";
	//else  std::cout << "interval gibt es nicht\n";
	//chec = removeInterval(itree,17,18,9);
	//if(chec) std::cout << "gefunden und gelöscht\n";
	//else  std::cout << "interval gibt es nicht\n";
	//chec = removeInterval(itree,0,25,10);
	//if(chec) std::cout << "gefunden und gelöscht\n";
	//else  std::cout << "interval gibt es nicht\n";
	//chec = removeInterval(itree,77,88,11);
	//if(chec) std::cout << "gefunden und gelöscht\n";
	//else  std::cout << "interval gibt es nicht\n";
	//chec = removeInterval(itree,-1,10,12);
	//if(chec) std::cout << "gefunden und gelöscht\n";
	//else  std::cout << "interval gibt es nicht\n";
	//
	//std::cout << "\n\njetzt müsste es wieder der alte baum sein:\n";

	//cout << itree.g;
	//for(int i = 0; i < length(itree.pm); ++i)
	//{
	//	std::cout << "\n"<< (itree.pm[i]).center << "\n";
	//	for(int j = 0; j < length((itree.pm[i]).list1); ++j)
	//		std::cout << leftBoundary((itree.pm[i]).list1[j])<<" ";
	//}
	//std::cout << "\n\n\n";

}


void Test_GraphIntervalTree() 
{
	srand((unsigned)time(NULL)); 

	testEasyIntervalTree<int>();

//	std::cout << "RandomCenter ";
	IntervalTreeTest<int,RandomCenter,StorePointsOnly>();

//	std::cout << "\nMidCenter ";
	IntervalTreeTest<double,MidCenter,StoreIntervals>();
	
//	std::cout << "\nComputeCenter ";
	IntervalTreeTest<int,ComputeCenter,StoreIntervals>();
	
//	std::cout << "\nRest ";
	IntervalTreeRestTest<int>();

	debug::verifyCheckpoints("projects/library/seqan/graph/graph_impl_interval_tree.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_impl_interval_types.h");


}

}

#endif
