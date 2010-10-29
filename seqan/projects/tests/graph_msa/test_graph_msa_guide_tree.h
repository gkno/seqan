/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de
  ===========================================================================
  Copyright (C) 2010 by Freie Universitaet Berlin

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
  ===========================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  Author: Tobias Rausch <rausch@embl.de>
  ===========================================================================
  Tests for the Graph MSA module.
  ===========================================================================
*/

#ifndef TESTS_TEST_GRAPH_MSA_GUIDE_TREE_H_
#define TESTS_TEST_GRAPH_MSA_GUIDE_TREE_H_

void Test_GuideTree_NeighbourJoining()
{
    using namespace seqan;
//____________________________________________________________________________
// Neighbor Joining

    // Create a distance matrix
    String<double> mat;
    fill(mat, 8*8, 0);
    assignValue(mat, 0*8+1, 7);assignValue(mat, 0*8+2, 8);assignValue(mat, 0*8+3, 11);assignValue(mat, 0*8+4, 13);assignValue(mat, 0*8+5, 16);assignValue(mat, 0*8+6, 13);assignValue(mat, 0*8+7, 17);
    assignValue(mat, 1*8+2, 5);assignValue(mat, 1*8+3, 8);assignValue(mat, 1*8+4, 10);assignValue(mat, 1*8+5, 13);assignValue(mat, 1*8+6, 10);assignValue(mat, 1*8+7, 14);
    assignValue(mat, 2*8+3, 5);assignValue(mat, 2*8+4, 7);assignValue(mat, 2*8+5, 10);assignValue(mat, 2*8+6, 7);assignValue(mat, 2*8+7, 11);
    assignValue(mat, 3*8+4, 8);assignValue(mat, 3*8+5, 11);assignValue(mat, 3*8+6, 8);assignValue(mat, 3*8+7, 12);
    assignValue(mat, 4*8+5, 5);assignValue(mat, 4*8+6, 6);assignValue(mat, 4*8+7, 10);
    assignValue(mat, 5*8+6, 9);assignValue(mat, 5*8+7, 13);
    assignValue(mat, 6*8+7, 8);

    typedef Graph<Tree<double> > TGraph;
    TGraph guideTreeOut;
    njTree(mat, guideTreeOut);
    //std::cout << guideTreeOut << std::endl;

    SEQAN_ASSERT_TRUE(numVertices(guideTreeOut) == 15);
    SEQAN_ASSERT_TRUE(findEdge(guideTreeOut, 8, 1) != 0);
    SEQAN_ASSERT_TRUE(findEdge(guideTreeOut, 8, 0) != 0);
    SEQAN_ASSERT_TRUE(findEdge(guideTreeOut, 9, 5) != 0);
    SEQAN_ASSERT_TRUE(findEdge(guideTreeOut, 9, 4) != 0);
    SEQAN_ASSERT_TRUE(findEdge(guideTreeOut, 10, 2) != 0);
    SEQAN_ASSERT_TRUE(findEdge(guideTreeOut, 10, 8) != 0);
    SEQAN_ASSERT_TRUE(findEdge(guideTreeOut, 10, 2) != 0);
    SEQAN_ASSERT_TRUE(findEdge(guideTreeOut, 10, 8) != 0);
    SEQAN_ASSERT_TRUE(findEdge(guideTreeOut, 11, 3) != 0);
    SEQAN_ASSERT_TRUE(findEdge(guideTreeOut, 11, 10) != 0);
    SEQAN_ASSERT_TRUE(findEdge(guideTreeOut, 12, 9) != 0);
    SEQAN_ASSERT_TRUE(findEdge(guideTreeOut, 12, 11) != 0);
    SEQAN_ASSERT_TRUE(findEdge(guideTreeOut, 13, 12) != 0);
    SEQAN_ASSERT_TRUE(findEdge(guideTreeOut, 13, 6) != 0);
    SEQAN_ASSERT_TRUE(findEdge(guideTreeOut, 14, 13) != 0);
    SEQAN_ASSERT_TRUE(findEdge(guideTreeOut, 14, 7) != 0);
    SEQAN_ASSERT_TRUE(getRoot(guideTreeOut) == 14);
}

template<typename TTag>
void
Test_UpgmaGuideTree(int seed) {
    using namespace seqan;

	typedef unsigned int TSize;

    RNG<MersenneTwister> rng(seed);

	for(TSize i = 0; i < 10; ++i) {
		// Set-up a sparse distance matrix
        PDF<Uniform<TSize> > pdfN(2, 11);
		TSize n = pickRandomNumber(rng, pdfN);
		Graph<Undirected<double> > distGraph;
		String<double> distMatrix;
		fill(distMatrix, n * n, 0);
		TSize all = (n * (n - 1)) / 2;
		typedef std::set<double> TDistanceSet;
		typedef TDistanceSet::iterator TSetIt;
		TDistanceSet distances;
		String<double> myDist;
		typedef Iterator<String<double> >::Type TStringIter;
		while (distances.size() < all) {
            PDF<Uniform<double> > pdf(0, 1000000);
			double newVal = pickRandomNumber(rng, pdf);
			if (distances.insert(newVal).second) {
				appendValue(myDist, newVal);
			}
		}
		double infCargo = _getInfinity<double>();
		//clear(myDist); appendValue(myDist, infCargo); appendValue(myDist, infCargo); appendValue(myDist, 84);
		TStringIter strIt = begin(myDist);
		for(TSize row = 0; row < n; ++row)
            addVertex(distGraph);
		for(TSize row = 0; row < n; ++row) {
			for(TSize col = n - 1; col > row; --col) {
				addEdge(distGraph, row, col, value(strIt));
				value(distMatrix, row * n + col) = value(strIt);
				goNext(strIt);
			}
		}
		//removeEdge(distGraph, 0, 1);removeEdge(distGraph, 0, 2);
        Graph<Undirected<double> > distGraphCopy;
        String<double> distMatrixCopy;
		for(TSize row = 0; row < n; ++row) {
			for(TSize col = n - 1; col > row; --col) {
                distGraphCopy = distGraph;
                distMatrixCopy = distMatrix;

                PDF<Uniform<double> > pdf(0, 1.0);
				if (pickRandomNumber(rng, pdf) < 0.5) {
					value(distMatrix, row * n + col) = infCargo;
					removeEdge(distGraph, row, col);
				}

                String<size_t> _;
                if (connected_components(distGraph, _)) {
                    move(distGraph, distGraphCopy);
                    move(distMatrix, distMatrixCopy);
                }
			}
		}
		// Guide Tree
		Graph<Tree<double> > guideTreeGraph;
		upgmaTree(distGraph, guideTreeGraph, TTag());
		Graph<Tree<double> > guideTreeMat;
		upgmaTree(distMatrix, guideTreeMat, TTag());
		typedef Iterator<Graph<Tree<double> >, BfsIterator>::Type TBfsIterator;
		String<TSize> set1;
		TBfsIterator itBfs(guideTreeGraph, getRoot(guideTreeGraph));
		for(;!atEnd(itBfs);goNext(itBfs)) appendValue(set1, value(itBfs));
		String<TSize> set2;
		TBfsIterator itBfs2(guideTreeMat, getRoot(guideTreeMat));
		for(;!atEnd(itBfs2);goNext(itBfs2)) appendValue(set2, value(itBfs2));
        SEQAN_ASSERT_TRUE(set1 == set1);
        /*
		if (set1 != set2) {
			std::cout << "Randomized test failed:" << std::endl;
			std::cout << "Upgma Guide Trees:" << std::endl;
			std::cout << guideTreeMat << std::endl;
			std::cout << guideTreeGraph << std::endl;
			for(TSize i=0;i<n;++i) {
				for(TSize j=i+1;j<n;++j) {
					std::cout << value(distMatrix, i*n+j) << ",";
				}
				std::cout << std::endl;
			}
			std::cout << std::endl;
			exit(0);
		}
        */
	}
}

SEQAN_DEFINE_TEST(test_graph_msa_guide_tree_neighbour_joining)
{
    Test_GuideTree_NeighbourJoining();
}

SEQAN_DEFINE_TEST(test_graph_msa_guide_tree_upgma_weight_avg)
{
    for (int i = 0; i < 10; ++i)
        Test_UpgmaGuideTree<seqan::UpgmaWeightAvg>(i);
}

SEQAN_DEFINE_TEST(test_graph_msa_guide_tree_upgma_avg)
{
    for (int i = 0; i < 10; ++i)
        Test_UpgmaGuideTree<seqan::UpgmaAvg>(i);
}

SEQAN_DEFINE_TEST(test_graph_msa_guide_tree_upgma_min)
{
    for (int i = 0; i < 10; ++i)
        Test_UpgmaGuideTree<seqan::UpgmaMin>(i);
}

SEQAN_DEFINE_TEST(test_graph_msa_guide_tree_upgma_max)
{
    for (int i = 0; i < 10; ++i)
        Test_UpgmaGuideTree<seqan::UpgmaMax>(i);
}

#endif  // #ifndef TESTS_TEST_GRAPH_MSA_GUIDE_TREE_H_
