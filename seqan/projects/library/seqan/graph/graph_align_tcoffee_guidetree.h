 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_GUIDETREE_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_GUIDETREE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// T-Coffee - Guide Tree
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSpec, typename TCargo, typename TSpec>
inline void
slowNjTree(String<double, TStringSpec>& mat, 
		   Graph<Tree<TCargo, TSpec> >& g) 
{
	SEQAN_CHECKPOINT
	
	typedef typename Size<String<double, TStringSpec> >::Type TSize;
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	
	TSize nseq = (TSize) std::sqrt((double)length(mat));

	//for(TSize i=0;i<nseq;++i) {
	//	for(TSize j=0;j<nseq;++j) {
	//		std::cout << getValue(mat, i*nseq+j) << ",";
	//	}
	//	std::cout << std::endl;
	//}

	// First initialization
	clearVertices(g);
	String<double> av;    // Average branch length to a combined node
	fill(av,nseq,0.0);

	String<TVertexDescriptor> connector;   // Nodes that need to be connected
	resize(connector, nseq);

	for(TSize i=0;i<nseq;++i) {
		addVertex(g);  // Add all the nodes that correspond to sequences
		assignValue(connector, i, i);
		assignValue(mat, i*nseq+i, 0.0);
	}

	// Main cycle
	double fnseqs=(double) nseq;
	for(TSize nc=0; nc<(nseq-3); ++nc) {
		double sumOfBranches = 0.0;

		// Copy upper triangle matrix to lower triangle
		for(TSize col=1; col<nseq; ++col) {
			for(TSize row=0; row<col; ++row) {
				assignValue(mat, col*nseq+row, getValue(mat, row*nseq+col));
				// Determine the sum of all branches
				sumOfBranches = sumOfBranches + getValue(mat, row*nseq+col);
			}
		}

		// Compute the sum of branch lengths for all possible pairs
		double tmin = 0.0;	
		TSize mini = 0;  // Next pair of seq i and j to join
		TSize minj = 0;
		for(TSize col=1; col<nseq; ++col)  {
			if (getValue(connector,col) != nilVertex) {
				for(TSize row=0; row<col; ++row) {
					if (getValue(connector,row) != nilVertex) {
						double diToAllOthers = 0.0;
						double djToAllOthers = 0.0;
						
						for(TSize i=0; i<nseq; ++i) {
							diToAllOthers += getValue(mat, i*nseq+row);
							djToAllOthers += getValue(mat, i*nseq+col);
						}

						double dij = getValue(mat, row*nseq+col);
						double total = diToAllOthers + djToAllOthers + (fnseqs - 2.0)*dij +2.0*(sumOfBranches - diToAllOthers - djToAllOthers);
						total /= (2.0*(fnseqs - 2.0));

						if ((tmin == 0) || (total < tmin)) {
							tmin = total;
							mini = row;
							minj = col;
						}
					}
				}
			}
		}

		// Print nodes that are about to be joined
		//std::cout << mini << std::endl;
		//std::cout << minj << std::endl;
		//std::cout << tmin << std::endl;
		//std::cout << std::endl;
		
		// Compute branch lengths
		double dMinIToOthers = 0.0;
		double dMinJToOthers = 0.0;
		for(TSize i=0; i<nseq; ++i) {
			dMinIToOthers += getValue(mat, i*nseq + mini);
			dMinJToOthers += getValue(mat, i*nseq + minj);
		}
		double dmin = getValue(mat, mini*nseq + minj);
		dMinIToOthers = dMinIToOthers / (fnseqs - 2.0);
		dMinJToOthers = dMinJToOthers / (fnseqs - 2.0);
		double iBranch = (dmin + dMinIToOthers - dMinJToOthers) * 0.5;
		double jBranch = dmin - iBranch;
		iBranch -= av[mini];
		jBranch -= av[minj];
		
		// Set negative branch length to zero
		if( fabs(iBranch) < 0.0001) iBranch = 0.0;
		if( fabs(jBranch) < 0.0001) jBranch = 0.0;
	
		// Print branch lengths
		//std::cout << iBranch << std::endl;
		//std::cout << jBranch << std::endl;
		//std::cout << std::endl;
		
		// Build tree
		TVertexDescriptor internalVertex = addVertex(g);
		addEdge(g, internalVertex, getValue(connector, mini), (TCargo) iBranch);
		addEdge(g, internalVertex, getValue(connector, minj), (TCargo) jBranch);

		// Remember the average branch length for the new combined node
		// Must be subtracted from all branches that include this node
		if(dmin <= 0.0) dmin = 0.000001;
		av[mini] = dmin * 0.5;


		// Re-initialisation
		// mini becomes the new combined node, minj is killed
		fnseqs = fnseqs - 1.0;
		assignValue(connector, minj, nilVertex);
		assignValue(connector, mini, internalVertex);

		for(TSize j=0; j<nseq; ++j) {
			if( getValue(connector, j) != nilVertex ) {
				double minIminJToOther = ( getValue(mat, mini*nseq+j) + getValue(mat, minj*nseq+j)) * 0.5;
				// Use upper triangle
				if((TSize) mini < j) assignValue(mat, mini*nseq+j, minIminJToOther);
				if((TSize) mini > j) assignValue(mat, j*nseq+mini, minIminJToOther);
			}
		}
		for(TSize j=0; j<nseq; ++j) {
			assignValue(mat, minj*nseq+j, 0.0);
			assignValue(mat, j*nseq+minj, 0.0);
		}
	}

	// Only three nodes left

	// Find the remaining nodes
	String<TSize> l;
	fill(l,3,0);
	TSize count = 0;
	for(TSize i=0; i<nseq; ++i) {
		if(getValue(connector, i) != nilVertex) {
			l[count] = i;
			++count;
		}
	}

	// Remaining nodes
	//std::cout << l[0] << std::endl;
	//std::cout << l[1] << std::endl;
	//std::cout << l[2] << std::endl;
	//std::cout << std::endl;

	String<double> branch;
	resize(branch, 3);
	branch[0] = (getValue(mat, l[0]*nseq+l[1]) + getValue(mat, l[0]*nseq+l[2]) - getValue(mat, l[1]*nseq+l[2])) * 0.5;
	branch[1] = (getValue(mat, l[1]*nseq+l[2]) + getValue(mat, l[0]*nseq+l[1]) - getValue(mat, l[0]*nseq+l[2])) * 0.5;
	branch[2] =  (getValue(mat, l[0]*nseq+l[2]) + getValue(mat, l[0]*nseq+l[1]) - getValue(mat, l[1]*nseq+l[2])) * 0.5;
    
	branch[0] -= av[l[0]];
	branch[1] -= av[l[1]];
	branch[2] -= av[l[2]];

	// Print branch lengths
	//std::cout << branch[0] << std::endl;
	//std::cout << branch[1] << std::endl;
	//std::cout << branch[2] << std::endl;
	//std::cout << std::endl;
    
	// Reset tiny negative and positive branch lengths to zero
	if( fabs(branch[0]) < 0.0001) branch[0] = 0.0;
	if( fabs(branch[1]) < 0.0001) branch[1] = 0.0;
	if( fabs(branch[2]) < 0.0001) branch[2] = 0.0;
    
	// Build tree
	TVertexDescriptor internalVertex = addVertex(g);
	addEdge(g, internalVertex, getValue(connector, l[0]), (TCargo) branch[0]);
	addEdge(g, internalVertex, getValue(connector, l[1]), (TCargo) branch[1]);
	TVertexDescriptor the_root = addVertex(g);
	addEdge(g, the_root, getValue(connector, l[2]), (TCargo) branch[2]);
	addEdge(g, the_root, internalVertex, (TCargo) 0);
	g.data_root = the_root;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TMatrix, typename TActive, typename TSize>
inline typename Value<TMatrix>::Type
_upgmaTreeMerge(TMatrix const& mat, 
				TActive const& active,
				TSize index_i,
				TSize index_j,
				TSize i,
				TSize nseq,
				UpgmaAvg) 
{
	SEQAN_CHECKPOINT
	typedef typename Value<TMatrix>::Type TValue;

	// Average
	return ((TValue) active[index_i] / (TValue) (active[index_i] + active[index_j])) * getValue(mat, index_i * nseq + i) + ((TValue) active[index_j] / (TValue) (active[index_i] + active[index_j])) * getValue(mat, index_j * nseq + i);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TMatrix, typename TActive, typename TSize>
inline typename Value<TMatrix>::Type
_upgmaTreeMerge(TMatrix const& mat, 
				TActive const& active,
				TSize index_i,
				TSize index_j,
				TSize i,
				TSize nseq,
				UpgmaMin) 
{
	SEQAN_CHECKPOINT
	typedef typename Value<TMatrix>::Type TValue;

	// Minimum
	TValue newDist = ((TValue) active[index_i] / (TValue) (active[index_i] + active[index_j])) * getValue(mat, index_i * nseq + i);
	TValue newDist2 = ((TValue) active[index_j] / (TValue) (active[index_i] + active[index_j])) * getValue(mat, index_j * nseq + i);
	if (newDist2 < newDist) return newDist2;
	else return newDist;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TMatrix, typename TActive, typename TSize>
inline typename Value<TMatrix>::Type
_upgmaTreeMerge(TMatrix const& mat, 
				TActive const& active,
				TSize index_i,
				TSize index_j,
				TSize i,
				TSize nseq,
				UpgmaMax) 
{
	SEQAN_CHECKPOINT
	typedef typename Value<TMatrix>::Type TValue;

	// Minimum
	TValue newDist = ((TValue) active[index_i] / (TValue) (active[index_i] + active[index_j])) * getValue(mat, index_i * nseq + i);
	TValue newDist2 = ((TValue) active[index_j] / (TValue) (active[index_i] + active[index_j])) * getValue(mat, index_j * nseq + i);
	if (newDist2 < newDist) return newDist;
	else return newDist2;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSpec, typename TCargo, typename TSpec, typename TTag>
inline void
upgmaTree(String<double, TStringSpec>& mat, 
		  Graph<Tree<TCargo, TSpec> >& g,
		  TTag) 
{
	SEQAN_CHECKPOINT
	
	typedef typename Size<String<double, TStringSpec> >::Type TSize;
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<String<double, TStringSpec> >::Type TValue;
	
	TSize nseq = (TSize) std::sqrt((TValue)length(mat));

	// First initialization
	TValue minVal = getValue(mat, 1*nseq + 0);
	TSize index_i = 0;
	TSize index_j = 1;
	clearVertices(g);
	// Which entries in the matrix are still active and how many members belong to this group
	String<unsigned int> active;
	fill(active, nseq, 1);
	// Vertex descriptor that represents that entry
	String<TVertexDescriptor> nodes;
	reserve(nodes, nseq);
	for(TSize i=0;i<nseq;++i) {
		for(TSize j=i+1;j<nseq;++j) {
			if (minVal > getValue(mat, i*nseq + j)) {
				minVal = getValue(mat, i*nseq + j);
				index_i = i;
				index_j = j;
			}
		}
		appendValue(nodes, addVertex(g));
	}

	// Property map for sum of weights for each node
	String<TValue> weights;
	fill(weights, nseq, 0);

	// Merge groups
	TSize m = nseq;
	while (m>1) {
		// Merge nodes
		TVertexDescriptor internalNode = addVertex(g);
		if (index_j < index_i) {
			TSize tmp = index_i;
			index_i = index_j;
			index_j = tmp;
		}


		//// Debug code
		//for(TSize i=0;i<nseq;++i) {
		//	if (!active[i]) continue;
		//	for(TSize j=0;j<nseq;++j) {
		//		if (!active[j]) continue;
		//		std::cout << getValue(mat, i*nseq+j) << ",";
		//	}
		//	std::cout << std::endl;
		//}
		//std::cout << minVal << ',' << index_i << ',' << index_j << ',' << std::endl;
		//std::cout << nodes[index_i] << ',' << nodes[index_j] << std::endl;
		//std::cout << std::endl;

		TCargo w = getValue(mat, index_i*nseq + index_j) / 2;
		addEdge(g, internalNode, nodes[index_i], w - getProperty(weights, nodes[index_i]));
		addEdge(g, internalNode, nodes[index_j], w - getProperty(weights, nodes[index_j]));
		appendValue(weights, w);		

		// Get the new distance values
		for(TSize i=0;i<nseq;++i) {
			if (!active[i]) continue;
			else if (i == index_i) continue;
			else if (i == index_j) continue;
			TValue newDist = _upgmaTreeMerge(mat, active, index_i, index_j, i, nseq, TTag());
			assignValue(mat, index_i*nseq + i, newDist);
			assignValue(mat, i*nseq + index_i, newDist);
		}
		// Inactivate one group, adjust the member count for the other one
		active[index_i] = active[index_i] + active[index_j];
		nodes[index_i] = internalNode;
		active[index_j] = 0;

		// Find new minimum
		bool first = true;
		for(TSize i=0;i<nseq;++i) {
			if (!active[i]) continue;
			for(TSize j=i+1;j<nseq;++j) {
				if (!active[j]) continue;
				if ((first) ||
					(minVal > getValue(mat, i*nseq + j))) 
				{
					first = false;
					minVal = getValue(mat, i*nseq + j);
					index_i = i;
					index_j = j;
				}
			}
		}
		--m;
	}
	g.data_root = numVertices(g) - 1;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSpec, typename TCargo, typename TSpec>
inline void
upgmaTree(String<double, TStringSpec>& mat, 
		  Graph<Tree<TCargo, TSpec> >& g) 
{
	SEQAN_CHECKPOINT
	upgmaTree(mat, g, UpgmaAvg());
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
