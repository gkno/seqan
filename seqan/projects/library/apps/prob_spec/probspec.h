 /*==========================================================================
                     PROBSPEC - Special local alignments (based on STELLAR)

 ============================================================================
  Copyright (C) 2011 by Knut Reinert

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your options) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ==========================================================================*/

#ifndef SEQAN_HEADER_PROBSPEC_H
#define SEQAN_HEADER_PROBSPEC_H

#include <iostream>
#include <seqan/index.h>
#include <seqan/seeds2.h>
#include "stellar_types.h"
#include "stellar_extension.h"

using namespace seqan;


template<typename TInfix, typename TSize, typename TEps, typename TAlign>
bool
_extendKmer(TInfix & a, TInfix & b, TSize minLength, TEps eps, TAlign & align) {
	typedef Seed<Simple> TSeed;
	typedef int TScore;
	
	TSeed seed(beginPosition(a), beginPosition(b), endPosition(a), endPosition(b));
	TSeed seedOld(seed);
	
	TScore penalty = static_cast<TScore>(ceil(-1/eps) + 1);
	Score<TScore> scoreMatrix(1, penalty, penalty);
	TScore scoreDropOff = -penalty * static_cast<TScore>(minLength * eps);
	
	
	ExtensionDirection direction = EXTEND_BOTH;
	extendSeed(seed, host(a), host(b), direction, scoreMatrix, scoreDropOff, GappedXDrop());
	
	return _bestExtension(a, b, seed, seedOld, unsigned(length(a)), 0u, scoreMatrix, direction, minLength, eps, align);
}




template<typename TSequence>
inline bool
_extendExactMatches(TSequence & database,
					StringSet<TSequence> & queries,
					StringSet<String<Align<TSequence > > >& result,
					ProbSpecOptions & options) {
	
	StringSet<String<Triple<unsigned,unsigned,unsigned> > > ematches;
	resize(ematches, length(queries));
	
	// built qGram index on database
	typedef Index<TSequence, IndexQGram<SimpleShape, OpenAddressing> > TQGramIndex;
    TQGramIndex qgramIndex(database);
	resize(indexShape(qgramIndex), options.lengthExact);
	
	std::cout << "Built qgram index for exact submatches" << ::std::endl;
	
	indexRequire(qgramIndex, QGramSADir());
	indexRequire(qgramIndex, QGramDir());
	
	typedef typename Iterator<TSequence>::Type TSequenceIt;
	
	for (unsigned int s=0; s < length(queries); s++) {
		//::std::cout << "Inspecting query " << s << " " << queries[s] << std::endl;
		
		TSequenceIt sit  = begin(queries[s]);
		TSequenceIt send = end(queries[s])-options.lengthExact;
		
		unsigned hv = hash(indexShape(qgramIndex),sit);
		typedef typename Infix<typename Fibre<TQGramIndex, QGramSA>::Type const>::Type TInfix;
		TInfix occs;
		
		int qp=0;
		for (; sit != send; ) { 
			occs = getOccurrences(qgramIndex, indexShape(qgramIndex));
			for (unsigned i = 0; i < length(occs); i++){
				//		::std::cout << s << " " << occs[i] << ::std::endl;
				Triple<unsigned, unsigned, unsigned> p(occs[i],s,qp);
				//		::std::cout << p << ::std::endl;
				appendValue(ematches[s], p);
			}
			sit++;
			qp++;
			hv = hashNext(indexShape(qgramIndex),sit);
			/*			
			 TSequence result;
			 unhash(result, hv, options.lengthExact);	
			 ::std::cout << result <<  ::std::endl;
			 */			
		}
	}
	
	// remove adjacent q-gram hits
	
	std::cout << "Filtering adjacent qgram hits" << ::std::endl;
	
	StringSet< String<Triple<unsigned,unsigned,unsigned> > > cematches;
	resize(cematches,length(ematches));
	
	
	for (unsigned s=0; s<length(ematches); s++) {
		unsigned cind = 0;		
		if( length(ematches[s]) > 0 ){
			std::sort(begin(ematches[s]), end(ematches[s]));
			resize(cematches[s],1);
			cematches[s][cind++] = ematches[s][0];
			//	::std::cout << "CMATCH "<< ematches[s][0] << ::std::endl;
			
		}
		if( length(ematches[s]) > 1 ){
			for (unsigned i=1; i<length(ematches[s]); ++i) {
				if( ematches[s][i].i2 != ematches[s][i-1].i2 || ematches[s][i].i1 != ematches[s][i-1].i1+1 ){
					appendValue(cematches[s],ematches[s][i]);
					//::std::cout << "CMATCH "<< ematches[s][i] << ::std::endl;
				}
			}
		}
		//		::std::cout << "size of cematches " << cind <<  " " << length(cematches[s]) << std::endl;
	}
	
	
	
	typedef typename Infix<TSequence>::Type TInfix;	
	TSequenceIt dbbegin = begin(database);
	
	Align<TSequence>  oldalign;
	resize(rows(oldalign), 2);
	
	std::cout << "Compute alignment extension and remove duplicate alignments" << ::std::endl;
	
	for (unsigned s=0; s<length(cematches); s++) {
		TSequenceIt qbegin = begin(queries[s]);
		
		// there are extensions to be made
		for (unsigned i=0; i<length(cematches[s]); i++) {
			unsigned dbpos = cematches[s][i].i1;
			unsigned qpos = cematches[s][i].i3;
			
			TInfix dbinfix = infix(database,dbbegin + dbpos,dbbegin + dbpos + options.lengthExact);
			TInfix qinfix = infix(queries[s],qbegin + qpos,qbegin + qpos + options.lengthExact);
			
			/*	::std::cout << "positions " << dbpos << " " << qpos << ::std::endl;	
			 ::std::cout << dbinfix << ::std::endl;
			 ::std::cout << qinfix << ::std::endl;
			 */
			Align<TSequence> align;
			resize(rows(align), 2);
			setSource(row(align, 0), host(dbinfix));
			setSource(row(align, 1), host(qinfix));
			
			
			if(_extendKmer(dbinfix, qinfix, options.minLength, options.epsilon, align)){			
				
				if(length(row(align,0))  > unsigned(options.minLength)){
					unsigned db1 = clippedBeginPosition(row(align,0));
					unsigned db2 = clippedEndPosition(row(align,0))-1;
					unsigned q1  = clippedBeginPosition(row(align,1));
					unsigned q2  = clippedEndPosition(row(align,1))-1;
					unsigned odb1 = clippedBeginPosition(row(oldalign,0));
					unsigned odb2 = clippedEndPosition(row(oldalign,0))-1;
					unsigned oq1  = clippedBeginPosition(row(oldalign,1));
					unsigned oq2  = clippedEndPosition(row(oldalign,1))-1;
					
					//		::std::cout << "Old Aligns Seq1[" << odb1 << ":" << odb2<< "]" << " and Seq2[" << oq1 << ":" <<  oq2 << "]";
					//			::std::cout << ::std::endl; 
					
					//	::std::cout << oldalign << ::std::endl;
					
					//		::std::cout << "Aligns Seq1[" << db1 << ":" << db2<< "]" << " and Seq2[" << q1 << ":" <<  q2 << "]";
					//		::std::cout << ::std::endl; 
					
					//		::std::cout << align << ::std::endl;
					
					if (odb1 != db1 || odb2 != db2 || oq1 != q1 || oq2 != q2 ) {
						appendValue(result[s],align);
						oldalign = align;
					}
					//		else	
					//			::std::cout << "double match " << ::std::endl;
					
				}
			}
		}
	}
	
	return true;
}


#endif
