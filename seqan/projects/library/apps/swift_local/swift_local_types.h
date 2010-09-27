 /*==========================================================================
                     SwiftLocal - Fast Local Alignment

 ============================================================================
  Copyright (C) 2010 by Birte Kehr

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

#include <seqan/align.h>
#include <seqan/store.h>

using namespace seqan;

struct SwiftLocalOptions {
	// i/o options
	CharString databaseFile;	// name of database file
	CharString queryFile;		// name of query file
	CharString outputFile;		// name of result file
	unsigned outputFormat;		// 1..gff
								// 2..??

	// main options
	unsigned qGram;				// length of the q-grams
	double epsilon;				// maximal error rate
	int minLength;				// minimal length of an epsilon-match
	double xDrop;				// maximal x-drop

	// more options
	bool reverse;				// compute also matches to reverse complemented database
	CharString fastOption;		// verification strategy: exact, bestLocal, bandedGlobal
	unsigned compactThresh;			// number of matches after which removal of overlaps and duplicates is started
	unsigned numMatches;		// maximal number of matches per query and database


	SwiftLocalOptions() {
		outputFile = "swift_local.gff";
		outputFormat = 1;

		qGram = 10;
		epsilon = 0.05;
		minLength = 100;
		xDrop = 5;

		reverse = false;
		fastOption = "exact";		// exact verification
		compactThresh = 500;
		numMatches = 50;
	}
}; 

template<typename _TSequence, typename _TId>
struct SwiftLocalMatch {
	typedef _TSequence							TSequence;
	typedef _TId								TId;
	typedef typename Position<TSequence>::Type	TPos;

	typedef Align<TSequence>					TAlign;
	typedef typename Row<TAlign>::Type			TRow;

	static const TId INVALID_ID;

	TPos begin1;
	TPos end1;
	TRow row1;

	TId id;
	TPos begin2;
	TPos end2;
	TRow row2;

	SwiftLocalMatch() {}

	template<typename TAlign, typename TId>
	SwiftLocalMatch(TAlign & _align, TId _id) {
		id = _id;
		row1 = row(_align, 0);
		row2 = row(_align, 1);

		begin1 = beginPosition(sourceSegment(row1));
		end1 = endPosition(sourceSegment(row1));

		begin2 = beginPosition(sourceSegment(row2));
		end2 = endPosition(sourceSegment(row2));
	}
};

template <typename TSequence, typename TId> 
const TId
SwiftLocalMatch<TSequence, TId>::INVALID_ID = "###########";


// to sort matches by position and remove overlapping matches
template <typename TMatch>
struct LessPos : public ::std::binary_function <TMatch, TMatch, bool> {		
	LessPos() {}
	
	inline int compare(TMatch const & a, TMatch const & b) const {
		// query number
		if ((a.id) < (b.id)) return -1;
		if ((a.id) > (b.id)) return 1;

		// database begin position
		typename TMatch::TPos aBegin1 = _min(a.begin1, a.end1);
		typename TMatch::TPos bBegin1 = _min(b.begin1, b.end1);
		if (aBegin1 < bBegin1) return -1;
		if (aBegin1 > bBegin1) return 1;

		// database end position
		typename TMatch::TPos aEnd1 = _max(a.begin1, a.end1);
		typename TMatch::TPos bEnd1 = _max(b.begin1, b.end1);
		if (aEnd1 < bEnd1) return -1;
		if (aEnd1 > bEnd1) return 1;

		// query begin position
		typename TMatch::TPos aBegin2 = _min(a.begin2, a.end2);
		typename TMatch::TPos bBegin2 = _min(b.begin2, b.end2);
		if (aBegin2 < bBegin2) return -1;
		if (aBegin2 > bBegin2) return 1;

		// query end position
		typename TMatch::TPos aEnd2 = _max(a.begin2, a.end2);
		typename TMatch::TPos bEnd2 = _max(b.begin2, b.end2);
		if (aEnd2 < bEnd2) return -1;
		if (aEnd2 > bEnd2) return 1;

		//// orientation
		//bool oa = a.begin1 < a.end1;
		//bool ob = b.begin1 < b.end1;
		//if (oa != ob) return oa;

		return 0;
	}
		
	inline bool operator() (TMatch const & a, TMatch const & b) const {
		return compare(a, b) == -1;
	}
};

// to sort matches by length
template <typename TMatch>
struct LessLength : public ::std::binary_function <TMatch, TMatch, bool> {		
	LessLength() {}

	inline int compare(TMatch const & a, TMatch const & b) const {
		typename TMatch::TPos aLength = abs((int)a.end1 - (int)a.begin1);
		typename TMatch::TPos bLength = abs((int)b.end1 - (int)b.begin1);

		if (a.id == TMatch::INVALID_ID) return 1;
		if (b.id == TMatch::INVALID_ID) return -1;

		if (aLength < bLength) return 1;
		if (aLength > bLength) return -1;

		return 0;
	}

	inline bool operator() (TMatch const & a, TMatch const & b) const {
		return compare(a, b) == -1;
	}
};

template<typename TSequence, typename TId, typename TRowNo, typename TSize>
inline bool
_isUpstream(SwiftLocalMatch<TSequence, TId> & match1, SwiftLocalMatch<TSequence, TId> & match2, TRowNo row, TSize minLength) {
SEQAN_CHECKPOINT 
	typedef typename SwiftLocalMatch<TSequence, TId>::TPos TPos;

	TPos e1, b2;
	if (row == 0) {
		e1 = match1.end1;
		b2 = match2.begin1;
	} else {
		e1 = match1.end2;
		b2 = match2.begin2;
	}

    if (e1 <= b2) return true;
    
	TPos b1, e2;
	if (row == 0) {
		e2 = match2.end1;
		b1 = match1.begin1;
	} else {
		e2 = match2.end2;
		b1 = match1.begin2;
	}

    if (b1 < b2 && (b2 - b1 >= minLength)) {
        if ((e1 < e2) && (e2 - e1 >= minLength)) return true;
    }
    
    return false;
}

template <typename TMatches, typename TFunctorLess>
inline void
sortMatches(TMatches & swiftLocalMatches, TFunctorLess const & less) {
	std::stable_sort(
		begin(swiftLocalMatches, Standard()), 
		end(swiftLocalMatches, Standard()), 
		less);
}

template <typename TSequence, typename TId>
inline typename Size<TSequence>::Type
length(SwiftLocalMatch<TSequence, TId> & match) {
	return _max(length(match.row1), length(match.row2));
}