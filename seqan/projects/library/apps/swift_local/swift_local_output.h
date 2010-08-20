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

#include <iostream>
#include <seqan/align.h>
//#include "swift_local_types.h"

using namespace seqan;

template<typename TAlign, typename TString>
void
_getCigarLine(TAlign const & align, TString & cigar, TString & mutations) { 
    typedef typename Size<typename Row<TAlign>::Type >::Type TSize;

    TSize dbPos = beginPosition(row(align, 0));
    TSize queryPos = beginPosition(row(align, 1));

    TSize dbEndPos = endPosition(row(align, 0));
    TSize queryEndPos = endPosition(row(align, 1));

    bool first = true;
    TSize readBasePos = queryPos;
    TSize readPos = 0;
	while (dbPos != dbEndPos && queryPos != queryEndPos) {
		int matched = 0;
		int inserted = 0;
		int deleted = 0;
		while (dbPos != dbEndPos && queryPos != queryEndPos &&
               !isGap(row(align, 0), dbPos) && !isGap(row(align, 1), queryPos)) {
            ++readPos;
			if (value(row(align, 0), dbPos) != value(row(align, 1), queryPos)) {
				if (first) first = false;
				else mutations << ",";
				mutations << readPos << value(source(row(align, 1)), readBasePos);
			}
			++readBasePos;
			++dbPos;
			++queryPos;
			++matched;
		}
		if (matched > 0) cigar << matched << "M" ;
		while (queryPos != queryEndPos && isGap(row(align, 1), queryPos)) {
			++dbPos;
			++queryPos;
			++deleted;
		}
		if (deleted > 0) cigar << deleted << "D";
		while (dbPos != dbEndPos && isGap(row(align, 0), dbPos)) {
			++dbPos;
			++queryPos;
			++readPos;
			if (first) first = false;
			else mutations << ",";
			mutations << readPos << value(source(row(align, 1)), readBasePos);
			++readBasePos;
			++inserted;
		}
		if (inserted > 0) cigar << inserted << "I";
	}
}

template<typename TAlign>
double
_calculateIdentity(TAlign const & align) {
    typedef typename Size<typename Row<TAlign>::Type >::Type TSize;
    TSize matches = 0;
    TSize len = _max(length(row(align, 0)), length(row(align, 1)));

    TSize pos0 = beginPosition(row(align, 0));
    TSize pos1 = beginPosition(row(align, 1));

    TSize end0 = endPosition(row(align, 0));
    TSize end1 = endPosition(row(align, 1));

    while ((pos0 < end0) && (pos1 < end1)) {
        if (!isGap(row(align, 0), pos0) && !isGap(row(align, 1), pos1)) {
            if (value(row(align, 0), pos0) == value(row(align , 1), pos1)) {
                ++matches;
            }
        }
        ++pos0;
        ++pos1;
    }

    return ((double)matches/(double)len)*100.0;
}

template<typename TId, typename TAlign, typename TFile>
void
_writeGffLine(TId const & databaseID,
              TId const & patternID,
              bool const databaseStrand,
              TAlign const & match,
              TFile & file) {
	typedef typename Row<TAlign>::Type TRow;
	TRow row0 = row(match, 0);
	TRow row1 = row(match, 1);
    
    for (typename Position<TId>::Type i = 0; i < length(databaseID) && value(databaseID, i) > 32; ++i) {
        file << value(databaseID, i);
    }

    file << "\tSwiftLocal";
    file << "\teps-matches";

    if (databaseStrand) {
        file << "\t" << 
			toSourcePosition(row0, beginPosition(row0)) + beginPosition(source(row0)) + 1;
        file << "\t" << 
			toSourcePosition(row0, endPosition(row0)) + beginPosition(source(row0));
    } else {
        file << "\t" << length(source(row0)) - 
            (toSourcePosition(row0, endPosition(row0)) + beginPosition(source(row0))) + 1;
        file << "\t" << length(source(row0)) - 
            (toSourcePosition(row0, beginPosition(row0)) + beginPosition(source(row0)));
    }

    file << "\t" << _calculateIdentity(match);

    file << "\t" << (databaseStrand ? '+' : '-');

    file << "\t.\t";
    for (typename Position<TId>::Type i = 0; i < length(patternID) && value(patternID, i) > 32; ++i) {
        file << value(patternID, i);
    }

	//file << ";seq2Length=" << length(source(row1));

    file << ";seq2Range=" << 
		toSourcePosition(row1, beginPosition(row1)) + beginPosition(source(row1)) + 1;
    file << "," << 
		toSourcePosition(row1, endPosition(row1)) + beginPosition(source(row1));

    std::stringstream cigar, mutations;
    _getCigarLine(match, cigar, mutations);
    file << ";cigar=" << cigar.str();
    file << ";mutations=" << mutations.str();
    file << "\n";
}

template<typename TAlign, typename TStrand, typename TFile>
void
_writeMatch(TAlign & match,
			TStrand databaseStrand,
			TFile & aliFile) {
	typedef typename Row<TAlign>::Type TRow;
	TRow row0 = row(match, 0);
	TRow row1 = row(match, 1);

	// write database positions
	if (databaseStrand) {
		aliFile << "< " <<
			toSourcePosition(row0, beginPosition(row0)) + beginPosition(source(row0));
		aliFile << " , " << 
			toSourcePosition(row0, endPosition(row0)) + beginPosition(source(row0));
	} else {
		aliFile << "< " << length(source(row0)) - 
			toSourcePosition(row0, beginPosition(row0)) + beginPosition(source(row0));
		aliFile << " , " << length(source(row0)) -
			toSourcePosition(row0, endPosition(row0)) + beginPosition(source(row0));
	}
	// write query positions
	aliFile << " >< " << 
		toSourcePosition(row1, beginPosition(row1)) + beginPosition(source(row1));
	aliFile << " , " << 
		toSourcePosition(row1, endPosition(row1)) + beginPosition(source(row1)) << " >\n";

	// write match
	aliFile << match;
}

template<typename TInfix, typename TQueryId, typename TNumber, typename TId, typename TIds, typename TFile>
void
_outputMatches(//StringSet<String<Align<TInfix> > > const & matches, 
			   StringSet<String<SwiftLocalMatch<TInfix, TQueryId> > > const & matches, 
			   TNumber const /*numSwiftHits*/,
			   TId const & databaseID,
			   bool const databaseStrand,
			   TIds const & ids,
			   TFile & file) {
	typedef SwiftLocalMatch<TInfix, TQueryId> TMatch;
	typedef typename Size<typename TMatch::TAlign>::Type TSize;

    //std::ofstream aliFile;
    //aliFile.open("swift_local.align"/*, std::ios::app*/);

    //TSize maxLength = 0;
    //TSize totalLength = 0;
    //TSize numMatches = 0;

    //aliFile << "Database sequence: " << databaseID;
    //if (!databaseStrand) aliFile << " complement\n";
    //else aliFile << "\n";

    for (TSize i = 0; i < length(matches); i++) {
        if (length(value(matches, i)) == 0) continue;
        //std::cout << "Pattern sequence: " << ids[i] << "\n";
        //aliFile << "Pattern sequence: " << ids[i] << "\n\n";
        for (TSize j = 0; j < length(value(matches, i)); j++) {
            TMatch m = value(value(matches, i), j);

            //TSize len = _max(length(row(m.align, 0)), length(row(m.align, 1)));
            //totalLength += len;
            //if(len > maxLength) maxLength = len;

            _writeGffLine(databaseID, ids[i], databaseStrand, m.align, file);
			//_writeMatch(m.align, databaseStrand, aliFile);
        }
        //numMatches += length(value(matches, i));
        //std::cout << "  # Eps-matches: " << length(value(matches, i)) << std::endl;
    }

	//if (numMatches > 0) {
	//	std::cout << "    Longest eps-match: " << maxLength << std::endl;
 //       std::cout << "    Avg match length : " << totalLength / numMatches << std::endl;
	//}
 //   std::cout << "    # SWIFT hits     : " << numSwiftHits << std::endl;
 //   std::cout << "    # Eps-matches    : " << numMatches << std::endl;

    //aliFile.close();
}