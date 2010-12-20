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
 Author: David Weese <david.weese@fu-berlin.de>
 ==========================================================================*/

#ifndef TESTS_INDEX_TEST_INDEX_CREATION_H
#define TESTS_INDEX_TEST_INDEX_CREATION_H

#define SEQAN_PROFILE
//#define SEQAN_DEBUG
//#define SEQAN_DEBUG_INDEX

//#define SEQAN_TEST
//#define SEQAN_TEST_SKEW3
//#define SEQAN_TEST_SKEW7


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

SEQAN_DEFINE_TEST(testIndexCreation)
{
		typedef String<char> TText;
		typedef String<unsigned> TArray;

		TText	text;
		TArray	sa;
		TArray	lcp;
		TArray	child, childExt;
		TText	bwt;

		const int runs = 2;					// conduct 10 test runs 
		const int maxSize = 20 * 1024 * 1024;	// max text size is 20 megabyte
		bool result = true;

		_proFloat timeDelta[12];
		_proFloat timeSum[12];
		for(int i = 0; i < 10; ++i)
			timeSum[i] = 0;
		__int64 textSum = 0;

		static const char* algNames[] = {
			"Skew3        ", 
			"Skew7        ", 
			"ManberMyers  ", 
			"LarssonSadake", 
			"SAQSort      ", 
			"SAQSortQGSR  ", 
			"Skew3Ext     ", 
			"Skew7Ext     ",
			"Kasai        ",
			"KasaiInPlace ",
			"KasaiExt     ",
			"ChildTab     ",
			"ChildTabExt  "
		};

		int TI;
		for(int i = 0; i < runs; ++i) {

			std::cout << "*** RUN " << i << " ***";
			
			int size = rand() % maxSize;
			TI = 0;

//___randomize_text___________________________________________________________

			resize(text,size);
/*			if (i < runs/2)
				randomize(text);
			else
*/				textRandomize(text);
/*			String<char,External<> > errorText;	// read in text causing an error
			open(errorText,"error.txt");
			text = errorText;
*/
/*			text = "MISSISSIPPI";
			size = length(text);
			std::cout << "text created (n=" << size << ")" << std::endl;
*/
			std::cout << "   textSize: " << length(text) << std::endl;

//___create_suffix_array______________________________________________________

			resize(sa, size);
/*			timeDelta[TI] = -SEQAN_PROGETTIME;
			createSuffixArray(sa, text, Skew3());
			timeDelta[TI++] += SEQAN_PROGETTIME;
			if (!isSuffixArray(sa, text)) {
				std::cout << "suffix array creation (internal Skew3) failed" << std::endl;
				result = false;
			}
*/			std::cout << "."; std::cout.flush();

			blank(sa);
/*			timeDelta[TI] = -SEQAN_PROGETTIME;
			createSuffixArray(sa, text, Skew7());
			timeDelta[TI++] += SEQAN_PROGETTIME;
			if (!isSuffixArray(sa, text)) {
				std::cout << "suffix array creation (internal Skew7) failed" << std::endl;
				result = false;
			}
*/			std::cout << "."; std::cout.flush();

			blank(sa);
/*			timeDelta[TI] = -SEQAN_PROGETTIME;
			createSuffixArray(sa, text, ManberMyers());
			timeDelta[TI++] += SEQAN_PROGETTIME;
			if (!isSuffixArray(sa, text)) {
				std::cout << "suffix array creation (internal ManberMyers) failed" << std::endl;
				result = false;
			}
*/			std::cout << "."; std::cout.flush();

			blank(sa);
/*			timeDelta[TI] = -SEQAN_PROGETTIME;
			createSuffixArrayExt(sa, text, LarssonSadakane());
			timeDelta[TI++] += SEQAN_PROGETTIME;
			if (!isSuffixArray(sa, text)) {
				std::cout << "suffix array creation (external LarssonSadakane) failed" << std::endl;
				result = false;
			}
*/			std::cout << "."; std::cout.flush();

			blank(sa);
			timeDelta[TI] = -SEQAN_PROGETTIME;
			createSuffixArray(sa, text, SAQSort());
			timeDelta[TI++] += SEQAN_PROGETTIME;
			if (!isSuffixArray(sa, text)) {
				std::cout << "suffix array creation (internal SAQSort) failed" << std::endl;
				result = false;
			}
			std::cout << "."; std::cout.flush();
/*
			blank(sa);
			timeDelta[TI] = -SEQAN_PROGETTIME;
			createSuffixArray(sa, text, QSQGSR(), 3);
			timeDelta[TI++] += SEQAN_PROGETTIME;
			if (!isSuffixArray(sa, text)) {
				std::cout << "suffix array creation (internal QSQGSR) failed" << std::endl;
				result = false;
			}
			std::cout << "."; std::cout.flush();
*/
			blank(sa);
/*			timeDelta[TI] = -SEQAN_PROGETTIME;
			createSuffixArrayExt(sa, text, Skew3());
			timeDelta[TI++] += SEQAN_PROGETTIME;
			if (!isSuffixArray(sa, text)) {
				std::cout << "suffix array creation (external Skew3) failed" << std::endl;
				result = false;
			}
			std::cout << "."; std::cout.flush();

			blank(sa);
			timeDelta[TI] = -SEQAN_PROGETTIME;
			createSuffixArrayExt(sa, text, Skew7());
			timeDelta[TI++] += SEQAN_PROGETTIME;
			if (!isSuffixArray(sa, text)) {
				std::cout << "suffix array creation (external Skew7) failed" << std::endl;
				result = false;
			}
			std::cout << "."; std::cout.flush();

//___create_lcp_table_________________________________________________________

			resize(lcp, size);
			timeDelta[TI] = -SEQAN_PROGETTIME;
			createLcpTable(lcp, text, sa, KasaiOriginal());
			timeDelta[TI++] += SEQAN_PROGETTIME;
			if (!isLCPTable(lcp, sa, text)) {
				std::cout << "suffix array creation (internal Kasai) failed" << std::endl;
				result = false;
			}
			std::cout << "."; std::cout.flush();

			blank(lcp);
			timeDelta[TI] = -SEQAN_PROGETTIME;
			createLcpTable(lcp, text, sa, Kasai());
			timeDelta[TI++] += SEQAN_PROGETTIME;
			if (!isLCPTable(lcp, sa, text)) {
				std::cout << "suffix array creation (internal in-place Kasai) failed" << std::endl;
				result = false;
			}
			std::cout << "."; std::cout.flush();

			blank(lcp);
			timeDelta[TI] = -SEQAN_PROGETTIME;
			createLcpTableExt(lcp, text, sa, Kasai());
			timeDelta[TI++] += SEQAN_PROGETTIME;
			if (!isLCPTable(lcp, sa, text)) {
				std::cout << "suffix array creation (external Kasai) failed" << std::endl;
				result = false;
			}
			std::cout << "."; std::cout.flush();

//___create_child_table_______________________________________________________

			resize(child, size);
			for(int i=0; i<size; ++i)
				child[i] = maxValue<unsigned>();
			timeDelta[TI] = -SEQAN_PROGETTIME;
			createChildTable(child, lcp);
			timeDelta[TI++] += SEQAN_PROGETTIME;
			std::cout << "."; std::cout.flush();

			unsigned undefs=0;
			for(int i=0; i<size; ++i)
				if (child[i] == maxValue<unsigned>()) ++undefs;
			if (undefs) ::std::cout << undefs << " undefined values";

			resize(childExt, size);
			timeDelta[TI] = -SEQAN_PROGETTIME;
			createChildTableExt(childExt, lcp);
			timeDelta[TI++] += SEQAN_PROGETTIME;
			std::cout << "."; std::cout.flush();

			if (!isEqual(child, childExt)) {
				std::cout << "child table creation failed" << std::endl;
				result = false;
			}
*/
//___update_performance_table_________________________________________________

			for(int i=0; i<TI; ++i) {
				timeSum[i] += timeDelta[i];
				textSum += length(text);
			}

			std::cout << " OK!" << std::endl;

		}
		std::cout << "*** TIME RESULTS (sec/MB) ***" << std::endl;
		for(int i=0; i<TI; ++i)
			std::cout << algNames[i] << " " << 1024.0*1024.0 * timeSum[i] / textSum << std::endl;
}

//////////////////////////////////////////////////////////////////////////////


} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
