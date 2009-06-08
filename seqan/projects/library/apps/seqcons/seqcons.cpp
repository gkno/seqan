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
==========================================================================*/

#define SEQAN_PROFILE

#include <seqan/basic.h>

// Profiling
#ifdef SEQAN_PROFILE
		SEQAN_PROTIMESTART(__myProfileTime); 
#endif

#include <seqan/consensus.h>

#include <iostream>
#include <fstream>


using namespace seqan;

//////////////////////////////////////////////////////////////////////////////////

void
printVersion() {
	::std::cerr << "**************************************************" << ::std::endl;
	::std::cerr << "* Consensus Computation                          *" << ::std::endl;
	::std::cerr << "*                                                *" << ::std::endl;
	::std::cerr << "* SeqCons                                        *" << ::std::endl;
	::std::cerr << "* Version: 0.204 (22. May 2009)                  *" << ::std::endl;
	::std::cerr << "**************************************************" << ::std::endl;
	::std::cerr << ::std::endl;
}

//////////////////////////////////////////////////////////////////////////////////


void 
printHelp() {
	::std::cerr <<  "Usage: ./seqcons -reads <FASTA File with Reads> [Options]\n" << ::std::endl;
	::std::cerr <<  "\nOptions\n" << ::std::endl;

	// Main options
	::std::cerr <<  "\nMain Options\n------------\n" << ::std::endl;
	::std::cerr <<  "-reads <FASTA File with Reads>\n" << ::std::endl;
	::std::cerr <<  "\tReads in FASTA format and approximate layout positions.\n\n" << ::std::endl;
	::std::cerr <<  "-afg <AMOS message file>\n" << ::std::endl;
	::std::cerr <<  "\tInput is read from an AMOS message file instead of a fasta file.\n\n" << ::std::endl;
	::std::cerr <<  "-matchlength <Length>\n" << ::std::endl;
	::std::cerr <<  "\tMinimum match-length for an overlap, default is 15.\n\n" << ::std::endl;
	::std::cerr <<  "-quality <Number>\n" << ::std::endl;
	::std::cerr <<  "\tMinimum quality of an overlap, default is 80 (for 80% identity).\n\n" << ::std::endl;
	::std::cerr <<  "-overlaps <Number>\n" << ::std::endl;
	::std::cerr <<  "\tNumber of overlaps that are computed per read, default is 3.\n\n" << ::std::endl;
	::std::cerr <<  "-bandwidth <Number>\n" << ::std::endl;
	::std::cerr <<  "\tSpecifies the bandwidth, default is 8.\n\n" << ::std::endl;
	::std::cerr <<  "-snp [majority | bayesian]\n" << ::std::endl;
	::std::cerr <<  "\tHow to call consensus bases, default is majority.\n\n" << ::std::endl;
	::std::cerr <<  "-outfile <Alignment Filename>\n" << ::std::endl;
	::std::cerr <<  "\tThe name of the output file, default is readAlign.txt.\n\n" << ::std::endl;
	::std::cerr <<  "-output [seqan | afg]\n" << ::std::endl;
	::std::cerr <<  "\tThe output format, default is seqan.\n\n" << ::std::endl;
	::std::cerr <<  "-h\n" << ::std::endl;
	::std::cerr <<  "\tThis help screen.\n\n" << ::std::endl;

	// Insert Sequencing
	::std::cerr <<  "\nInsert Sequencing\n------------\n" << ::std::endl;	
	::std::cerr <<  "-window <Window-Size>\n" << ::std::endl;
	::std::cerr <<  "\tIf this parameter is given all overlaps within a given window are computed (no banded alignment).\n\n" << ::std::endl;

	// Examples
	::std::cerr <<  "\n\n\nExamples\n" << ::std::endl;
	::std::cerr <<  "\nMulti-Read Alignment:\n" << ::std::endl;
	::std::cerr <<  "\t./seqcons -reads reads.fasta\n" << ::std::endl;
	::std::cerr <<  "\t./seqcons -reads reads.fasta -matchlength 15 -quality 80 -outfile align.txt\n" << ::std::endl;
	::std::cerr <<  "\nInsert Sequencing (Reads of length 35):\n" << ::std::endl;
	::std::cerr <<  "\t./seqcons -reads reads.fasta -window 200 -overlaps 100\n" << ::std::endl;
	::std::cerr <<  "\nInsert Sequencing (Reads of length 200):\n" << ::std::endl;
	::std::cerr <<  "\t./seqcons -reads reads.fasta -window 300 -overlaps 100\n" << ::std::endl;
	::std::cerr <<  "\nInsert Sequencing (Reads of length 800):\n" << ::std::endl;
	::std::cerr <<  "\t./seqcons -reads reads.fasta -window 1000 -overlaps 100\n" << ::std::endl;
}

//////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////

template <typename TAlignedRead, typename TContigId>
struct _SimpleLess : 
	public ::std::binary_function<TAlignedRead, TContigId, bool> 
{
	inline bool 
	operator()(const TAlignedRead& a1, const TAlignedRead& a2) const 
	{
		return a1.contigId < a2.contigId;
	}

	bool operator()(const TAlignedRead& a1, TContigId contigId) const {
		return(a1.contigId < contigId);
	}

	bool operator()(TContigId contigId, const TAlignedRead& a2) const {
		return(contigId < a2.contigId);
	} 
};


//////////////////////////////////////////////////////////////////////////////



template <typename TSpec, typename TConfig, typename TPos, typename TGapAnchor, typename TSpecAlign, typename TBeginClr, typename TEndClr>
inline void
getClrRange(FragmentStore<TSpec, TConfig> const& fragStore,
			AlignedReadStoreElement<TPos, TGapAnchor, TSpecAlign> const& alignEl,
			TBeginClr& begClr,		// Out-parameter: left / begin position of the clear range
			TEndClr& endClr)		// Out-parameter: right / end position of the clear range
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename Iterator<String<TGapAnchor>, Standard>::Type TGapIter;
	TSize lenRead = length((value(fragStore.readStore, alignEl.readId)).seq);

	TGapIter itGap = begin(alignEl.gaps, Standard() );
	TGapIter itGapEnd = end(alignEl.gaps, Standard() );
	
	// Any gaps or clipped characters?
	if (itGap == itGapEnd) {
		begClr = 0;
		endClr = lenRead;
		return;
	}

	// Begin clear range
	if (itGap->gapPos == 0) begClr = itGap->seqPos;
	else begClr = 0;

	TSize lenGaps = length(alignEl.gaps);
	if ((value(alignEl.gaps, lenGaps - 1)).seqPos != lenRead) {
		endClr = lenRead;
	} else {
		int diff = 0;
		if (lenGaps > 1) diff = (value(alignEl.gaps, lenGaps - 2)).gapPos - (value(alignEl.gaps, lenGaps - 2)).seqPos;
		int newDiff = (value(alignEl.gaps, lenGaps - 1)).gapPos - (value(alignEl.gaps, lenGaps - 1)).seqPos;
		if (newDiff < diff) {
			endClr = lenRead - (diff - newDiff);
		} else {
			endClr = lenRead;
		}
	}

	// For reverse reads adapt clear ranges
	if (alignEl.beginPos > alignEl.endPos) {
		TBeginClr tmp = begClr;
		begClr = lenRead - endClr;
		endClr = lenRead - tmp;
	}
}


//////////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TStrSpec, typename TPosPair, typename TStringSpec, typename TSpec, typename TConfig, typename TId>
inline void 
getContigReads(StringSet<TValue, Owner<TStrSpec> >& strSet,
			   String<TPosPair, TStringSpec>& startEndPos,
			   FragmentStore<TSpec, TConfig> const& fragStore,
			   TId const contigId)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename TFragmentStore::TReadPos TReadPos;

	// All fragment store element types
	typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadStoreElement;
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;

	// Sort aligned reads according to contig id
	sortAlignedReads(fragStore.alignedReadStore, SortContigId());

	// Retrieve all reads, limit them to the clear range and if required reverse complement them
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
	TAlignIter alignIt = ::std::lower_bound(begin(fragStore.alignedReadStore, Standard()), end(fragStore.alignedReadStore, Standard()), contigId, _SimpleLess<TAlignedElement, TId>());
	TAlignIter alignItEnd = end(fragStore.alignedReadStore, Standard() );
	for(TSize i = 0;alignIt != alignItEnd; goNext(alignIt), ++i) {
		if (alignIt->contigId != contigId) break;
		resize(strSet, i + 1);
		TSize offset = _min(alignIt->beginPos, alignIt->endPos);
		TReadPos begClr = 0;
		TReadPos endClr = 0;
		getClrRange(fragStore, value(alignIt), begClr, endClr);
		value(strSet, i) = infix((value(fragStore.readStore, alignIt->readId)).seq, begClr, endClr);
		TSize lenRead = length(value(strSet, i));
		if (alignIt->beginPos < alignIt->endPos) {
			appendValue(startEndPos, TPosPair(offset, offset + lenRead));
		} else {
			reverseComplementInPlace(value(strSet, i));
			appendValue(startEndPos, TPosPair(offset + lenRead, offset));
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TConfig, typename TPosTriple, typename TAlignMatrix, typename TGappedCons, typename TId>
inline void 
updateContigReads(FragmentStore<TSpec, TConfig>& fragStore,
				  TPosTriple& readBegEndRowPos,
				  TAlignMatrix& alignmentMatrix,
				  TGappedCons& gappedCons,
				  TId const contigId)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename TFragmentStore::TReadPos TReadPos;
	char gapChar = gapValue<char>();
	TSize len = length(alignmentMatrix) / length(readBegEndRowPos);

	// Fragment store typedefs
	typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigStoreElement;
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;
	
	// Update the contig
	TContigStoreElement& contigEl = value(fragStore.contigStore, contigId);
	clear(contigEl.gaps);
	clear(contigEl.seq);

	// Create the gap anchors
	typedef typename Iterator<String<char> >::Type TStringIter;
	TStringIter seqIt = begin(gappedCons);
	TStringIter seqItEnd = end(gappedCons);
	typedef typename TFragmentStore::TReadPos TPos;
	typedef typename TFragmentStore::TContigGapAnchor TContigGapAnchor;
	TPos ungappedPos = 0;
	TPos gappedPos = 0;
	bool gapOpen = false;
	for(;seqIt != seqItEnd; goNext(seqIt), ++gappedPos) {
		if (value(seqIt) == gapChar) gapOpen = true;				
		else {
			if (gapOpen) {
				appendValue(contigEl.gaps, TContigGapAnchor(ungappedPos, gappedPos));
				gapOpen = false;
			}
			Dna5Q letter = value(seqIt);
			assignQualityValue(letter, 'D');
			appendValue(contigEl.seq, letter);
			++ungappedPos;
		}
	}
	if (gapOpen) appendValue(contigEl.gaps, TContigGapAnchor(ungappedPos, gappedPos));


	// Update all aligned reads
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
	TAlignIter alignIt = ::std::lower_bound(begin(fragStore.alignedReadStore, Standard()), end(fragStore.alignedReadStore, Standard()), contigId, _SimpleLess<TAlignedElement, TId>());
	TAlignIter alignItEnd = end(fragStore.alignedReadStore, Standard() );
	for(TSize i = 0;alignIt != alignItEnd; goNext(alignIt), ++i) {
		if (alignIt->contigId != contigId) break;
		TSize lenRead = length((value(fragStore.readStore, alignIt->readId)).seq);
		TReadPos begClr = 0;
		TReadPos endClr = 0;
		getClrRange(fragStore, value(alignIt), begClr, endClr);
		clear(alignIt->gaps);
		ungappedPos = begClr;
		if (alignIt->beginPos > alignIt->endPos) ungappedPos = lenRead - endClr;
		if (ungappedPos != 0) appendValue(alignIt->gaps, TContigGapAnchor(ungappedPos, 0));
		gappedPos = 0;
		gapOpen = false;
		for(TSize column = (readBegEndRowPos[i]).i1; column<(readBegEndRowPos[i]).i2; ++column, ++gappedPos) {
			if (value(alignmentMatrix, (readBegEndRowPos[i]).i3 * len + column) == gapChar) gapOpen = true;				
			else {
				if (gapOpen) {
					appendValue(alignIt->gaps, TContigGapAnchor(ungappedPos, gappedPos));
					gapOpen = false;
				}
				++ungappedPos;
			}
		}
		if (gapOpen) appendValue(alignIt->gaps, TContigGapAnchor(ungappedPos, gappedPos));
		if (alignIt->beginPos < alignIt->endPos) {
			if (endClr != lenRead) {
				appendValue(alignIt->gaps, TContigGapAnchor(lenRead, lenRead + (gappedPos - ungappedPos) - (lenRead - endClr)));
			}
		} else {
			if (begClr != 0) {
				appendValue(alignIt->gaps, TContigGapAnchor(lenRead, lenRead + (gappedPos - ungappedPos) - begClr));
			}
		}

		// Set new begin and end position
		if (alignIt->beginPos < alignIt->endPos) {
			alignIt->beginPos = (readBegEndRowPos[i]).i1;
			alignIt->endPos = (readBegEndRowPos[i]).i2;
		} else {
			alignIt->beginPos = (readBegEndRowPos[i]).i2;
			alignIt->endPos = (readBegEndRowPos[i]).i1;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char *argv[]) {
	
	typedef unsigned int TSize;

	// Version
#ifdef SEQAN_PROFILE
	printVersion();
#endif

	// At least two arguments
	if (argc < 2) {	printHelp(); return 1; }

	// Set default consensus options
	ConsensusOptions consOpt;
#ifdef CELERA_OFFSET
	consOpt.bandwidth = 15;
	consOpt.overlaps = 5;
#else
	consOpt.bandwidth = 8;
	consOpt.overlaps = 3;
#endif

	// Command line parsing
	for(int arg = 1; arg < argc; ++arg) {
		if (argv[arg][0] == '-') {
			if (strcmp(argv[arg], "-reads") == 0) {
				if (arg + 1 < argc) {
					++arg;
					consOpt.readsfile = argv[arg];
				}
			}
			else if (strcmp(argv[arg], "-afg") == 0) {
				if (arg + 1 < argc) {
					++arg;
					consOpt.afgfile = argv[arg];
				}
			}
			else if (strcmp(argv[arg], "-bandwidth") == 0) {
				if (arg + 1 < argc) {
					++arg;
					::std::istringstream istr(argv[arg]);
					istr >> consOpt.bandwidth;
					if (istr.fail()) { printHelp(); exit(1); }
				}
			}
			else if (strcmp(argv[arg], "-overlaps") == 0) {
				if (arg + 1 < argc) {
					++arg;
					::std::istringstream istr(argv[arg]);
					istr >> consOpt.overlaps;
					if (istr.fail()) { printHelp(); exit(1); }
				}
			}
			else if (strcmp(argv[arg], "-matchlength") == 0) {
				if (arg + 1 < argc) {
					++arg;
					::std::istringstream istr(argv[arg]);
					istr >> consOpt.matchlength;
					if (istr.fail()) { printHelp(); exit(1); }
				}
			}
			else if (strcmp(argv[arg], "-quality") == 0) {
				if (arg + 1 < argc) {
					++arg;
					::std::istringstream istr(argv[arg]);
					istr >> consOpt.quality;
					if (istr.fail()) { printHelp(); exit(1); }
				}
			}
			else if (strcmp(argv[arg], "-window") == 0) {
				if (arg + 1 < argc) {
					++arg;
					::std::istringstream istr(argv[arg]);
					istr >> consOpt.window;
					if (istr.fail()) { printHelp(); exit(1); }
				}
			}
			else if (strcmp(argv[arg], "-snp") == 0) {
				if (arg + 1 < argc) {
					++arg;
					if (strcmp(argv[arg], "bayesian") == 0) {
						consOpt.snp = 1;
					} 
					else if (strcmp(argv[arg], "majority") == 0) {
						consOpt.snp = 0;
					}
				}
			}
			else if (strcmp(argv[arg], "-output") == 0) {
				if (arg + 1 < argc) {
					++arg;
					if (strcmp(argv[arg], "afg") == 0) {
						consOpt.output = 1;
					} 
					else if (strcmp(argv[arg], "seqan") == 0) {
						consOpt.output = 0;
					}
				}
			}
			else if (strcmp(argv[arg], "-outfile") == 0) {
				if (arg + 1 < argc) {
					++arg;
					consOpt.outfile = argv[arg];
				}
			}
			else if (strcmp(argv[arg], "-convert") == 0) {
				if (arg + 1 < argc) {
					++arg;
					if (strcmp(argv[arg], "afg") == 0) {
						consOpt.convert = 1;
					} 
					else if (strcmp(argv[arg], "frg") == 0) {
						consOpt.convert = 2;
					}
					else if (strcmp(argv[arg], "cgb") == 0) {
						consOpt.convert = 3;
					}
				}
			}
			else if (strcmp(argv[arg], "-moveToFront") == 0) {
				if (arg + 1 < argc) {
					++arg;
					if (strcmp(argv[arg], "false") == 0) {
						consOpt.moveToFront = false;
					} 
					else if (strcmp(argv[arg], "true") == 0) {
						consOpt.moveToFront = true;
					}
				}
			}
			else if (strcmp(argv[arg], "-source") == 0) {
				if (arg + 1 < argc) {
					++arg;
					consOpt.source = argv[arg];
				}
			}
			else if ((strcmp(argv[arg], "-h") == 0) || (strcmp(argv[arg], "-help") == 0) || (strcmp(argv[arg], "-?") == 0)) {
				printHelp(); 
				return 0; 
			}
		}
	}


	// Fragment store
	typedef FragmentStore<> TFragmentStore;
	TFragmentStore fragStore;

	// Load the reads and layout positions
	TSize numberOfContigs = 0;
	if (!consOpt.readsfile.empty()) {
		// Load simple read file
		std::fstream strmReads;
		strmReads.open(consOpt.readsfile.c_str(), ::std::ios_base::in | ::std::ios_base::binary);
		bool success = _convertSimpleReadFile(strmReads, fragStore, consOpt.readsfile, consOpt.moveToFront);
		strmReads.close();
		if (!success) { printHelp(); exit(1); }
		numberOfContigs = 1;
	} else if (!consOpt.afgfile.empty()) {
		// Load Amos message file
		std::fstream strmReads;
		strmReads.open(consOpt.afgfile.c_str(), ::std::ios_base::in | ::std::ios_base::binary);
		read(strmReads, fragStore, Amos());	
		strmReads.close();
		numberOfContigs = length(fragStore.contigStore);
	} else {
		printHelp();
		exit(1);
	}

	// Just convert the input file
	if (consOpt.convert != 0) {
		if (consOpt.convert == 1) {
			std::fstream strmWrite;
			strmWrite.open(consOpt.outfile.c_str(), ::std::ios_base::out | ::std::ios_base::trunc);
			write(strmWrite, fragStore, Amos());	
			strmWrite.close();
		} else if (consOpt.convert == 2) {
			//ToDo: Celera output 

			//std::fstream strmWrite;
			//strmWrite.open(consOpt.outfile.c_str(), ::std::ios_base::out | ::std::ios_base::trunc);
			//write(strmWrite,readSt,frgSt,libSt,ctgSt,CeleraFrg());	
			//strmWrite.close();
		} else if (consOpt.convert == 3) {
			//ToDo: Celera output 

			//std::fstream strmWrite;
			//strmWrite.open(consOpt.outfile.c_str(), ::std::ios_base::out | ::std::ios_base::trunc);
			//write(strmWrite,readSt,frgSt,libSt,ctgSt,CeleraCgb());	
			//strmWrite.close();
		}
		return 0;
	}

	// Iterate over all contigs
	for(TSize currentContig = 0; currentContig < numberOfContigs; ++currentContig) {
	
// Profiling
#ifdef SEQAN_PROFILE
		SEQAN_PROTIMEUPDATE(__myProfileTime); 
#endif

		// Import all reads of the given contig
		typedef TFragmentStore::TReadSeq TReadSeq;
		typedef Id<TFragmentStore>::Type TId;
		StringSet<TReadSeq, Owner<> > readSet;
		String<Pair<TSize, TSize> > begEndPos;

		// No just get begin and end pointers
		getContigReads(readSet, begEndPos, fragStore, currentContig);
		TSize nseq = length(readSet);
		for(TSize i = 0; i<nseq; ++i) {
			std::cout << value(readSet, i) << std::endl;
			std::cout << value(begEndPos, i).i1 << ',' << value(begEndPos, i).i2 << std::endl;
		}
		if (nseq == 0) continue;
#ifdef SEQAN_PROFILE
		::std::cout << "Import sequences done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << ::std::endl;
#endif

		// Align the reads
		Graph<Alignment<StringSet<TReadSeq, Dependent<> >, void, WithoutEdgeId> > gOut(readSet);
		consensusAlignment(gOut, begEndPos, consOpt);

		// Build the read alignment matrix
		TSize alignDepth;
		String<char> alignmentMatrix;
		String<Triple<unsigned int, unsigned int, unsigned int> > readBegEndRowPos;
		multireadAlignment(gOut, alignmentMatrix, readBegEndRowPos, alignDepth);
		clear(gOut);
#ifdef SEQAN_PROFILE
		std::cout << "Multi-read Alignment done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

		// Call the consensus
		// ToDo: Qualtiy-based consensus calling
		String<unsigned int> coverage;
		String<char> gappedConsensus;
		String<Dna> consensusSequence;
		if (consOpt.snp == 0) consensusCalling(alignmentMatrix, consensusSequence, gappedConsensus, coverage, alignDepth, Majority_Vote() );
		else consensusCalling(alignmentMatrix, consensusSequence, gappedConsensus, coverage, alignDepth, Bayesian() );
#ifdef SEQAN_PROFILE		
		std::cout << "Consensus done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

		// Output of aligned reads
		if (consOpt.output == 0) {
			std::fstream strm;
			if (currentContig == 0) strm.open(consOpt.outfile.c_str(), ::std::ios_base::out | ::std::ios_base::trunc);
			else strm.open(consOpt.outfile.c_str(), ::std::ios_base::out | ::std::ios_base::app);
			write(strm, readSet, alignmentMatrix, begEndPos, readBegEndRowPos, gappedConsensus, coverage, FastaReadFormat());
			strm.close();

			//// Debug code for CA
			//mtRandInit();
			//String<char> fileTmp1 = "tmp1";
			//String<char> fileTmp2 = "tmp2";
			//for(int i = 0; i<10; ++i) {
			//	int file = (mtRand() % 20) + 65;
			//	appendValue(fileTmp1, char(file));
			//	appendValue(fileTmp2, char(file));
			//}
			//std::fstream strm3;
			//strm3.open(toCString(fileTmp2), std::ios_base::out | std::ios_base::trunc);
			//for(int i = 0;i<(int) length(origStrSet); ++i) {
			//	std::stringstream name;
			//	name << value(begEndPos, i).i1 << "," << value(begEndPos, i).i2;
			//	String<char> myTitle = name.str();
			//	write(strm3, origStrSet[i], myTitle, Fasta());			
			//	if (value(begEndPos, i).i1 > value(begEndPos, i).i2) reverseComplementInPlace(origStrSet[i]);
			//}
			//strm3.close();
			//std::fstream strm2;
			//strm2.open(toCString(fileTmp1), std::ios_base::out | std::ios_base::trunc);
			//write(strm2, seqSet, alignmentMatrix, begEndPos, readBegEndRowPos, gappedConsensus, coverage, FastaReadFormat());
			//strm2.close();
		} 
		else if (consOpt.output == 1) {
			updateContigReads(fragStore, readBegEndRowPos, alignmentMatrix, gappedConsensus, currentContig);
		}
#ifdef SEQAN_PROFILE
		std::cout << "Output done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
	}
	
	// Write the AMOS message file
	if (consOpt.output == 1) {
		// Write Amos
		std::fstream strmWrite;
		strmWrite.open(consOpt.outfile.c_str(), ::std::ios_base::out | ::std::ios_base::trunc);
		write(strmWrite, fragStore, Amos());	
		strmWrite.close();
	}

	return 0;
}
