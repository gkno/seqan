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

template <typename TName, typename TKey>
inline void
__findThis(TName const& searchText,
		   TKey const& keyText,
		   TName& result) 
{
	typedef typename Iterator<TName>::Type TNameIter;
	TName key = keyText;
	TNameIter nameIt = begin(searchText);
	TNameIter nameItEnd = end(searchText);
	for(;nameIt != nameItEnd;++nameIt) {
		TNameIter keyIt = begin(key);
		TNameIter keyItEnd = end(key);
		TNameIter localNameIt = nameIt;
		while((localNameIt != nameItEnd) && (keyIt != keyItEnd) && (value(localNameIt) == value(keyIt))) {
			++keyIt;
			++localNameIt;
		}
		if (keyIt == keyItEnd) {
			while ((value(localNameIt) != ']') && (value(localNameIt) != ',')) {
				appendValue(result, value(localNameIt));
				++localNameIt;
			}
			break;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TOptions, typename TAlphabet, typename TSpec, typename TFragmentStore, typename TLibraryStore, typename TContigStore>
inline bool
convertSimpleReadFile(ReadStore<TAlphabet, TSpec>& readSt,
					  TFragmentStore& frgSt,
					  TLibraryStore& libSt,
					  TContigStore& ctgSt,
					  TOptions& consOpt,
					  bool moveToFront)
{
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef String<char> TName;



	// Convert reads
	typedef String<TAlphabet> TSequence;
	StringSet<String<TAlphabet>, Owner<> > origStrSet;
	String<TName> names;
	String<char>  filePath;
	filePath = consOpt.readsfile;
	TSize count = _loadSequences(filePath, origStrSet, names);
	if (count == 0) return false;
	else {
		for(TSize i = 0; i<length(origStrSet); ++i) {
			// Get the begin and end position
			TSize begRead;
			if (readSt.data_pos_count == 0) begRead = 0;
			else begRead = (value(readSt.data_begin_end, readSt.data_pos_count - 1)).i2;
			TSize endRead = begRead + length(value(origStrSet, i));
			appendValue(readSt.data_begin_end, Pair<TSize,TSize>(begRead, endRead));

			// Save the read
			typedef typename Iterator<String<TAlphabet> >::Type TSeqIter;
			TSeqIter itSeq = begin(value(origStrSet, i));
			TSeqIter itSeqEnd = end(value(origStrSet, i));
			for(;itSeq != itSeqEnd; ++itSeq) {
				appendValue(readSt.data_reads, *itSeq);
				appendValue(readSt.data_qualities, char(48+60)); // Dummy quality value
			}
			
			// Get the fragment id
			TSize fragId = i;
			TName valKey;
			__findThis(value(names, i), "fragId=", valKey);
			if (!empty(valKey)) fragId = _stringToNumber<TSize>(valKey);
			appendValue(readSt.data_frg_id, fragId);
			
			// Build eid string
			std::stringstream input;
			clear(valKey);
			__findThis(value(names, i), "id=", valKey);
			if (empty(valKey)) input << "R" << i;
			else input << valKey;

			// Possibly append repeatId
			clear(valKey);
			__findThis(value(names, i), "repeatId=", valKey);
			if (!empty(valKey)) input << "-" << valKey;
			String<char> tmp(input.str().c_str());
			appendValue(readSt.data_names, tmp);

			// Append a clear range
			appendValue(readSt.data_clr, Pair<TSize, TSize>(0, length(value(origStrSet, i))));
			++readSt.data_pos_count;
		}
	}




	// Write minimal contig

	// Create a contig name
	String<char> tmp = "C0";
	appendValue(ctgSt.data_names, tmp);

	// Create a string of gapped reads
	appendValue(ctgSt.data_reads, String<GappedRead<> >());
	TSize lastLetter = 0;
	for(TSize i = 0; i<length(names); ++i) {
		// Add the original read id
		GappedRead<> gapRead;
		gapRead.data_source = i;

		// Add the offset and clear range
		typedef typename Iterator<TName>::Type TNameIter;
		TNameIter nameIt = begin(value(names, i));
		TNameIter nameItEnd = end(value(names, i));
		bool before = true;
		String<char> inf1; String<char> inf2;
		for(;nameIt != nameItEnd; goNext(nameIt)) {
			if (value(nameIt) == ',') before = false;
			else {
				if (before) appendValue(inf1, value(nameIt));
				else appendValue(inf2, value(nameIt));
			}
		}
		TSize posI = _stringToNumber<TSize>(inf1);
		TSize posJ = _stringToNumber<TSize>(inf2);
		if (posI < posJ) {
			if (posJ > lastLetter) lastLetter = posJ;
			gapRead.data_offset = posI;
			gapRead.data_clr = Pair<TSize, TSize>(0, length(value(origStrSet, i)));
		} else {
			if (posI > lastLetter) lastLetter = posI;
			gapRead.data_offset = posJ;
			gapRead.data_clr = Pair<TSize, TSize>(length(value(origStrSet, i)), 0);
		}
		appendValue(value(ctgSt.data_reads, ctgSt.data_pos_count), gapRead);
	}

	// Create a contig sequence	
	StringSet<String<char>, Owner<> > seqSet;
	String<TName> seqNames;
	clear(filePath);
	if (!consOpt.source.empty()) {
		filePath = consOpt.source;
	} else {
		filePath = consOpt.readsfile;
		appendValue(filePath, 'S');
	}
	count = _loadSequences(filePath, seqSet, seqNames);
	if (count != 0) {
		for(TSize l = 0;l < length(value(seqSet,0)); ++l) {
			appendValue(ctgSt.data_contig, value(value(seqSet, 0), l));
			appendValue(ctgSt.data_quality, char(48+60)); // Dummy quality value
		}
		appendValue(ctgSt.data_begin_end, Pair<TSize,TSize>(0, length(value(seqSet,0))));
	} else {
		for(TSize l = 0;l < lastLetter; ++l) {
			appendValue(ctgSt.data_contig, '-'); // Dummy letter
			appendValue(ctgSt.data_quality, char(48+60)); // Dummy quality value
		}
		appendValue(ctgSt.data_begin_end, Pair<TSize, TSize>(0, lastLetter));
	}
	++ctgSt.data_pos_count; // Added just one contig



	// Convert Libraries
	clear(filePath);
	filePath = consOpt.readsfile;
	appendValue(filePath, 'L');
	StringSet<String<char>, Owner<> > libSet;
	String<TName> libIds;
	count = _loadSequences(filePath, libSet, libIds);
	if (count != 0) {
		for(TSize i = 0; i<length(libSet); ++i) {
			// Create the library name
			std::stringstream input;
			input << "L" << i;
			String<char> tmp(input.str().c_str());
			appendValue(libSt.data_names, tmp);

			// Get the mean and std
			typedef typename Iterator<String<char> >::Type TStrIter;
			TStrIter nameIt = begin(value(libSet, i));
			TStrIter nameItEnd = end(value(libSet, i));
			bool before = true;
			String<char> inf1; String<char> inf2;
			for(;nameIt != nameItEnd; goNext(nameIt)) {
				if (value(nameIt) == ',') before = false;
				else {
					if (before) appendValue(inf1, value(nameIt));
					else appendValue(inf2, value(nameIt));
				}
			}
			TSize mean = _stringToNumber<TSize>(inf1);
			TSize std = _stringToNumber<TSize>(inf2);
			appendValue(libSt.data_mean, mean);
			appendValue(libSt.data_std, std);
			++libSt.data_pos_count;
		}
	} else {
		appendValue(libSt.data_names, "L0");
		appendValue(libSt.data_mean, 0);
		appendValue(libSt.data_std, 0);
		++libSt.data_pos_count;
	}



	// Convert Fragments
	filePath = consOpt.readsfile;
	appendValue(filePath, 'F');
	StringSet<String<char>, Owner<> > fragSet;
	String<TName> fragIds;
	count = _loadSequences(filePath, fragSet, fragIds);
	if (count != 0) {
		for(TSize i = 0;i<length(fragSet); ++i) {
			// Create the fragment name
			std::stringstream input;
			input << "F" << i;
			String<char> tmp(input.str().c_str());
			appendValue(frgSt.data_names, tmp);
			
			// Find the library id
			TSize libId = i;
			TName valKey;
			__findThis(value(fragIds, i), "libId=", valKey);
			if (!empty(valKey)) libId = _stringToNumber<TSize>(valKey);
			appendValue(frgSt.data_lib_id, libId);
			
			// Get rds1 and rds2
			typedef typename Iterator<String<char> >::Type TStrIter;
			TStrIter nameIt = begin(value(fragSet, i));
			TStrIter nameItEnd = end(value(fragSet, i));
			bool before = true;
			String<char> inf1; String<char> inf2;
			for(;nameIt != nameItEnd; goNext(nameIt)) {
				if (value(nameIt) == ',') before = false;
				else {
					if (before) appendValue(inf1, value(nameIt));
					else appendValue(inf2, value(nameIt));
				}
			}
			TSize rds1 = _stringToNumber<TSize>(inf1);
			TSize rds2 = _stringToNumber<TSize>(inf2);
			if (rds1 < rds2) appendValue(frgSt.data_rds, Pair<TSize, TSize>(rds1, rds2));
			else appendValue(frgSt.data_rds, Pair<TSize, TSize>(rds2, rds1));
			++frgSt.data_pos_count;
		}
	} else {
		for(TSize i = 0; i<readSt.data_pos_count; ++i) {
			std::stringstream input;
			input << "F" << i;
			String<char> tmp(input.str().c_str());
			appendValue(frgSt.data_names, tmp);
			appendValue(frgSt.data_lib_id, 0);
			appendValue(frgSt.data_rds, Pair<TSize, TSize>(0, 0));
			++frgSt.data_pos_count;
		}
	}



	// Find smallest offset and reset all contig offsets
	if (moveToFront) {
		for(TSize i = 0; i<length(ctgSt); ++i) {
			TSize smallestOffset = length(ctgSt.data_contig);
			String<GappedRead<> >& gappedReads = value(ctgSt.data_reads, i);
			for(TSize k = 0; k<length(gappedReads); ++k) {
				if ((value(gappedReads, k)).data_offset < smallestOffset) smallestOffset = (value(gappedReads, k)).data_offset;
			}
			if (smallestOffset > 0) {
				for(TSize k = 0; k<length(gappedReads); ++k) {
					(value(gappedReads, k)).data_offset -= smallestOffset;
				}
			}
		}
	}
	return true;
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

	// Load simulation file, afg file or celera consensus file
	ReadStore<> readSt;
	FrgStore<> frgSt;
	LibStore<> libSt;
	CtgStore<> ctgSt;

	TSize numberOfContigs = 0;
	if (!consOpt.readsfile.empty()) {
		bool success = convertSimpleReadFile(readSt, frgSt, libSt, ctgSt, consOpt, consOpt.moveToFront);
		if (!success) { printHelp(); exit(1); }
		numberOfContigs = 1;
	} else if (!consOpt.afgfile.empty()) {
		// Read Amos
		std::fstream strmReads;
		strmReads.open(consOpt.afgfile.c_str(), ::std::ios_base::in | ::std::ios_base::binary);
		read(strmReads, readSt, frgSt, libSt, ctgSt, Amos());	
		strmReads.close();
		numberOfContigs = length(ctgSt);
	} else {
		printHelp();
		exit(1);
	}

	// Just convert the input file
	if (consOpt.convert != 0) {
		if (consOpt.convert == 1) {
			std::fstream strmWrite;
			strmWrite.open(consOpt.outfile.c_str(), ::std::ios_base::out | ::std::ios_base::trunc);
			write(strmWrite,readSt,frgSt,libSt,ctgSt,Amos());	
			strmWrite.close();
		} else if (consOpt.convert == 2) {
			std::fstream strmWrite;
			strmWrite.open(consOpt.outfile.c_str(), ::std::ios_base::out | ::std::ios_base::trunc);
			write(strmWrite,readSt,frgSt,libSt,ctgSt,CeleraFrg());	
			strmWrite.close();
		} else if (consOpt.convert == 3) {
			std::fstream strmWrite;
			strmWrite.open(consOpt.outfile.c_str(), ::std::ios_base::out | ::std::ios_base::trunc);
			write(strmWrite,readSt,frgSt,libSt,ctgSt,CeleraCgb());	
			strmWrite.close();
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
		typedef String<Dna> TSequence;
		typedef Id<TSequence>::Type TId;
		StringSet<TSequence, Owner<> > origStrSet;
		String<Pair<TSize, TSize> > begEndPos;
		loadReadsClr(readSt, ctgSt, currentContig, origStrSet, begEndPos);
		TSize nseq = length(origStrSet);
		if (nseq == 0) continue;
#ifdef SEQAN_PROFILE
		::std::cout << "Import sequences done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << ::std::endl;
#endif

		// Align the reads
		Graph<Alignment<StringSet<TSequence, Dependent<> >, void, WithoutEdgeId> > gOut(origStrSet);
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
			write(strm, origStrSet, alignmentMatrix, begEndPos, readBegEndRowPos, gappedConsensus, coverage, FastaReadFormat());
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
			TSize len = length(gappedConsensus);
			TSize begContig = (value(ctgSt.data_begin_end, ctgSt.data_pos_count - 1)).i2;
			TSize endContig = begContig + length(gappedConsensus);
			for(TSize i = 0; i<len; ++i) appendValue(ctgSt.data_contig, gappedConsensus[i]);
			for(TSize i = 0; i<len; ++i) appendValue(ctgSt.data_quality, 'D');
			value(ctgSt.data_begin_end, currentContig) = Pair<TSize, TSize>(begContig, endContig);
			for(TSize i = 0; i < length(readBegEndRowPos); ++i) {
				String<TSize> gaps;
				TSize letterCount = 0;
				TSize gapCount = 0;
				for(TSize column = (readBegEndRowPos[i]).i1; column<(readBegEndRowPos[i]).i2; ++column) {
					if (value(alignmentMatrix, (readBegEndRowPos[i]).i3 * len + column) == '-') {
						++gapCount;
						appendValue(gaps, letterCount);
					} else ++letterCount;
				}
				value(value(ctgSt.data_reads, currentContig), i).data_gap = gaps;
				value(value(ctgSt.data_reads, currentContig), i).data_offset = (readBegEndRowPos[i]).i1;
			}
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
		write(strmWrite,readSt,frgSt,libSt,ctgSt,Amos());	
		strmWrite.close();
	}

	return 0;
}
