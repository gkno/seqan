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

#include <seqan/basic.h>
#include <seqan/graph_align.h>
#include "../seqan_tcoffee/rna_alphabet.h"

#include <iostream>
#include <fstream>


using namespace seqan;


//////////////////////////////////////////////////////////////////////////////////

void
printVersion() {
	::std::cerr << "*********************************************" << ::std::endl;
	::std::cerr << "* Pairwise alignment                        *" << ::std::endl;
	::std::cerr << "*                                           *" << ::std::endl;
	::std::cerr << "* PairAlign                                 *" << ::std::endl;
	::std::cerr << "* Version: 1.0   (07. July 2009)            *" << ::std::endl;
	::std::cerr << "*********************************************" << ::std::endl;
	::std::cerr << ::std::endl;
}

//////////////////////////////////////////////////////////////////////////////////

void 
printHelp() {
	::std::cerr <<  "Usage: pair_align -seq <FASTA Sequence File> [Options]\n" << ::std::endl;
	// Main options
	::std::cerr << "\nMain Options\n------------" << ::std::endl;
	::std::cerr <<  "-seq <FASTA Sequence File>" << ::std::endl;
	::std::cerr <<  "/* File with two sequences in FASTA format. */" << ::std::endl;
	::std::cerr <<  "-alphabet [protein | dna | rna]" << ::std::endl;
	::std::cerr <<  "/* Sequence alphabet, default is protein. */" << ::std::endl;
	::std::cerr <<  "-h" << ::std::endl;
	::std::cerr <<  "/* This help screen. */" << ::std::endl;
	::std::cerr <<  "-method [nw | gotoh | sw]" << ::std::endl;
	::std::cerr <<  "/* Alignment method, default is gotoh. */" << ::std::endl;
	::std::cerr <<  "-outfile <Alignment Filename>" << ::std::endl;
	::std::cerr <<  "/* Name of the output file, default is out.fasta. */" << ::std::endl;
	::std::cerr <<  "-output [fasta | msf]" << ::std::endl;
	::std::cerr <<  "/* Output format, default is fasta. */" << ::std::endl;
	// Scoring
	::std::cerr <<  "\nScoring Options\n------------" << ::std::endl;
	::std::cerr <<  "-gop <Number>" << ::std::endl;
	::std::cerr <<  "/* Gap open penalty, default is -11. */" << ::std::endl;
	::std::cerr <<  "-gex <Number>" << ::std::endl;
	::std::cerr <<  "/* Gap extension penalty, default is -1. */" << ::std::endl;
	::std::cerr <<  "-matrix <Score-Matrix File>" << ::std::endl;
	::std::cerr <<  "/* Specifies a score-matrix file, default is Blosum62. */" << ::std::endl;
	::std::cerr <<  "-msc <Number>" << ::std::endl;
	::std::cerr <<  "/* Match score for Dna / Rna alphabet, default is 5. */" << ::std::endl;
	::std::cerr <<  "-mmsc <Number>" << ::std::endl;
	::std::cerr <<  "/* Mismatch penalty for Dna / Rna alphabet, default is -4. */" << ::std::endl;
	// Banded alignment
	::std::cerr <<  "\nBanded Alignment Options (only nw and gotoh)\n------------" << ::std::endl;
	::std::cerr <<  "-low <Number>" << ::std::endl;
	::std::cerr <<  "/* Lower diagonal (has to be greater than the negative length of the second string). */" << ::std::endl;
	::std::cerr <<  "-high <Number>" << ::std::endl;
	::std::cerr <<  "/* Upper diagonal (has to be smaller than the length of the first string). */" << ::std::endl;
	// Alignment configuration
	::std::cerr <<  "\nAlignment configuration\n------------" << ::std::endl;
	::std::cerr <<  "-config [true | false], [true | false], [true | false], [true | false]" << ::std::endl;
	::std::cerr <<  "/* Firt value indicates whether the first row of the dynamic programming matrix is initialized with 0's. Second row indicates whether the left side is initialized with 0's. Third value indicates whether the maximum is searched in the whole last column. Fourth value indicates whether the maximum is searched in the whole last row. Default is false, false, false, false. */" << ::std::endl;
	

	::std::cerr <<  "\n\n\nExamples\n------------" << ::std::endl;
	::std::cerr <<  "\nProtein Alignment:" << ::std::endl;
	::std::cerr <<  "\t./pair_align -seq seq.fasta" << ::std::endl;
	::std::cerr <<  "\t./pair_align -seq seq.fasta -method nw" << ::std::endl;
	::std::cerr <<  "\nDna Alignment:" << ::std::endl;
	::std::cerr <<  "\t./pair_align -seq seq.fasta -alphabet dna" << ::std::endl;
	::std::cerr <<  "\nLocal Alignment:" << ::std::endl;
	::std::cerr <<  "\t./pair_align -seq seq.fasta -method sw" << ::std::endl;

}


//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TAlignConfig, typename TScore, typename TSeqFile, typename TMethod, typename TDiag, typename TOutputFormat, typename TOutfile>
inline void
pairwise_align(TScore const& sc,
			   TSeqFile& seqfile,
			   TMethod method,
			   TDiag low,
			   TDiag high,
			   bool banded,
			   TOutputFormat outputFormat,
			   TOutfile& outfile) 
{
	// Load the 2 sequences
	typedef String<TAlphabet> TSequence;
	StringSet<TSequence, Owner<> > sequenceSet;
	StringSet<String<char> > sequenceNames;
	_loadSequences(seqfile, sequenceSet, sequenceNames);

	// Fix low and high diagonal.
	low = _max(low, -1 * (int) length(sequenceSet[1]));
	high = _min(high, (int) length(sequenceSet[0]));

	// Align the sequences
	Graph<Alignment<StringSet<TSequence, Dependent<> >, void, WithoutEdgeId> > gAlign(sequenceSet);
	
	int aliScore = 0;
	// Banded alignment?
	if (!banded) {
		if (method == 0) aliScore = globalAlignment(gAlign, sc, TAlignConfig(), NeedlemanWunsch());
		else if (method == 1) aliScore = globalAlignment(gAlign, sc, TAlignConfig(), Gotoh());
		else if (method == 2) aliScore = localAlignment(gAlign, sc, SmithWaterman());
	} else {
		if (method == 0) aliScore = globalAlignment(gAlign, sc, TAlignConfig(), low, high, BandedNeedlemanWunsch());
		else if (method == 1) aliScore = globalAlignment(gAlign, sc, TAlignConfig(), low, high, BandedGotoh());
	}

	
	//// Debug options
	//std::cout << "Seq1: " << sequenceSet[0] << std::endl;
	//std::cout << "Seq2: " << sequenceSet[1] << std::endl;
	//std::cout << "Scoring parameters:" << std::endl;
	//std::cout << "*Gap opening: " << scoreGapOpen(sc) << std::endl;
	//std::cout << "*Gap extension: " << scoreGapExtend(sc) << std::endl;
	//std::cout << "*Scoring matrix: " << std::endl;
	//unsigned int alphSize = ValueSize<TAlphabet>::VALUE;
	//std::cout << "   ";
	//for(unsigned int col = 0; col<alphSize; ++col) std::cout << TAlphabet(col) << ',';
	//std::cout << std::endl;
	//for(unsigned int row = 0; row<alphSize; ++row) {
	//	for(unsigned int col = 0; col<alphSize; ++col) {
	//		if (col == 0) std::cout << TAlphabet(row) << ": ";
	//		std::cout << score(sc, TAlphabet(row), TAlphabet(col));
	//		if (col < alphSize - 1) std::cout << ',';
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;
	//if (method == 0) std::cout << "Method: NeedlemanWunsch" << std::endl;
	//if (method == 1) std::cout << "Method: Gotoh" << std::endl;
	//if (method == 2) std::cout << "Method: SmithWaterman" << std::endl;
	//std::cout << "Output file: " << outfile << std::endl;
	//if (outputFormat == 0) std::cout << "Output format: Fasta" << std::endl;
	//if (outputFormat == 1) std::cout << "Output format: Msf" << std::endl;
	//if (banded) std::cout << "Alignment bands: " << low << ',' << high << std::endl;
	//if ((method == 0) || (method == 1)) std::cout << "AlignConfig<" << __myInitTop(TAlignConfig()) << ',' << __myInitLeft(TAlignConfig()) << ',' << __myInitRight(TAlignConfig()) << ',' << __myInitBottom(TAlignConfig()) << '>' << std::endl;
	//std::cout << "Alignment: " << std::endl;
	//std::cout << gAlign << std::endl;
	
	// Alignment output
	std::cout << "Alignment score: " << aliScore << std::endl;
	if (outputFormat == 0) {
		std::fstream strm;
		strm.open(outfile.c_str(), ::std::ios_base::out | ::std::ios_base::trunc);
		write(strm, gAlign, sequenceNames, FastaFormat());
		strm.close();
	} else if (outputFormat == 1) {
		std::fstream strm;
		strm.open(outfile.c_str(), ::std::ios_base::out | ::std::ios_base::trunc);
		write(strm, gAlign, sequenceNames, MsfFormat());
		strm.close();
	}
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TScore, typename TSc>
inline void
_setMatchScore(TScore&, TSc) {
	// No operation
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TScore, typename TSc>
inline void
_setMismatchScore(TScore&, TSc) {
	// No operation
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TSc>
inline void
_setMatchScore(Score<int, Simple>& sc, TSc msc) {
	sc.data_match = msc;
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TSc>
inline void
_setMismatchScore(Score<int, Simple>& sc, TSc mmsc) {
	sc.data_mismatch = mmsc;
}


//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TScore>
inline void
_initAlignParams(int argc, const char *argv[], TScore& sc) {
	// Set default options
	sc.data_gap_open = -11;
	sc.data_gap_extend = -1;
	_setMatchScore(sc, 5);
	_setMismatchScore(sc, -4);
	::std::string outfile = "out.fasta";
	::std::string seqfile;
	String<bool> bools;
	int low = 0;
	int high = 0;
	unsigned int method = 1;
	unsigned int outputFormat = 0;
	bool banded = false;
	
	// Command line parsing
	for(int arg = 1; arg < argc; ++arg) {
		if (argv[arg][0] == '-') {
			if (strcmp(argv[arg], "-seq") == 0) {
				if (arg + 1 < argc) {
					++arg;
					seqfile = argv[arg];
				}
			}
			else if (strcmp(argv[arg], "-outfile") == 0) {
				if (arg + 1 < argc) {
					++arg;
					outfile = argv[arg];
				}
			}
			else if (strcmp(argv[arg], "-method") == 0) {
				if (arg + 1 < argc) {
					++arg;
					::std::string meth = argv[arg];
					if (meth == "nw") method = 0;
					else if (meth == "gotoh") method = 1;
					else if (meth == "sw") method = 2;
				}
			}
			else if (strcmp(argv[arg], "-output") == 0) {
				if (arg + 1 < argc) {
					++arg;
					::std::string output = argv[arg];
					if (output == "fasta") outputFormat = 0;
					else if (output == "msf") outputFormat = 1;
				}
			}
			else if (strcmp(argv[arg], "-gop") == 0) {
				if (arg + 1 < argc) {
					++arg;
					::std::istringstream istr(argv[arg]);
					istr >> sc.data_gap_open;
					if (istr.fail()) { printHelp(); exit(1); }
				}
			}
			else if (strcmp(argv[arg], "-gex") == 0) {
				if (arg + 1 < argc) {
					++arg;
					::std::istringstream istr(argv[arg]);
					istr >> sc.data_gap_extend;
					if (istr.fail()) { printHelp(); exit(1); }
				}
			}
			else if (strcmp(argv[arg], "-msc") == 0) {
				if (arg + 1 < argc) {
					++arg;
					int msc;
					::std::istringstream istr(argv[arg]);
					istr >> msc;
					if (istr.fail()) { printHelp(); exit(1); }
					_setMatchScore(sc, msc);
				}
			}
			else if (strcmp(argv[arg], "-mmsc") == 0) {
				if (arg + 1 < argc) {
					++arg;
					int mmsc;
					::std::istringstream istr(argv[arg]);
					istr >> mmsc;
					if (istr.fail()) { printHelp(); exit(1); }
					_setMismatchScore(sc, mmsc);
				}
			}
			else if (strcmp(argv[arg], "-low") == 0) {
				if (arg + 1 < argc) {
					++arg;
					::std::istringstream istr(argv[arg]);
					istr >> low;
					if (istr.fail()) { printHelp(); exit(1); }
					banded = true;
				}
			}
			else if (strcmp(argv[arg], "-high") == 0) {
				if (arg + 1 < argc) {
					++arg;
					::std::istringstream istr(argv[arg]);
					istr >> high;
					if (istr.fail()) { printHelp(); exit(1); }
					banded = true;
				}
			}
			else if (strcmp(argv[arg], "-config") == 0) {
				while ((arg+1 < argc) && (argv[arg+1][0] != '-')) {
					++arg;
					::std::string boolval = argv[arg];
					::std::string::size_type pos = boolval.find(',');
					if (pos != ::std::string::npos) boolval.erase(pos);
					if (boolval == "true") appendValue(bools, true);
					else appendValue(bools, false);
				}
			}
		}
	}
	if (seqfile.empty()) { printHelp(); exit(1); }
	if (low > high) banded = false;
	if (length(bools) != 4) pairwise_align<TAlphabet, AlignConfig<false, false, false, false> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
	else {
		if ((bools[0]) && (bools[1]) && (bools[2]) && (bools[3])) pairwise_align<TAlphabet, AlignConfig<true, true, true, true> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if ((bools[0]) && (bools[1]) && (bools[2]) && (!bools[3])) pairwise_align<TAlphabet, AlignConfig<true, true, true, false> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if ((bools[0]) && (bools[1]) && (!bools[2]) && (bools[3])) pairwise_align<TAlphabet, AlignConfig<true, true, false, true> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if ((bools[0]) && (bools[1]) && (!bools[2]) && (!bools[3])) pairwise_align<TAlphabet, AlignConfig<true, true, false, false> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if ((bools[0]) && (!bools[1]) && (bools[2]) && (bools[3])) pairwise_align<TAlphabet, AlignConfig<true, false, true, true> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if ((bools[0]) && (!bools[1]) && (bools[2]) && (!bools[3])) pairwise_align<TAlphabet, AlignConfig<true, false, true, false> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if ((bools[0]) && (!bools[1]) && (!bools[2]) && (bools[3])) pairwise_align<TAlphabet, AlignConfig<true, false, false, true> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if ((bools[0]) && (!bools[1]) && (!bools[2]) && (!bools[3])) pairwise_align<TAlphabet, AlignConfig<true, false, false, false> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if ((!bools[0]) && (bools[1]) && (bools[2]) && (bools[3])) pairwise_align<TAlphabet, AlignConfig<false, true, true, true> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if ((!bools[0]) && (bools[1]) && (bools[2]) && (!bools[3])) pairwise_align<TAlphabet, AlignConfig<false, true, true, false> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if ((!bools[0]) && (bools[1]) && (!bools[2]) && (bools[3])) pairwise_align<TAlphabet, AlignConfig<false, true, false, true> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if ((!bools[0]) && (bools[1]) && (!bools[2]) && (!bools[3])) pairwise_align<TAlphabet, AlignConfig<false, true, false, false> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if ((!bools[0]) && (!bools[1]) && (bools[2]) && (bools[3])) pairwise_align<TAlphabet, AlignConfig<false, false, true, true> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if ((!bools[0]) && (!bools[1]) && (bools[2]) && (!bools[3])) pairwise_align<TAlphabet, AlignConfig<false, false, true, false> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if ((!bools[0]) && (!bools[1]) && (!bools[2]) && (bools[3])) pairwise_align<TAlphabet, AlignConfig<false, false, false, true> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if ((!bools[0]) && (!bools[1]) && (!bools[2]) && (!bools[3])) pairwise_align<TAlphabet, AlignConfig<false, false, false, false> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
	}
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TMatrixFile>
inline void
_initScoreMatrix(int argc, const char *argv[], TMatrixFile& matrix, Dna5 const) {
	if (!empty(matrix)) {
		Score<int, ScoreMatrix<> > sc;
		loadScoreMatrix(sc, matrix);
		_initAlignParams<Dna5>(argc, argv, sc);
	} else {
		Score<int> sc;
		_initAlignParams<Dna5>(argc, argv, sc);
	}
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TMatrixFile>
inline void
_initScoreMatrix(int argc, const char *argv[], TMatrixFile& matrix, Rna5 const) {
	if (!empty(matrix)) {
		Score<int, ScoreMatrix<> > sc;
		loadScoreMatrix(sc, matrix);
		_initAlignParams<Rna5>(argc, argv, sc);
	} else {
		Score<int> sc;
		_initAlignParams<Rna5>(argc, argv, sc);
	}
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TMatrixFile>
inline void
_initScoreMatrix(int argc, const char *argv[], TMatrixFile& matrix, AminoAcid const) {
	if (!empty(matrix)) {
		Score<int, ScoreMatrix<> > sc;
		loadScoreMatrix(sc, matrix);
		_initAlignParams<AminoAcid>(argc, argv, sc);
	} else {
		Blosum62 sc;
		_initAlignParams<AminoAcid>(argc, argv, sc);
	}
}

//////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char *argv[]) {

	// Version
	printVersion();

	// At least two arguments
	if (argc < 2) {	printHelp(); return 1; }

	// Basic command line options
	String<char> alphabet = "protein";
	String<char> matrix;
	for(int arg = 1; arg < argc; ++arg) {
		if (argv[arg][0] == '-') {
			// Dna, Rna or AminoAcid alignment
			if (strcmp(argv[arg], "-alphabet") == 0) {
				if (arg + 1 < argc) {
					++arg;
					alphabet = argv[arg];
				}
			}
			else if (strcmp(argv[arg], "-matrix") == 0) {
				if (arg + 1 < argc) {
					++arg;
					matrix = argv[arg];
				}
			}
			else if ((strcmp(argv[arg], "-h") == 0) || (strcmp(argv[arg], "-help") == 0) || (strcmp(argv[arg], "-?") == 0)) {
				printHelp(); 
				return 0; 
			}
		}
	}

	// Initialize scoring matrices
	if (alphabet == "dna") _initScoreMatrix(argc, argv, matrix, Dna5());
	else if (alphabet == "rna") _initScoreMatrix(argc, argv, matrix, Rna5());
	else _initScoreMatrix(argc, argv, matrix, AminoAcid());

	return 0;
}
