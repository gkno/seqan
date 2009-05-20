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

#include <seqan/graph_msa.h>
#include "rna_alphabet.h"

using namespace seqan;


//////////////////////////////////////////////////////////////////////////////////

void
printVersion() {
	::std::cerr << "*********************************************" << ::std::endl;
	::std::cerr << "* Segment-based multiple sequence alignment *" << ::std::endl;
	::std::cerr << "*                                           *" << ::std::endl;
	::std::cerr << "* SeqAn::T-Coffee                           *" << ::std::endl;
	::std::cerr << "* Version: 1.101 (20. May 2009)             *" << ::std::endl;
	::std::cerr << "*********************************************" << ::std::endl;
	::std::cerr << ::std::endl;
}

//////////////////////////////////////////////////////////////////////////////////

void 
printHelp() {
	::std::cerr <<  "Usage: seqan_tcoffee -seq <FASTA Sequence File> [Options]\n" << ::std::endl;
	::std::cerr << "\nOptions\n" << ::std::endl;
	// Main options
	::std::cerr << "\nMain Options\n------------\n" << ::std::endl;
	::std::cerr <<  "-seq <FASTA Sequence File>\n" << ::std::endl;
	::std::cerr <<  "\tFile with multiple sequences in FASTA format.\n\n" << ::std::endl;
	::std::cerr <<  "-alphabet [protein | dna | rna]\n" << ::std::endl;
	::std::cerr <<  "\tSequence alphabet, default is protein.\n\n" << ::std::endl;
	::std::cerr <<  "-h\n" << ::std::endl;
	::std::cerr <<  "\tThis help screen.\n\n" << ::std::endl;
	::std::cerr <<  "-outfile <Alignment Filename>\n" << ::std::endl;
	::std::cerr <<  "\tName of the output file, default is out.fasta.\n\n" << ::std::endl;
	::std::cerr <<  "-output [fasta | msf]\n" << ::std::endl;
	::std::cerr <<  "\tOutput format, default is fasta.\n\n" << ::std::endl;
	// Segment match generation options
	::std::cerr <<  "\nSegment-Match Generation Options\n------------\n" << ::std::endl;
	::std::cerr <<  "-method [global], [local], [overlap], [lcs]\n" << ::std::endl;
	::std::cerr <<  "\tMethods to generate segment matches, default is global, local.\n\n" << ::std::endl;
	::std::cerr <<  "-blast <BLAST matches>, <BLAST matches>, ...\n" << ::std::endl;
	::std::cerr <<  "\tFiles with gapless segment matches in BLAST tabular format (-m 8 -g F).\n\n" << ::std::endl;
	::std::cerr <<  "-mummer <MUMmer matches>, <MUMmer matches>, ...\n" << ::std::endl;
	::std::cerr <<  "\tFiles with gapless segment matches in MUMmer format.\n\n" << ::std::endl;
	::std::cerr <<  "-aln <FASTA Alignment File>, <FASTA Alignment File>, ...\n" << ::std::endl;
	::std::cerr <<  "\tSegment-matches are extracted from given alignments (Meta-Alignment).\n\n" << ::std::endl;
	::std::cerr <<  "-lib <T-Coffee library>, <T-Coffee library>, ...\n" << ::std::endl;
	::std::cerr <<  "\tSegment-matches are extracted from T-Coffee libraries.\n\n" << ::std::endl;
	// Scoring
	::std::cerr <<  "\nScoring Options\n------------\n" << ::std::endl;
	::std::cerr <<  "-gop <Number>\n" << ::std::endl;
	::std::cerr <<  "\tGap open penalty, default is -11.\n\n" << ::std::endl;
	::std::cerr <<  "-gex <Number>\n" << ::std::endl;
	::std::cerr <<  "\tGap extension penalty, default is -1.\n\n" << ::std::endl;
	::std::cerr <<  "-matrix <Score-Matrix File>\n" << ::std::endl;
	::std::cerr <<  "\tSpecifies a score-matrix file, default is Blosum62.\n\n" << ::std::endl;
	::std::cerr <<  "-msc <Number>\n" << ::std::endl;
	::std::cerr <<  "\tMatch score for Dna / Rna alphabet, default is 5.\n\n" << ::std::endl;
	::std::cerr <<  "-mmsc <Number>\n" << ::std::endl;
	::std::cerr <<  "\tMismatch penalty for Dna / Rna alphabet, default is -4.\n\n" << ::std::endl;
	::std::cerr <<  "-rescore [ true | false ]\n" << ::std::endl;
	::std::cerr <<  "\tRe-score all segment-matches, default is true.\n\n" << ::std::endl;
	// Guide Tree
	::std::cerr <<  "Guide Tree Options\n------------\n" << ::std::endl;
	::std::cerr <<  "-usetree <Newick Guide Tree>\n" << ::std::endl;
	::std::cerr <<  "\tA guide tree in newick format.\n\n" << ::std::endl;
	// Alignment evaluation
	::std::cerr <<  "Alignment Evaluation Options\n------------\n" << ::std::endl;
	::std::cerr <<  "-infile <FASTA Alignment File>\n" << ::std::endl;
	::std::cerr <<  "\tEvaluate the given multiple alignment.\n\n" << ::std::endl;


	::std::cerr <<  "\n\n\nExamples\n" << ::std::endl;
	::std::cerr <<  "\nProtein Alignment:\n" << ::std::endl;
	::std::cerr <<  "\t./seqan_tcoffee -seq seq.fasta\n" << ::std::endl;
	::std::cerr <<  "\t./seqan_tcoffee -seq seq.fasta -method global, local\n" << ::std::endl;
	::std::cerr <<  "\nDna Alignment:\n" << ::std::endl;
	::std::cerr <<  "\t./seqan_tcoffee -seq seq.fasta -alphabet dna\n" << ::std::endl;
	::std::cerr <<  "\nGenome Alignment:\n" << ::std::endl;
	::std::cerr <<  "\t./seqan_tcoffee -seq seq.fasta -method lcs -alphabet dna\n" << ::std::endl;
	::std::cerr <<  "\t./seqan_tcoffee -seq seq.fasta -mummer my.mums -blast my1, my2 -alphabet dna\n" << ::std::endl;
	::std::cerr <<  "\nMeta-Alignment:\n" << ::std::endl;
	::std::cerr <<  "\t./seqan_tcoffee -seq seq.fasta -aln sub1.fasta, sub2.fasta\n" << ::std::endl;
	::std::cerr <<  "\nAlignment Evaluation:\n" << ::std::endl;
	::std::cerr <<  "\t./seqan_tcoffee -infile align.fasta\n" << ::std::endl;
	::std::cerr <<  "\t./seqan_tcoffee -infile align.fasta -gop -5 -gex -4 -alphabet dna\n" << ::std::endl;
	::std::cerr <<  "\nAlignment from a T-Coffee Library:\n" << ::std::endl;
	::std::cerr <<  "\t./seqan_tcoffee -lib tcoffee.tc_lib\n" << ::std::endl;
}


//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TScore>
inline void
customizedMsaAlignment(MsaOptions<TAlphabet, TScore> const& msaOpt) {
	typedef String<TAlphabet> TSequence;
	StringSet<TSequence, Owner<> > sequenceSet;
	StringSet<String<char> > sequenceNames;

	// Read the sequences from a T-Coffee library file or a regular sequence file in FASTA format
	if ((msaOpt.seqfile.empty()) && (!empty(msaOpt.libfiles))) {
		std::fstream strm;
		strm.open(value(msaOpt.libfiles, 0).c_str(), ::std::ios_base::in | ::std::ios_base::binary);
		read(strm, sequenceSet, sequenceNames, TCoffeeLib());	
		strm.close();			
	} else if (!msaOpt.seqfile.empty()) {
		_loadSequences(msaOpt.seqfile, sequenceSet, sequenceNames);
	} else {
		::std::cerr << "No input sequences!" << ::std::endl;
		exit(1);
	}
	

#ifdef SEQAN_PROFILE
	std::cout << "Number of sequences: " << length(sequenceSet) << std::endl;
	std::cout << "Import of sequences done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

	// Alignment of the sequences
	Graph<Alignment<StringSet<TSequence, Dependent<> >, void, WithoutEdgeId> > gAlign;
	
	// MSA
	globalMsaAlignment(gAlign, sequenceSet, sequenceNames, msaOpt);
		
	// Alignment output
	if (msaOpt.outputFormat == 0) {
		std::fstream strm;
		strm.open(msaOpt.outfile.c_str(), ::std::ios_base::out | ::std::ios_base::trunc);
		write(strm, gAlign, sequenceNames, FastaFormat());
		strm.close();
	} else if (msaOpt.outputFormat == 1) {
		std::fstream strm;
		strm.open(msaOpt.outfile.c_str(), ::std::ios_base::out | ::std::ios_base::trunc);
		write(strm, gAlign, sequenceNames, MsfFormat());
		strm.close();
	}

#ifdef SEQAN_PROFILE
	std::cout << "Output done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TScore, typename TSc>
inline void
_setMatchScore(MsaOptions<TAlphabet, TScore>&, TSc) {
	// No operation
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TScore, typename TSc>
inline void
_setMismatchScore(MsaOptions<TAlphabet, TScore>&, TSc) {
	// No operation
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSc>
inline void
_setMatchScore(MsaOptions<TAlphabet, Score<int, Simple> >& msaOpt, TSc msc) {
	msaOpt.sc.data_match = msc;
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSc>
inline void
_setMismatchScore(MsaOptions<TAlphabet, Score<int, Simple> >& msaOpt, TSc mmsc) {
	msaOpt.sc.data_mismatch = mmsc;
}


//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TScore>
inline void
_initMsaParams(int argc, const char *argv[], MsaOptions<TAlphabet, TScore>& msaOpt) {
	// Set default options
	msaOpt.sc.data_gap_open = -11;
	msaOpt.sc.data_gap_extend = -1;
	_setMatchScore(msaOpt, 5);
	_setMismatchScore(msaOpt, -4);
	msaOpt.outfile = "out.fasta";
	
	// Command line parsing
	for(int arg = 1; arg < argc; ++arg) {
		if (argv[arg][0] == '-') {
			if (strcmp(argv[arg], "-seq") == 0) {
				if (arg + 1 < argc) {
					++arg;
					msaOpt.seqfile = argv[arg];
				}
			}
			else if (strcmp(argv[arg], "-aln") == 0) {
				while ((arg+1 < argc) && (argv[arg+1][0] != '-')) {
					++arg;
					::std::string aln = argv[arg];
					::std::string::size_type pos = aln.find(',');
					if (pos != ::std::string::npos) aln.erase(pos);
					appendValue(msaOpt.alnfiles, aln);
				}
			}
			else if (strcmp(argv[arg], "-lib") == 0) {
				while ((arg+1 < argc) && (argv[arg+1][0] != '-')) {
					++arg;
					::std::string lib = argv[arg];
					::std::string::size_type pos = lib.find(',');
					if (pos != ::std::string::npos) lib.erase(pos);
					appendValue(msaOpt.libfiles, lib);
				}
			}
			else if (strcmp(argv[arg], "-blast") == 0) {
				while ((arg+1 < argc) && (argv[arg+1][0] != '-')) {
					++arg;
					::std::string blast = argv[arg];
					::std::string::size_type pos = blast.find(',');
					if (pos != ::std::string::npos) blast.erase(pos);
					appendValue(msaOpt.blastfiles, blast);
				}
			}
			else if (strcmp(argv[arg], "-mummer") == 0) {
				while ((arg+1 < argc) && (argv[arg+1][0] != '-')) {
					++arg;
					::std::string mummer = argv[arg];
					::std::string::size_type pos = mummer.find(',');
					if (pos != ::std::string::npos) mummer.erase(pos);
					appendValue(msaOpt.mummerfiles, mummer);
				}
			}
			else if (strcmp(argv[arg], "-usetree") == 0) {
				if (arg + 1 < argc) {
					++arg;
					msaOpt.treefile = argv[arg];
				}
			}
			else if (strcmp(argv[arg], "-outfile") == 0) {
				if (arg + 1 < argc) {
					++arg;
					assign(msaOpt.outfile, argv[arg]);
				}
			}
			else if (strcmp(argv[arg], "-method") == 0) {
				while ((arg+1 < argc) && (argv[arg+1][0] != '-')) {
					++arg;
					::std::string method = argv[arg];
					::std::string::size_type pos = method.find(',');
					if (pos != ::std::string::npos) method.erase(pos);
					if (method == "global") appendValue(msaOpt.method, 0);
					else if (method == "local") appendValue(msaOpt.method, 1);
					else if (method == "overlap") appendValue(msaOpt.method, 2);
					else if (method == "lcs") appendValue(msaOpt.method, 3);
				}
			}
			else if (strcmp(argv[arg], "-output") == 0) {
				if (arg + 1 < argc) {
					++arg;
					::std::string output = argv[arg];
					if (output == "fasta") msaOpt.outputFormat = 0;
					else if (output == "msf") msaOpt.outputFormat = 1;
				}
			}
			else if (strcmp(argv[arg], "-infile") == 0) {
				if (arg + 1 < argc) {
					++arg;
					msaOpt.infile = argv[arg];
				}
			}
			else if (strcmp(argv[arg], "-gop") == 0) {
				if (arg + 1 < argc) {
					++arg;
					::std::istringstream istr(argv[arg]);
					istr >> msaOpt.sc.data_gap_open;
					if (istr.fail()) { printHelp(); exit(1); }
				}
			}
			else if (strcmp(argv[arg], "-gex") == 0) {
				if (arg + 1 < argc) {
					++arg;
					::std::istringstream istr(argv[arg]);
					istr >> msaOpt.sc.data_gap_extend;
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
					_setMatchScore(msaOpt, msc);
				}
			}
			else if (strcmp(argv[arg], "-mmsc") == 0) {
				if (arg + 1 < argc) {
					++arg;
					int mmsc;
					::std::istringstream istr(argv[arg]);
					istr >> mmsc;
					if (istr.fail()) { printHelp(); exit(1); }
					_setMismatchScore(msaOpt, mmsc);
				}
			}
			else if (strcmp(argv[arg], "-rescore") == 0) {
				if (arg + 1 < argc) {
					++arg;
					if (strcmp(argv[arg], "false") == 0) {
						msaOpt.rescore = false;
					} else {
						msaOpt.rescore = true;
					}
				}
			}
		}
	}

	// Check if any segment-match generation procedure is selected, otherwise set the default
	if ((empty(msaOpt.blastfiles)) && (empty(msaOpt.mummerfiles)) && (empty(msaOpt.libfiles)) && (empty(msaOpt.alnfiles)) && (empty(msaOpt.method))) {
		appendValue(msaOpt.method, 0);
		appendValue(msaOpt.method, 1);
	}

	// Evaluation mode?
	if (!empty(msaOpt.infile)) {
#ifdef SEQAN_PROFILE
		::std::cout << "Alignment evaluation" << ::std::endl;
#endif
		evaluateAlignment(msaOpt);
	} else { // or alignment mode?
#ifdef SEQAN_PROFILE
		::std::cout << "Multiple Sequence Alignment" << ::std::endl;
#endif
		customizedMsaAlignment(msaOpt);
	}
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TMatrixFile>
inline void
_initScoreMatrix(int argc, const char *argv[], TMatrixFile& matrix, Dna5 const) {
	if (!empty(matrix)) {
		MsaOptions<Dna5, Score<int, ScoreMatrix<> > > msaOpt;
		loadScoreMatrix(msaOpt.sc, matrix);
		_initMsaParams(argc, argv, msaOpt);
	} else {
		MsaOptions<Dna5, Score<int> > msaOpt;
		_initMsaParams(argc, argv, msaOpt);
	}
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TMatrixFile>
inline void
_initScoreMatrix(int argc, const char *argv[], TMatrixFile& matrix, Rna5 const) {
	if (!empty(matrix)) {
		MsaOptions<Rna5, Score<int, ScoreMatrix<> > > msaOpt;
		loadScoreMatrix(msaOpt.sc, matrix);
		_initMsaParams(argc, argv, msaOpt);
	} else {
		MsaOptions<Rna5, Score<int> > msaOpt;
		_initMsaParams(argc, argv, msaOpt);
	}
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TMatrixFile>
inline void
_initScoreMatrix(int argc, const char *argv[], TMatrixFile& matrix, AminoAcid const) {
	if (!empty(matrix)) {
		MsaOptions<AminoAcid, Score<int, ScoreMatrix<> > > msaOpt;
		loadScoreMatrix(msaOpt.sc, matrix);
		_initMsaParams(argc, argv, msaOpt);
	} else {
		MsaOptions<AminoAcid, Blosum62> msaOpt;
		_initMsaParams(argc, argv, msaOpt);
	}
}

//////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char *argv[]) {

	// Version
#ifdef SEQAN_PROFILE
	printVersion();
#endif

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
					assign(alphabet, argv[arg]);
				}
			}
			else if (strcmp(argv[arg], "-matrix") == 0) {
				if (arg + 1 < argc) {
					++arg;
					assign(matrix, argv[arg]);
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
