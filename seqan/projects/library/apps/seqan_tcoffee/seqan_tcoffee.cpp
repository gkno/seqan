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
#include "rna_alphabet.h"
#include "seqan_tcoffee.h"

using namespace seqan;

int main(int argc, const char *argv[]) {
	//////////////////////////////////////////////////////////////////////////////
	// Command line parsing
	//////////////////////////////////////////////////////////////////////////////
	
	// Set the keys
	typedef String<char> TKey;
	typedef String<char> TValue;
	typedef Size<TKey>::Type TSize;
	ConfigOptions<TKey, TValue> cfgOpt;
	TKey keys[] = {"seq", "alphabet", "aln", "lib","blast","mummer", "usetree","outfile","method","output","gop", "gex", "matrix", "infile", "msc", "mmsc", "rescore"};
	assignKeys(cfgOpt, keys, 17);

	// Set default options
	assign(cfgOpt, "rescore", "true");
	assign(cfgOpt, "output", "fasta");
	assign(cfgOpt, "outfile", "out.fasta");
	assign(cfgOpt, "alphabet", "protein");
	assign(cfgOpt, "gop", "11");
	assign(cfgOpt, "gex", "1");
	assign(cfgOpt, "msc", "5");
	assign(cfgOpt, "mmsc", "-4");

	// Help Message
	String<char> helpMsg;
	append(helpMsg, "Usage: seqan_tcoffee -seq <FASTA Sequence File> [Options]\n");
	append(helpMsg, "\nOptions\n");
	// Main options
	append(helpMsg, "\nMain Options\n------------\n");
	append(helpMsg, "-seq <FASTA Sequence File>\n");
	append(helpMsg, "\tFile with multiple sequences in FASTA format.\n\n");
	append(helpMsg, "-alphabet [protein | dna | rna]\n");
	append(helpMsg, "\tSequence alphabet, default is protein.\n\n");
	append(helpMsg, "-outfile <Alignment Filename>\n");
	append(helpMsg, "\tName of the output file, default is out.fasta.\n\n");
	append(helpMsg, "-output [fasta | msf]\n");
	append(helpMsg, "\tOutput format, default is fasta.\n\n");
	// Segment match generation options
	append(helpMsg, "\nSegment-Match Generation Options\n------------\n");
	append(helpMsg, "-method [global], [local], [overlap], [lcs]\n");
	append(helpMsg, "\tMethods to generate segment matches, default is global, local.\n\n");
	append(helpMsg, "-blast <BLAST matches>, <BLAST matches>, ...\n");
	append(helpMsg, "\tFiles with gapless segment matches in BLAST tabular format (-m 8 -g F).\n\n");
	append(helpMsg, "-mummer <MUMmer matches>, <MUMmer matches>, ...\n");
	append(helpMsg, "\tFiles with gapless segment matches in MUMmer format.\n\n");
	append(helpMsg, "-aln <FASTA Alignment File>, <FASTA Alignment File>, ...\n");
	append(helpMsg, "\tSegment-matches are extracted from given alignments (Meta-Alignment).\n\n");
	append(helpMsg, "-lib <T-Coffee library>, <T-Coffee library>, ...\n");
	append(helpMsg, "\tSegment-matches are extracted from T-Coffee libraries.\n\n");
	// Scoring
	append(helpMsg, "\nScoring Options\n------------\n");
	append(helpMsg, "-gop <Number>\n");
	append(helpMsg, "\tGap open penalty, default is 11.\n\n");
	append(helpMsg, "-gex <Number>\n");
	append(helpMsg, "\tGap extension penalty, default is 1.\n\n");
	append(helpMsg, "-matrix <Score-Matrix File>\n");
	append(helpMsg, "\tSpecifies a score-matrix file, default is Blosum62.\n\n");
	append(helpMsg, "-msc <Number>\n");
	append(helpMsg, "\tMatch score for Dna / Rna alphabet, default is 5.\n\n");
	append(helpMsg, "-mmsc <Number>\n");
	append(helpMsg, "\tMismatch penalty for Dna / Rna alphabet, default is 4.\n\n");
	append(helpMsg, "-rescore [ true | false ]\n");
	append(helpMsg, "\tRe-score all segment-matches, default is true.\n\n");
	// Guide Tree
	append(helpMsg, "Guide Tree Options\n------------\n");
	append(helpMsg, "-usetree <Newick Guide Tree>\n");
	append(helpMsg, "\tA guide tree in newick format.\n\n");
	// Alignment evaluation
	append(helpMsg, "Alignment Evaluation Options\n------------\n");
	append(helpMsg, "-infile <FASTA Alignment File>\n");
	append(helpMsg, "\tEvaluate the given multiple alignment.\n\n");


	append(helpMsg, "\n\n\nExamples\n");
	append(helpMsg, "\nProtein Alignment:\n");
	append(helpMsg, "\t./seqan_tcoffee -seq seq.fasta\n");
	append(helpMsg, "\t./seqan_tcoffee -seq seq.fasta -method global, local\n");
	append(helpMsg, "\nDna Alignment:\n");
	append(helpMsg, "\t./seqan_tcoffee -seq seq.fasta -alphabet dna\n");
	append(helpMsg, "\nGenome Alignment:\n");
	append(helpMsg, "\t./seqan_tcoffee -seq seq.fasta -method lcs -alphabet dna\n");
	append(helpMsg, "\t./seqan_tcoffee -seq seq.fasta -mummer my.mums -blast my1, my2 -alphabet dna\n");
	append(helpMsg, "\nMeta-Alignment:\n");
	append(helpMsg, "\t./seqan_tcoffee -seq seq.fasta -aln sub1.fasta, sub2.fasta\n");
	append(helpMsg, "\nAlignment Evaluation:\n");
	append(helpMsg, "\t./seqan_tcoffee -infile align.fasta\n");
	append(helpMsg, "\t./seqan_tcoffee -infile align.fasta -gop 5 -gex 4 -alphabet dna\n");
	append(helpMsg, "\nAlignment from a T-Coffee Library:\n");
	append(helpMsg, "\t./seqan_tcoffee -lib tcoffee.tc_lib\n");
	assignHelp(cfgOpt, helpMsg);

#ifdef SEQAN_PROFILE
	std::cout << "*********************************************" << std::endl;
	std::cout << "* Segment-based multiple sequence alignment *" << std::endl;
	std::cout << "*                                           *" << std::endl;
	std::cout << "* SeqAn::T-Coffee                           *" << std::endl;
	std::cout << "* Version: 1.1 (10. September 2008)         *" << std::endl;
	std::cout << "*********************************************" << std::endl;
	std::cout << std::endl;
#endif

	if (argc < 2) {	std::cerr << valueHelp(cfgOpt) << std::endl; return -1; }
	if (!parseCmdLine(argc, argv, cfgOpt)) return -1;

	// Check if any segment-match generation procedure is selected, otherwise set the default
	if ((!length(value(cfgOpt, "blast"))) && (!length(value(cfgOpt, "mummer"))) && (!length(value(cfgOpt, "method"))) && (!length(value(cfgOpt, "aln"))) && (!length(value(cfgOpt, "lib")))) {
		assign(cfgOpt, "method", "global, local");
	}

	//////////////////////////////////////////////////////////////////////////////
	// Evaluation mode?
	//////////////////////////////////////////////////////////////////////////////

	if (length(value(cfgOpt, "infile"))) {
		if ((value(cfgOpt, "alphabet") == "dna")) {
#ifdef SEQAN_PROFILE
			std::cout << "Alignment evaluation for Dna sequences" << std::endl;
#endif
			return evaluateAlignment(cfgOpt, Dna5() );
		} else if ((value(cfgOpt, "alphabet") == "rna")) {
#ifdef SEQAN_PROFILE	
			std::cout << "Alignment evaluation for Rna sequences" << std::endl;
#endif
			return evaluateAlignment(cfgOpt, Rna5() );
		} else {
#ifdef SEQAN_PROFILE
			std::cout << "Alignment evaluation for AminoAcid sequences" << std::endl;
#endif
			return evaluateAlignment(cfgOpt, AminoAcid() );
		}
	}

	//////////////////////////////////////////////////////////////////////////////
	// Alignment of Dna, Rna or Amino Acid Sequences
	//////////////////////////////////////////////////////////////////////////////

	if ((value(cfgOpt, "alphabet") == "dna")) {
#ifdef SEQAN_PROFILE
		std::cout << "Multiple Sequence Alignment of Dna sequences" << std::endl;
#endif
		return customizedMsaAlignment(cfgOpt, Dna5() );
	} else if ((value(cfgOpt, "alphabet") == "protein")) {
#ifdef SEQAN_PROFILE		
		std::cout << "Multiple Sequence Alignment of AminoAcid sequences" << std::endl;
#endif
		return customizedMsaAlignment(cfgOpt, AminoAcid() );
	} else if ((value(cfgOpt, "alphabet") == "rna")) {
#ifdef SEQAN_PROFILE		
		std::cout << "Multiple Sequence Alignment of Rna sequences" << std::endl;
#endif
		return customizedMsaAlignment(cfgOpt, Rna5() );
	} else {	std::cerr << valueHelp(cfgOpt) << std::endl; return -1; }
}
