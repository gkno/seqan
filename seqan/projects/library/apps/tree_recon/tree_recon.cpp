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
#include <seqan/graph_msa.h>
#include <seqan/modifier.h>
#include <seqan/misc/misc_cmdparser.h>

#include <iostream>
#include <fstream>


using namespace seqan;

//////////////////////////////////////////////////////////////////////////////////

inline void
_addVersion(CommandLineParser& parser) {
	::std::string rev = "$Revision: 4566 $";
	addVersionLine(parser, "Version 1.0 (15. July 2009) Revision: " + rev.substr(11, 4) + "");
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TMat, typename TNames>
inline void 
_readPhylipMatrix(TFile& file,
				  TMat& matrix,
				  TNames& names)
{
	typedef typename Value<TFile>::Type TValue;
	typedef typename Value<TMat>::Type TDistance;
	typedef typename Size<TMat>::Type TSize;
	typedef typename Value<TNames>::Type TName;
	typedef typename Iterator<TMat, Standard>::Type TMatIter;

	// Parse the file and convert the internal ids
	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);
	while (!_streamEOF(file)) {
		if (_streamEOF(file)) break;
		_parse_skipWhitespace(file, c);
		TSize nseq = _parse_readNumber(file, c);
		_parse_skipWhitespace(file, c);
		resize(matrix, nseq * nseq);
		resize(names, nseq);
		TMatIter it = begin(matrix, Standard());
		for(TSize row = 0; row<nseq; ++row) {
			names[row] = "label = \"";
			_parse_readIdentifier(file, names[row], c);
			append(names[row], '"');
			_parse_skipWhitespace(file, c);
			for(TSize col = 0; col<nseq; ++col, ++it) {
				*it = _parse_readDouble(file, c);
				_parse_skipWhitespace(file, c);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char *argv[]) {

	// Command line parsing
	CommandLineParser parser;
	_addVersion(parser);
	
	addTitleLine(parser, "***************************************");
	addTitleLine(parser, "* Tree reconstrucion - TreeRecon      *");
	addTitleLine(parser, "* (c) Copyright 2009 by Tobias Rausch *");
	addTitleLine(parser, "***************************************");

	addUsageLine(parser, "-m <Phylip distance matrix> [Options]");

	addSection(parser, "Main Options:");
	addOption(parser, addArgumentText(CommandLineOption("m", "matrix", "file with distance matrix", OptionType::String), "<Phylip distance matrix>"));
	addOption(parser, addArgumentText(CommandLineOption("b", "build", "tree building method", OptionType::String, "nj"), "[nj, min, max, avg, wavg]"));
	addHelpLine(parser, "nj = Neighbor-joining");
	addHelpLine(parser, "min = UPGMA single linkage");
	addHelpLine(parser, "max = UPGMA complete linkage");
	addHelpLine(parser, "avg = UPGMA average linkage");
	addHelpLine(parser, "wavg = UPGMA weighted average linkage");
	addHelpLine(parser, "/*Neighbor-joining creates an");
	addHelpLine(parser, "  unrooted tree. We root that tree");
	addHelpLine(parser, "  at the last joined pair.*/");
	addOption(parser, addArgumentText(CommandLineOption("o", "outfile", "output filename", OptionType::String, "tree.dot"), "<Filename>"));
		
	if (argc == 1)
	{
		shortHelp(parser, std::cerr);	// print short help and exit
		return 0;
	}

	if (!parse(parser, argc, argv, ::std::cerr)) return 1;
	if (isSetLong(parser, "help") || isSetLong(parser, "version")) return 0;	// print help or version and exit

	// Tree reconstruction
	typedef double TDistanceValue;
	typedef String<char> TName;
	typedef Size<TName>::Type TSize;

	// Read the options	
	::std::string infile;
	getOptionValueLong(parser, "matrix", infile);
	::std::string outfile;
	getOptionValueLong(parser, "outfile", outfile);
	TSize build = 0;
	String<char> meth;
	getOptionValueLong(parser, "build", meth);
	if (meth == "nj") build = 0;
	else if (meth == "min") build = 1;
	else if (meth == "max") build = 2;
	else if (meth == "avg") build = 3;
	else if (meth == "wavg") build = 4;

	// Read the distance matrix
	String<TName> names;
	String<TDistanceValue> matrix;
	FILE* strmMat = fopen(infile.c_str(), "rb");
	_readPhylipMatrix(strmMat, matrix, names);	
	fclose(strmMat);

	// Create the tree
	Graph<Tree<TDistanceValue> > tree;
	if (build == 0) njTree(matrix, tree);
	else if (build == 1) upgmaTree(matrix, tree, UpgmaMin());
	else if (build == 2) upgmaTree(matrix, tree, UpgmaMax());
	else if (build == 3) upgmaTree(matrix, tree, UpgmaAvg());
	else if (build == 4) upgmaTree(matrix, tree, UpgmaWeightAvg());
	
	// Append emty names for internal vertices
	TSize nameLen = length(names);
	resize(names, numVertices(tree));
	for(;nameLen < length(names); ++nameLen) {
		names[nameLen] = "label = \"\"";
	}

	// Write the result
	FILE* strmDot;
	strmDot = fopen(outfile.c_str(), "w");
	write(strmDot, tree, names, DotDrawing());
	fclose(strmDot);

	return 0;
}
