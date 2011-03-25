 /*==========================================================================
                  HSA - Hierarchical Segment-based Alignment

 ============================================================================
  Copyright (C) 2011 by Birte Kehr

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

#define SEQAN_PROFILE

#include <iostream>
#include <fstream>
#include <seqan/graph_msa.h>
#include <seqan/misc/misc_cmdparser.h>
#include "hsa.h"

using namespace seqan;


///////////////////////////////////////////////////////////////////////////////
// Parses a string of comma separated file names and appends them to options.fileNames
template<typename TCharString1, typename TCharString2>
void parseSequenceFileNames(TCharString1 & files, String<TCharString2> & str) {
	typedef typename Iterator<TCharString1>::Type TIterator;
	TIterator it = begin(files);
	TIterator itEnd = end(files);

	TCharString2 fileName;
	while (it != itEnd) {
		if (*it == ',') {
			if (length(fileName) > 0) {
				appendValue(str, fileName);
			}
			clear(fileName);
		}
		else {
			appendValue(fileName, *it);
		}
		++it;
	}
	if (length(fileName) > 0) {
		appendValue(str, fileName);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
// Computes and outputs some statistics of the multiple alignment
template<typename TSequence, typename TSize>
void _printAlignmentStatistics(Graph<Alignment<TSequence> > & g, TSize totalSeqLength) {
	typedef typename Iterator<Graph<Alignment<TSequence> >, VertexIterator>::Type TIterator;

	// compute number of characters per vertex degree
	String<TSize> matchedChars, matchedVertices;
	resize(matchedChars, length(stringSet(g)), 0);
	resize(matchedVertices, length(stringSet(g)), 0);
	for (TIterator itV(g); !atEnd(itV); goNext(itV)) {
		TSize deg = outDegree(g, *itV);
		matchedVertices[deg]++;
		matchedChars[deg] += fragmentLength(g, *itV);
	}

	std::cout << "# edges: " << numEdges(g) << std::endl;
	std::cout << "# segments: " << numVertices(g) << std::endl;
	std::cout << "Average segment length: " << totalSeqLength/(double)numVertices(g) << std::endl;
	std::cout << "Aligned with at least one other sequence: ";
	std::cout << (totalSeqLength-matchedChars[0]) * 100 / (double)totalSeqLength << "% of characters" << std::endl;
	std::cout << "Average length of aligned segments: ";
	std::cout << (totalSeqLength-matchedChars[0]) / (double)(numVertices(g)-matchedVertices[0]) << std::endl;
	std::cout << std::endl;

	for (TSize i = 0; i < length(matchedChars); ++i) {
		std::cout << "Aligned with " << i << " other sequences:" << std::endl;

		std::cout << "  #chars: ";
		std::cout << matchedChars[i]*100/(double)totalSeqLength << "% (" << matchedChars[i] << ")" << std::endl;

		std::cout << "  #vertices: ";
		std::cout << matchedVertices[i]*100/(double)numVertices(g) << "% (" << matchedVertices[i] << ")" << std::endl;

		std::cout << "  avg segment length: ";
		std::cout << matchedChars[i]/(double)matchedVertices[i] << std::endl;

		std::cout << std::endl;
	}
}


/////////////////////////////////////////////////////////////////////////////////////////
// Computes and outputs a pairwise distance matrix of a multiple alignment g
//   based on the sequence fraction in matched alignment graph vertices.
template<typename TSequence, typename TDistance>
void _computeDistanceMatrix(Graph<Alignment<TSequence> > & g,
							String<TDistance> & matrix) {
	typedef typename Size<String<TDistance> >::Type TSize;
	typedef Graph<Alignment<TSequence> > TAlignmentGraph;
	typedef typename Iterator<TAlignmentGraph, EdgeIterator>::Type TEdgeIt;
	typedef typename Id<TAlignmentGraph>::Type TId;

	TSize numSeq = length(stringSet(g));
	resize(matrix, numSeq*numSeq);

	// init counters for pairwisely matched characters
	String<String<TSize> > chars, edges;
	resize(chars, numSeq);
	resize(edges, numSeq);

	for (TSize i = 0; i < numSeq; ++i) {
		resize(chars[i], numSeq, 0);
		resize(edges[i], numSeq, 0);
	}

	// count pairwisely matched characters
	for (TEdgeIt itE(g); !atEnd(itE); goNext(itE)) {
		TId sourceId = sequenceId(g, sourceVertex(itE));
		TId targetId = sequenceId(g, targetVertex(itE));
		chars[sourceId][targetId] += fragmentLength(g, sourceVertex(itE)) + fragmentLength(g, targetVertex(itE));
		chars[targetId][sourceId] += fragmentLength(g, targetVertex(itE)) + fragmentLength(g, sourceVertex(itE));
		edges[sourceId][targetId]++;
		edges[targetId][sourceId]++;
	}

	// Output heading
	std::cout << "Number of pairwisely mapped characters:" << std::endl;
	std::cout << "    |\t";
	for (TSize j = 0; j < numSeq; ++j) {
		std::cout << j << "\t";
	}
	std::cout << std::endl;
	for (TSize j = 0; j <= numSeq; ++j) std::cout << "--------";
	std::cout << std::endl;

	// fill and output distance matrix
	for (TSize i = 0; i < numSeq; ++i) {
		std::cout << i << "   |\t";
		for (TSize j = 0; j < numSeq; ++j) {
			std::cout << chars[i][j] << "\t";
		}
		
		std::cout << "   |\t";
		for (TSize j = 0; j < numSeq; ++j) {
			TDistance totalSeqLen = length(stringSet(g)[i]) + length(stringSet(g)[j]);
			matrix[numSeq*i+j] = 100.0 - chars[i][j] / totalSeqLen;
			matrix[numSeq*j+i] = 100.0 - chars[i][j] / totalSeqLen;

			std::cout << chars[i][j] * 100.0 / totalSeqLen << "%\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// Writes parsed command line options to screen
template<typename TOptions>
void
_writeParams(TOptions & options) {
	typedef typename Size<CharString>::Type TSize;
	std::cout << "Sequence files: " << options.fileNames[0] << std::endl;
	for (TSize i = 1; i < length(options.fileNames); ++i) {
		std::cout << "                " << options.fileNames[i] << std::endl;
	}
	std::cout << std::endl;
	std::cout << "Guide tree    : " << (options.globalGuideTree?"global":"local") << std::endl;
	std::cout << "Fixed matches : " << (options.fixedHigherLevelMatches?"true":"false") << std::endl;
	std::cout << "Comparison    : " << (options.anchoredPairwiseComparison?"reduced":"full") << std::endl;
	std::cout << std::endl;
	std::cout << "Hier. levels  : " << options.recursions << std::endl;
	std::cout << "Initial minLen: " << options.initialMinLength << std::endl;
	std::cout << "Delta minLen  : " << options.deltaMinLength << std::endl;
	std::cout << "Initial eps   : " << options.initialEpsilon << std::endl;
	std::cout << "Delta eps     : " << options.deltaEpsilon << std::endl;
	std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// Parses options from command line parser and writes them into options object
template<typename TParser, typename TOptions>
bool
_parseOptions(TParser & parser, TOptions & options) {

	CharString sequenceFiles;
	getOptionValueShort(parser, 's', sequenceFiles);
	parseSequenceFileNames(sequenceFiles, options.fileNames);

    if (isSetShort(parser, 'o')) getOptionValueShort(parser, 'o', options.outputFile);
    if (isSetShort(parser, 'p')) getOptionValueShort(parser, 'p', options.phylogeny);
    if (isSetShort(parser, 'v')) getOptionValueShort(parser, 'v', options.verbose);

    if (isSetShort(parser, 'r')) getOptionValueShort(parser, 'r', options.recursions);
	if (isSetShort(parser, 'l')) getOptionValueShort(parser, 'l', options.initialMinLength);
	if (isSetShort(parser, 'e')) getOptionValueShort(parser, 'e', options.initialEpsilon);
	if (isSetShort(parser, "dl")) getOptionValueShort(parser, "dl", options.deltaMinLength);
	if (isSetShort(parser, "de")) getOptionValueShort(parser, "de", options.deltaEpsilon);

    if (isSetShort(parser, "gt")) getOptionValueShort(parser, "gt", options.globalGuideTree);
    if (isSetShort(parser, "f")) getOptionValueShort(parser, "f", options.fixedHigherLevelMatches);
    if (isSetShort(parser, "a")) getOptionValueShort(parser, "a", options.anchoredPairwiseComparison);

	return 1;
}

///////////////////////////////////////////////////////////////////////////////
// Set-Up of Command Line Parser
template<typename TParser>
void
_setParser(TParser & parser) {
    addTitleLine(parser, "****************************************");
	addTitleLine(parser, "* Hierarchical Segment-based Alignment *");
	addTitleLine(parser, "* (c) Copyright 2010 by Birte Kehr     *");
	addTitleLine(parser, "****************************************");

	addUsageLine(parser, "-s <FASTA sequence file>,...,<FASTA sequence file> [Options]");

	addLine(parser, "");
    addLine(parser, "Short description will follow soon.");

	addSection(parser, "File I/O Options:");
    addOption(parser, CommandLineOption('s', "seqs", "Fasta files containing sequences",
              (OptionType::String | OptionType::Mandatory)));
	addHelpLine(parser, "File names separated by comma.");
	addOption(parser, CommandLineOption('o', "outFile", "Output filename", OptionType::String, "ReSeAl.dot"));
	addOption(parser, CommandLineOption('p', "phylo", "Compute and output distance matrix and UPGMA tree from final alignment", OptionType::Bool, "false"));
	addOption(parser, CommandLineOption('v', "verbose", "Verbosity mode", OptionType::Bool, "false"));

	addSection(parser, "Parameter Options:");
	addOption(parser, CommandLineOption('r', "recursions", "Number of recursion steps", OptionType::Int, 3));
	addOption(parser, CommandLineOption('l', "minLength", "Initial minimum match length", OptionType::Int, 100));
	addOption(parser, CommandLineOption("dl", "deltaLength", "Step size for minimal match length (decreasing)", OptionType::Int, 30));
	addOption(parser, CommandLineOption('e', "epsilon", "(Initial) error rate", OptionType::Double, "0.1"));
	addOption(parser, CommandLineOption("de", "deltaEps", "Step size for error rate (increasing)", OptionType::Double, "0.0"));

	addSection(parser, "Algorithm Options:");
	addOption(parser, CommandLineOption("gt", "globalTree", "Use guide tree of full sequences in all recursions", OptionType::Bool, "false"));
	addOption(parser, CommandLineOption("f", "fixed", "Guarantee that previous level matches are contained in final alignment", OptionType::Bool, "false"));
	addOption(parser, CommandLineOption("a", "anchored", "Anchored pairwise comparison", OptionType::Bool, "false"));
}

int main (int argc, const char *argv[]) {
	// arguments: "Z:\GenomeData\NC_001405_short.fa" "Z:\GenomeData\NC_001460_short.fa" or
	//            "Z:\GenomeData\testSeq1.fa" "Z:\GenomeData\testSeq2.fa"
	typedef String<Dna5> TSequence;
	typedef StringSet<TSequence> TSequenceSet;

	typedef const void * TId;
    typedef Infix<TSequence>::Type TInfix;
	typedef std::map<TId, TInfix> TInfixMap;

	typedef StringSet<TSequence, Dependent<> > TDepSequenceSet;
	typedef Graph<Alignment<TDepSequenceSet> > TAlignmentGraph;

	typedef Size<TSequence>::Type TSize;

	// command line parsing
	CommandLineParser parser("stellar");

	_setParser(parser);
	if (!parse(parser, argc, argv)) {
		if (isSetShort(parser, 'h')) return 0; 
		shortHelp(parser, std::cerr);
		return 1;
	}

	MyOptions<TSequence> options = MyOptions<TSequence>();
	if (!_parseOptions(parser, options)) {
		return 1;
	}

	// output header
	_title(parser, std::cout);
	std::cout << std::endl;

	// output parameters
	_writeParams(options);

	TInfixMap segments;
	TSequenceSet seqs;
	resize(seqs, length(options.fileNames));
	TSize totalSeqLength = 0;

	// read sequences
	for (unsigned i = 0; i < length(options.fileNames); ++i) {
		std::fstream fstream;
		fstream.open(toCString(options.fileNames[i]), std::ios_base::in | std::ios_base::binary);
		if (fstream.is_open()) read(fstream, seqs[i], Fasta());
		else return 1;
		fstream.close();

		SEQAN_ASSERT_EQ(segments.count(id(seqs[i])), 0u);
		segments[id(seqs[i])] = infix(seqs[i], 0, length(seqs[i]));
		if (options.verbose) {
			std::cout << id(seqs[i]) << ": " << options.fileNames[i] << std::endl;
		}
		totalSeqLength += length(seqs[i]);
	}
	if (options.verbose) std::cout << std::endl << std::cout << std::endl;

    SEQAN_PROTIMESTART(timeRecSegmAlign);
	unsigned iterations = 1;
	
	// ----------- main algorithm call --------------
	TAlignmentGraph g(seqs);
	for (unsigned t = 0; t < iterations; t++) {
		clearEdges(g); clearVertices(g);
		recurseSegmentAlignment(segments, g, options, 0u, g);
	}

	double runningTime = SEQAN_PROTIMEDIFF(timeRecSegmAlign)/(double)iterations;

	//std::cout << g << std::endl;

	// output alignment graph to file in dot-format
	std::fstream fstream;
	fstream.open(toCString(options.outputFile), std::ios_base::out);
	if (!fstream.is_open()) {
		std::cerr << "Could not open output file: " << options.outputFile << std::endl;
	} else {
		write(fstream, g, DotDrawing());
	}
	fstream.close();

	// screen output
	std::cout << "Running time: " << runningTime << "s" << std::endl;
	_printAlignmentStatistics(g, totalSeqLength);

	if (options.phylogeny) {
		String<double> matrix;
		_computeDistanceMatrix(g, matrix);

		Graph<Tree<double> > tree;
		upgmaTree(matrix, tree);

		std::cout << "Phylogeny: " << std::endl;
		std::cout << tree << std::endl;
	}

	return 0;
}
