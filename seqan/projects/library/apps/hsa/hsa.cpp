#define SEQAN_PROFILE

#include <iostream>
#include <fstream>
#include <seqan/graph_msa.h>
#include <seqan/misc/misc_cmdparser.h>
#include "hsa.h"

using namespace seqan;

//template<typename TAlignmentGraph, typename TSegmentSet, typename TSegment, typename TSize, typename TPositions>
//void determineFirst(TAlignmentGraph & g,
//					TSegmentSet & segs,
//					TSegment & segment,
//					TSize seqId,
//					TPositions & first)
//{
//	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;
//	typedef typename Iterator<TAlignmentGraph, OutEdgeIterator>::Type TOutEdgeIterator;
//
//	int first_i = value( first, seqId );
//	if( first_i < (int)beginPosition( segment ) )
//	{
//		value( first, seqId ) = beginPosition( segment );
//
//		TVertexDescriptor source = findVertex( g, seqId, first_i +1 );
//		while( source != getNil<TVertexDescriptor>() && fragmentBegin(g,source) < beginPosition(segment) )
//		{
//			if( outDegree(g,source) > 0 )
//			{
//				TOutEdgeIterator eIt(g,source);
//				for( ; !atEnd(eIt); goNext(eIt) )
//				{
//					// find next unmatched segment on target sequence of the edge
//					TVertexDescriptor target = targetVertex( eIt );
//					TSegment newSegment;
//					if( fragmentBegin( g, target ) + fragmentLength( g, target ) <
//						length( value( stringSet(g), sequenceId( g, target ) ) ) )
//					{
//						target = findVertex( g, sequenceId( g, target ), 
//							fragmentBegin( g, target ) + fragmentLength( g, target ) );
//						newSegment = infix( value( stringSet(g), sequenceId(g,target) ),
//							fragmentBegin(g,target),
//							fragmentBegin(g,target) + fragmentLength( g, target ));
//					}
//					else
//					{
//						newSegment = infix( value( stringSet(g), sequenceId(g,target) ),
//							fragmentBegin(g,target) + fragmentLength( g, target ),
//							fragmentBegin(g,target) + fragmentLength( g, target ));
//					}
//					TSize newSeqId = sequenceId(g,target);
//					determineFirst( g, segs, newSegment, newSeqId, first );
//				}
//			}
//			source = findVertex( g, seqId, fragmentBegin(g,source) + fragmentLength(g,source) );
//		}
//	}
//}
//
//template<typename TAlignmentGraph, typename TSegmentSet, typename TSegment, typename TSize, typename TPositions>
//void determineLast(TAlignmentGraph & g,
//					TSegmentSet & segs,
//					TSegment & segment,
//					TSize seqId,
//					TPositions & last)
//{
//	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;
//	typedef typename Iterator<TAlignmentGraph, OutEdgeIterator>::Type TOutEdgeIterator;
//
//	TSize last_i = value( last, seqId );
//	if( last_i > endPosition( segment ) )
//	{
//		value( last, seqId ) = endPosition( segment );
//
//		TVertexDescriptor source = findVertex( g, seqId, endPosition(segment) );
//		while( source != getNil<TVertexDescriptor>() && fragmentBegin(g,source) < last_i ) 
//		{
//			if( outDegree(g,source) > 0 )
//			{
//				TOutEdgeIterator eIt(g,source);
//				for( ; !atEnd(eIt); goNext(eIt) )
//				{
//					// find closest leading and unmatched segment on target sequence of the edge
//					TVertexDescriptor target = targetVertex( eIt );
//					TSegment newSegment;
//					if( fragmentBegin( g, target ) > 0 )
//					{
//						target = findVertex( g, sequenceId( g, target ), 
//							fragmentBegin( g, target ) -1 );
//						newSegment = infix( value( stringSet(g), sequenceId(g,target) ), 
//							fragmentBegin(g,target), 
//							fragmentBegin(g,target) + fragmentLength( g, target ) );
//					}
//					else
//					{
//						newSegment = infix( value( stringSet(g), sequenceId(g,target) ), 
//							fragmentBegin(g,target), 
//							fragmentBegin(g,target) );
//					}
//					TSize newSeqId = sequenceId(g,target);
//					determineLast( g, segs, newSegment, newSeqId, last );
//				}
//			}
//			source = findVertex( g, seqId, fragmentBegin( g, source ) + fragmentLength( g, source ) );
//		}
//	}
//}
//
//template<typename TAlignmentGraph, typename TSegmentSet, typename TSize>
//TSize segSetId(TAlignmentGraph & g,
//			   TSegmentSet & segs,
//			   TSize hostId,
//			   TSize segBegin,
//			   TSize segEnd)
//{
//	typedef typename Iterator<TSegmentSet,Rooted>::Type TIter;
//	TIter iter = begin(segs);
//	TIter itEnd = end(segs);
//	
//	for( ; iter != itEnd; iter++ )
//	{
//		if( id(host( *iter )) == id( value(stringSet(g),hostId) ) &&
//			(TSize)beginPosition( *iter ) == segBegin &&
//			(TSize)endPosition( *iter ) == segEnd)
//		{
//			return position( iter );
//		}
//	}
//	//cout << "schlecht" << endl;
//	return 0;
//}
//
//template<typename TSegmentSet, typename TSegPairList, typename TAlignmentGraph, typename TSize>
//void findUnalignedRegions(TSegmentSet & segs,
//						  TSegPairList & segPairList,
//						  TAlignmentGraph & g,
//						  TSegmentSet & segInfixes,
//						  TSegPairList & newPairList,
//						  TSize minLength )
//{
//	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;
//	typedef typename Iterator<TAlignmentGraph, VertexIterator>::Type TVertexIterator;
//	typedef typename Value<TSegPairList>::Type TSegPair;
//
//	// determine the unmatched segments and append them to segInfixes
//	TVertexIterator vertexIter(g);
//	for( ; !atEnd(vertexIter); goNext(vertexIter) )
//	{
//		if( outDegree( g, *vertexIter ) == 0 )
//		{
//			TSize seqId = sequenceId( g, *vertexIter );
//			TSize segBeg = fragmentBegin( g, *vertexIter );
//			TSize segLength = fragmentLength( g, *vertexIter );
//
//			if( segLength >= minLength ) 
//			{
//				appendValue( segInfixes, infix( value( stringSet(g), seqId ), segBeg, segBeg+segLength ) );
//				//cout << length(segInfixes)-1 << ": " << infix( value( stringSet(g), seqId ), segBeg, segBeg+segLength ) << endl;
//			}
//		}
//	}
//
//	// determine pairs of segments from segInfixes that shall be compared
//	for( TSize pair = 0; pair < length( segPairList ); pair++ )
//	{
//		typedef typename Value<TSegmentSet>::Type TSegment;
//		TSize segmentId1 = value( segPairList, pair ).first;
//		TSize segmentId2 = value( segPairList, pair ).second;
//		TSize seqId1;
//		for( TSize s = 0; s < length( stringSet(g) ); s++ )
//		{
//			if( id( value( stringSet(g), s ) ) == id( segs[ segmentId1 ] ) )
//			{
//				seqId1 = s;
//			}
//		}
//		TSize seqId2;
//		for( TSize s = 0; s < length( stringSet(g) ); s++ )
//		{
//			if( id( value( stringSet(g), s ) ) == id( segs[ segmentId2 ] ) )
//			{
//				seqId2 = s;
//			}
//		}
//		TSegment seg1 = value( segs, segmentId1 );
//		TSegment seg2 = value( segs, segmentId2 );
//
//		// iterate over unmatched segments (s1) in seg1
//		TSize pos1 = beginPosition(seg1);
//		while(pos1 < endPosition(seg1))
//		{
//			TVertexDescriptor v = findVertex(g,seqId1,pos1);
//			pos1 += fragmentLength(g,v);
//			
//			TSize begin = fragmentBegin( g, v );
//			TSize len = fragmentLength( g, v );
//
//			if( outDegree( g, v ) == 0 && len >= minLength ) 
//			{
//				TSegment s1 = infix( host(seg1), begin, begin+len );
//				TSize s1Id = segSetId( g, segInfixes, seqId1, begin, begin+len );
//			
//				// positions in the string that mark the current relevant regions
//				typedef String<TSize> TPositions;
//				TPositions first, last;
//
//				// initialize first and last
//				resize( first, length(segs) );
//				resize( last, length(segs) );
//				for( int i = 0; i < length( stringSet(g) ); i++ )
//				{
//					value( first, i ) = -1;
//					value( last, i ) = length( value(stringSet(g),i) );
//				}
//
//				determineFirst( g, segs, s1, seqId1, first );
//				determineLast( g, segs, s1, seqId1, last ); // bis hierher unabhaengig von zweiter Sequenz, koennte man also schon frueher machen...
//					
//				//cout << "s1=" << s1 << "  begin(s1)=" << beginPosition(s1) << "  end(s1)=" << endPosition(s1) << endl;
//				//cout << "first(2)=" << value( first, seqId2 ) << " last(2)=" << value( last, seqId2 ) << endl;
//
//				// find all unmatched segments between first and last
//				TSize pos2 = value( first, seqId2 );
//				if(pos2 < 0) pos2 = 0;
//				while( pos2 < value( last, seqId2 ))
//				{
//					TVertexDescriptor v2 = findVertex(g,seqId2,pos2);
//					pos2 += fragmentLength(g,v2);
//
//					TSize begin2 = fragmentBegin(g,v2);// + beginPosition(seg2);
//					TSize len2 =  fragmentLength(g,v2);
//
//					if( outDegree( g, v2 ) == 0 && len2 >= minLength )
//					{
//						TSize s2Id = segSetId( g, segInfixes, seqId2, begin2, begin2+len2 );
//	
//						appendValue( newPairList, TSegPair(s1Id,s2Id) );
//						//cout << "appended to pair comparison list: " << s1Id << "," << s2Id << endl;
//					}
//				}
//				clear(last);
//				clear(first);
//			}
//		}
//	}
//}

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


///////////////////////////////////////////////////////////////////////////////
// Parses options from command line parser and writes them into options object
template<typename TParser, typename TOptions>
bool
_parseOptions(TParser & parser, TOptions & options) {

	CharString sequenceFiles;
	getOptionValueShort(parser, 's', sequenceFiles);
	parseSequenceFileNames(sequenceFiles, options.fileNames);

    if (isSetShort(parser, 'o')) getOptionValueShort(parser, 'o', options.outputFile);

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

	addSection(parser, "Parameter Options:");
	addOption(parser, CommandLineOption('r', "recursions", "Number of recursion steps", OptionType::Int, 3));
	addOption(parser, CommandLineOption('l', "minLength", "Initial minimum match length", OptionType::Int, 100));
	addOption(parser, CommandLineOption("dl", "deltaLength", "Step size for minimal match length (decreasing)", OptionType::Int, 30));
	addOption(parser, CommandLineOption('e', "epsilon", "(Initial) error rate", OptionType::Double, "0.1"));
	addOption(parser, CommandLineOption("de", "deltaEps", "Step size for error rate (increasing)", OptionType::Double, "0.0"));

	addSection(parser, "Algorithm Options:");
	addOption(parser, CommandLineOption("gt", "globalTree", "Use guide tree of full sequences in all recursions", OptionType::Bool, "false"));
	addOption(parser, CommandLineOption("f", "fixed", "Guarantee that higher level matches are contained in final alignment", OptionType::Bool, "false"));
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
		totalSeqLength += length(seqs[i]);
	}

    SEQAN_PROTIMESTART(timeRecSegmAlign);
	unsigned iterations = 1;
	
	TAlignmentGraph g(seqs);
	for (unsigned t = 0; t < iterations; t++) {
		clearEdges(g); clearVertices(g);
		recurseSegmentAlignment(segments, g, options, 0u, g);
	}

	double runningTime = SEQAN_PROTIMEDIFF(timeRecSegmAlign)/(double)iterations;

	//std::cout << g << std::endl;

	std::fstream fstream;
	fstream.open(toCString(options.outputFile), std::ios_base::out);
	if (!fstream.is_open()) {
		std::cerr << "Could not open output file: " << options.outputFile << std::endl;
	} else {
		write(fstream, g, DotDrawing());
	}
	fstream.close();

	std::cout << "running time: " << runningTime << "s" << std::endl;
	std::cout << "# edges: " << numEdges(g) << std::endl;
	std::cout << "# segments: " << numVertices(g) << std::endl;
	std::cout << "average segment length: " << totalSeqLength/(double)numVertices(g) << std::endl;

	Iterator<TAlignmentGraph, VertexIterator>::Type itV(g);

	unsigned int matchedChars = 0, matchedVertices = 0;
	for (; !atEnd(itV); goNext(itV)) {
		if (outDegree(g, *itV) > 0) {
			matchedVertices++;
			matchedChars += fragmentLength(g, *itV);
		}
	}
	std::cout << "matched characters: " << matchedChars * 100 / (double)totalSeqLength << "%" << std::endl;
	std::cout << "average length of matched segments: " << matchedChars / (double)matchedVertices << std::endl;
	std::cout << "average length of unmatched segments: ";
	std::cout << (totalSeqLength-matchedChars) / (double)(numVertices(g)-matchedVertices) << std::endl;

	std::cout << length(segments);
	std::cout << " " << totalSeqLength / (double)length(segments);
	std::cout << " " << options.initialMinLength;
	std::cout << " " << options.initialMinLength - options.deltaMinLength * (options.recursions-1);
	std::cout << " " << options.deltaMinLength;
	std::cout << " " << matchedChars * 100 / (double)totalSeqLength;
	std::cout << " " << totalSeqLength / (double)numVertices(g);
	std::cout << " " << matchedChars / (double)matchedVertices;
	std::cout << " " << (totalSeqLength-matchedChars) / (double)(numVertices(g)-matchedVertices);
	std::cout << " " << runningTime;
	std::cout << " " << options.recursions;
	std::cout << std::endl;

	return 0;
}
