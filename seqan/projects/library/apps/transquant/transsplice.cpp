//#define SEQAN_ENABLE_DEBUG 1
#include <fstream>
#include <iostream>
#include <seqan/store.h>
#include <seqan/misc/misc_cmdparser.h>

using namespace seqan;

// forward and reverse strands are separated,
// they have no common subexons
#define SEPARATE_STRANDS
#define COMPACT_NODE_IDS

#ifdef COMPACT_NODE_IDS
#define LARGE_NODE_ID "LargeNodeId"
#else
#define LARGE_NODE_ID "NodeId"
#endif

template <typename TValue>
inline TValue absInt(TValue t)
{
	return (t >= 0)? t: -t;
}

template <typename TSequence>
void removeEqualElements(TSequence &seq)
{
	typedef typename Iterator<TSequence, Standard>::Type TIter;

	TIter itBegin = begin(seq, Standard());
	TIter itEnd = end(seq, Standard());
	TIter src = itBegin;
	TIter dst = itBegin;
	for (; src != itEnd; ++src)
	{
		if (*dst != *src)
		{
			++dst;
			*dst = *src;
		}
	}
	if (itBegin != itEnd)
		resize(seq, dst - itBegin + 1);
}


// Step 3a: Get Exon Boundaries
//
// Go through all annotated exons and collect their boundaries
// group according to strand and orientation.
// Afterwards sort boundaries and remove duplicated elements

template <typename TExonBoundString, typename TFragmentStore>
void getExonBoundaries(TExonBoundString &exonBounds, TFragmentStore &store)
{
	typedef typename TFragmentStore::TAnnotationStore TAnnotationStore;
	typedef typename Iterator<TFragmentStore, AnnotationTree<> >::Type TAnnoTreeIter;
	typedef typename Value<TAnnotationStore>::Type TAnnotation;
	
	clear(exonBounds);
	TAnnoTreeIter dfsIt = begin(store, AnnotationTree<>());
		
	while (!atEnd(dfsIt))
	{
		TAnnotation &anno = getAnnotation(dfsIt);
		if (anno.typeId == TFragmentStore::ANNO_EXON)
		{
			// the reverse strands are interleaved with the forward strands
			if (anno.contigId == TAnnotation::INVALID_ID)
			{
				std::cerr << "Exon annotation has no contig (parentId=" << getParentName(dfsIt) << ")." << std::endl;
				continue;
			}
			unsigned contigNum = anno.contigId;
#ifdef SEPARATE_STRANDS
			contigNum *= 2;
			if (anno.beginPos > anno.endPos) ++contigNum;
#endif
			
			if (length(exonBounds) <= contigNum)
				resize(exonBounds, contigNum + 1);
			appendValue(exonBounds[contigNum], anno.beginPos);
			appendValue(exonBounds[contigNum], anno.endPos);
			goNextRight(dfsIt);
		} 
		else
			goNext(dfsIt);
	}
	for (unsigned i = 0; i < length(exonBounds); ++i)
	{
		std::sort(begin(exonBounds[i], Standard()), end(exonBounds[i], Standard()));
		removeEqualElements(exonBounds[i]);
	}
}


// Step 3b: Refine Exon Boundaries
//
// For every exon interval get all contained boundaries and
// divide exon into subexons using these boundaries.
// Insert subexons as children of the exon (even if there is only one subexon).
// NodeIds are set according to the boundary rank in the 
// sorted list of boundaries (the total rank is determined by a partial sum).

template <typename TExonBoundString, typename TFragmentStore, typename TNodeIds>
void refineExonBoundaries(TExonBoundString &exonBounds, TFragmentStore &store, TNodeIds &nodeIds)
{
	typedef typename Value<TExonBoundString>::Type TContigExonBounds;
	typedef typename Iterator<TContigExonBounds, Standard>::Type TContigExonBoundIterator;
	typedef typename TFragmentStore::TAnnotationStore TAnnotationStore;
	typedef typename  Iterator<TFragmentStore, AnnotationTree<> >::Type TAnnoTreeIter;
	typedef typename  Value<TAnnotationStore>::Type TAnnotation;
	typedef typename TFragmentStore::TContigPos TContigPos;
	typedef typename  Id<TAnnotation>::Type TId;
		
	String<int> partialSum;
	for (unsigned i = 0, sum = 0; i < length(exonBounds); ++i)
	{
		appendValue(partialSum, sum);
		sum += length(exonBounds[i]);
	}
	
	TAnnoTreeIter dfsIt = begin(store, AnnotationTree<>());
	while (!atEnd(dfsIt))
	{
		TAnnotation &anno = getAnnotation(dfsIt);
		if (anno.typeId == TFragmentStore::ANNO_EXON)
		{
			// the reverse strands are interleaved with the forward strands
			if (anno.contigId == TAnnotation::INVALID_ID)
			{
				std::cerr << "Exon annotation has no contig (parentId=" << getParentName(dfsIt) << ")." << std::endl;
				continue;
			}
			
			unsigned contigNum = anno.contigId;
#ifdef SEPARATE_STRANDS
			contigNum *= 2;
			if (anno.beginPos > anno.endPos) ++contigNum;
#endif

			TContigExonBoundIterator boundsBegin = begin(exonBounds[contigNum], Standard());
			TContigExonBoundIterator boundsEnd = end(exonBounds[contigNum], Standard());
			
			// determine exon orientation and order boundaries
			bool forward = true;
			TContigPos beginPos = anno.beginPos;
			TContigPos endPos = anno.endPos;
			if (beginPos > endPos)
			{
				forward = false;
				beginPos = endPos;
				endPos = anno.beginPos;
			}

			TContigExonBoundIterator first = std::upper_bound(boundsBegin, boundsEnd, beginPos);
			TContigExonBoundIterator last = std::lower_bound(boundsBegin, boundsEnd, endPos);
			
			TContigPos subExonBeginPos = beginPos;
			for (TContigExonBoundIterator it = first; it <= last; ++it)
			{
				std::stringstream tmp;
				TContigPos subExonEndPos = *it;
				TAnnoTreeIter childIt = createRightChild(dfsIt);
				clearValues(childIt);
				int nodeId;
				if (forward)
				{	// positive strand
					getAnnotation(childIt).beginPos = subExonBeginPos;
					getAnnotation(childIt).endPos = subExonEndPos;
					nodeId = (it - boundsBegin) + partialSum[contigNum];
				} else
				{	// negative strand
					getAnnotation(childIt).endPos = subExonBeginPos;
					getAnnotation(childIt).beginPos = subExonEndPos;
					nodeId = (boundsEnd - it - 1) + partialSum[contigNum];
				}
				SEQAN_ASSERT_LT(beginPos, endPos);
				SEQAN_ASSERT_LEQ(beginPos, subExonBeginPos);
				SEQAN_ASSERT_LEQ(subExonEndPos, endPos);
				tmp << nodeId;
#ifdef COMPACT_NODE_IDS
				append(nodeIds, nodeId);
#endif
				setType(childIt, "subexon");
				assignValueByKey(childIt, LARGE_NODE_ID, tmp.str());
				subExonBeginPos = subExonEndPos;
			}
			goNextRight(dfsIt);
		}
		else
			goNext(dfsIt);
	}
}

// Step 4: Calculate Orderings
//
// For every exon interval get all contained boundaries and
// divide exon into subexons using these boundaries.
// Fill ordering array with subexon ids for every transcript.
// NodeIds are set according to the boundary rank in the sorted list of boundaries.

template <typename TContigOrderings, typename TFragmentStore>
void calculateOrderings(TContigOrderings &orderings, TFragmentStore &store)
{
	typedef typename TFragmentStore::TAnnotationStore TAnnotationStore;
	typedef typename Iterator<TFragmentStore, AnnotationTree<> >::Type TAnnoTreeIter;
	typedef typename Value<TAnnotationStore>::Type TAnnotation;
	
	unsigned subExonType = 0;
	_storeAppendType(store, subExonType, "subexon");
	String<char, CStyle> str;
	Dna5String tmp;
	
	TAnnoTreeIter dfsIt = begin(store, AnnotationTree<>());
	while (!atEnd(dfsIt))
	{
		TAnnotation &anno = getAnnotation(dfsIt);
		if (anno.typeId == subExonType)
		{
			unsigned mRNAId = value(nodeUp(nodeUp(dfsIt)));
			
			if (length(orderings) <= mRNAId)
				resize(orderings, mRNAId + 1);
			
			TAnnotation &anno = getAnnotation(dfsIt);
			
			unsigned i = 0;
			if (anno.beginPos < anno.endPos)
			{
				for (; i < length(orderings[mRNAId]); ++i)
				{
					TAnnotation &anno2 = store.annotationStore[orderings[mRNAId][i]];
					if (anno.beginPos < anno2.beginPos) break;
				}
			} 
			else 
			{
				for (; i < length(orderings[mRNAId]); ++i)
				{
					TAnnotation &anno2 = store.annotationStore[orderings[mRNAId][i]];
					if (anno.beginPos > anno2.beginPos) break;
				}
			}
			
			insertValue(orderings[mRNAId], i, value(dfsIt));
			goNextRight(dfsIt);
		} 
		else
			goNext(dfsIt);
	}
}

template <typename TStream, typename TContigOrderings, typename TFragmentStore>
void writeTranscripts(TStream &target, TContigOrderings &orderings, TFragmentStore &store)
{
	CharString id;
	CharString tmp, transcript;
	typedef typename TFragmentStore::TAnnotationStore TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type TAnnotation;

	for (unsigned i = 0; i < length(orderings); ++i)
		if (!empty(orderings[i]))
		{
			std::stringstream sstream; 
			_streamWrite(sstream, "Locus_");
			_streamWrite(sstream, store.annotationNameStore[store.annotationStore[i].parentId]);
			_streamWrite(sstream, "_Transcript_");
			_streamWrite(sstream, store.annotationNameStore[i]);
			_streamWrite(sstream, "_Confidence_99.99");
			id = sstream.str();

			clear(transcript);
			for (unsigned j = 0; j < length(orderings[i]); ++j)
			{
				TAnnotation &anno = store.annotationStore[orderings[i][j]];
				if (_max(anno.beginPos, anno.endPos) < length(store.contigStore[anno.contigId].seq))
				{
					if (anno.beginPos < anno.endPos)
						append(transcript, infix(store.contigStore[anno.contigId].seq, anno.beginPos, anno.endPos));
					else
					{
						tmp = infix(store.contigStore[anno.contigId].seq, anno.endPos, anno.beginPos);
						reverseComplementInPlace(tmp);
						toLowerInPlace(tmp);
						append(transcript, tmp);
					}
				} 
				else 
				{
					std::cerr << "Genomic sequence missing for " << store.contigNameStore[anno.contigId];
					std::cerr << " from " << anno.beginPos << " to " << anno.endPos << std::endl;
				}

			}
			write(target, transcript, id, Fasta());			
		}
}

template <typename TStream, typename TOrderings, typename TFragmentStore, typename TNodeIds>
void writeOrderings(TStream &target, TOrderings &orderings, TFragmentStore &store, TNodeIds &nodeIds)
{
	typedef typename TFragmentStore::TAnnotationStore TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type TAnnotation;

#ifdef COMPACT_NODE_IDS
		std::sort(begin(nodeIds, Standard()), end(nodeIds, Standard()));	
#endif

	CharString nodeId;
	String<char, CStyle> idStr;
	for (unsigned i = 0; i < length(orderings); ++i)
		if (!empty(orderings[i]))
		{
			_streamWrite(target, ">Locus_");
			_streamWrite(target, store.annotationNameStore[store.annotationStore[i].parentId]);
			_streamWrite(target, "_Transcript_");
			_streamWrite(target, store.annotationNameStore[i]);
			_streamWrite(target, "_Confidence_99.99\n");
			int endPos = 0;
			int lastId = -1;
			__int64 lastPos = -1;
			bool silent = false;
			for (unsigned j = 0; j < length(orderings[i]); ++j)
			{
				TAnnotation &anno = store.annotationStore[orderings[i][j]];
				if (j > 0)
					_streamWrite(target, "-(0)->");
				if (anno.beginPos > anno.endPos)
					_streamPut(target, '-');
				bool hasNodeId = annotationGetValueByKey(store, anno, LARGE_NODE_ID, idStr);
				SEQAN_ASSERT_EQ(hasNodeId, true);
				int nodeId = atoi(idStr);
#ifdef COMPACT_NODE_IDS
				typename Iterator<TNodeIds, Standard>::Type it = std::lower_bound(begin(nodeIds, Standard()), end(nodeIds, Standard()), nodeId);
				nodeId = it - begin(nodeIds, Standard()) + 1;
				std::stringstream tmp;
				tmp << nodeId;
				annotationAssignValueByKey(store, anno, "NodeId", tmp.str());
#endif
				
				if (absInt(nodeId) <= lastId && !silent)
				{
					std::cerr << "nodeId is not increasing in "
						<< "Locus_" << store.annotationNameStore[store.annotationStore[i].parentId] 
						<< "_Transcript_" << store.annotationNameStore[i]
						<< "_Confidence_99.99" << std::endl;
					for (unsigned k = 0; k < length(orderings[i]); ++k)
						std::cerr << '\t' << orderings[i][k];
					std::cerr << std::endl;
					silent = true;
				}
				
				if (((anno.beginPos < anno.endPos && (__int64)anno.beginPos <= lastPos) ||
					(anno.beginPos > anno.endPos && (__int64)anno.beginPos >= lastPos)) && !silent && lastId != -1)
				{
					std::cerr << "beginPos is not monotonic in "
						<< "Locus_" << store.annotationNameStore[store.annotationStore[i].parentId] 
						<< "_Transcript_" << store.annotationNameStore[i]
						<< "_Confidence_99.99" << std::endl;
					for (unsigned k = 0; k < length(orderings[i]); ++k)
						std::cerr << '\t' << store.annotationStore[orderings[i][k]].beginPos << '-' << store.annotationStore[orderings[i][k]].endPos;
					std::cerr << std::endl;
					silent = true;
				}
				_streamPutInt(target, nodeId);
				_streamPut(target, ':');
				endPos += absInt((__int64)anno.endPos - (__int64)anno.beginPos);
				_streamPutInt(target, endPos);
				lastId = absInt(nodeId);
				lastPos = anno.beginPos;
			}
			_streamPut(target, '\n');
		}
}

int main(int argc, char const * argv[])
{
	typedef	FragmentStore<> TFragmentStore;
	typedef String<int> TContigExonBounds;
	typedef String<TContigExonBounds> TExonBounds;
	typedef Id<TFragmentStore>::Type TId;

	//////////////////////////////////////////////////////////////////////////////
	// Define options
	CommandLineParser	parser;
	addTitleLine(parser, "*****************************************");
	addTitleLine(parser, "***        Transcript Splicer         ***");
	addTitleLine(parser, "*** (c) Copyright 2010 by David Weese ***");
	addTitleLine(parser, "*****************************************");
	addUsageLine(parser, "[OPTION]... <genome.fa> <annotation.gff> <output-prefix>");
	addUsageLine(parser, "[OPTION]... <genome.fa> <annotation.gtf> <output-prefix>");
	addUsageLine(parser, "[OPTION]... <genome.fa> <knownGene.txt> <knownIsoforms.txt> <output-prefix>");
	requiredArguments(parser, 3);
	
	if (argc == 1)
	{
		shortHelp(parser, std::cerr);	// print short help and exit
		return 0;
	}
	
	if (!parse(parser, argc, argv, std::cerr)) return 0;

	double t0 = sysTime();
	TFragmentStore store;

	CharString prefix = getArgumentValue(parser, argumentCount(parser) - 1);
	if (empty(prefix))
		prefix = ".";
	CharString orderingFileName(prefix);
	CharString transcriptsFileName(prefix);
	CharString refinedFileName(prefix);
	append(orderingFileName, "/ordering.fa");
	append(transcriptsFileName, "/transcripts.fa");
	append(refinedFileName, "/refined.gff");
		
	
	//////////////////////////////////////////////////////////////////////////////
	// Step 1: Reading the genome
	loadContigs(store, toCString(getArgumentValue(parser, 0)), true);
	double t1 = sysTime();		std::cout << "Reading the genome took "<< t1-t0 << " seconds." << std::endl;

	//////////////////////////////////////////////////////////////////////////////
	// Step 2: Reading the annotation
	if (argumentCount(parser) == 3) // GTF/GFF?
	{
		std::ifstream annotationFile(toCString(getArgumentValue(parser, 1)), ::std::ios_base::in | ::std::ios_base::binary);
		if (!annotationFile.is_open())
		{
			std::cerr << "Could not open " << getArgumentValue(parser, 1) << std::endl;
			return 1;
		}
		read(annotationFile, store, GFF());
	} 
	else 
	{
		std::ifstream knownGenes(toCString(getArgumentValue(parser, 1)), ::std::ios_base::in | ::std::ios_base::binary);
		std::ifstream knownIsoforms(toCString(getArgumentValue(parser, 2)), ::std::ios_base::in | ::std::ios_base::binary);
		if (!knownGenes.is_open())
		{
			std::cerr << "Could not open " << getArgumentValue(parser, 1) << std::endl;
			return 1;
		}
		if (!knownIsoforms.is_open())
		{
			std::cerr << "Could not open " << getArgumentValue(parser, 2) << std::endl;
			return 1;
		}
		read(knownGenes, store, UCSC());
		read(knownIsoforms, store, UCSC());
	}
	double t2 = sysTime();		std::cout << "Reading the annotation took "<< t2-t1 << " seconds." << std::endl;	

	//////////////////////////////////////////////////////////////////////////////
	// Step 3: Get exon boundaries and refine into subexons
	TExonBounds exonBounds;
	String<int> nodeIds;
	getExonBoundaries(exonBounds, store);
								std::cout << "Collected all exon boundaries." << std::endl;
	refineExonBoundaries(exonBounds, store, nodeIds);
	double t3 = sysTime();		std::cout << "Refining exons took "<< t3-t2 << " seconds." << std::endl;	

	//////////////////////////////////////////////////////////////////////////////
	// Step 4: Calculate orderings
	StringSet<String<int> > orderings;
	calculateOrderings(orderings, store);
	double t4 = sysTime();		std::cout << "Calculating orderings took "<< t4-t3 << " seconds." << std::endl;

	//////////////////////////////////////////////////////////////////////////////
	// Step 5: Write ordering file
	std::ofstream orderingsFile(toCString(orderingFileName), ::std::ios_base::out | ::std::ios_base::binary);
	writeOrderings(orderingsFile, orderings, store, nodeIds);
	orderingsFile.close();
	double t5 = sysTime();		std::cout << "Writing contig orderings took "<< t5-t4 << " seconds." << std::endl;

	//////////////////////////////////////////////////////////////////////////////
	// Step 6: Write transcripts
	std::ofstream transcriptsFile(toCString(transcriptsFileName), ::std::ios_base::out | ::std::ios_base::binary);
	writeTranscripts(transcriptsFile, orderings, store);
	transcriptsFile.close();
	double t6 = sysTime();		std::cout << "Writing spliced transcipts took "<< t6-t5 << " seconds." << std::endl;

	//////////////////////////////////////////////////////////////////////////////
	// Step 7: Write the annotation enhanced by subexons
	std::ofstream file2(toCString(refinedFileName), ::std::ios_base::out | ::std::ios_base::binary);
	write(file2, store, GFF());
	file2.close();
	double t7 = sysTime();		std::cout << "Writing the enhanced annotation took "<< t7-t6 << " seconds." << std::endl;

	return 0;
}
