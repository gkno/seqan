///A tutorial about global alignments.
#include <fstream>
#include <iostream>
#include <seqan/store.h>
#include <seqan/store/store_io_gff.h>
#include <seqan/misc/misc_cmdparser.h>

using namespace seqan;

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
			unsigned contigNum = anno.contigId * 2;
			if (anno.beginPos > anno.endPos) ++contigNum;
			
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

template <typename TExonBoundString, typename TSubExonIds, typename TFragmentStore>
void refineExonBoundaries(TExonBoundString &exonBounds, TSubExonIds &subExonIds, TFragmentStore &store)
{
	typedef typename Value<TExonBoundString>::Type TContigExonBounds;
	typedef typename Iterator<TContigExonBounds, Standard>::Type TContigExonBoundIterator;
	typedef typename TFragmentStore::TAnnotationStore TAnnotationStore;
	typedef typename  Iterator<TFragmentStore, AnnotationTree<> >::Type TAnnoTreeIter;
	typedef typename  Value<TAnnotationStore>::Type TAnnotation;
	typedef typename TFragmentStore::TContigPos TContigPos;
	typedef typename  Id<TAnnotation>::Type TId;

	unsigned nextSubExonId = 1;
	resize(subExonIds, length(exonBounds));
	for (unsigned i = 0; i < length(exonBounds); ++i)
		fill(subExonIds[i], length(exonBounds[i]), TAnnotation::INVALID_ID);
	
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
			
			unsigned contigNum = anno.contigId * 2;
			if (anno.beginPos > anno.endPos) ++contigNum;

			TContigExonBoundIterator boundsBegin = begin(exonBounds[contigNum], Standard());
			TContigExonBoundIterator boundsEnd = end(exonBounds[contigNum], Standard());
			
			TContigPos beginPos = anno.beginPos;
			TContigPos endPos = anno.endPos;
			if (beginPos > endPos)
			{
				beginPos = endPos;
				endPos = anno.beginPos;
			}

			TContigExonBoundIterator first = std::upper_bound(boundsBegin, boundsEnd, beginPos);
			TContigExonBoundIterator last = std::lower_bound(boundsBegin, boundsEnd, endPos);
			
			TContigPos subExonBeginPos = beginPos;
			for (TContigExonBoundIterator it = first; it <= last; ++it)
			{
				TContigPos subExonEndPos = *it;
				TAnnoTreeIter childIt = createRightChild(dfsIt);
				clearValues(childIt);
				if ((contigNum & 1) == 0)
				{	// positive strand
					getAnnotation(childIt).beginPos = subExonBeginPos;
					getAnnotation(childIt).endPos = subExonEndPos;
				} else
				{	// negative strand
					getAnnotation(childIt).endPos = subExonBeginPos;
					getAnnotation(childIt).beginPos = subExonEndPos;
				}
				if (subExonIds[contigNum][it - boundsBegin] == TAnnotation::INVALID_ID)
					subExonIds[contigNum][it - boundsBegin] = nextSubExonId++;
				setType(childIt, "subexon");
				std::stringstream tmp;
				tmp << subExonIds[contigNum][it - boundsBegin];
				assignValueByKey(childIt, "NodeId", tmp.str());
				subExonBeginPos = subExonEndPos;
			}
			goNextRight(dfsIt);
		}
		else
			goNext(dfsIt);
	}
}

template <typename TTranscriptSet, typename TContigOrderings, typename TFragmentStore>
void spliceTranscripts(TTranscriptSet &transcripts, TContigOrderings &orderings, TFragmentStore &store)
{
	typedef typename TFragmentStore::TAnnotationStore TAnnotationStore;
	typedef typename Iterator<TFragmentStore, AnnotationTree<> >::Type TAnnoTreeIter;
	typedef typename Value<TAnnotationStore>::Type TAnnotation;
	
	unsigned subExonId = 0;
	_storeAppendType(store, subExonId, "subexon");
	String<char, CStyle> str;
	Dna5String tmp;
	
	TAnnoTreeIter dfsIt = begin(store, AnnotationTree<>());
	while (!atEnd(dfsIt))
	{
		TAnnotation &anno = getAnnotation(dfsIt);
		if (anno.typeId == subExonId)
		{
			unsigned mRNAId = value(nodeUp(nodeUp(dfsIt)));
			
			if (length(transcripts) <= mRNAId)
			{
				resize(transcripts, mRNAId + 1);
				resize(orderings, mRNAId + 1);
			}
			
			TAnnotation &anno = getAnnotation(dfsIt);
			
			if (_max(anno.beginPos, anno.endPos) < length(store.contigStore[anno.contigId].seq))
			{
				unsigned i = 0;
				__int64 len = 0;
				if (anno.beginPos < anno.endPos)
				{
					for (; i < length(orderings[mRNAId]); ++i)
					{
						TAnnotation &anno2 = store.annotationStore[orderings[mRNAId][i]];
						if (anno.beginPos < anno2.beginPos) break;
						len += absInt((__int64)anno2.endPos - (__int64)anno2.beginPos);
					}
					replace(transcripts[mRNAId], len, len, infix(store.contigStore[anno.contigId].seq, anno.beginPos, anno.endPos));
				} else {
					for (; i < length(orderings[mRNAId]); ++i)
					{
						TAnnotation &anno2 = store.annotationStore[orderings[mRNAId][i]];
						if (anno.beginPos > anno2.beginPos) break;
						len += absInt((__int64)anno2.endPos - (__int64)anno2.beginPos);
					}
					tmp = infix(store.contigStore[anno.contigId].seq, anno.endPos, anno.beginPos);
					reverseComplementInPlace(tmp);
					replace(transcripts[mRNAId], len, len, tmp);
				}
				insertValue(orderings[mRNAId], i, value(dfsIt));
			}
//			else
//				std::cerr << "Error: Exon (" << getUniqueName(nodeUp(dfsIt)) << ") is annotated outside the contig sequence (" << store.contigNameStore[anno.contigId] << ")." << std::endl;
			goNextRight(dfsIt);
		} 
		else
			goNext(dfsIt);
	}
}

template <typename TStream, typename TTranscriptSet, typename TFragmentStore>
void writeTranscripts(TStream &target, TTranscriptSet &transcripts, TFragmentStore &store)
{
	CharString id;
	for (unsigned i = 0; i < length(transcripts); ++i)
		if (!empty(transcripts[i]))
		{
			std::stringstream sstream; 
			_streamWrite(sstream, "Locus_");
			_streamWrite(sstream, store.annotationNameStore[store.annotationStore[i].parentId]);
			_streamWrite(sstream, "_Transcript_");
			_streamWrite(sstream, store.annotationNameStore[i]);
			_streamWrite(sstream, "_Confidence_99.99");
			id = sstream.str();
			write(target, transcripts[i], id, Fasta());
		}
}

template <typename TStream, typename TOrderings, typename TFragmentStore>
void writeOrderings(TStream &target, TOrderings &orderings, TFragmentStore &store)
{
	typedef typename TFragmentStore::TAnnotationStore TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type TAnnotation;

	CharString nodeId;
	for (unsigned i = 0; i < length(orderings); ++i)
		if (!empty(orderings[i]))
		{
			_streamWrite(target, ">Locus_");
			_streamWrite(target, store.annotationNameStore[store.annotationStore[i].parentId]);
			_streamWrite(target, "_Transcript_");
			_streamWrite(target, store.annotationNameStore[i]);
			_streamWrite(target, "_Confidence_99.99\n");
			int endPos = 0;
			for (unsigned j = 0; j < length(orderings[i]); ++j)
			{
				TAnnotation &anno = store.annotationStore[orderings[i][j]];
				if (j > 0)
					_streamWrite(target, "-(0)->");
				if (anno.beginPos > anno.endPos)
					_streamPut(target, '-');
				_streamWrite(target, annotationGetValueByKey(store, anno, "NodeId"));
				_streamPut(target, ':');
				endPos += absInt((__int64)anno.endPos - (__int64)anno.beginPos);
				_streamPutInt(target, endPos);
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
	addUsageLine(parser, "[OPTION]... <genome> <annotation> <output-prefix>");
	requiredArguments(parser, 3);
	
	if (argc == 1)
	{
		shortHelp(parser, std::cerr);	// print short help and exit
		return 0;
	}
	
	if (!parse(parser, argc, argv, std::cerr)) return 0;

	double t0 = sysTime();
	TFragmentStore store;
	
	//////////////////////////////////////////////////////////////////////////////
	// Step 1: Reading the genome
	loadContigs(store, toCString(getArgumentValue(parser, 0)), true);
	double t1 = sysTime();		std::cout << "Reading the genome took "<< t1-t0 << " seconds." << std::endl;

	//////////////////////////////////////////////////////////////////////////////
	// Step 2: Reading the annotation
	std::ifstream annotationFile(toCString(getArgumentValue(parser, 1)), ::std::ios_base::in | ::std::ios_base::binary);
	read(annotationFile, store, GFF());
	double t2 = sysTime();		std::cout << "Reading the annotation took "<< t2-t1 << " seconds." << std::endl;	

	//////////////////////////////////////////////////////////////////////////////
	// Step 3: Get exon boundaries and refine into subexons
	TExonBounds exonBounds;
	getExonBoundaries(exonBounds, store);
	String<String<TId> > subExonIds;
	refineExonBoundaries(exonBounds, subExonIds, store);
	double t3 = sysTime();		std::cout << "Refining exons took "<< t3-t2 << " seconds." << std::endl;	

	//////////////////////////////////////////////////////////////////////////////
	// Step 4: Transcribe
	StringSet<Dna5String> transcripts;
	StringSet<String<int> > orderings;
	spliceTranscripts(transcripts, orderings, store);
	double t4 = sysTime();		std::cout << "Splicing transcipts took "<< t4-t3 << " seconds." << std::endl;

	//////////////////////////////////////////////////////////////////////////////
	// Step 5: Write transcripts
	CharString prefix = getArgumentValue(parser, 2);
	if (empty(prefix))
		prefix = ".";
	CharString transcriptsFileName(prefix);
	CharString orderingFileName(prefix);
	CharString refinedFileName(prefix);
	append(transcriptsFileName, "/transcripts.fa");
	append(orderingFileName, "/ordering.fa");
	append(refinedFileName, "/refined.gff");
	std::ofstream transcriptsFile(toCString(transcriptsFileName), ::std::ios_base::out | ::std::ios_base::binary);
	writeTranscripts(transcriptsFile, transcripts, store);
	transcriptsFile.close();
	double t5 = sysTime();		std::cout << "Writing spliced transcipts took "<< t5-t4 << " seconds." << std::endl;

	//////////////////////////////////////////////////////////////////////////////
	// Step 6: Write ordering file
	std::ofstream orderingsFile(toCString(orderingFileName), ::std::ios_base::out | ::std::ios_base::binary);
	writeOrderings(orderingsFile, orderings, store);
	orderingsFile.close();
	double t6 = sysTime();		std::cout << "Writing contig orderings took "<< t6-t5 << " seconds." << std::endl;

	//////////////////////////////////////////////////////////////////////////////
	// Step 7: Write the annotation enhanced by subexons
	std::ofstream file2(toCString(refinedFileName), ::std::ios_base::out | ::std::ios_base::binary);
	write(file2, store, GFF());
	file2.close();
	double t7 = sysTime();		std::cout << "Writing the enhanced annotation took "<< t7-t6 << " seconds." << std::endl;

	return 0;
}
