 /*==========================================================================
             RazerS - Fast Read Mapping with Controlled Loss Rate
                   http://www.seqan.de/projects/razers.html

 ============================================================================
  Copyright (C) 2008 by David Weese

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

#ifndef SEQAN_HEADER_RAZERS_MATEPAIRS_H
#define SEQAN_HEADER_RAZERS_MATEPAIRS_H

namespace SEQAN_NAMESPACE_MAIN
{

// We expect mate-pairs to be stored together in one read string.
// Pair i has mates at positions 2*i and 2*i+1 in the read string.



template <typename TValue, typename TSpec = Alloc<> >
class Dequeue
{
public:
	typedef String<TValue, TSpec>						TString;
	typedef typename Iterator<TString, Standard>::Type	TIter;

	String<TValue, TSpec> data_string;

	TIter data_begin;	// string begin
	TIter data_end;		// string end

	TIter data_front;	// front fifo character
	TIter data_back;	// back fifo character
	bool data_empty;	// fifo is empty

//____________________________________________________________________________

public:
	inline Dequeue()
	{
		clear(*this);
	}
};


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline bool
empty(Dequeue<TValue, TSpec> const &me)
{
	return me.empty;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline void
clear(Dequeue<TValue, TSpec> &me)
{
	clear(me.data_string);
	me.data_begin = begin(me.data_string, Standard());
	me.data_end = end(me.data_string, Standard());

	me.data_front = me.data_back = me.data_begin;
	me.data_empty = true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline TValue &
front(Dequeue<TValue, TSpec> &me)
{
	return *me.data_front;
}

template <typename TValue, typename TSpec>
inline TValue const &
front(Dequeue<TValue, TSpec> const &me)
{
	return *me.data_front;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline TValue &
back(Dequeue<TValue, TSpec> &me)
{
	return *me.data_back;
}

template <typename TValue, typename TSpec>
inline TValue const &
back(Dequeue<TValue, TSpec> const &me)
{
	return *me.data_back;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline bool
popFront(Dequeue<TValue, TSpec> &me)
{
	if (me.data_empty) return false;

	if (me.data_front == me.data_back)
		me.data_empty = true;

	if (++me.data_front == me.data_end)
		me.data_front = me.data_begin;

	return true;
}

template <typename TValue, typename TSpec>
inline bool
popBack(Dequeue<TValue, TSpec> &me)
{
	if (me.data_empty) return false;

	if (me.data_front == me.data_back)
		me.data_empty = true;

	if (me.data_back == me.data_begin)
		me.data_back = me.data_end;
	--me.data_back;

	return true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline void
pushFront(Dequeue<TValue, TSpec> &me, TValue const & _value)
{
	if (me.data_empty) 
		me.data_empty = false;
	else 
	{
		if (me.data_front == me.data_begin)
			me.data_front = me.data_end;
		--me.data_front;
	}

	assign(*me.data_front, _value);
}

template <typename TValue, typename TSpec>
inline void
pushBack(Dequeue<TValue, TSpec> &me, TValue const & _value)
{
	if (me.data_empty) 
		me.data_empty = false;
	else 
	{
		if (++me.data_back == me.data_end)
			me.data_back = me.data_begin;
	}

	assign(*me.data_back, _value);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline typename Size<Dequeue<TValue, TSpec> >::Type
length(Dequeue<TValue, TSpec> const &me)
{
	if (empty(me)) return 0;

	if (me.data_front < me.data_back)
		return (me.data_back - me.data_front) + 1;
	else
		return (me.data_end - me.data_begin) - (me.data_front - me.data_back) + 1;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TSize_, typename TExpand>
inline typename Size<Dequeue<TValue, TSpec> >::Type
reserve(Dequeue<TValue, TSpec> &me, TSize_ new_capacity, Tag<TExpand> const tag)
{
	typedef typename Size<Dequeue<TValue, TSpec> >::Type TSize;
	
	TSize len = length(me);
	if (len < new_capacity && capacity(me.data_string) != new_capacity)
	{
		TSize pos_front = me.data_front - me.data_begin;
		TSize pos_back  = me.data_back  - me.data_begin;
		TSize new_freeSpace = new_capacity - len;

		if (pos_front < pos_back)
		{
			TSize freeSpace = capacity(me.data_string) - len;
			if (new_freeSpace > freeSpace)
				reserve(me.data_string, new_capacity, tag);
			else
			{
				freeSpace -= new_freeSpace;	// reduce the free space by <freeSpace>
				if (freeSpace <= pos_front)
					_clearSpace(me, freeSpace, 0, pos_front + 1, tag);
				else
				{
					freeSpace -= pos_front;
					_clearSpace(me, 0, 0, pos_front + 1, tag);
					_clearSpace(me, freeSpace, pos_back + 1, capacity(me.data_string), tag);
				}
			}
		} else
			_clearSpace(me, new_freeSpace, pos_back, pos_front + 1, tag);

		me.data_begin = begin(me.data_string, Standard());
		me.data_end = end(me.data_string, Standard());
		me.data_front = me.data_begin + pos_front;
		me.data_back = me.data_begin + pos_back;
	}
}

//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta sequences
template <typename TReadSet, typename TNameSet, typename TRazerSOptions>
bool loadReads(
	TReadSet &reads,			// resulting mate sequences
	TNameSet &fastaIDs,			// resulting mate ids
	const char *fileNameL,		// left mates file
	const char *fileNameR,		// right mates file
	TRazerSOptions &options)
{
	bool countN = !(options.matchN || options.outputFormat == 1);

	MultiFasta leftMates;
	MultiFasta rightMates;

	if (!open(leftMates.concat, fileNameL, OPEN_RDONLY)) return false;
	if (!open(rightMates.concat, fileNameR, OPEN_RDONLY)) return false;

	AutoSeqFormat formatL;
	guessFormat(leftMates.concat, formatL);
	split(leftMates, formatL);

	AutoSeqFormat formatR;
	guessFormat(rightMates.concat, formatR);
	split(rightMates, formatR);

	unsigned seqCount = length(leftMates);
	if (seqCount != length(rightMates))
	if (options._debugLevel > 1) 
	{
		::std::cerr << "Numbers of mates differ: " << seqCount << "(left) != " << length(rightMates) << "(right).\n";
		return false;
	}


#ifndef RAZERS_CONCATREADS
	resize(reads, 2*seqCount, Exact());
#endif
	if (options.readNaming == 0)
		resize(fastaIDs, 2*seqCount, Exact());
	
	Dna5String seq[2];
	unsigned kickoutcount = 0;
	for (unsigned i = 0; i < seqCount; ++i) 
	{
		if (options.readNaming == 0)
		{
			assignSeqId(fastaIDs[2*i], leftMates[i], formatL);		// read Fasta id
			assignSeqId(fastaIDs[2*i+1], rightMates[i], formatR);	// read Fasta id
		}
		assignSeq(seq[0], leftMates[i], formatL);					// read Read sequence
		assignSeq(seq[1], rightMates[i], formatR);					// read Read sequence

		if (countN)
		{
			for (int j = 0; j < 2; ++j)
			{
				int allowedNs = (int)(options.errorRate * length(seq[j]));
				for (unsigned k = 0; k < length(seq[j]); ++k)
					if (getValue(seq[j], k) == 'N')
						if (allowedNs-- == 0)
						{
							clear(seq[0]);
							clear(seq[1]);
							++kickoutcount;
							break;
						}
			}
		}
#ifdef RAZERS_CONCATREADS
		appendValue(reads, seq[0], Generous());
		appendValue(reads, seq[1], Generous());
#else
		assign(reads[2*i], seq[0], Exact());
		assign(reads[2*i+1], seq[1], Exact());
#endif
	}
#ifdef RAZERS_CONCATREADS
	reserve(reads.concat, length(reads.concat), Exact());
#endif

	if (options._debugLevel > 1 && kickoutcount > 0) 
		::std::cerr << "Ignoring " << kickoutcount << " low quality mate-pairs.\n";
	return (seqCount > 0);
}

#ifndef RAZERS_PARALLEL
//////////////////////////////////////////////////////////////////////////////
// Find read matches in one genome sequence
template <
	typename TMatches, 
	typename TGenome,
	typename TReadIndex, 
	typename TSwiftSpec, 
	typename TVerifier,
	typename TSpec >
void mapMatePairReads(
	TMatches &matches,				// resulting matches
	TGenome &genome,				// genome ...
	unsigned gseqNo,				// ... and its sequence number
	Pattern<TReadIndex, Swift<TSwiftSpec> > &swiftPatternL,
	Pattern<TReadIndex, Swift<TSwiftSpec> > &swiftPatternR,
	TVerifier &forwardPatternsL,
	TVerifier &forwardPatternsR,
	char orientation,				// q-gram index of reads
	RazerSOptions<TSpec> &options)
{
	typedef typename Fibre<TReadIndex, Fibre_Text>::Type	TReadSet;
	typedef typename Size<TGenome>::Type					TSize;
	typedef typename Position<TGenome>::Type				TGPos;
	typedef typename Value<TMatches>::Type					TMatch;
	typedef typename Infix<TGenome>::Type					TGenomeInf;

	// FILTRATION
	typedef Finder<TGenome, Swift<TSwiftSpec> >				TSwiftFinderL;
	typedef Finder<TGenomeInf, Swift<TSwiftSpec> >			TSwiftFinderR;
	typedef Pattern<TReadIndex, Swift<TSwiftSpec> >			TSwiftPattern;

	// iterate all genomic sequences
	if (options._debugLevel >= 1)
	{
		::std::cerr << ::std::endl << "Process genome seq #" << gseqNo;
		if (orientation == 'F')
			::std::cerr << "[fwd]";
		else
			::std::cerr << "[rev]";
	}

	// Mate-Pair parameters
	TGPos firstPosR = options.libraryLength;
	if (firstPosR >= options.libraryError)
		firstPosR -= options.libraryError;
	else
		firstPosR = 0;

	TSize distanceCutOff = options.libraryLength + options.libraryError;
	

	// exit if contig is shorter than library size
	if (length(genome) < firstPosR)
		return;

	TGenomeInf genomeInf = infix(genome, firstPosR, length(genome));
	TReadSet &readSet = host(host(swiftPatternL));
	TSwiftFinderL swiftFinderL(genome, options.repeatLength, 1);
	TSwiftFinderR swiftFinderR(genomeInf, options.repeatLength, 1);

	typedef Pair<unsigned, Pair<TGPos> >	TDequeueValue;
	Dequeue<TDequeueValue>					fifo;		// stores left-mate potential matches
	String<unsigned>						potMatches;	// counts pot. matches in fifo

	resize(potMatches, length(host(swiftPatternL)), Exact());

	TMatch m = { 0, 0, 0, 0, 0, 0 };	// to supress uninitialized warnings
	// iterate all verification regions returned by SWIFT
	while (find(swiftFinderR, swiftPatternR, options.errorRate, options._debugLevel)) 
	{
		unsigned rseqNo = (*swiftFinderR.curHit).ndlSeqNo;

		// remove left mates pot. matches too distant to have a right mate
		while (front(fifo).i2.i1 + distanceCutOff < beginPosition(swiftFinderR))
			pop(fifo);

	}
}
#endif


//////////////////////////////////////////////////////////////////////////////
// Find read matches in many genome sequences (import from Fasta)
template <
	typename TMatches, 
	typename TReadSet_, 
	typename TSpec, 
	typename TShape,
	typename TSwiftSpec >
int mapMatePairReads(
	TMatches &				matches,
	StringSet<CharString> &	genomeFileNameList,
	StringSet<CharString> &	genomeNames,	// genome names, taken from the Fasta file
	::std::map<unsigned,::std::pair< ::std::string,unsigned> > & gnoToFileMap,
	TReadSet_ const &		readSet,
	RazerSOptions<TSpec> &	options,
	TShape const &			shape,
	Swift<TSwiftSpec> const)
{
	typedef typename Value<TReadSet_>::Type				TRead;
	typedef StringSet<TRead, Dependent<> >				TReadSet;
	typedef Index<TReadSet, Index_QGram<TShape> >		TIndex;			// q-gram index
	typedef Pattern<TIndex, Swift<TSwiftSpec> >			TSwiftPattern;	// filter
	typedef Pattern<TRead, MyersUkkonen>				TMyersPattern;	// verifier

	// split mate-pairs over two indices
	TReadSet readSetL, readSetR;
	unsigned readCount = length(readSet) / 2;
	reserve(readSetL, readCount, Exact());
	reserve(readSetR, readCount, Exact());

	for (unsigned i = 0; i < readCount; ++i)
	{
		appendValue(readSetL, readSet[2*i]);
		appendValue(readSetR, readSet[2*i+1]);
	}

	// configure q-gram index
	TIndex swiftIndexL(readSetL, shape);
	TIndex swiftIndexR(readSetR, shape);
	cargo(swiftIndexL).abundanceCut = options.abundanceCut;
	cargo(swiftIndexR).abundanceCut = options.abundanceCut;
	cargo(swiftIndexL)._debugLevel = options._debugLevel;
	cargo(swiftIndexR)._debugLevel = options._debugLevel;

	// configure Swift
	TSwiftPattern swiftPatternL(swiftIndexL);
	TSwiftPattern swiftPatternR(swiftIndexR);
	swiftPatternL.params.minThreshold = options.threshold;
	swiftPatternR.params.minThreshold = options.threshold;
	swiftPatternL.params.tabooLength = options.tabooLength;
	swiftPatternR.params.tabooLength = options.tabooLength;

	// init edit distance verifiers
	String<TMyersPattern> forwardPatternsL;
	String<TMyersPattern> forwardPatternsR;
	options.compMask[4] = (options.matchN)? 15: 0;
	if (!options.hammingOnly)
	{
		resize(forwardPatternsL, readCount, Exact());
		resize(forwardPatternsR, readCount, Exact());
		for(unsigned i = 0; i < readCount; ++i)
		{
			setHost(forwardPatternsL[i], readSetL[i]);
			setHost(forwardPatternsR[i], readSetR[i]);
			_patternMatchNOfPattern(forwardPatternsL[i], options.matchN);
			_patternMatchNOfPattern(forwardPatternsR[i], options.matchN);
			_patternMatchNOfFinder(forwardPatternsL[i], options.matchN);
			_patternMatchNOfFinder(forwardPatternsR[i], options.matchN);
		}
	}

#ifdef RAZERS_MASK_READS
	// init read mask
	clear(options.readMask);
	fill(options.readMask, (readCount + options.WORD_SIZE - 1) / options.WORD_SIZE, (unsigned long)-1);
#endif

	// clear stats
	options.FP = 0;
	options.TP = 0;
	options.timeMapReads = 0;
	options.timeDumpResults = 0;

	unsigned filecount = 0;
	unsigned numFiles = length(genomeFileNameList);
	unsigned gseqNo = 0;

	// open genome files, one by one	
	while (filecount < numFiles)
	{
		// open genome file	
		::std::ifstream file;
		file.open(toCString(genomeFileNameList[filecount]), ::std::ios_base::in | ::std::ios_base::binary);
		if (!file.is_open())
			return RAZERS_GENOME_FAILED;

		// remove the directory prefix of current genome file
		::std::string genomeFile(toCString(genomeFileNameList[filecount]));
		size_t lastPos = genomeFile.find_last_of('/') + 1;
		if (lastPos == genomeFile.npos) lastPos = genomeFile.find_last_of('\\') + 1;
		if (lastPos == genomeFile.npos) lastPos = 0;
		::std::string genomeName = genomeFile.substr(lastPos);
		

		CharString	id;
		Dna5String	genome;
		unsigned gseqNoWithinFile = 0;
		// iterate over genome sequences
		SEQAN_PROTIMESTART(find_time);
		for(; !_streamEOF(file); ++gseqNo)
		{
			if (options.genomeNaming == 0)
			{
				//readID(file, id, Fasta());			// read Fasta id
				readShortID(file, id, Fasta());			// read Fasta id up to first whitespace
				appendValue(genomeNames, id, Generous());
			}
			read(file, genome, Fasta());			// read Fasta sequence
			
			gnoToFileMap.insert(::std::make_pair<unsigned,::std::pair< ::std::string,unsigned> >(gseqNo,::std::make_pair< ::std::string,unsigned>(genomeName,gseqNoWithinFile)));
			
			if (options.forward)
				findReads(matches, genome, gseqNo, swiftPatternL, swiftPatternR, forwardPatternsL, forwardPatternsR, 'F', options);

			if (options.reverse)
			{
				reverseComplementInPlace(genome);
				findReads(matches, genome, gseqNo, swiftPatternL, swiftPatternR, forwardPatternsL, forwardPatternsR, 'R', options);
			}
			++gseqNoWithinFile;

		}
		options.timeMapReads += SEQAN_PROTIMEDIFF(find_time);
		file.close();
		++filecount;
	}

	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Finding reads took               \t" << options.timeMapReads << " seconds" << ::std::endl;
	if (options._debugLevel >= 2) {
		::std::cerr << ::std::endl;
		::std::cerr << "___FILTRATION_STATS____" << ::std::endl;
		::std::cerr << "Swift FP: " << options.FP << ::std::endl;
		::std::cerr << "Swift TP: " << options.TP << ::std::endl;
	}
	return 0;
}



}

#endif
