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

// We require mate-pairs to be stored together in one read string.
// Pair i has mates at positions 2*i and 2*i+1 in the read string.



template <typename TValue, typename TSpec = Alloc<> >
class Dequeue
{
public:
	typedef String<TValue, TSpec>						TString;
	typedef typename Iterator<TString, Standard>::Type	TIter;

	String<TValue, TSpec> data_string;

	TIter data_begin;	// string beginning
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
// Iterators
//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Iterator.param.T.type:Spec.Dequeue

template<typename TValue, typename TSpec>
struct Iterator<Dequeue<TValue, TSpec>, Standard> 
{
	typedef Iter<Dequeue<TValue, TSpec>, PositionIterator> Type;
};

template<typename TValue, typename TSpec>
struct Iterator<Dequeue<TValue, TSpec> const, Standard> 
{
	typedef Iter<Dequeue<TValue, TSpec> const, PositionIterator> Type;
};

template<typename TValue, typename TSpec>
struct Iterator<Dequeue<TValue, TSpec>, Rooted> 
{
	typedef Iter<Dequeue<TValue, TSpec>, PositionIterator> Type;
};

template<typename TValue, typename TSpec>
struct Iterator<Dequeue<TValue, TSpec> const, Rooted> 
{
	typedef Iter<Dequeue<TValue, TSpec> const, PositionIterator> Type;
};


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline bool
empty(Dequeue<TValue, TSpec> const &me)
{
	return me.data_empty;
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

template <typename TValue, typename TSpec, typename TPos>
inline TValue &
value(Dequeue<TValue, TSpec> &me, TPos pos)
{
	typedef typename Size<Dequeue<TValue, TSpec> >::Type TSize;
	TSize wrap = length(me) - (me.data_front - me.data_begin);
	
	if (pos < wrap)
		return value(me.data_front + pos);
	else
		return value(me.data_front + (pos - wrap));
}

template <typename TValue, typename TSpec, typename TPos>
inline TValue const &
value(Dequeue<TValue, TSpec> const &me, TPos pos)
{
	typedef typename Size<Dequeue<TValue, TSpec> >::Type TSize;
	TSize wrap = length(me) - (me.data_front - me.data_begin);
	
	if (pos < wrap)
		return value(me.data_front + pos);
	else
		return value(me.data_front + (pos - wrap));
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
	else
	{
		if (++me.data_front == me.data_end)
			me.data_front = me.data_begin;
	}

	return true;
}

template <typename TValue, typename TSpec>
inline bool
popBack(Dequeue<TValue, TSpec> &me)
{
	if (me.data_empty) return false;

	if (me.data_front == me.data_back)
		me.data_empty = true;
	else
	{
		if (me.data_back == me.data_begin)
			me.data_back = me.data_end;
		--me.data_back;
	}

	return true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline void
pushFront(Dequeue<TValue, TSpec> &me, TValue const & _value)
{
	typedef typename Dequeue<TValue, TSpec>::TIter TIter;

	if (me.data_empty) 
	{
		if (me.data_begin == me.data_end)
			reserve(me, computeGenerousCapacity(me, length(me.data_string) + 1), Generous());
		me.data_empty = false;
	}
	else 
	{
		TIter new_front = me.data_front;
		if (new_front == me.data_begin)
			new_front = me.data_end;
		--new_front;

		if (new_front == me.data_back)
		{
			reserve(me, computeGenerousCapacity(me, length(me.data_string) + 1), Generous());

			if (me.data_front == me.data_begin)
				me.data_front = me.data_end;
			--me.data_front;
		} else
			me.data_front = new_front;
	}
	assign(*me.data_front, _value);
}

template <typename TValue, typename TSpec>
inline void
pushBack(Dequeue<TValue, TSpec> &me, TValue const & _value)
{
	typedef typename Dequeue<TValue, TSpec>::TIter TIter;

	if (me.data_empty) 
	{
		if (me.data_begin == me.data_end)
			reserve(me, computeGenerousCapacity(me, length(me.data_string) + 1), Generous());
		me.data_empty = false;
	}
	else 
	{
		TIter new_back = me.data_back;
		if (++new_back == me.data_end)
			new_back = me.data_begin;

		if (new_back == me.data_front)
		{
			reserve(me, computeGenerousCapacity(me, length(me.data_string) + 1), Generous());
			// in this case reserve adds new space behind data_back
			++me.data_back;
		} else
			me.data_back = new_back;
	}
	assign(*me.data_back, _value);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline typename Size<Dequeue<TValue, TSpec> >::Type
length(Dequeue<TValue, TSpec> const &me)
{
	if (empty(me)) return 0;

	if (me.data_front <= me.data_back)
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
	::std::cout << "resize to "<<new_capacity<<::std::endl;
	TSize len = length(me);
	if (len < new_capacity && length(me.data_string) != new_capacity)
	{
		TSize pos_front = me.data_front - me.data_begin;
		TSize pos_back  = me.data_back  - me.data_begin;
		TSize new_freeSpace = new_capacity - len;

		if (pos_front <= pos_back)
		{
			// |empty|data|empty|
			// 0
			TSize freeSpace = length(me.data_string) - len;
			if (new_freeSpace > freeSpace)
				resize(me.data_string, new_capacity, tag);
			else
			{
				freeSpace -= new_freeSpace;	// reduce the free space by <freeSpace>
				if (pos_front >= freeSpace)
				{
					resizeSpace(me.data_string, pos_front - freeSpace, (TSize)0, pos_front, tag);
					pos_back -= freeSpace;
					pos_front -= freeSpace;
				}
				else
				{
					freeSpace -= pos_front;
					resizeSpace(me.data_string, length(me.data_string) - freeSpace, pos_back + 1, length(me.data_string), tag);
					resizeSpace(me.data_string, (TSize)0, (TSize)0, pos_front, tag);
					pos_back -= pos_front;
					pos_front = 0;
				}
			}
		}
		else
		{
			// |data|empty|data|
			// 0
			resizeSpace(me.data_string, new_freeSpace, pos_back + 1, pos_front, tag);
			pos_front += new_freeSpace;
		}

		me.data_begin = begin(me.data_string, Standard());
		me.data_end = end(me.data_string, Standard());
		me.data_front = me.data_begin + pos_front;
		me.data_back = me.data_begin + pos_back;
	}
	return length(me.data_string);
}

//////////////////////////////////////////////////////////////////////////////
// Definitions

typedef StringSet<TRead const, Dependent<> >	TMPReadSet;

template <typename TShape>
struct Cargo< Index<TMPReadSet const, TShape> > {
	typedef struct {
		double		abundanceCut;
		int			_debugLevel;
	} Type;
};

#ifdef RAZERS_PRUNE_QGRAM_INDEX

//////////////////////////////////////////////////////////////////////////////
// Repeat masker
template <typename TShape>
inline bool _qgramDisableBuckets(Index<TMPReadSet, Index_QGram<TShape> > &index) 
{
	typedef Index<TReadSet, Index_QGram<TShape>	>		TReadIndex;
	typedef typename Fibre<TReadIndex, QGram_Dir>::Type	TDir;
	typedef typename Iterator<TDir, Standard>::Type		TDirIterator;
	typedef typename Value<TDir>::Type					TSize;
	
	TDir &dir    = indexDir(index);
	bool result  = false;
	unsigned counter = 0;
	TSize thresh = (TSize)(length(index) * cargo(index).abundanceCut);
	if (thresh < 100) thresh = 100;
	
	TDirIterator it = begin(dir, Standard());
	TDirIterator itEnd = end(dir, Standard());
	for (; it != itEnd; ++it)
		if (*it > thresh) 
		{
			*it = (TSize)-1;
			result = true;
			++counter;
		}
	
	if (counter > 0 && cargo(index)._debugLevel >= 1)
		::std::cerr << "Removed " << counter << " k-mers" << ::std::endl;
	
	return result;
}

#endif

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
	CharString qual[2];
	String<Dna5Q> hybridSeq[2];
	
	unsigned kickoutcount = 0;
	for(unsigned i = 0; i < seqCount; ++i) 
	{
		if (options.readNaming == 0)
		{
			assignSeqId(fastaIDs[2*i], leftMates[i], formatL);		// read left Fasta id
			assignSeqId(fastaIDs[2*i+1], rightMates[i], formatR);	// read right Fasta id
			append(fastaIDs[2*i], "/L");
			append(fastaIDs[2*i+1], "/R");
		}
		
		assignSeq(seq[0], leftMates[i], formatL);					// read left Read sequence
		assignSeq(seq[1], rightMates[i], formatR);					// read right Read sequence
		reverseComplementInPlace(seq[1]);

		assignQual(qual[0], leftMates[i], formatL);					// read left ascii quality values  
		assignQual(qual[1], rightMates[i], formatR);				// read right ascii quality values  
		reverseInPlace(qual[1]);
		
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
		
		for (int j = 0; j < 2; ++j)
		{
			resize(hybridSeq[j], length(seq[j]));
			unsigned p = 0;
		
			// store dna and quality together
			for (; p < length(qual[j]) && p < length(seq[j]); ++p)
				hybridSeq[j][p] = (unsigned int) (
					(((ordValue(qual[j][p]) <= 64)? ordValue(qual[j][p]) - 33: 31) << 3) 
					| ordValue(seq[j][p]));
		
			// fill non-existent qualities with q40
			for (; p < length(seq[j]); ++p)
				hybridSeq[j][p] = (unsigned int) ((80 << 3) | ordValue(seq[j][p]));
		
			/*		std::cout << "read = " << (Dna5)((unsigned char)seq[0]& (unsigned char)0x07)<< (Dna5)((unsigned char)seq[1]& (unsigned char)0x07)<< "... ";
			 unsigned char check = seq[0];
			 unsigned intQual = (check>>3);
			 std::cout << "qual = " <<  (int)intQual << ::std::endl;*/
			
			if (options.trimLength > 0 && length(hybridSeq[j]) > (unsigned)options.trimLength)
				resize(hybridSeq[j], options.trimLength);
		}
#ifdef RAZERS_CONCATREADS
		appendValue(reads, hybridSeq[0], Generous());
		appendValue(reads, hybridSeq[1], Generous());
#else
		assign(reads[2*i], hybridSeq[0], Exact());
		assign(reads[2*i+1], hybridSeq[1], Exact());
#endif
	}
#ifdef RAZERS_CONCATREADS
	reserve(reads.concat, length(reads.concat), Exact());
#endif

	if (options._debugLevel > 1 && kickoutcount > 0) 
		::std::cerr << "Ignoring " << kickoutcount << " low quality mate-pairs.\n";
	return (seqCount > 0);
}

	template <typename TReadMatch>
	struct LessPairErrors : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
			if ((a.rseqNo >> 1) < (b.rseqNo >> 1)) return true;
			if ((a.rseqNo >> 1) > (b.rseqNo >> 1)) return false;

			// quality
#ifdef RAZERS_MATEPAIRS
			if (a.pairScore > b.pairScore) return true;
			if (a.pairScore < b.pairScore) return false;
			return a.pairId < b.pairId;
#else
			return a.editDist < b.editDist;
#endif
		}
	};
	


//////////////////////////////////////////////////////////////////////////////
// Remove low quality matches
template < typename TMatches, typename TCounts, typename TSpec >
void compactPairMatches(TMatches &matches, TCounts & /*cnts*/, RazerSOptions<TSpec> &options)
{
	typedef typename Value<TMatches>::Type					TMatch;
	typedef typename Iterator<TMatches, Standard>::Type		TIterator;
	
	unsigned readNo = -1;
	unsigned hitCount = 0;
	unsigned hitCountCutOff = options.maxHits;
	int scoreDistCutOff = InfimumValue<int>::VALUE;

	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());
	TIterator dit = it;
	TIterator ditBeg = it;

	// sort 
	::std::sort(it, itEnd, LessPairErrors<TMatch>());

	for (; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-') continue;
		if (readNo == ((*it).rseqNo >> 1))
		{ 
			if ((*it).pairScore <= scoreDistCutOff) continue;
			if (++hitCount >= hitCountCutOff)
			{
#ifdef RAZERS_MASK_READS
				if (hitCount == hitCountCutOff)
				{
					if (options.purgeAmbiguous)
					{
						dit = ditBeg;
						if (options._debugLevel >= 2)
							::std::cerr << "(read #" << readNo << " disabled)";
						options.readMask[readNo / options.WORD_SIZE] &= ~(1ul << (readNo % options.WORD_SIZE));
					} else
						if ((*it).editDist == 0)
						{
							if (options._debugLevel >= 2)
								::std::cerr << "(read #" << readNo << " disabled)";
							options.readMask[readNo / options.WORD_SIZE] &= ~(1ul << (readNo % options.WORD_SIZE));
						}
				}
#endif
				continue;
			}
		}
		else
		{
			readNo = (*it).rseqNo >> 1;
			hitCount = 0;
			if (options.distanceRange > 0)
				scoreDistCutOff = (*it).pairScore - options.distanceRange;
			ditBeg = dit;
		}
		*dit = *it;	++dit; ++it;
		*dit = *it;	++dit;
	}
	resize(matches, dit - begin(matches, Standard()));
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
	typename TCounts,
	typename TSpec >
void mapMatePairReads(
	TMatches &matches,				// resulting matches
	TGenome &genome,				// genome ...
	unsigned gseqNo,				// ... and its sequence number
	Pattern<TReadIndex, Swift<TSwiftSpec> > &swiftPatternL,
	Pattern<TReadIndex, Swift<TSwiftSpec> > &swiftPatternR,
	TVerifier &forwardPatternsL,
	TVerifier &forwardPatternsR,
	TCounts & cnts,
	char orientation,				// q-gram index of reads
	RazerSOptions<TSpec> &options)
{
	typedef typename Fibre<TReadIndex, Fibre_Text>::Type	TReadSet;
	typedef typename Size<TGenome>::Type					TSize;
	typedef typename Position<TGenome>::Type				TGPos;
	typedef typename _MakeSigned<TGPos>::Type				TSignedGPos;
	typedef typename Value<TMatches>::Type					TMatch;
	typedef typename Infix<TGenome>::Type					TGenomeInf;

	// FILTRATION
	typedef Finder<TGenome, Swift<TSwiftSpec> >				TSwiftFinderL;
	typedef Finder<TGenomeInf, Swift<TSwiftSpec> >			TSwiftFinderR;
	typedef Pattern<TReadIndex, Swift<TSwiftSpec> >			TSwiftPattern;

	// MATE-PAIR FILTRATION
//	typedef Pair<unsigned, Pair<TGPos> >					TDequeueValue;
	typedef TMatch											TDequeueValue;
	typedef Dequeue<TDequeueValue>							TDequeue;
	typedef typename TDequeue::TIter						TDequeueIterator;

	const unsigned NOT_VERIFIED = 1u << (8*sizeof(unsigned)-1);

	// iterate all genomic sequences
	if (options._debugLevel >= 1)
	{
		::std::cerr << ::std::endl << "Process genome seq #" << gseqNo;
		if (orientation == 'F')
			::std::cerr << "[fwd]";
		else
			::std::cerr << "[rev]";
	}

	TReadSet &readSetL = host(host(swiftPatternL));
	TReadSet &readSetR = host(host(swiftPatternR));

	if (empty(readSetL))
		return;

	// distance <= libLen + libErr + 2*(parWidth-readLen) - shapeLen
	// distance >= libLen - libErr - 2*parWidth + shapeLen
	TSize readLength = length(readSetL[0]);
	TSignedGPos maxDistance = options.libraryLength + options.libraryError - 2 * (int)readLength - (int)length(indexShape(host(swiftPatternL)));
	TSignedGPos minDistance = options.libraryLength - options.libraryError + (int)length(indexShape(host(swiftPatternL)));
	TGPos scanShift = (minDistance < 0)? 0: minDistance;

	// exit if contig is shorter than library size
	if (length(genome) <= scanShift)
		return;

	TGenomeInf genomeInf = infix(genome, scanShift, length(genome));
	TSwiftFinderL swiftFinderL(genome, options.repeatLength, 1);
	TSwiftFinderR swiftFinderR(genomeInf, options.repeatLength, 1);

	TDequeue fifo;						// stores left-mate potential matches
	String<TGPos> lastSeen;				// last position the left-mate was seen
	Pair<TGPos> gPair;

	fill(lastSeen, length(host(swiftPatternL)), ~(TGPos)maxDistance, Exact());

	TSize gLength = length(genome);
	TMatch mL = {	// to supress uninitialized warnings
		0, 0, 0, 0,
#ifdef RAZERS_MATEPAIRS
		0, 0, 0,
#endif
		0,
#ifdef RAZERS_DIRECT_MAQ_MAPPING
		0, 0,
#endif
		0
	};
	TMatch mR = mL;	// to supress uninitialized warnings
	mL.gseqNo = gseqNo;
	mR.gseqNo = gseqNo;
	mL.orientation = orientation;
	mR.orientation = (orientation == 'F')? 'R': 'F';

	// iterate all verification regions returned by SWIFT
	while (find(swiftFinderR, swiftPatternR, options.errorRate, options._debugLevel)) 
	{
		unsigned rseqNo = swiftPatternR.curSeqNo;
#ifdef RAZERS_MASK_READS
		if ((options.readMask[rseqNo / options.WORD_SIZE] & (1ul << (rseqNo % options.WORD_SIZE))) == 0)
			continue;
#endif

		TGPos rEndPos = endPosition(swiftFinderR) + scanShift;
		TGPos doubleParWidth = 2 * (*swiftFinderR.curHit).bucketWidth;

		// remove out-of-window left mates from fifo
		while (!empty(fifo) && front(fifo).gEnd + maxDistance + (TSignedGPos)doubleParWidth < (TSignedGPos)rEndPos)
			popFront(fifo);

		// add within-window left mates to fifo
		while (empty(fifo) || back(fifo).gEnd + minDistance < (TSignedGPos)(rEndPos + doubleParWidth))
		{
			if (find(swiftFinderL, swiftPatternL, options.errorRate, false))
			{
				mL.rseqNo = swiftPatternL.curSeqNo | NOT_VERIFIED;
				gPair = positionRange(swiftFinderL);
				mL.gBegin = gPair.i1;
				mL.gEnd = gPair.i2;
				lastSeen[swiftPatternL.curSeqNo] = mL.gEnd;
				pushBack(fifo, mL);
			} else
				break;
		}

		if (lastSeen[rseqNo] + maxDistance + doubleParWidth >= rEndPos)
		{
			// both mates have potential matches within window
			
			if (empty(fifo))
				continue;

			TDequeueIterator it = fifo.data_front;
			// iterate over fifo, if not empty
			// ignore the last element (is not within correct distance)
			do
			{
				// search left mate
				if (((*it).rseqNo & ~NOT_VERIFIED) == rseqNo)
				{
					// verify left mate (equal seqNo), if not done already
					if ((*it).rseqNo & NOT_VERIFIED)
					{
						if (matchVerify(
								*it, infix(genome, (*it).gBegin, (*it).gEnd), 
								rseqNo, readSetL, forwardPatternsL, 
								options, TSwiftSpec()))
							(*it).rseqNo &= ~NOT_VERIFIED;		// has been verified positively
						else
							(*it).rseqNo = ~NOT_VERIFIED;		// has been verified negatively
					}

					// verify right mate, if left mate matches
					if ((*it).rseqNo == rseqNo)
						if (matchVerify(
								mR, range(swiftFinderR, genomeInf),
								rseqNo, readSetR, forwardPatternsR,
								options, TSwiftSpec()))
						{
							// distance between left mate beginning and right mate end
							__int64 dist = (__int64)mR.gEnd - (__int64)(*it).gBegin;
							if (dist <= options.libraryLength + options.libraryError &&
								options.libraryLength <= dist + options.libraryError)
							{
								mL = *it;

								// transform mate readNo to global readNo
								mL.rseqNo = (rseqNo << 1);
								mR.rseqNo = (rseqNo << 1) + 1;

								// transform coordinates to the forward strand
								if (orientation == 'R') 
								{
									TSize temp = mL.gBegin;
									mL.gBegin = gLength - mL.gEnd;
									mL.gEnd = gLength - temp;
									temp = mR.gBegin;
									mR.gBegin = gLength - mR.gEnd;
									mR.gEnd = gLength - temp;
									dist = -dist;
								}

								// set a unique pair id
								mL.pairId = mR.pairId = options.nextMatePairId;
								if (++options.nextMatePairId == 0)
									options.nextMatePairId = 1;

								// score the whole match pair
								mL.pairScore = mR.pairScore = 0 - mL.editDist - mR.editDist;

								// relative positions
								mL.mateDelta = dist;
								mR.mateDelta = -dist;

								// both mates match with correct library size
/*								std::cout << "found " << rseqNo << " on " << orientation << gseqNo;
								std::cout << " dist:" << dist;
								if (orientation=='F')
									std::cout << " \t_" << mL.gBegin+1 << "_" << mR.gEnd;
								else
									std::cout << " \t_" << mR.gBegin+1 << "_" << mL.gEnd;
	//							std::cout << " L_" << (*it).gBegin << "_" << (*it).gEnd << "_" << (*it).editDist;
	//							std::cout << " R_" << mR.gBegin << "_" << mR.gEnd << "_" << mR.editDist;
								std::cout << std::endl;
*/
								if (!options.spec.DONT_DUMP_RESULTS)
								{
									appendValue(matches, mL, Generous());
									appendValue(matches, mR, Generous());
									if (length(matches) > options.compactThresh)
									{
										typename Size<TMatches>::Type oldSize = length(matches);
	//									maskDuplicates(matches);	// overlapping parallelograms cause duplicates
										compactPairMatches(matches, cnts, options);
										options.compactThresh += (options.compactThresh >> 1);
										if (options._debugLevel >= 2)
											::std::cerr << '(' << oldSize - length(matches) << " matches removed)";
									}
								}
								++options.TP;
							}
							else
								++options.FP;
						}
				}

				if (it == fifo.data_back)
					break;

				if (++it == fifo.data_end)
					it = fifo.data_begin;

			} while (true);
		}
	}
}
#endif


//////////////////////////////////////////////////////////////////////////////
// Find read matches in many genome sequences (import from Fasta)
template <
	typename TMatches, 
	typename TReadSet_, 
	typename TCounts,
	typename TSpec, 
	typename TShape,
	typename TSwiftSpec >
int mapMatePairReads(
	TMatches &				matches,
	StringSet<CharString> &	genomeFileNameList,
	StringSet<CharString> &	genomeNames,	// genome names, taken from the Fasta file
	::std::map<unsigned,::std::pair< ::std::string,unsigned> > & gnoToFileMap,
	TReadSet_ const &		readSet,
	TCounts &				cnts,
	RazerSOptions<TSpec> &	options,
	TShape const &			shape,
	Swift<TSwiftSpec> const)
{
	typedef typename Value<TReadSet_ const>::Type		TRead;
//	typedef typename Infix<TRead const>::Type			TReadInfix;
#ifdef RAZERS_CONCATREADS
	typedef StringSet<TRead, Owner<ConcatDirect<> > >	TReadSet;
#else
	typedef StringSet<TRead, Dependent<> >				TReadSet;
#endif
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
#ifdef RAZERS_CONCATREADS
		appendValue(readSetL, readSet[2*i], Generous());
		appendValue(readSetR, readSet[2*i+1], Generous());
#else
		assign(readSetL[i], readSet[2*i]);
		assign(readSetR[i], readSet[2*i+1]);
#endif
	}

	// configure q-gram index
	TIndex swiftIndexL(readSetL, shape);
	TIndex swiftIndexR(readSetR, shape);
/*	cargo(swiftIndexL).abundanceCut = options.abundanceCut;
	cargo(swiftIndexR).abundanceCut = options.abundanceCut;
	cargo(swiftIndexL)._debugLevel = options._debugLevel;
	cargo(swiftIndexR)._debugLevel = options._debugLevel;
*/
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
				mapMatePairReads(matches, genome, gseqNo, swiftPatternL, swiftPatternR, forwardPatternsL, forwardPatternsR, cnts, 'F', options);

			if (options.reverse)
			{
				reverseComplementInPlace(genome);
				mapMatePairReads(matches, genome, gseqNo, swiftPatternL, swiftPatternR, forwardPatternsL, forwardPatternsR, cnts, 'R', options);
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
