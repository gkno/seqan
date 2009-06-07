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

 ============================================================================
  $Id: graph_consensus_readstore.h 1811 2008-03-31 15:38:54Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_CONSENSUS_READSTORE_H
#define SEQAN_HEADER_GRAPH_CONSENSUS_READSTORE_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Read Store
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet = Dna, typename TSpec = Default>
class ReadStore;

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec>
class ReadStore
{
	public:
		typedef typename Size<ReadStore>::Type TSize;
		String<TAlphabet, External<> > data_reads;
		//String<TAlphabet, Alloc<> > data_reads;
		String<char, External<> > data_qualities;
		//String<char, Alloc<> > data_qualities;
		StringSet<String<char, External<> >, Owner<ConcatDirect<> > > data_names;
		//StringSet<String<char, Alloc<> > > data_names;
		String<Pair<TSize, TSize> > data_begin_end;
		String<Pair<TSize, TSize> > data_clr;
		String<TSize> data_frg_id;
		TSize data_pos_count;
		
	public:
		ReadStore() : data_pos_count(0)
		{
			SEQAN_CHECKPOINT
			clear(data_reads);
			clear(data_qualities);
			clear(data_names);
			clear(data_begin_end);
			clear(data_clr);
			clear(data_frg_id);
		}

		~ReadStore() 
		{
			SEQAN_CHECKPOINT
		}


	private:
		ReadStore(ReadStore const & _other)
		{
			SEQAN_CHECKPOINT
			data_reads = _other.data_reads;
			data_qualities = _other.data_qualities;
			data_names = _other.data_names;
			data_begin_end = _other.data_begin_end;
			data_frg_id = _other.data_frg_id;
			data_clr = _other.data_clr;
			data_pos_count = _other.data_pos_count;
		}

		ReadStore const& 
		operator = (ReadStore const& _other) 
		{
			SEQAN_CHECKPOINT
			if (this == &_other) return *this;
			data_reads = _other.data_reads;
			data_qualities = _other.data_qualities;
			data_names = _other.data_names;
			data_begin_end = _other.data_begin_end;
			data_frg_id = _other.data_frg_id;
			data_clr = _other.data_clr;
			data_pos_count = _other.data_pos_count;
			return *this;
		}
};
	


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline void 
reverseComplementInPlace(String<Dna5Q, TSpec>& str1)
{
	typedef typename Size<String<Dna5Q, TSpec> >::Type TSize;
	TSize pos1 = 0;
	TSize pos2 = length(str1)-1;
	for(;((pos1 < length(str1)) && (pos1<=pos2)); ++pos1, --pos2) {
		if ((str1[pos1] == 'A') && (str1[pos2] == 'A')) {
			Dna5Q c = 'T';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'T';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'A') && (str1[pos2] == 'C')) {
			Dna5Q c = 'G';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'T';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'A') && (str1[pos2] == 'G')) {
			Dna5Q c = 'C';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'T';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'A') && (str1[pos2] == 'T')) {
			Dna5Q c = 'A';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'T';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'A') && (str1[pos2] == 'N')) {
			Dna5Q c = 'N';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'T';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'C') && (str1[pos2] == 'A')) {
			Dna5Q c = 'T';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'G';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'C') && (str1[pos2] == 'C')) {
			Dna5Q c = 'G';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'G';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'C') && (str1[pos2] == 'G')) {
			Dna5Q c = 'C';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'G';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'C') && (str1[pos2] == 'T')) {
			Dna5Q c = 'A';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'G';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'C') && (str1[pos2] == 'N')) {
			Dna5Q c = 'N';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'G';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'G') && (str1[pos2] == 'A')) {
			Dna5Q c = 'T';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'C';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'G') && (str1[pos2] == 'C')) {
			Dna5Q c = 'G';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'C';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'G') && (str1[pos2] == 'G')) {
			Dna5Q c = 'C';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'C';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'G') && (str1[pos2] == 'T')) {
			Dna5Q c = 'A';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'C';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'G') && (str1[pos2] == 'N')) {
			Dna5Q c = 'N';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'C';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		}  else if ((str1[pos1] == 'T') && (str1[pos2] == 'A')) {
			Dna5Q c = 'T';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'A';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'T') && (str1[pos2] == 'C')) {
			Dna5Q c = 'G';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'A';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'T') && (str1[pos2] == 'G')) {
			Dna5Q c = 'C';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'A';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'T') && (str1[pos2] == 'T')) {
			Dna5Q c = 'A';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'A';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'T') && (str1[pos2] == 'N')) {
			Dna5Q c = 'N';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'A';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'N') && (str1[pos2] == 'A')) {
			Dna5Q c = 'T';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'N';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'N') && (str1[pos2] == 'C')) {
			Dna5Q c = 'G';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'N';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'N') && (str1[pos2] == 'G')) {
			Dna5Q c = 'C';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'N';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'N') && (str1[pos2] == 'T')) {
			Dna5Q c = 'A';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'N';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'N') && (str1[pos2] == 'N')) {
			Dna5Q c = 'N';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'N';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		}
	}
}


template<typename TSpec>
inline void 
reverseComplementInPlace(String<Dna5Q, TSpec> const& str1)
{
	reverseComplementInPlace(const_cast<String<Dna5Q, TSpec>&>(str1));
}





//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TSize, typename TRead>
inline void 
loadRead(ReadStore<TAlphabet, TSpec>& readSt, 
		 TSize index,
		 TRead& seq) 
{
	SEQAN_CHECKPOINT
	clear(seq);
	seq = infix(readSt.data_reads, (value(readSt.data_begin_end, index)).i1, (value(readSt.data_begin_end, index)).i2);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TSpec2, typename TSize, typename TRead>
inline void 
loadRead(ReadStore<TAlphabet, TSpec>& readSt,
		 String<TAlphabet, TSpec2>& readString, 
		 TSize index,
		 TRead& seq) 
{
	SEQAN_CHECKPOINT
	clear(seq);
	seq = infix(readString, (value(readSt.data_begin_end, index)).i1, (value(readSt.data_begin_end, index)).i2);
}




//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TSize, typename TRead>
inline void 
loadQuality(ReadStore<TAlphabet, TSpec>& readSt,
			TSize index,
			TRead& seq) 
{
	SEQAN_CHECKPOINT
	clear(seq);
	seq = infix(readSt.data_qualities, (value(readSt.data_begin_end, index)).i1, (value(readSt.data_begin_end, index)).i2);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TSpec2, typename TSize, typename TRead>
inline void 
loadQuality(ReadStore<TAlphabet, TSpec>& readSt,
			String<char, TSpec2>& qualityString, 
			TSize index,
			TRead& seq) 
{
	SEQAN_CHECKPOINT
	clear(seq);
	seq = infix(qualityString, (value(readSt.data_begin_end, index)).i1, (value(readSt.data_begin_end, index)).i2);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TSize, typename TName>
inline void 
loadExtId(ReadStore<TAlphabet, TSpec>& readSt,
		  TSize index,
		  TName& name) 
{
	SEQAN_CHECKPOINT
	name = value(readSt.data_names, index);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TSize>
inline typename Id<ReadStore<TAlphabet, TSpec> >::Type
loadFrgId(ReadStore<TAlphabet, TSpec>& readSt, 
		  TSize index) 
{
	SEQAN_CHECKPOINT
	return value(readSt.data_frg_id, index);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec>
inline typename Size<ReadStore<TAlphabet, TSpec> >::Type
length(ReadStore<TAlphabet, TSpec>& readSt) 
{
	SEQAN_CHECKPOINT
	return readSt.data_pos_count;
}



//////////////////////////////////////////////////////////////////////////////


template<typename TFile, typename TAlphabet, typename TSpec, typename TFragmentStore, typename TLibraryStore, typename TContigStore>
inline void 
write(TFile& target,
	  ReadStore<TAlphabet, TSpec>& readSt,
	  TFragmentStore&,
	  TLibraryStore&,
	  TContigStore&,
	  CeleraFrg) 
{
	SEQAN_CHECKPOINT
	typedef ReadStore<TAlphabet, TSpec> TReadStore;
	typedef typename Id<TReadStore>::Type TId;
	typedef typename Size<TReadStore>::Type TSize;
	typedef typename Value<TFile>::Type TValue;
	typedef String<char> TCharString;

	// Write Reads
	for(TSize i = 0; i<length(readSt); ++i) {
		_streamWrite(target,"{FRG\n");
		_streamWrite(target,"act:");
		_streamPut(target, 'A');
		_streamPut(target, '\n');
		_streamWrite(target,"acc:");
		_streamPutInt(target, i + 1);
		_streamPut(target, '\n');
		_streamWrite(target,"typ:");
		_streamPut(target, 'R');
		_streamPut(target, '\n');
		_streamWrite(target,"src:\n");
		_streamWrite(target, value(readSt.data_names, i));
		_streamWrite(target, "\n.\n");
		_streamWrite(target,"etm:");
		_streamPut(target, '0');
		_streamPut(target, '\n');
		_streamWrite(target,"seq:\n");
		String<TAlphabet> read;
		loadRead(readSt, i, read);
		for(TSize k = 0;k<length(read); k+=70) {
			TSize endK = k + 70;
			if (endK > length(read)) endK = length(read);
			_streamWrite(target, infix(read, k, endK));
			_streamPut(target, '\n');
		}
		_streamWrite(target, ".\n");
		_streamWrite(target,"qlt:\n");
		String<char> qlt;
		loadQuality(readSt, i, qlt);
		for(TSize k = 0;k<length(qlt); k+=70) {
			TSize endK = k + 70;
			if (endK > length(qlt)) endK = length(qlt);
			_streamWrite(target, infix(qlt, k, endK));
			_streamPut(target, '\n');
		}
		_streamWrite(target, ".\n");
		_streamWrite(target,"clr:");
		_streamPutInt(target, value(readSt.data_clr, i).i1);
		_streamPut(target, ',');
		_streamPutInt(target, value(readSt.data_clr, i).i2);
		_streamPut(target, '\n');
		_streamWrite(target,"}\n");
	}
}


//////////////////////////////////////////////////////////////////////////////


template<typename TFile, typename TAlphabet, typename TSpec, typename TFragmentStore, typename TLibraryStore, typename TContigStore>
inline void 
write(TFile& target,
	  ReadStore<TAlphabet, TSpec>& readSt,
	  TFragmentStore&,
	  TLibraryStore&,
	  TContigStore& ctgSt,
	  CeleraCgb) 
{
	SEQAN_CHECKPOINT
	typedef ReadStore<TAlphabet, TSpec> TReadStore;
	typedef typename Id<TReadStore>::Type TId;
	typedef typename Size<TReadStore>::Type TSize;
	typedef typename Value<TFile>::Type TValue;
	typedef String<char> TCharString;
	String<GappedRead<> >& origGappedReads = value(ctgSt.data_reads, 0);
	String<char> contig;
	loadContig(ctgSt, 0, contig);

	// Find smallest offset
	TSize offsetLeft = length(contig);
	for(TSize k = 0; k<length(origGappedReads); ++k) if ((value(origGappedReads, k)).data_offset < offsetLeft) offsetLeft = (value(origGappedReads, k)).data_offset;

	// Sort the reads
	typedef std::set<std::pair<TSize, TSize> > TOffsetIndexPair;
	TOffsetIndexPair offsetIndexSet;
	for(TSize k = 0; k<length(origGappedReads); ++k) {
		TSize clr1 = (value(origGappedReads, k)).data_clr.i1;
		TSize clr2 = (value(origGappedReads, k)).data_clr.i2;
		clr1 += ((value(origGappedReads, k)).data_offset - offsetLeft);
		clr2 += ((value(origGappedReads, k)).data_offset - offsetLeft);
		if (clr1 > clr2) { TSize tmp = clr1; clr1 = clr2; clr2 = tmp; }
		offsetIndexSet.insert(std::make_pair(clr1, k));
	}
	String<GappedRead<> > gappedReads;
	for(typename TOffsetIndexPair::const_iterator itPos = offsetIndexSet.begin(); itPos != offsetIndexSet.end(); ++itPos) {
		appendValue(gappedReads, value(origGappedReads, (*itPos).second));
	}

	//// Write IAF record for all reads
	//for(TSize k = 0; k<length(gappedReads); ++k) {
	//	_streamWrite(target,"{IAF\n");
	//	_streamWrite(target,"acc:");
	//	_streamPutInt(target, (value(gappedReads, k)).data_source + 1);
	//	_streamPut(target, '\n');
	//	_streamWrite(target,"typ:");
	//	_streamPut(target, 'R');
	//	_streamPut(target, '\n');
	//	_streamWrite(target,"chi:0\ncha:0\nclr:-1,-1\nmst:U\n}\n");
	//}

	// Write Header
	_streamWrite(target,"{IUM\nacc:0\nsrc:\ngen> @@ [0,0]\n.\ncov:0.000\nsta:X\nfur:X\nabp:0\nbbp:0\n");
	_streamWrite(target,"len:");
	_streamPutInt(target, length(contig));
	_streamPut(target, '\n');
	_streamWrite(target,"cns:\n.\nqlt:\n.\nfor:0\n");
	_streamWrite(target,"nfr:");
	_streamPutInt(target, length(readSt));
	_streamPut(target, '\n');
	
	// Write gapped reads
	for(TSize k = 0; k<length(gappedReads); ++k) {
		_streamWrite(target,"{IMP\n");
		_streamWrite(target,"typ:");
		_streamPut(target, 'R');
		_streamPut(target, '\n');
		_streamWrite(target,"mid:");
		_streamPutInt(target, (value(gappedReads, k)).data_source + 1);
		_streamPut(target, '\n');
		TSize clr1 = (value(gappedReads, k)).data_clr.i1;
		TSize clr2 = (value(gappedReads, k)).data_clr.i2;
		clr1 += ((value(gappedReads, k)).data_offset - offsetLeft);
		clr2 += ((value(gappedReads, k)).data_offset - offsetLeft);
		_streamWrite(target,"con:");
		//TSize orig1 = clr1;
		//TSize orig2 = clr2;
		//if (orig1 > orig2) {TSize tmp = orig1; orig1 = orig2; orig2 = tmp; }
		//TSize best = 0;
		//TSize bestDist = 0;
		//for(TSize other = 0; other<length(gappedReads); ++other) {
		//	if (other == k) continue;
		//	TSize thisRead1 = (value(gappedReads, other)).data_clr.i1;
		//	TSize thisRead2 = (value(gappedReads, other)).data_clr.i2;
		//	thisRead1 += ((value(gappedReads, other)).data_offset - offsetLeft);
		//	thisRead2 += ((value(gappedReads, other)).data_offset - offsetLeft);
		//	if (thisRead1 > thisRead2) {TSize tmp = thisRead1; thisRead1 = thisRead2; thisRead2 = tmp; }
		//	if ((orig1 > thisRead1) && (orig2 < thisRead2)) {
		//		if ((best == 0) ||
		//			(bestDist > (thisRead2 - thisRead1))) {
		//			bestDist = (thisRead2 - thisRead1);
		//			best = other + 1;
		//		}
		//	}
		//}

		//_streamPutInt(target, best);
		_streamPut(target, '0');
		_streamPut(target, '\n');
		_streamWrite(target,"pos:");
		_streamPutInt(target, clr1);
		_streamPut(target, ',');
		_streamPutInt(target, clr2);
		_streamPut(target, '\n');
		_streamWrite(target,"dln:0\n");
		_streamWrite(target,"del:\n");
		_streamWrite(target,"}\n");
	}
	_streamWrite(target,"}\n");
}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
