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

#ifndef SEQAN_MMAP_FASTA_H
#define SEQAN_MMAP_FASTA_H

namespace SEQAN_NAMESPACE_MAIN
{

	// define memory mapped stringset
	typedef StringSet<String<char, MMap<> >, Owner<ConcatDirect<> > >	MultiFasta;


	template <typename TValue>
	inline bool
	_isLineBreak(TValue value)
	{
		return (value == '\n' || value == '\r');
	}

	template <typename TIterator>
	inline bool
	_seekLineBreak(TIterator &it, TIterator itEnd)
	{
		while (!_isLineBreak(*it))
			if (++it == itEnd) return false;
		return true;
	}

	template <typename TIterator>
	inline bool
	_seekNonLineBreak(TIterator &it, TIterator itEnd)
	{
		while (_isLineBreak(*it))
			if (++it == itEnd) return false;
		return true;
	}


//////////////////////////////////////////////////////////////////////////////
// File Formats - Fasta
//////////////////////////////////////////////////////////////////////////////

	// split stringset into single Fasta sequences
	template < typename TValue, typename TConfig, typename TDelimiter >
	inline void
	split(
		StringSet<String<TValue, MMap<TConfig> >, Owner<ConcatDirect<TDelimiter> > > &me, 
		Fasta)
	{
		typedef String<TValue, MMap<TConfig> >						TString;
		typedef StringSet<TString, ConcatDirect<TDelimiter> >		TStringSet;
		typedef typename Iterator<TString const, Standard>::Type	TIterator;

		clear(me.limits);

		TIterator itBeg = begin(me.concat, Standard());
		TIterator itEnd = end(me.concat, Standard());
		bool newLine = true;
		for (TIterator it = itBeg; it != itEnd; ++it)
		{
			TValue c = *it;
			if (newLine && c == '>')
				appendValue(me.limits, it - itBeg, Generous());
			newLine = (c == '\n' || c == '\r');
		}
		if (empty(me.limits))
			appendValue(me.limits, 0);
		appendValue(me.limits, itEnd - itBeg);
	}

	template <typename TSeq, typename TFastaSeq>
	inline void
	assignSeq(
		TSeq & dst,
		TFastaSeq const & fasta,
		Fasta)
	{
		typedef typename Iterator<TFastaSeq const, Standard>::Type	TIterator;
		typedef typename Iterator<TSeq, Standard>::Type				TDstIterator;

		TIterator it = begin(fasta, Standard());
		TIterator itEnd = end(fasta, Standard());

		clear(dst);
		
		// skip Fasta id
		if (it == itEnd) return;
		if (*it == '>')
		{
			if (!_seekLineBreak(it, itEnd)) return;
			if (!_seekNonLineBreak(it, itEnd)) return;
		}

		// copy sequence
		resize(dst, itEnd - it);		
		TDstIterator dit = begin(dst, Standard());
		for (; it != itEnd; ++it)
			if (!_isLineBreak(*it))
			{
				*dit = *it;
				++dit;
			}
		resize(dst, dit - begin(dst, Standard()));
	}

	template <typename TSeq, typename TFastaSeq>
	inline void
	assignSeqId(
		TSeq & dst,
		TFastaSeq const & fasta,
		Fasta)
	{
		typedef typename Iterator<TFastaSeq const, Standard>::Type	TIterator;

		TIterator itBeg = begin(fasta, Standard());
		TIterator itEnd = end(fasta, Standard());
		TIterator it = itBeg;
		
		clear(dst);
		if (it == itEnd) return;
		if (*it == '>')
		{
			_seekLineBreak(it, itEnd);
			assign(dst, infix(fasta, 1, it - itBeg));
		}
	}

//////////////////////////////////////////////////////////////////////////////
// File Formats - Fastq (Fasta extension for quality values)
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format.tag.Fastq:
	FASTQ file format for sequences.
*/
struct TagFastq_;
typedef Tag<TagFastq_> const Fastq;

	// split stringset into single Fasta sequences
	template < typename TValue, typename TConfig, typename TDelimiter >
	inline void
	split(
		StringSet<String<TValue, MMap<TConfig> >, Owner<ConcatDirect<TDelimiter> > > &me, 
		Fastq)
	{
		typedef String<TValue, MMap<TConfig> >						TString;
		typedef StringSet<TString, ConcatDirect<TDelimiter> >		TStringSet;
		typedef typename Iterator<TString const, Standard>::Type	TIterator;

		clear(me.limits);

		TIterator itBeg = begin(me.concat, Standard());
		TIterator itEnd = end(me.concat, Standard());
		bool newLine = true;
		for (TIterator it = itBeg; it != itEnd; ++it)
		{
			TValue c = *it;
			if (newLine && c == '@')
				appendValue(me.limits, it - itBeg, Generous());
			newLine = (c == '\n' || c == '\r');
		}
		if (empty(me.limits))
			appendValue(me.limits, 0);
		appendValue(me.limits, itEnd - itBeg);
	}

	template <typename TSeq, typename TFastaSeq>
	inline void
	assignSeq(
		TSeq & dst,
		TFastaSeq const & fasta,
		Fastq)
	{
		typedef typename Iterator<TFastaSeq const, Standard>::Type	TIterator;
		typedef typename Iterator<TSeq, Standard>::Type				TDstIterator;

		TIterator it = begin(fasta, Standard());
		TIterator itEnd = end(fasta, Standard());

		clear(dst);
		
		// skip Fasta id
		if (it == itEnd) return;
		if (*it == '@')
		{
			if (!_seekLineBreak(it, itEnd)) return;
			if (!_seekNonLineBreak(it, itEnd)) return;
		}

		// copy sequence
		resize(dst, itEnd - it);		
		TDstIterator dit = begin(dst, Standard());
		for (; it != itEnd; ++it) 
		{
			if (_isLineBreak(*it))
			{
				if (!_seekNonLineBreak(it, itEnd)) break;
				if (*it == '+') break;
			}
			*dit = *it;
			++dit;
		}
		resize(dst, dit - begin(dst, Standard()));
	}

	template <typename TSeq, typename TFastaSeq>
	inline void
	assignSeqId(
		TSeq & dst,
		TFastaSeq const & fasta,
		Fastq)
	{
		typedef typename Iterator<TFastaSeq const, Standard>::Type	TIterator;

		TIterator itBeg = begin(fasta, Standard());
		TIterator itEnd = end(fasta, Standard());
		TIterator it = itBeg;
		
		clear(dst);
		if (it == itEnd) return;
		if (*it == '@')
		{
			_seekLineBreak(it, itEnd);
			assign(dst, infix(fasta, 1, it - itBeg));
		}
	}

	template <typename TSeq, typename TFastaSeq>
	inline void
	assignQual(
		TSeq & dst,
		TFastaSeq const & fasta,
		Fastq)
	{
		typedef typename Iterator<TFastaSeq const, Standard>::Type	TIterator;
		typedef typename Iterator<TSeq, Standard>::Type				TDstIterator;

		TIterator it = begin(fasta, Standard());
		TIterator itEnd = end(fasta, Standard());

		clear(dst);
		
		if (it == itEnd) return;
		if (*it == '@')
		{
			// seek quality id
			do {
				if (!_seekLineBreak(it, itEnd)) return;
				if (!_seekNonLineBreak(it, itEnd)) return;
			} while (*it != '+');

			// skip quality id
			if (!_seekLineBreak(it, itEnd)) return;
			if (!_seekNonLineBreak(it, itEnd)) return;

			// copy sequence
			resize(dst, itEnd - it);		
			TDstIterator dit = begin(dst, Standard());
			for (; it != itEnd; ++it)
				if (!_isLineBreak(*it))
				{
					*dit = *it;
					++dit;
				}
			resize(dst, dit - begin(dst, Standard()));
		}
	}

	template <typename TSeq, typename TFastaSeq>
	inline void
	assignQualId(
		TSeq & dst,
		TFastaSeq const & fasta,
		Fastq)
	{
		typedef typename Iterator<TFastaSeq const, Standard>::Type	TIterator;

		TIterator itBeg = begin(fasta, Standard());
		TIterator itEnd = end(fasta, Standard());
		TIterator it1 = itBeg;
		
		clear(dst);
		if (it1 == itEnd) return;
		if (*it1 == '@')
		{
			do {
				if (!_seekLineBreak(it1, itEnd)) return;
				if (!_seekNonLineBreak(it1, itEnd)) return;
			} while (*it1 != '+');
			TIterator it2 = it1;
			_seekLineBreak(it2, itEnd);
			assign(dst, infix(fasta, (it1 - itBeg) + 1, it2 - itBeg));
		}
	}

}

#endif
