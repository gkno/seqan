#ifndef SEQAN_HEADER_FIND_MULTI_H
#define SEQAN_HEADER_FIND_MULTI_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

struct _MultipatternFinder;
typedef Tag<_MultipatternFinder> MultipatternFinder;
	
//____________________________________________________________________________

template <typename THaystack>
class Finder<THaystack, MultipatternFinder>
{
//____________________________________________________________________________
private:
	unsigned int data_pattern;

public:
	Finder():
		data_pattern(0)
	{
SEQAN_CHECKPOINT
	}

	Finder(Finder const & other_):
		data_pattern(other_.data_pattern)
	{
SEQAN_CHECKPOINT
	}

	~Finder()
	{
SEQAN_CHECKPOINT
	}
//____________________________________________________________________________

	Finder & 
	operator = (Finder const & other_)
	{
SEQAN_CHECKPOINT
		data_pattern = other_.data_pattern;
		return *this;
	}
//____________________________________________________________________________

	friend inline unsigned int &
	needle(Finder & me)
	{
SEQAN_CHECKPOINT
		return me.data_pattern;
	}
	friend inline unsigned int const &
	needle(Finder const & me)
	{
SEQAN_CHECKPOINT
		return me.data_pattern;
	}
//____________________________________________________________________________

	friend inline void
	setNeedle(Finder & me, unsigned int const needleIndex_)
	{
SEQAN_CHECKPOINT
		me.data_pattern = needleIndex_;
	}

//____________________________________________________________________________

	friend inline void
	init(Finder & me)
	{
SEQAN_CHECKPOINT
		me.data_pattern = 0;
	}
//____________________________________________________________________________


//////////////////////////////////////////////////////////////////////////////

template <typename THaystack, typename TNeedle>
friend inline bool
find(Finder & me,
	 THaystack & hstk,
	 TNeedle const & ndl)
{
SEQAN_CHECKPOINT
	while ( needle(me) < length(ndl) )
	{
		Finder<THaystack, Horspool> horspool(ndl[needle(me)]);
		bool found = find(horspool, hstk, ndl[needle(me)]);
		if (found)
		{
			return true;
		}
		setPosition(hstk, 0);
		++needle(me);
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////////
/*
template <typename THaystack, typename TNeedle>
bool
findNext(Finder & me,
		 THaystack & hstk,
		 TNeedle const & ndl)
{
SEQAN_CHECKPOINT
	++hstk;
	return find(me, hstk, ndl);
}*/

//////////////////////////////////////////////////////////////////////////////

};

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
