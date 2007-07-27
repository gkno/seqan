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
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_GAPS_SKIPLIST_H
#define SEQAN_HEADER_GAPS_SKIPLIST_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////


//??? TODO: UNDER CONSTRUCTION


//////////////////////////////////////////////////////////////////////////////

/*
unsigned char DEPTH [] = {
	1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1, 5, //31
	1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1, 6, //32
	1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1, 5,	//31
	1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1, 6} //32
*/


//////////////////////////////////////////////////////////////////////////////
// war als Blatt-Datenstruktur fuer die Skipliste gedacht
// noch ungetestet

struct MiniGaps
{
//____________________________________________________________________________

	enum 
	{
		SIZE = 32;
		LINK_OFFSET = 0x100 - SIZE;
	};
	
	unsigned char data_values[SIZE];
	int data_sum;
	int * data_ext;
	unsigned char data_values_length;
	unsigned char data_ext_length;
	unsigned char data_free;
//____________________________________________________________________________

	MiniGaps():
		data_sum(0),
		data_ext(0),
		data_values_length(0),
		data_ext_length(0),
		data_free(0)
	{
	}
	MiniGaps(MiniGaps const & other_):
		data_sum(other_.data_sum),
		data_values_length(other_.data_values_length),
		data_ext_length(other_.data_ext_length),
		data_free(other_.data_free)
	{
		//copy values array
		arrayCopyForward(other_.data_values, other_.data_values + other_.data_values_length, data_values);

		//copy ext array
		if (other_.data_ext)
		{
			size_t size = _computeSize4Length(*this, other_.data_ext_length);
			allocate(*this, data_ext, size);
			arrayCopyForward(other_.data_ext, other_.data_ext + other_.data_ext_length, data_ext);
		}
		else
		{
			data_ext = 0;
		}
	}

	MiniGaps const &
	operator = (MiniGaps const & other_)
	{
		//copy members
		data_sum = other_.data_sum,
		data_values_length = other_.data_values_length,
		data_free = other_.data_free;

		//copy values array
		arrayCopyForward(other_.data_values, other_.data_values + other_.data_values_length, data_values);

		//copy ext array
		size_t old_size = _computeSize4Length(*this, data_ext_length);
		size_t new_size = _computeSize4Length(*this, other_.data_ext_length);

		bool do_copy = other_.data_ext;
		bool do_allocate = do_copy && (old_size != new_size);
		bool do_deallocate = data_ext && (need_new || !do_copy);

		if (do_deallocate)
		{
			deallocate(*this, data_ext, old_size);
		}
		if (do_allocate)
		{
			allocate(*this, data_ext, new_size);
		}
		if (do_copy)
		{
			arrayCopyForward(other_.data_ext, other_.data_ext + other_.data_ext_length, data_ext);
		}

		data_ext_length = other_.data_ext_length,
	}
		
	~MiniGaps()
	{
		if (data_ext)
		{
			//delete ext array
			size_t size = _computeSize4Length(*this, data_ext_length);
			deallocate(*this, data_ext, size);
			data_ext = 0;
		}
	}
//____________________________________________________________________________

};

//////////////////////////////////////////////////////////////////////////////

inline size_t
_computeSize4Length(MiniGaps const &,
					size_t length)
{
	return (length | 0x01) + 1;
}

//////////////////////////////////////////////////////////////////////////////

inline unsigned int
getValue(MiniGaps const & me, 
		 unsigned int pos)
{
	SEQAN_ASSERT(pos < me.data_values_length)

	unsigned int ret = me.data_values[pos];
	if (ret >= MiniGaps::LINK_OFFSET)
	{
		SEQAN_ASSERT(me.data_ext)
		SEQAN_ASSERT(ret - MiniGaps::LINK_OFFSET < me.data_ext_length);

		ret = me.data_ext[ret - MiniGaps::LINK_OFFSET];
	}
	return ret;
}

//////////////////////////////////////////////////////////////////////////////

inline void
assignValue(MiniGaps & me, 
			unsigned int pos,
			unsigned int val)
{
	SEQAN_ASSERT(pos < me.data_values_length)

	unsigned char old_val = me.data_values[pos];
	me.data_sum = me.data_sum - old_val + val;

	bool old_is_large = (old_val >= MiniGaps::LINK_OFFSET);
	bool new_is_large = (val >= MiniGaps::LINK_OFFSET);

	if (new_is_large)
	{
		if (old_is_large)
		{//use current ext
			me.data_ext[old_val] = val;
		}
		else
		{//get new ext
			if (me.data_free)
			{//use recycled
				//store next in free_list
				unsigned char next_free = me.data_ext[me.data_free];

				//link in values_array
				me.data_ext[me.data_free] = val;
				me.data_values[pos] = me.data_free;

				//unlink from free_list
				me.data_free = next_free;
			}
			else
			{//create ext
				size_t old_size = _computeSize4Length(*this, me.data_ext_length);
				size_t new_size = _computeSize4Length(*this, me.data_ext_length + 1);

				if (old_size != new_size)
				{//expand buffer
					unsigned int * new_data_ext;
            		allocate(*this, new_data_ext, new_size);
					arrayCopyForward(me.data_ext, me.data_ext + me.data_ext_length, new_data_ext);
					deallocate(*this, me.data_ext, old_size);
					me.data_ext = new_data_ext;
				}

				//link in values_array
				me.data_ext[me.data_ext_length] = val;
				me.data_values[pos] = me.data_ext_length;
				++me.data_ext_length;
			}
		}
	}
	else
	{
		if (old_is_large)
		{//dont need ext anymore
			me.data_ext[old_val] = me.data_free;
			me.data_free = old_val;
		}

		//set value in values_array
		me.data_values[pos] = val;
	}
	
}
//////////////////////////////////////////////////////////////////////////////

inline void
insertValue(MiniGaps & me, 
			unsigned int pos,
			unsigned int val)
{
	SEQAN_ASSERT(pos <= me.data_values_length)
	SEQAN_ASSERT(me.data_values_length < MiniGaps::SIZE)

	if (pos < me.data_values_length)
	{//make room
		arrayCopyBackward(me.data_values + pos, me.data_values + me.data_values_length, me.data_values + pos + 1);
	}

	++me.data_ext_length;
	assignValue(me, pos, val);

	me.data_sum += val;
}

//////////////////////////////////////////////////////////////////////////////

inline void
removeValue(MiniGaps & me, 
			unsigned int pos)
{
	SEQAN_ASSERT(pos < me.data_values_length)

	unsigned char old_val = me.data_values[pos];
	if (old_val >= MiniGaps::LINK_OFFSET)
	{//dont need ext anymore
		me.data_ext[old_val] = me.data_free;
		me.data_free = old_val;
	}

	if (pos < me.data_values_length - 1)
	{//remove room
		arrayCopyForward(me.data_values + pos + 1, me.data_values + me.data_values_length, me.data_values + pos);
	}

	--me.data_ext_length;
	me.data_sum -= old_val;
}

//////////////////////////////////////////////////////////////////////////////

//note: function returns rest value in val
inline unsigned int
findEntry(MiniGaps & me, 
		  unsigned int & val)
{
	if (val >= me.data_sum)
	{
		val -= me.data_sum;
		return me.data_values_length;
	}

	unsigned short * it = me.data_values;
	while (true)
	{
		unsigned int it_val = *it;
		if (it_val >= MiniGaps::LINK_OFFSET)
		{
			it_val = me.data_ext[it_val];
		}
		if (it_val >= val)
		{//stop search
			return it - me.data_values;
		}

		val -= it_val;
		++it;
	}
}
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
}// namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////
#endif //#ifndef SEQAN_HEADER_...
