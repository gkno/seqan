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
  $Id: $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SUMLIST_H
#define SEQAN_HEADER_SUMLIST_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
class SumList;

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
struct Value< SumList<TValue, TSpec> >
{
	typedef TValue Type;
};


//////////////////////////////////////////////////////////////////////////////
// MiniSumList
//////////////////////////////////////////////////////////////////////////////

template <unsigned short SIZE = 0x0040, typename TSpec = Default>
struct MiniSumList;

template <typename TValue, unsigned short SIZE, typename TSpec>
class SumList<TValue, MiniSumList<SIZE, TSpec> >
{
public:
	typedef SumList<TValue, MiniSumList<SIZE, TSpec> > TSumList;
	typedef typename Size<TSumList>::Type TSize;

	unsigned char data_ [SIZE];
	TSize data_length;  //number of elements in list
	TSize data_size;  //number of bytes used in data_
	TValue data_sum; //sum of all list numbers

	SumList()
		: data_length(0)
		, data_size(0)
		, data_sum(TValue())
	{}

	SumList(SumList const & other)
		: data_length(other.data_length)
		, data_size(other.data_size)
		, data_sum(other.data_sum)
	{
		arrayCopyForward(other.data_, other.data_ + SIZE, data_);
	}

	~SumList() 
	{}

	SumList const &
	operator = (SumList const & other)
	{
		assign(*this, other);
	}


//____________________________________________________________________________

	struct Entry
	{
		enum
		{
			LIMIT_0 = 1 << 6,
			LIMIT_1 = 1 << 14,
			LIMIT_2 = 1 << 30
		};

		static const unsigned char SIZES [4]; 

		union _union
		{
			struct struct_0
			{
				unsigned char select: 2;
				unsigned char value_0: 6;
			} _0;
			struct struct_1
			{
				unsigned short select_1: 2;
				unsigned short value_1: 14;
			} _1;
			struct struct_2
			{
				unsigned int select_2: 2;
				unsigned int value_2: 30;
			} _2;
			unsigned char value_3[sizeof(TValue)+1];
		} data;

		Entry() {}
		Entry(TValue val) { this->assignValue(val); }
		~Entry() {}

		inline TValue getValue()
		{
			switch(data._0.select)
			{
			case 0: return data._0.value_0;
			case 1: return data._1.value_1;
			case 2: return data._2.value_2;
			default: return * reinterpret_cast<TValue *>(data.value_3 + 1);
			}
		}

		inline int size()
		{
			return SIZES[data._0.select];
		}

		inline void assignValue(TValue val)
		{
			if (val < LIMIT_1)
			{
				if (val < LIMIT_0)
				{
					data._0.select = 0;
					data._0.value_0 = val;
				}
				else
				{
					data._0.select = 1;
					data._1.value_1 = val;
				}
			}
			else
			{
				if (val < LIMIT_2)
				{
					data._0.select = 2;
					data._2.value_2 = val;
				}
				else
				{
					data._0.select = 3;
					*reinterpret_cast<TValue *>(data.value_3 + 1) = val;
				}
			}
		}
	};

//____________________________________________________________________________
};


template <typename TValue, unsigned short SIZE, typename TSpec>
const unsigned char SumList<TValue, MiniSumList<SIZE, TSpec> >::Entry::SIZES [4] = {1, 2, 4, 1 + sizeof(TValue)};


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned short SIZE, typename TSpec>
struct Size< SumList<TValue, MiniSumList<SIZE, TSpec> > >
{
	typedef unsigned short Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned short SIZE, typename TSpec>
struct Position< SumList<TValue, MiniSumList<SIZE, TSpec> > >
{
	typedef unsigned short Type;
};

//////////////////////////////////////////////////////////////////////////////

//update current entry
//returns false on overflow, otherwise true
template <typename TValue, unsigned short SIZE, typename TSpec, typename TPosition, typename TValue2>
inline bool _MiniSumList_updateValue(SumList<TValue, MiniSumList<SIZE, TSpec> > & me,
									 TPosition byte_pos,
									 TValue2 new_value)
{
SEQAN_ASSERT(byte_pos < me.data_size)

	typedef typename SumList<TValue, MiniSumList<SIZE, TSpec> >::Entry TEntry;

	TEntry new_entr(new_value);
	int new_entry_size = new_entr.size();

	TEntry & old_entr = * reinterpret_cast<TEntry *>(me.data_ + byte_pos);
	int old_entry_size = old_entr.size();

	int new_size = me.data_size + new_entry_size - old_entry_size;
	if (new_size > SIZE) return false; //not enough space

	//update size and sum
	me.data_size = new_size;
	me.data_sum += new_value;
	me.data_sum -= old_entr.getValue();

	if (new_entry_size != old_entry_size)
	{// make room
		arrayCopy(me.data_ + byte_pos + old_entry_size, me.data_ + me.data_size, me.data_ + byte_pos + new_entry_size);
	}

	//set value
	arrayCopy(new_entr.data.value_3, new_entr.data.value_3 + new_entry_size, me.data_ + byte_pos);
	return true;	
}

//////////////////////////////////////////////////////////////////////////////

//appends value
//returns false on overflow, otherwise true
template <typename TValue, unsigned short SIZE, typename TSpec, typename TValue2>
inline bool appendValue(SumList<TValue, MiniSumList<SIZE, TSpec> > & me,
						TValue2 new_value)
{
	typedef typename SumList<TValue, MiniSumList<SIZE, TSpec> >::Entry TEntry;

	TEntry new_entr(new_value);
	int new_entry_size = new_entr.size();

	int new_size = me.data_size + new_entry_size;
	if (new_size > SIZE) return false; //not enough space

	//set value
	arrayCopy(new_entr.data.value_3, new_entr.data.value_3 + new_entry_size, me.data_ + me.data_size);

	//update size and sum
	me.data_size = new_size;
	++me.data_length;
	me.data_sum += new_value;

	return true;	
}

//////////////////////////////////////////////////////////////////////////////

//update current entry and inserts two more entries behind the current
//returns false on overflow, otherwise true
template <typename TValue, unsigned short SIZE, typename TSpec, typename TPosition, typename TValue2, typename TValue3, typename TValue4 >
inline bool _MiniSumList_updateAndInsert(SumList<TValue, MiniSumList<SIZE, TSpec> > & me,
										 TPosition byte_pos,
										 TValue2 update_value,
										 TValue3 insert_value1,
										 TValue4 insert_value2)
{
SEQAN_ASSERT(byte_pos < me.data_size)

	typedef typename SumList<TValue, MiniSumList<SIZE, TSpec> >::Entry TEntry;

	TEntry update_entr(update_value);
	int update_entr_size = update_entr.size();

	TEntry insert_entr1(insert_value1);
	int insert_entr1_size = insert_entr1.size();

	TEntry insert_entr2(insert_value2);
	int insert_entr2_size = insert_entr2.size();

	int new_entry_size = update_entr_size + insert_entr1_size + insert_entr2_size;

	TEntry & old_entr = * reinterpret_cast<TEntry *>(me.data_ + byte_pos);
	int old_entry_size = old_entr.size();

	int new_size = me.data_size + new_entry_size - old_entry_size;
	if (new_size > SIZE) return false;

	//update size and sum
	me.data_size = new_size;
	me.data_length += 2;
	me.data_sum += update_value + insert_value1 + insert_value2;
	me.data_sum -= old_entr.getValue();

	if (new_entry_size != old_entry_size)
	{// make room
		arrayCopy(me.data_ + byte_pos + old_entry_size, me.data_ + me.data_size, me.data_ + byte_pos + new_entry_size);
	}

	//set values
	arrayCopy(update_entr.data.value_3, update_entr.data.value_3 + update_entr_size, me.data_ + byte_pos);
	arrayCopy(insert_entr1.data.value_3, insert_entr1.data.value_3 + insert_entr1_size, me.data_ + byte_pos + update_entr_size);
	arrayCopy(insert_entr2.data.value_3, insert_entr2.data.value_3 + insert_entr2_size, me.data_ + byte_pos + update_entr_size + insert_entr1_size);
	return true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned short SIZE, typename TSpec>
inline void
assign(SumList<TValue, MiniSumList<SIZE, TSpec> > & target,
	   SumList<TValue, MiniSumList<SIZE, TSpec> > & source)
{
	target.data_length = source.data_length;
	arrayCopyForward(source.data_, source.data_ + SIZE, target.data_);
	target.data_sum = source.data_sum;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned short SIZE, typename TSpec>
inline typename Size< SumList<TValue, MiniSumList<SIZE, TSpec> > >::Type
length(SumList<TValue, MiniSumList<SIZE, TSpec> > & me)
{
	return me.data_length;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned short SIZE, typename TSpec>
inline typename Value< SumList<TValue, MiniSumList<SIZE, TSpec> > >::Type
sum(SumList<TValue, MiniSumList<SIZE, TSpec> > & me)
{
	return me.data_sum;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned short SIZE, typename TSpec>
inline void
clear(SumList<TValue, MiniSumList<SIZE, TSpec> > & me)
{
	me.data_length = 0;
	me.data_size = 0;
	me.data_sum = 0;
}

//////////////////////////////////////////////////////////////////////////////

//moves the right half of $me$ into $right$
template <typename TValue, unsigned short SIZE, typename TSpec>
inline void
splitSumList(SumList<TValue, MiniSumList<SIZE, TSpec> > & me,
			 SumList<TValue, MiniSumList<SIZE, TSpec> > & right)
{
	typedef SumList<TValue, MiniSumList<SIZE, TSpec> > TSumList;
	typedef typename Iterator<TSumList>::Type TIterator;
	typedef typename Value<TSumList>::Type TValue2;

	clear(right);

	TIterator it(me);
	if (atEnd(it)) return; //nothing to split

	TValue2 sum = getValue(it);

	while (true)
	{
		if (it.data_bytepos >= SIZE/2) break;
		goNext(it);
		if (atEnd(it)) return; //nothing to split
		sum += getValue(it);
	}

	arrayCopyForward(me.data_ + it.data_bytepos, me.data_ + me.data_size, right.data_);
	right.data_length = me.data_length - it.data_position;
	right.data_size = me.data_size - it.data_bytepos;
	right.data_sum = me.data_sum - sum;

	me.data_length = it.data_position;
	me.data_size = it.data_bytepos;
	me.data_sum = sum;
}

//////////////////////////////////////////////////////////////////////////////
// Iterator for MiniSumList

struct MiniSumListIterator;

template <typename TValue, unsigned short SIZE, typename TSpec>
class Iter<SumList<TValue, MiniSumList<SIZE, TSpec> >, MiniSumListIterator>
{
public:
	typedef SumList<TValue, MiniSumList<SIZE, TSpec> > TContainer;
	typedef typename Size<TContainer>::Type TContainerSize;
	typedef typename Position<TContainer>::Type TContainerPosition;

	TContainer & data_container;
	TContainerSize data_bytepos;  //position of entry in container.data_
	TContainerPosition data_position;

	Iter(TContainer & cont, TContainerSize bytepos = 0, TContainerPosition pos = 0)
		: data_container(cont)
		, data_bytepos(bytepos)
		, data_position(pos)
	{
	}
	Iter(Iter const & other)
		: data_container(other.data_container)
		, data_bytepos(other.bytepos)
		, data_position(other.data_position)
	{
	}
	~Iter()
	{
	}
	inline Iter const & 
	operator = (Iter const & other)
	{
		data_container = other.data_container;
		data_bytepos = other.bytepos;
		data_position = other.data_position;
	}
};

template <typename TValue, unsigned short SIZE, typename TSpec, typename TIteratorSpec>
struct Iterator< SumList<TValue, MiniSumList<SIZE, TSpec> >, TIteratorSpec>
{
	typedef Iter< SumList<TValue, MiniSumList<SIZE, TSpec> >, MiniSumListIterator> Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned short SIZE, typename TSpec>
struct Value< Iter< SumList<TValue, MiniSumList<SIZE, TSpec> >, MiniSumListIterator > >:
	Value< SumList<TValue, MiniSumList<SIZE, TSpec> > >
{
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned short SIZE, typename TSpec>
struct Size< Iter< SumList<TValue, MiniSumList<SIZE, TSpec> >, MiniSumListIterator > >:
	Size< SumList<TValue, MiniSumList<SIZE, TSpec> > >
{
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned short SIZE, typename TSpec>
struct Position< Iter< SumList<TValue, MiniSumList<SIZE, TSpec> >, MiniSumListIterator > >:
	Position< SumList<TValue, MiniSumList<SIZE, TSpec> > >
{
};


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned short SIZE, typename TSpec>
inline typename Position< Iter< SumList<TValue, MiniSumList<SIZE, TSpec> >, MiniSumListIterator > >::Type
position(Iter< SumList<TValue, MiniSumList<SIZE, TSpec> >, MiniSumListIterator > & it)
{
	return it.data_position;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned short SIZE, typename TSpec>
inline SumList<TValue, MiniSumList<SIZE, TSpec> > &
container(Iter< SumList<TValue, MiniSumList<SIZE, TSpec> >, MiniSumListIterator > & it)
{
	return * it.data_container;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned short SIZE, typename TSpec>
inline void
goNext(Iter< SumList<TValue, MiniSumList<SIZE, TSpec> >, MiniSumListIterator > & it)
{
	typedef typename SumList<TValue, MiniSumList<SIZE, TSpec> >::Entry TEntry;
	TEntry & entr = * reinterpret_cast<TEntry *>(it.data_container.data_ + it.data_bytepos);
	it.data_bytepos += entr.size();
	++it.data_position;
}

//////////////////////////////////////////////////////////////////////////////

//returns false on overflow, true otherwise
template <typename TValue, unsigned short SIZE, typename TSpec, typename TValue2>
inline bool
assignValue(Iter< SumList<TValue, MiniSumList<SIZE, TSpec> >, MiniSumListIterator > & it,
			TValue2 val)
{
	return _MiniSumList_updateValue(it.data_container, it.data_bytepos, val);
}

//////////////////////////////////////////////////////////////////////////////

//returns false on overflow, true otherwise
template <typename TValue, unsigned short SIZE, typename TSpec, typename TValue2, typename TValue3, typename TValue4>
inline bool
_assignAndInsert(Iter< SumList<TValue, MiniSumList<SIZE, TSpec> >, MiniSumListIterator > & it,
				 TValue2 val,
				 TValue3 insert1,
				 TValue4 insert2)
{
	return _MiniSumList_updateAndInsert(it.data_container, it.data_bytepos, val, insert1, insert2);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned short SIZE, typename TSpec>
inline TValue
getValue(Iter< SumList<TValue, MiniSumList<SIZE, TSpec> >, MiniSumListIterator > & it)
{
	typedef typename SumList<TValue, MiniSumList<SIZE, TSpec> >::Entry TEntry;
	TEntry & entr = * reinterpret_cast<TEntry *>(it.data_container.data_ + it.data_bytepos);
	return entr.getValue();
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned short SIZE, typename TSpec>
inline bool
atEnd(Iter< SumList<TValue, MiniSumList<SIZE, TSpec> >, MiniSumListIterator > & it)
{
	return it.data_bytepos >= it.data_container.data_size;
}

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

//____________________________________________________________________________


//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
