// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Andres Gogol-Döring <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Implementation of the Holder class and specializations.
// ==========================================================================

#ifndef SEQAN_BASIC_BASIC_HOLDER_H_
#define SEQAN_BASIC_BASIC_HOLDER_H_

// TODO(holtgrew): What about const holders?
// TODO(holtgrew): Are holders on pointers used anywhere?
// TODO(holtgrew): Split into holder_base.h, holder_simple.h, and holder_tristate.h?

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename T> struct IsSimple;

// Used in functions section.
// TODO(holtgrew): This will go away, right?
template <typename TValue> inline size_t length(TValue const * me);
template <typename TValue> inline size_t length(TValue * me);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Holder
// ----------------------------------------------------------------------------

/**
.Class.Holder:
..cat:Basic
..summary:Manages relationship to another object.
..signature:Holder<TValue, TSpec>
..param.TValue:Type of the managed object.
...metafunction:Metafunction.Value
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:$Tristate$
..remarks.text:The main purpose of this class is to facilitate the handling of
member objects. If we want class $A$ to be dependent on or the owner of another object of class $B$, 
then we add a data member of type $Holder<B>$ to $A$. 
$Holder$ offers some useful access functions and stores the kind of relationship between $A$ and $B$.
..include:seqan/basic.h

.Memfunc.Holder:
..class:Class.Holder
..summary:Constructor
..signature:Holder<TValue, TSpec>()
..signature:Holder<TValue, TSpec>(holder)
..signature:Holder<TValue, TSpec>(value)
..param.holder:Another holder object.
..param.value:An object of type $TValue$.
..remarks.text:
The default constructor creates a holder that is in state 'empty'.
If a $value$ is passed to the constructor, the holder will be in state 'dependent'.
 */

// Default for class holder.
struct Tristate_;
typedef Tag<Tristate_> Tristate;

template <typename TValue, typename TSpec = Tristate>
struct Holder;

// ----------------------------------------------------------------------------
// Specialization Tristate Holder
// ----------------------------------------------------------------------------

/**
.Spec.Tristate Holder
..cat:Holders
..summary:Holder that can be empty, dependent, or owner.
..signature:Holder<TValue, Tristate>
..param.TValue:Type of the managed object.
..general:Class.Holder
..remarks.text:A tristate holder $A$ that holds an object $B$ has one of the following states:
..remarks.text:- owner: $A$ is the owner of $B$. If $A$ is destroyed, $B$ will be destroyed automatically.
..remarks.text:- dependent: $A$ depends on $B$. $B$ should not be destroyed as long as $A$ is used.
..remarks.text:- empty: there is currently no object reference stored in the holder $A$.
..remarks.text:The state of the holder can be determined by @Function.empty@ and @Function.dependent@.
..remarks.text:If a holder object is in owner state when destructed, the owned object is destructed as well.
..include:seqan/basic.h
 */

template <typename TValue>
struct Holder<TValue, Tristate>
{
	enum EHolderState
	{
		EMPTY = 0,
		OWNER = 1,
		DEPENDENT = ~0
	};

	typedef typename Value<Holder>::Type THostValue;

    // ------------------------------------------------------------------------
    // Members
    // ------------------------------------------------------------------------

	THostValue * data_value;
	EHolderState data_state;

    // ------------------------------------------------------------------------
    // Constructors; Destructor
    // ------------------------------------------------------------------------

	Holder() : data_state(EMPTY)
	{
        SEQAN_CHECKPOINT;
	}

	Holder(Holder const & source_) : data_state(EMPTY)
	{
        SEQAN_CHECKPOINT;
		assign(*this, source_);
	}

	Holder(THostValue & value_) : data_state(EMPTY)
	{
        SEQAN_CHECKPOINT;
		setValue(*this, value_);
	}

	Holder(THostValue const & value_) : data_state(EMPTY)
	{
        SEQAN_CHECKPOINT;
		assignValue(*this, value_);
	}

	~Holder()
	{
        SEQAN_CHECKPOINT;
		clear(*this);
	}

    // ------------------------------------------------------------------------
    // Assignment Operators;  Must be defined in class.
    // ------------------------------------------------------------------------

	inline Holder const &
	operator=(Holder const & source_)
	{
        SEQAN_CHECKPOINT;
		assign(*this, source_);
		return *this;
	}

	inline Holder const &
	operator=(THostValue const & value_)
	{
        SEQAN_CHECKPOINT;
		assignValue(*this, value_);
		return *this;
	}

    // ------------------------------------------------------------------------
    // Conversion Operators;  Must be defined in class.
    // ------------------------------------------------------------------------

	inline operator THostValue()
	{
        SEQAN_CHECKPOINT;
		return _dataValue(*this);
	}
};

// TODO(holtgrew): Why does this not reliably work without copying where it can?
template <typename TValue>
struct Holder<TValue const, Tristate>
{
	enum EHolderState
	{
		EMPTY = 0,
		OWNER = 1,
		DEPENDENT = ~0
	};

	typedef typename Value<Holder>::Type THostValue;

    // ------------------------------------------------------------------------
    // Members
    // ------------------------------------------------------------------------

	THostValue * data_value;
	EHolderState data_state;

    // ------------------------------------------------------------------------
    // Constructors; Destructor
    // ------------------------------------------------------------------------

	Holder() : data_state(EMPTY)
	{
        SEQAN_CHECKPOINT;
	}

	Holder(Holder const & source_) : data_state(EMPTY)
	{
        SEQAN_CHECKPOINT;
		assign(*this, source_);
	}

	Holder(TValue & value_) : data_state(EMPTY)
	{
        SEQAN_CHECKPOINT;
		setValue(*this, value_);
	}

	Holder(TValue const & value_): data_state(EMPTY)
	{
        SEQAN_CHECKPOINT;
		setValue(*this, value_);
	}

	~Holder()
	{
        SEQAN_CHECKPOINT;
		clear(*this);
	}

    // ------------------------------------------------------------------------
    // Assignment Operators;  Has to be defined in class.
    // ------------------------------------------------------------------------

	inline Holder &
	operator=(Holder const & source_)
	{
        SEQAN_CHECKPOINT;
		assign(*this, source_);
		return *this;
	}

	inline Holder &
	operator=(THostValue const & value_)
	{
        SEQAN_CHECKPOINT;
		assignValue(*this, value_);
		return *this;
	}

    // ------------------------------------------------------------------------
    // Conversion Operators;  Has to be defined in class.
    // ------------------------------------------------------------------------

	inline operator THostValue()
	{
        SEQAN_CHECKPOINT;
		return _dataValue(*this);
	}
};

template <typename TValue>
struct Holder<TValue *, Tristate>
{

	enum EHolderState
	{
		EMPTY = 0,
		OWNER = 1,
		DEPENDENT = ~0
	};

	typedef typename Value<Holder>::Type THostValue;

    // ------------------------------------------------------------------------
    // Members
    // ------------------------------------------------------------------------

	THostValue data_value;
	EHolderState data_state;

    // ------------------------------------------------------------------------
    // Constructor; Destructors
    // ------------------------------------------------------------------------

	Holder() : data_state(EMPTY)
	{
        SEQAN_CHECKPOINT;
	}

	Holder(Holder const & source_) : data_state(EMPTY)
	{
        SEQAN_CHECKPOINT;
		assign(*this, source_);
	}

	Holder(TValue * value_) : data_state(EMPTY)
	{
        SEQAN_CHECKPOINT;
		setValue(*this, value_);
	}

	~Holder()
	{
        SEQAN_CHECKPOINT;
		clear(*this);
	}

    // ------------------------------------------------------------------------
    // Assignment Operators;  Have to be defined in class.
    // ------------------------------------------------------------------------

	inline Holder const &
	operator = (Holder const & source_)
	{
        SEQAN_CHECKPOINT;
		assign(*this, source_);
		return *this;
	}

	inline Holder const &
	operator = (THostValue value_)
	{
        SEQAN_CHECKPOINT;
		setValue(*this, value_);
		return *this;
	}

    // ------------------------------------------------------------------------
    // Conversion Operators;  Have to be defined in class.
    // ------------------------------------------------------------------------

	inline operator THostValue()
	{
        SEQAN_CHECKPOINT;
		return _dataValue(*this);
	}
};

template <typename TValue>
struct Holder<TValue * const, Tristate>
{
	enum EHolderState
	{
		EMPTY = 0,
		OWNER = 1,
		DEPENDENT = ~0
	};

	typedef typename Value<Holder>::Type THostValue;

    // ------------------------------------------------------------------------
    // Members
    // ------------------------------------------------------------------------

	THostValue data_value;
	EHolderState data_state;

    // ------------------------------------------------------------------------
    // Constructors; Destructor
    // ------------------------------------------------------------------------

	Holder() : data_state(EMPTY)
	{
        SEQAN_CHECKPOINT;
	}

	Holder(Holder const & source_) : data_state(EMPTY)
	{
        SEQAN_CHECKPOINT;
		assign(*this, source_);
	}

	Holder(TValue * value_) : data_state(EMPTY)
	{
        SEQAN_CHECKPOINT;
		setValue(*this, value_);
	}

	~Holder()
	{
        SEQAN_CHECKPOINT;
		clear(*this);
	}

    // ------------------------------------------------------------------------
    // Assignment Operators;  Have to be defined in class.
    // ------------------------------------------------------------------------

	inline Holder &
	operator=(Holder const & source_)
	{
        SEQAN_CHECKPOINT;
		assign(*this, source_);
		return *this;
	}

	inline Holder &
	operator=(THostValue value_)
	{
        SEQAN_CHECKPOINT;
		setValue(*this, value_);
		return *this;
	}

    // ------------------------------------------------------------------------
    // Conversion Operators;  Have to be defined in class.
    // ------------------------------------------------------------------------

	inline operator THostValue()
	{
        SEQAN_CHECKPOINT;
		return _dataValue(*this);
	}
};

// ----------------------------------------------------------------------------
// Specialization Simple Holder
// ----------------------------------------------------------------------------

#ifdef PLATFORM_WINDOWS_VS
// Disable warning C4521 locally (multiple copy constructors).
#pragma warning( push )
#pragma warning( disable: 4521 )
#endif  // PLATFORM_WINDOWS_VS

// TODO(holtgrew): Use struct Simple_; typedef Tag<Simple_> Simple?
struct Simple;

/**
.Spec.Simple Holder
..cat:Holders
..summary:Simple copying wrapper without any additional state.
..signature:Holder<TValue, Simle>
..param.TValue:Type of the managed object.
..general:Class.Holder
..remarks.text:This holder stores a copy of the value.
..include:seqan/basic.h
 */

template <typename TValue>
struct Holder<TValue, Simple>
{
	typedef typename Value<Holder>::Type THolderValue;
	typedef typename Parameter_<THolderValue>::Type THolderParameter;

    // ------------------------------------------------------------------------
    // Members
    // ------------------------------------------------------------------------

	mutable typename RemoveConst_<THolderValue>::Type data_value;

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

	Holder() {}

	Holder(Holder & source_) : data_value(source_.data_value)
	{
        SEQAN_CHECKPOINT;
	}

	Holder(Holder const & source_) : data_value(source_.data_value)
	{
        SEQAN_CHECKPOINT;
	}

	template <typename TSource>
	Holder(TSource & value_) : data_value(value_)
	{
        SEQAN_CHECKPOINT;
	}

	template <typename TSource>
	Holder(TSource const & value_) : data_value(value_)
	{
        SEQAN_CHECKPOINT;
	}

    // ------------------------------------------------------------------------
    // Assignment Operators;  Have to be defined in class.
    // ------------------------------------------------------------------------

	Holder &
	operator=(Holder const & source_)
	{
        SEQAN_CHECKPOINT;
		data_value = source_.data_value;
		return *this;
	}

	Holder &
	operator=(THolderValue const & value_)
	{
        SEQAN_CHECKPOINT;
		data_value = value_;
		return *this;
	}

    // ------------------------------------------------------------------------
    // Conversion Operators;  Have to be defined in class.
    // ------------------------------------------------------------------------

	operator THolderParameter()
	{
        SEQAN_CHECKPOINT;
		return *data_value;
	}
};

#ifdef PLATFORM_WINDOWS_VS
// Set old warning C4521 state again (multiple copy constructors).
#pragma warning( pop )
#endif  // PLATFORM_WINDOWS_VS

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

///.Metafunction.Value.param.T.type:Class.Holder

template <typename TValue, typename TSpec>
struct Value<Holder<TValue, TSpec> >
{
	typedef TValue Type;
};

template <typename TValue, typename TSpec>
struct Value<Holder<TValue, TSpec> const>
{
	typedef TValue Type;
};

// TODO(holtgrew): What about holders on pointers?
template <typename TValue, typename TSpec>
struct Value<Holder<TValue * const, TSpec> >
{
	typedef TValue * Type;
};

// ----------------------------------------------------------------------------
// Metafunction Spec
// ----------------------------------------------------------------------------

///.Metafunction.Spec.param.T.type:Class.Holder

template <typename TValue, typename TSpec>
struct Spec<Holder<TValue, TSpec> >
{
	typedef TSpec Type;
};

template <typename TValue, typename TSpec>
struct Spec<Holder<TValue, TSpec> const>
{
	typedef TSpec Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

///.Metafunction.Reference.param.T.type:Class.Holder

template <typename TValue, typename TSpec>
struct Reference<Holder<TValue, TSpec> >
{
	typedef typename Value<Holder<TValue, TSpec> >::Type & Type;
};

template <typename TValue, typename TSpec>
struct Reference< Holder<TValue, TSpec> const>
{
	typedef typename Value<Holder<TValue, TSpec> const>::Type & Type;
};

template <typename TValue, typename TSpec>
struct Reference<Holder<TValue *, TSpec> const>
{
	typedef typename Value<Holder<TValue *, TSpec> const>::Type const & Type;
};

// ============================================================================
// Functions
// ============================================================================

// ****************************************************************************
// Specialization Tristate Holder
// ****************************************************************************

// ----------------------------------------------------------------------------
// Function _dataValue()
// ----------------------------------------------------------------------------

template <typename TValue>
inline typename Reference<Holder<TValue, Tristate> >::Type
_dataValue(Holder<TValue, Tristate> & me)
{
	return * me.data_value;
}

template <typename TValue>
inline typename Reference<Holder<TValue, Tristate> const>::Type
_dataValue(Holder<TValue, Tristate> const & me)
{
	return * me.data_value;
}

template <typename TValue>
inline typename Reference<Holder<TValue *, Tristate> >::Type
_dataValue(Holder<TValue *, Tristate> & me)
{
	return me.data_value;
}

template <typename TValue>
inline typename Reference<Holder<TValue *, Tristate> const>::Type
_dataValue(Holder<TValue *, Tristate> const & me)
{
	return me.data_value;
}

template <typename TValue>
inline typename Reference<Holder<TValue * const, Tristate> >::Type
_dataValue(Holder<TValue * const, Tristate> & me)
{
	return me.data_value;
}

template <typename TValue>
inline typename Reference<Holder<TValue * const, Tristate> const>::Type
_dataValue(Holder<TValue * const, Tristate> const & me)
{
	return me.data_value;
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

///.Function.empty.param.object.type:Class.Holder

template <typename TValue>
inline bool
empty(Holder<TValue, Tristate> const & me)
{
SEQAN_CHECKPOINT
	return (me.data_state == Holder<TValue, Tristate>::EMPTY);
}

// ----------------------------------------------------------------------------
// Function dependent()
// ----------------------------------------------------------------------------

///.Function.dependent.param.object.type:Class.Holder

template <typename TValue>
inline bool
dependent(Holder<TValue, Tristate> const & me)
{
    SEQAN_CHECKPOINT;
	return me.data_state == Holder<TValue, Tristate>::DEPENDENT;
}

/// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename THolder, typename TValue>
inline void
_holderDeallocate(THolder & me, TValue const & data)
{	
	valueDestruct(& data);
	deallocate(me, & data, 1);
}

template <typename THolder, typename TValue>
inline void
_holderDeallocate(THolder & me, TValue * data, True)				// is a pointer to an *array* of objects
{
	size_t len = length(data)+1;
	arrayDestruct(data, data+len);
	deallocate(me, data, len);
}

template <typename THolder, typename TValue>
inline void
_holderDeallocate(THolder & /*me*/, TValue * /*data*/, False)				// is a pointer to *one* object
{
}

template <typename THolder, typename TValue>
inline void
_holderDeallocate(THolder & me, TValue * data)
{
	return _holderDeallocate(me, data, IsSimple<TValue>());			// try to distinguish between a pointer to one/array of object(s)
}

///.Function.clear.param.object.type:Class.Holder
///.Function.clear.remarks.text:If $clear$ is applied on a @Class.Holder@ object,
///the state of this object is set to 'empty'.

template <typename TValue>
inline void
clear(Holder<TValue, Tristate> & me)
{
	switch (me.data_state) {
        case Holder<TValue, Tristate>::EMPTY:
            break;
        case Holder<TValue, Tristate>::DEPENDENT:
            SEQAN_CHECKPOINT;
            me.data_state = Holder<TValue, Tristate>::EMPTY;
            break;
        default:  // case Holder<TValue, TSpec>::OWNER
            SEQAN_CHECKPOINT;
            _holderDeallocate(me, _dataValue(me));
            me.data_state = Holder<TValue, Tristate>::EMPTY;
            break;
	}
}

// ----------------------------------------------------------------------------
// Function create()
// ----------------------------------------------------------------------------

template <typename THolder, typename TValue>
inline typename Value<THolder, 0>::Type *
_holderAllocateObject(THolder & me, TValue const & data)
{	
	typename Value<THolder>::Type * ret;
	allocate(me, ret, 1);
	valueConstruct(ret, data);
	return ret;
}

template <typename THolder, typename TValue>
inline typename Value<THolder, 0>::Type
_holderAllocatePointer(THolder & me, TValue * data, True)			// is a pointer to an *array* of objects
{
	typename Value<THolder>::Type ret;
	size_t len = length(data)+1;
	allocate(me, ret, len);
	arrayConstructCopy(data, data + len, ret);
	return ret;
}

template <typename THolder, typename TValue>
inline typename Value<THolder, 0>::Type
_holderAllocatePointer(THolder & /*me*/, TValue * data, False)			// is a pointer to *one* object
{
	return data;
}

template <typename THolder, typename TValue>
inline typename Value<THolder, 0>::Type
_holderAllocatePointer(THolder & me, TValue * data)
{
	return _holderAllocatePointer(me, data, IsSimple<TValue>());	// try to distinguish between a pointer to one/array of object(s)
}

/**
.Function.create:
..summary:Makes an object to owner of its content.
..cat:Dependent Objects
..signature:create(holder [, object])
..param.holder:A holder object.
...type:Class.Holder
..param.object:Object from which a copy is made and stored in $holder$. (optional)
...type:Metafunction.Value.Value<Holder>::Type
..remarks.text:After this operation, $holder$ will be in state 'owner'.
If $object$ is specified, $holder$ will hold a copy of $object$ at the end of this function.
If $object$ is not specified, the action depends on the former state of $holder$:
..remarks.text:- If the state of $holder$ was 'empty', a new object is default constructed and stored into $holder$.
..remarks.text:- If the state of $holder$ was 'dependent', a copy of the former object is made and stored into $holder$. 
..remarks.text:- If the state of $holder$ was already 'owner', nothing happens.
..see:Class.Holder
..include:seqan/basic.h
*/

template <typename TValue>
inline void
create(Holder<TValue, Tristate> & me)
{
	typedef Holder<TValue, Tristate> THolder;

	switch (me.data_state)
	{
        case Holder<TValue, Tristate>::EMPTY:
            SEQAN_CHECKPOINT;
            allocate(me, me.data_value, 1);
            valueConstruct(me.data_value);
            me.data_state = THolder::OWNER;
            break;
            
        case THolder::DEPENDENT:
            SEQAN_CHECKPOINT;
            create(me, _dataValue(me));
            break;
        default:;
	}
}


template <typename TValue>
inline void
create(Holder<TValue *, Tristate> & me)
{
	typedef Holder<TValue *, Tristate> THolder;

	switch (me.data_state) {
        case Holder<TValue *, Tristate>::EMPTY:
            SEQAN_CHECKPOINT;
            valueConstruct(& me.data_value);
            me.data_state = THolder::OWNER;
            break;
        case THolder::DEPENDENT:
            SEQAN_CHECKPOINT;
            create(me, _dataValue(me));
            break;
        default:;
	}
}

template <typename TValue>
inline void
create(Holder<TValue * const, Tristate> & me)
{
	typedef Holder<TValue *, Tristate> THolder;

	switch (me.data_state) {
        case Holder<TValue *, Tristate>::EMPTY:
            SEQAN_CHECKPOINT;
            valueConstruct(& me.data_value);
            me.data_state = THolder::OWNER;
            break;
        case THolder::DEPENDENT:
            SEQAN_CHECKPOINT;
            create(me, _dataValue(me));
            break;
        default:;
	}
}

template <typename TValue, typename TValue2>
inline void
create(Holder<TValue, Tristate> & me,
	   TValue2 & value_)
{
    SEQAN_CHECKPOINT;

	if (me.data_state == Holder<TValue, Tristate>::OWNER) {
		assign(_dataValue(me), value_);
		return;
	}

	clear(me);
	me.data_value = _holderAllocateObject(me, value_);
	me.data_state = Holder<TValue, Tristate>::OWNER;
}

template <typename TValue, typename TValue2>
inline void
create(Holder<TValue const, Tristate> & me,
	   TValue2 & value_)
{
    SEQAN_CHECKPOINT;

	clear(me);
	me.data_value = _holderAllocateObject(me, value_);
	me.data_state = Holder<TValue const, Tristate>::OWNER;
}

template <typename TValue, typename TValue2>
inline void
create(Holder<TValue *, Tristate> & me,
	   TValue2 & value_)
{
SEQAN_CHECKPOINT

	clear(me);
	me.data_value = _holderAllocatePointer(me, value_);
	me.data_state = Holder<TValue *, Tristate>::OWNER;
}

template <typename TValue, typename TValue2>
inline void
create(Holder<TValue * const, Tristate> & me,
	   TValue2 & value_)
{
SEQAN_CHECKPOINT

	clear(me);
	me.data_value = _holderAllocatePointer(me, value_);
	me.data_state = Holder<TValue *, Tristate>::OWNER;
}


//////////////////////////////////////////////////////////////////////////////
/**
.Function.detach:
..summary:Makes an object independent from other objects.
..cat:Dependent Objects
..signature:detach(object)
..param.object:An object.
...type:Class.Holder
..remarks:
After this function, $object$ does not depends from any other entity outside of $object$,
like a @Function.source@ or a @Function.host@, and @Function.dependent.dependent(object)@ returns $false$ 
..see:Function.source
..see:Function.host
..see:Function.createSource
..see:Function.create
..include:seqan/basic.h
*/

template <typename TValue>
inline void
detach(Holder<TValue, Tristate> & me)
{
    SEQAN_CHECKPOINT;
	create(me);
}

// ----------------------------------------------------------------------------
// Function setValue()
// ----------------------------------------------------------------------------

/**
.Function.setValue:
..cat:Content Manipulation
..summary:Makes holder dependent.
..signature:setValue(holder, object)
..param.holder:A holder object.
...type:Class.Holder
..param.object:Object from which $holder$ will be dependent.
...type:Metafunction.Value.Value<Holder>::Type
..remarks.text:After this operation, $holder$ will be dependent in state 'dependent'.
..see:Class.Holder
..include:seqan/basic.h
*/

template <typename TValue>
inline void
setValue(Holder<TValue, Tristate> & me,
		 TValue & value_)
{
    SEQAN_CHECKPOINT;
	clear(me);
	me.data_value = & value_;
	me.data_state = Holder<TValue, Tristate>::DEPENDENT;
}

template <typename TValue>
inline void
setValue(Holder<TValue const, Tristate> & me,
		 TValue & value_)
{
    SEQAN_CHECKPOINT;
	clear(me);
	me.data_value = & value_;
	me.data_state = Holder<TValue const, Tristate>::DEPENDENT;
}

// ----------------------------------------------------------------------------
// Function setValue()
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
setValue(Holder<TValue *, Tristate> & me,
		 TValue * & value_)
{
    SEQAN_CHECKPOINT;
	clear(me);
	me.data_value = value_;
	me.data_state = Holder<TValue *, Tristate>::DEPENDENT;
}

template <typename TValue>
inline void
setValue(Holder<TValue *, Tristate> & me,
		 TValue * const & value_)
{
    SEQAN_CHECKPOINT;
	clear(me);
	me.data_value = value_;
	me.data_state = Holder<TValue *, Tristate>::DEPENDENT;
}

template <typename TValue>
inline void
setValue(Holder<TValue * const, Tristate> & me,
		 TValue * & value_)
{
    SEQAN_CHECKPOINT;
	clear(me);
	me.data_value = value_;
	me.data_state = Holder<TValue *, Tristate>::DEPENDENT;
}

template <typename TValue>
inline void
setValue(Holder<TValue * const, Tristate> & me,
		 TValue * const & value_)
{
    SEQAN_CHECKPOINT;
	clear(me);
	me.data_value = value_;
	me.data_state = Holder<TValue *, Tristate>::DEPENDENT;
}

template <typename TValue, size_t I>
inline void
setValue(Holder<TValue *, Tristate> & me,
		 TValue (& value_)[I])
{
    SEQAN_CHECKPOINT;
	clear(me);
	me.data_value = value_;
	me.data_state = Holder<TValue *, Tristate>::DEPENDENT;
}

template <typename TValue, size_t I>
inline void
setValue(Holder<TValue *, Tristate> & me,
		 TValue const (& value_)[I])
{
    SEQAN_CHECKPOINT;
	clear(me);
	me.data_value = value_;
	me.data_state = Holder<TValue *, Tristate>::DEPENDENT;
}

template <typename TValue, size_t I>
inline void
setValue(Holder<TValue * const, Tristate> & me,
		 TValue (& value_)[I])
{
    SEQAN_CHECKPOINT;
	clear(me);
	me.data_value = value_;
	me.data_state = Holder<TValue *, Tristate>::DEPENDENT;
}

template <typename TValue, size_t I>
inline void
setValue(Holder<TValue * const, Tristate> & me,
		 TValue const (& value_)[I])
{
    SEQAN_CHECKPOINT;
	clear(me);
	me.data_value = value_;
	me.data_state = Holder<TValue *, Tristate>::DEPENDENT;
}

template <typename TValue, typename TValue2>
inline void
setValue(Holder<TValue, Tristate> & me,
		 TValue2 & value_)
{
    SEQAN_CHECKPOINT;
	set(value(me), value_);
}

template <typename TValue, typename TValue2>
inline void
setValue(Holder<TValue, Tristate> & me,
		 TValue2 const & value_)
{
    SEQAN_CHECKPOINT;
	set(value(me), value_);
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

///.Function.value.param.object.type:Class.Holder

template <typename TValue>
inline typename Reference<Holder<TValue, Tristate> >::Type
value(Holder<TValue, Tristate> & me)
{
    SEQAN_CHECKPOINT;
	typedef Holder<TValue, Tristate> THolder;

	if (empty(me))
	{
		create(me);
	}

	return _dataValue(me);
}

template <typename TValue>
inline typename Reference<Holder<TValue, Tristate> const>::Type
value(Holder<TValue, Tristate> const & me)
{
    SEQAN_CHECKPOINT;
	SEQAN_ASSERT_NOT(empty(me));

	return _dataValue(me);
}

// ----------------------------------------------------------------------------
// Function assignValue()
// ----------------------------------------------------------------------------

///.Function.assignValue.param.object.type:Class.Holder

template <typename TValue, typename TSource>
inline void
assignValue(Holder<TValue, Tristate> & me,
			TSource const & value_)
{
    SEQAN_CHECKPOINT;
	typedef typename Value<Holder<TValue, Tristate> >::Type THostValue;
	if (empty(me)) {
		create(me, value_);
	} else {
		assign(_dataValue(me), value_);
	}
}

// ----------------------------------------------------------------------------
// Function moveValue()
// ----------------------------------------------------------------------------

///.Function.moveValue.param.object.type:Class.Holder

template <typename TValue, typename TSource>
inline void
moveValue(Holder<TValue, Tristate> & me,
		  TSource const & value_)
{
    SEQAN_CHECKPOINT;
	if (empty(me)) {
		create(me, value_);
	} else {
		move(_dataValue(me), value_);
	}
}

// ----------------------------------------------------------------------------
// Function assign()
// ----------------------------------------------------------------------------

///.Function.assign.param.target.type:Class.Holder
///.Function.assign.param.source.type:Class.Holder

template <typename TValue>
inline void
assign(Holder<TValue, Tristate> & target_,
	   Holder<TValue, Tristate> const & source_)
{
    SEQAN_CHECKPOINT;
	switch(source_.data_state) {
        case Holder<TValue, Tristate>::EMPTY:
            clear(target_);
            break;
            
        case Holder<TValue, Tristate>::OWNER:
            assignValue(target_, value(source_));
            break;
            
        default:  // case Holder<TValue, Tristate>::DEPENDENT
            setValue(target_, value(source_));
            break;
	}
}

// ****************************************************************************
// Specialization Simple Holder
// ****************************************************************************

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

template <typename TValue>
inline bool
empty(Holder<TValue, Simple> const & /*me*/)
{
    SEQAN_CHECKPOINT;
	return false;
}

// ----------------------------------------------------------------------------
// Function dependent()
// ----------------------------------------------------------------------------

template <typename TValue>
inline bool
dependent(Holder<TValue, Simple> const & /*me*/)
{
    SEQAN_CHECKPOINT;
	return false;
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
clear(Holder<TValue, Simple> & /*me*/)
{
    SEQAN_CHECKPOINT;
}

// ----------------------------------------------------------------------------
// Function create()
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
create(Holder<TValue, Simple> & /*me*/)
{
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): Should be create(me.data_value), right?
}

template <typename TValue>
inline void
create(Holder<TValue, Simple> & me,
	   TValue const & value_)
{
    SEQAN_CHECKPOINT;
	me.data_value = value_;
}

// ----------------------------------------------------------------------------
// Function detach()
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
detach(Holder<TValue, Simple> & /*me*/)
{
    SEQAN_CHECKPOINT;
}

// ----------------------------------------------------------------------------
// Function setValue()
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
setValue(Holder<TValue, Simple> & me,
		 TValue const & value_)
{
    SEQAN_CHECKPOINT;
	me.data_value = value_;
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TValue>
inline typename Reference<Holder<TValue, Simple> >::Type
value(Holder<TValue, Simple> & me)
{
    SEQAN_CHECKPOINT;
	return me.data_value;
}

template <typename TValue>
inline typename Reference<Holder<TValue, Simple> const>::Type
value(Holder<TValue, Simple> const & me)
{
    SEQAN_CHECKPOINT;
	return me.data_value;
}

// ----------------------------------------------------------------------------
// Function assignValue()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSource>
inline void
assignValue(Holder<TValue, Simple> & me,
			TSource const & value_)
{
    SEQAN_CHECKPOINT;
	assignValue(me.data_value, value_);
}

// ----------------------------------------------------------------------------
// Function moveValue()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSource>
inline void
moveValue(Holder<TValue, Simple> & me,
		  TSource const & value_)
{
    SEQAN_CHECKPOINT;
	move(me.data_value, value_);
}

// ----------------------------------------------------------------------------
// Function assign()
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
assign(Holder<TValue, Simple> & target_,
	   Holder<TValue, Simple> const & source_)
{
    SEQAN_CHECKPOINT;
	assignValue(target_, source_);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_BASIC_HOLDER_H_
