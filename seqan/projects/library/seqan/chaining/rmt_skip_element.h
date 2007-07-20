#ifndef SEQAN_RMT_SKIP_ELEMENT_H
#define SEQAN_RMT_SKIP_ELEMENT_H

namespace seqan
{

// Modifications of the struct SkipElement for use in a SkipListStatic< True, RT< MaxTag > >

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct _RangeMaxCargo
	{

		SkipList< TObject, TModus, TSpec, TStructuring > * _assocStruct;
		TObject * _maxObj;

		_RangeMaxCargo()
			: _assocStruct( NULL )
			, _maxObj( NULL )
		{}

	};

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	SkipList< TObject, TModus, TSpec, TStructuring > *
	_getAssoc( _RangeMaxCargo< TObject, TModus, TSpec, TStructuring > & me )
	{
		return me._assocStruct;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void
	_setAssoc( _RangeMaxCargo< TObject, TModus, TSpec, TStructuring > & me,
				SkipList< TObject, TModus, TSpec, TStructuring > * list )
	{
		me._assocStruct = list;
	}

	
		// handling for max objects

	template< typename TObject, typename TSpec, typename TStructuring > inline
	TObject *
	_getMaxObject( _RangeMaxCargo< TObject, SkipListStatic, TSpec, TStructuring > & me )
	{
		return me._maxObj;
	}

	template< typename TObject, typename TSpec, typename TStructuring > inline
	void
	_setMaxObject( _RangeMaxCargo< TObject, SkipListStatic, TSpec, TStructuring > & me,
				TObject * obj )
	{
		me._maxObj = obj;
	}

	template< typename TObject,typename TSpec, typename TStructuring > inline
	TObject *
	_getMaxObject( SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * me )
	{
		return me->_cargo._maxObj;
	}

	template< typename TObject, typename TSpec, typename TStructuring > inline
	void
	_setMaxObject( SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * me,
					TObject * maxObj )
	{
		me->_cargo._maxObj = maxObj;
	}
	
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct Cargo< SkipElement< TObject, TModus, RT< MaxTree< TSpec > >, TStructuring > >
	{
		typedef _RangeMaxCargo< TObject, TModus, RT< MaxTree< TSpec > >, TStructuring > Type;
	};

		// get the score value of the related object
	template< typename TObject, typename TSpec, typename TStructuring > inline
	typename Weight< TObject >::Type
	weight( SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * me )
	{
		return weight( *me->_cargo._maxObj );
	}

		// get the chain score value of the related object
	template< typename TObject, typename TSpec, typename TStructuring > inline
	typename Weight< TObject >::Type
	priority( SkipElement< TObject, SkipListStatic, RT< MaxTree< TSpec > >, TStructuring > * me )
	{
		return priority( *me->_cargo._maxObj );
	}

} // namespace SEQAN_NAMESPACE_SKIPLIST
#endif // SEQAN_RMT_SKIP_ELEMENT_H
