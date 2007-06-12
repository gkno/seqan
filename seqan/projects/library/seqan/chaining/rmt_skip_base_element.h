#ifndef SEQAN_RMT_SKIP_BASE_ELEMENT_H
#define SEQAN_RMT_SKIP_BASE_ELEMENT_H


namespace seqan
{

		// get the score value of the related object
	template< typename TObject, typename TSpec, typename TStructuring > inline
	typename Weight< TObject >::Type
	weight( SkipBaseElement< TObject, Static, RT< MaxTree< TSpec > >, TStructuring > * me )
	{
		return weight( *getObject( me ) );
	}

		// get the chain score value of the related object
	template< typename TObject, typename TSpec, typename TStructuring > inline
	typename Weight< TObject >::Type
	priority( SkipBaseElement< TObject, Static, RT< MaxTree< TSpec > >, TStructuring > * me )
	{
		return priority( *getObject( me ) );
	}

} // namespace SEQAN_NAMESPACE_SKIPLIST
#endif //SEQAN_RMT_SKIP_BASE_ELEMENT_H


