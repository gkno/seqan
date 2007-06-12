#ifndef _RT_MAX_SKIP_BASE_ELEMENT_H
#define _RT_MAX_SKIP_BASE_ELEMENT_H


namespace seqan
{


//___________________________ struct SkipBaseElement< TObject, TModus, RT< TSpec >, Deferred > > _______________________
// Adaption of the struct SkipBaseElement for use in a RangeTree
// Instead saving it's own theKey and value, it has a pointer to the related 
// RTEntry object, which stores the theKey and value
// To get the correct theKey, SkipBaseElement needs information about the dimension
// of the RangeTree it is part of.


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void 
	dump( SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		typename Size< TObject >::Type dim = dimension( *getObject( &me ) );
		typename Size< TObject >::Type counter = 0;
		while( counter < dim )
		{
			if( key( &getObject( &me ), counter ) == infimumValue< typename Key< TObject >::Type >( ) )
					std::cout << std::left << "L";
			else
				std::cout << key( getObject( &me ), counter );
			++counter;
			std::cout << " ";
		}
		std::cout << std::endl;
	}


	
} // namespace SEQAN_NAMESPACE_SKIPLIST
#endif //_RT_SKIP_BASE_ELEMENT_H
