#ifndef SEQAN_HEADER_TEST_CHAIN
#define SEQAN_HEADER_TEST_CHAIN


namespace seqan{

	template< typename TContainer, typename TBorder >
	void
	testChain( TContainer & chain )
	{
		typename Iterator< TContainer >::Type it = begin( chain );
		typename Iterator< TContainer >::Type next = begin( chain );
		goNext( next );
		bool b = true;
		while( next != end( chain ) )
		{	
			for( typename Size< typename Position< TContainer >::Type >::Type dim = 0; dim < dimension( *it ); ++dim )
			{
				b = b && ( rightPosition( *it, dim ) < leftPosition( *next, dim ) );
			}
			if( !b )
			{
				std::cout << std::endl;
				dump( *it );
				dump( *next );
			}
			goNext( it );
			goNext( next );
			
		}
	}


}

#endif


