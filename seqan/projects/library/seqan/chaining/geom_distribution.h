/*	
*	Random bit generator
*	
*	from Numerical Recipes in C
*	
*/


#ifndef SEQAN_HEADER_RAND_GEOM
#define SEQAN_HEADER_RAND_GEOM


#define SEQAN_RG_IB1 1
#define SEQAN_RG_IB2 2
#define SEQAN_RG_IB5 16
#define SEQAN_RG_IB18 131072
#define SEQAN_RG_MASK ( SEQAN_RG_IB1 + SEQAN_RG_IB2 + SEQAN_RG_IB5 )

namespace seqan {

	template< typename T > inline
	T
	_geomRand( )
	{
		static unsigned long seed = rand();
		T value = 0;
		while( true )
		{
			if( ( seed & SEQAN_RG_IB18 ) ){
				seed = ( ( seed ^ SEQAN_RG_MASK ) << 1 ) | SEQAN_RG_IB1;
				++value;
			}
			else {
				seed <<= 1;
				break;
			}
		}
		return value;
	}

}

#endif
