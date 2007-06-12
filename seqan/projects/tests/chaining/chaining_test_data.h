#ifndef SEQAN_HEADER_FRAGMENT_DATA_H
#define SEQAN_HEADER_FRAGMENT_DATA_H

#include <math.h>

namespace seqan
{

	template< typename TContainer >
	void
	_generateRandomFrags( TContainer & dest,
							int num,
							int min,
							int max,
							int minwidth,
							int maxwidth,
							int dim )
	{
		typedef typename Value< TContainer >::Type FragType;
		typename Key< FragType >::Type * left_pos;
		typename Key< FragType >::Type * right_pos;
		reserve( dest, num );
		double d_dim = static_cast< double >( dim );
		for( int i = 0; i < num; ++i )
		{
			
			double width_sum = 0;
			left_pos = new typename Key< FragType >::Type[ dim ];
			right_pos = new typename Key< FragType >::Type[ dim ];
			for( int d = 0; d < dim; ++d )
			{
				left_pos[ d ] = ( rand() % ( max - min ) ) + min;
				int width = ( rand() % ( maxwidth - minwidth ) )+ minwidth;
				width_sum += width;
				right_pos[ d ] = left_pos[ d ] + width;			
				
			}
			FragType frag( left_pos, right_pos, dim, static_cast< typename Weight< FragType >::Type >(exp(log(width_sum)/d_dim)) * 100 );

			delete[] left_pos;
			delete[] right_pos;

			appendValue( dest, frag );
		}
	}


}

#endif
