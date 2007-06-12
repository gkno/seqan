#ifndef SEQAN_HEADER_FILE_CHAINER_H
#define SEQAN_HEADER_FILE_CHAINER_H

#include <stdio.h>
//#include <seqan/file.h>

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File Formats - FragFile
//////////////////////////////////////////////////////////////////////////////


/**
.Tag.File Format.value.FragFile:
	CHAINER file format for sequences.
*/
struct TagFragFile_;
typedef Tag<TagFragFile_> const FragFile;

struct TagChainStatistics_;
typedef Tag<TagChainStatistics_> const ChainStatistics;

struct TagFragGnuplot_;
typedef Tag<TagFragGnuplot_> const FragGnuplot;

struct TagFragChaFile_;
typedef Tag<TagFragChaFile_> const FragChaFile;

/////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
// read
//////////////////////////////////////////////////////////////////////////////

template < typename TFile, typename TData >
void
read( TFile & file,
		 TData & data,
		 FragFile )
{
SEQAN_CHECKPOINT


	//determine length and end key
	typename Key< TFile >::Type begin_pos = 0;
	typename Size< TData >::Type count = 0;
	typename Size< TData >::Type max_count = 256;
	
	typename Iterator< TData >::Type data_it = end( data );
	typedef typename Value< TData >::Type FragType;
	typename Size< FragType >::Type dim = 0;
	float weight = 0;
	char * line;
	allocate( line, line, max_count );
	char const frag_pattern[] = "%d,%d]";
	
	file.getline( line, max_count );

	char c;
	c = line[ 0 ];

	reserve( data, 300000 );

	if (c == '>')
	{
		sscanf( line, ">CHA %d", &dim );
	}
	else
	{
		return;
	}

	typename Size< FragType >::Type dim_counter = 0;
	
	typename Key< FragType >::Type * left_pos = new typename Key< FragType >::Type[dim];
	typename Key< FragType >::Type * right_pos = new typename Key< FragType >::Type[dim];
	
	while( true )
	{
		file.getline( line, max_count );

		c = line[ 0 ];
		
		if (c == '#')
		{
			sscanf( line, "#%f", &weight );
						
			dim_counter = 0;
			file.getline( line, max_count );
			char * pToken = strtok(line, "[");
			if( pToken )
			{
				sscanf( pToken, frag_pattern, left_pos + dim_counter , right_pos + dim_counter );
				++dim_counter;
				while ( (pToken = strtok(NULL, "[")) && dim_counter < dim )
				{
					sscanf( pToken, frag_pattern, left_pos + dim_counter , right_pos + dim_counter );
					++dim_counter;
				}
			}
			FragType newFrag( left_pos, right_pos, dim, weight );
			appendValue( data, newFrag );
		}
		else
			return;	
	}
}


	template < typename TFile, typename TData >
	void
	writeFragments( TFile & file,
			 TData & data,
			 FragChaFile )
	{
		file << ">CHA " << dimension( value ( begin( data ) ) ) << std::endl;
		typename Iterator< TData >::Type fragIt = begin( data );
		while( fragIt != end( data ) )
		{
			file<< "# "<< static_cast< double >( weight( value (fragIt ) ) ) << std::endl;
			for( typename Size< typename Value< TData >::Type >::Type counter = 0; counter < dimension( value( fragIt ) ); ++counter )
			{
				file << "[" << leftPosition( value( fragIt ), counter ) << "," << rightPosition( value( fragIt ), counter ) << "] ";

			}
			file << std::endl;
			goNext( fragIt );
		}
	}

	template < typename TFile, typename TData >
	void
	writeFragments( TFile & file,
			 TData & data,
			 FragFile )
	{
		SEQAN_CHECK( lenght( data ) != 0 )
		file << ">CHA " << dimension( value ( begin( data ) ) ) << std::endl;
		typename Iterator< TData >::Type fragIt = begin( data );
		int num = 0;
		while( fragIt != end( data ) )
		{
			file<< num << "\t";
			for( typename Size< typename Value< TData >::Type >::Type counter = 0; counter < dimension( value( fragIt ) ); ++counter )
			{
				file	<< leftPosition( value( fragIt ), counter ) << "\t"
						<< rightPosition( value( fragIt ), counter ) << "\t";

			}
			file << weight( value( fragIt ) ) << '\n';
			goNext( fragIt );
			++num;
		}
	}


	template < typename TFile, typename TData >
	void
	writeFragments( TFile & file,
					 TData & data,
					 ChainStatistics )
	{
		file << "# format: id s1 e1 ... sk ek score " << std::endl;
		file << "# id: fragment id  " << std::endl;
		file << "# s1: fragment's start position in genome 1 " << std::endl;
		file << "# e1: fragment's end position in genome 1 " << std::endl;
		file << "# sk: fragment's start position in genome k " << std::endl;
		file << "# ek: fragment's end position in genome k " << std::endl;
		file << "# score: score of the chain ending at this fragment  " << std::endl;
		file << "# num of genomes: " << dimension( value ( begin( data ) ) ) << std::endl;
		file << "# total num. of genomes including origin and terminus: " << length( data ) << std::endl;
		file << "# global chain score: 1111 " << std::endl;
		typename Size< typename Value< TData >::Type >::Type dim = dimension( value ( begin( data ) ) );
		int num = 0;
		typename Iterator< TData >::Type fragIt = begin( data );
		while( fragIt != end( data ) )
		{
			file<< num << "\t";
			for( typename Size< typename Value< TData >::Type >::Type counter = 0; counter < dim; ++counter )
			{
				file	<< leftPosition( value( fragIt ), counter ) << "\t"
						<< rightPosition( value( fragIt ), counter ) << "\t";

			}
			file << weight( value( fragIt ) ) << '\n';
			goNext( fragIt );
			++num;
		}
	}

	template < typename TFile, typename TData >
	void
	writeFragments( TFile & file,
					 TData & data,
					 FragGnuplot )
	{
		typename Iterator< TData >::Type fragIt = begin( data );
		while( fragIt != end( data ) )
		{
			file << weight( value( fragIt ) ) <<' '
			   << leftPosition( value( fragIt ), 0 ) <<' '
			   << rightPosition( value( fragIt ), 0 ) <<' '
			   << leftPosition( value( fragIt ), 1 ) <<' '
			   << rightPosition( value( fragIt ), 1 ) <<' '
			   << '\n'
			  ;
			goNext( fragIt );
		}
		return;
	}  


}

#endif //#ifndef SEQAN_HEADER_...

