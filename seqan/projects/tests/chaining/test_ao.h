#ifndef SEQAN_HEADER_TEST_AO_CHAINING
#define SEQAN_HEADER_TEST_AO_CHAINING

#include "test_chain.h"

namespace seqan{

	template< typename TSource, typename TCount >
	void
	testChainerG0( TSource & frag_set,
					TCount count )
	{

		String< Fragment< int > > chain;
		reserve( chain, count );

		std::cout << "chaining.." << std::endl;
		compute_chain( frag_set, chain, Score< int, Zero >(), Chainer(), Complete() );
		std::cout << "testing.." << std::endl;

		testChain< String< Fragment< int > >, int >( chain );

		std::fstream gnu_file;
		gnu_file.open( "g0.data", std::ios::out );
		writeFragments( gnu_file, frag_set, FragGnuplot() );
		gnu_file.close();

		std::fstream chain_file;
		chain_file.open( "g0.data.cha.dbf", std::ios::out );
		writeFragments( chain_file, chain, ChainStatistics() );
		chain_file.close();

	}

	template< typename TSource, typename TCount >
	void
	testChainerG1( TSource & frag_set,
					TCount count )
	{

		String< Fragment< int > > chain;
		reserve( chain, count );

		compute_chain( frag_set, chain, Score< int, Manhattan >(), Chainer(), Complete() );

		testChain< String< Fragment< int > >, int >( chain );

		std::fstream gnu_file;
		gnu_file.open( "g1.data", std::ios::out );
		writeFragments( gnu_file, frag_set, FragGnuplot() );
		gnu_file.close();

		std::fstream chain_file;
		chain_file.open( "g1.data.cha.dbf", std::ios::out );
		writeFragments( chain_file, chain, ChainStatistics() );
		chain_file.close();

	}


	template< typename TSource, typename TCount >
	void
	testChainerGSoP( TSource & frag_set,
						TCount count )
	{

		String< Fragment< int > > chain;
		reserve( chain, count );

		compute_chain( frag_set, chain, Score< int, ChainSoP >(), Chainer(), Complete() );

		testChain< String< Fragment< int > >, int >( chain );

		std::fstream gnu_file;
		gnu_file.open( "gSoP.data", std::ios::out );
		writeFragments( gnu_file, frag_set, FragGnuplot() );
		gnu_file.close();

		std::fstream chain_file;
		chain_file.open( "gSoP.data.cha.dbf", std::ios::out );
		writeFragments( chain_file, chain, ChainStatistics() );
		chain_file.close();
	}
}

#endif
