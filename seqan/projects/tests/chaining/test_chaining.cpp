
#include <vector>

//#define SEQAN_DEBUG
//#define _SEQAN_CHAIN_DEBUG
//
#define CHAINGG0

#include <seqan/chaining.h>

#include "file_chainer.h"
#include <iostream>
#include <fstream>
#include "chaining_test_data.h"
#include "test_chain.h"
#include "test_ao.h"

#include "clock.h"
#include <time.h>

#include "chain_benchmark.h"

#include <seqan/sequence.h>


using namespace seqan;


int main(int argc, char * argv[])
{
	//int numFrags = 1000;
	//int maxCoord = 2 * numFrags;

	//String< Fragment< int > > frag_set;
	//_generateRandomFrags( frag_set, numFrags, 1, maxCoord, 2, 5, 3 );

	srand( time( 0 ) );

	//std::fstream frag_file2;
	//frag_file2.open( "fragmentfile.mat", std::ios::in );
	//read( frag_file2, frag_set, FragFile() );

	//std::cout << "compute chain g0" << std::endl;

	//String< Fragment< int > > chain;
	//reserve( chain, length( frag_set ) );
	//compute_chain( frag_set, chain, Score< int, Zero >(), Chainer(), Complete() );

	//std::fstream gnu_file;
	//gnu_file.open( "g0.data", std::ios::out );
	//writeFragments( gnu_file, frag_set, FragGnuplot() );
	//gnu_file.close();

	//std::fstream chain_file;
	//chain_file.open( "g0.data.cha.dbf", std::ios::out );
	//writeFragments( chain_file, chain, ChainStatistics() );
	//chain_file.close();

	//std::cout << "compute chain g1" << std::endl;

	//String< Fragment< int > > chain_1;
	//reserve( chain_1, length( frag_set ) );
	//compute_chain( frag_set, chain_1, Score< int, Manhattan >(), Chainer(), Complete() );

	//std::fstream gnu_file_1;
	//gnu_file_1.open( "g1.data", std::ios::out );
	//writeFragments( gnu_file_1, frag_set, FragGnuplot() );
	//gnu_file_1.close();

	//std::fstream chain_file_1;
	//chain_file_1.open( "g1.data.cha.dbf", std::ios::out );
	//writeFragments( chain_file_1, chain_1, ChainStatistics() );
	//chain_file_1.close();

	//std::cout << "compute chain gSoP" << std::endl;

	//String< Fragment< int > > chain_sop;
	//reserve( chain_sop, length( frag_set ) );
	//compute_chain( frag_set, chain_sop, Score< int, ChainSoP >(), Chainer(), Complete() );

	//std::fstream gnu_file_sop;
	//gnu_file_sop.open( "gSoP.data", std::ios::out );
	//writeFragments( gnu_file_sop, frag_set, FragGnuplot() );
	//gnu_file_sop.close();

	//std::fstream chain_file_sop;
	//chain_file_sop.open( "gSoP.data.cha.dbf", std::ios::out );
	//writeFragments( chain_file_sop, chain_sop, ChainStatistics() );
	//chain_file_sop.close();

	//testChainerG0( frag_set, numFrags );

	//testChainerG1( frag_set, numFrags );
	//testChainerGSoP( frag_set, numFrags );
	//std::cout<< std::endl;

	timetest_chaining();

	//testGusfield();

	return 0;
}

