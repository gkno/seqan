#ifndef _SEQAN_RT_BENCHMARK
#define _SEQAN_RT_BENCHMARK


namespace seqan{
	
	int dimensions[] = { 3, 3, 4, 5, 10 };

	template< typename TSource, typename TStructuring >
	void
	benchChainerG0( TSource & frag_set,
					TSource & chain,
					TStructuring structuring,
					Clock & c )
	{
		c.Start();
		volatile int i = compute_chain( frag_set, chain, Score< int, Zero >(), Chainer(), structuring );
		c.Stop();
		std::cout << i << std::endl;
	}

	template< typename TSource, typename TStructuring >
	void
	benchChainerG1( TSource & frag_set,
				 TSource & chain,
				 TStructuring structuring,
					Clock & c )
	{
		c.Start();
		volatile int i = compute_chain( frag_set, chain, Score< int, Manhattan >(), Chainer(), structuring );
		c.Stop();
		std::cout << i << std::endl;
	}


	template< typename TSource, typename TStructuring >
	void
	benchChainerGSoP( TSource & frag_set,
						 TSource & chain,
						 TStructuring structuring,
							Clock & c )
	{
		c.Start();
		volatile int i = compute_chain( frag_set, chain, Score< int, ChainSoP >(), Chainer(), structuring );
		c.Stop();
		std::cout << i << std::endl;
	}


	void chain_benchmarkG0( int elements, 
							int dim, 
							int iterations_bm, 
							std::fstream & def_file,
							std::fstream & semi_file,
							std::fstream & compl_file )
//-------------------------------------- benchmarking, comparison to stl map ----------------------------------
	{
		
		
		std::ostringstream s;
		Clock c1;
		Clock c2;
		Clock c3;
		
		for( int times = 0; times < iterations_bm; ++ times )
		{
			//struct tm * act_time_format;
			//time_t act_time_simple;
			//time( &act_time_simple );
			//act_time_format = localtime( &act_time_simple );

			//char * time_string = new char[10];
			//sprintf( time_string, "%2i%2i_%2i%2i%2i", act_time_format->tm_mday,
			//										act_time_format->tm_mon,
			//										act_time_format->tm_hour,
			//										act_time_format->tm_min,
			//										act_time_format->tm_sec );
			//std::cout << time_string << std::endl;
			//std::ostringstream file_name;
			//file_name <<  "GSoP" << time_string;

			//std::fstream frag_file;
			//std::ostringstream frag_file_name;
			//frag_file_name << elements << "_fragments.cha";
			//frag_file.open( frag_file_name.str().c_str(), std::ios::out );

			//std::fstream data_file;
			//std::ostringstream data_file_name;
			//data_file_name << elements << ".data";
			//data_file.open( data_file_name.str().c_str(), std::ios::out );

			//std::fstream def_chain_file;
			//std::ostringstream def_file_name;
			//def_file_name << "defg0.data.cha.dbf";
			//def_chain_file.open( def_file_name.str().c_str(), std::ios::out );

			//std::fstream semi_chain_file;
			//std::ostringstream semi_file_name;
			//semi_file_name << "semig0.data.cha.dbf";
			//semi_chain_file.open( semi_file_name.str().c_str(), std::ios::out );

			//std::fstream compl_chain_file;
			//std::ostringstream compl_file_name;
			//compl_file_name << elements << "_compl.data.cha.dbf";
			//compl_chain_file.open( compl_file_name.str().c_str(), std::ios::out );


			String< Fragment<int > > data;
			reserve( data, elements );
			
			_generateRandomFrags( data, elements, 1, 3 * elements, 10, 20, dim );
		
			//writeFragments( frag_file, data, FragFile() );
			//writeFragments( data_file, data, FragGnuplot() );

			String< Fragment<int> > chain1;
			reserve( chain1, elements );
			//benchChainerG0( data, chain1, Deferred(), c1 );

			//writeFragments( def_chain_file, chain1, ChainStatistics() );

			String< Fragment<int> > chain2;
			reserve( chain2, elements );
			benchChainerG0( data, chain2, SemiDeferred(), c2 );

			//writeFragments( semi_chain_file, chain2, ChainStatistics() );

			String< Fragment<int> > chain3;
			reserve( chain3, elements );
			benchChainerG0( data, chain3, Complete(), c3 );

			//writeFragments( compl_chain_file, chain3, ChainStatistics() );

			//data_file.flush();
			//data_file.close();
			//frag_file.flush();
			//frag_file.close();
			//def_chain_file.flush();
			//def_chain_file.close();
			//semi_chain_file.flush();
			//semi_chain_file.close();
			//compl_chain_file.flush();
			//compl_chain_file.close();
		}
		
		s << c1.Usr()/iterations_bm << "\t";

		def_file << s.str();
		def_file.flush();
		s.str("");
		s.clear();

		s << c2.Usr()/iterations_bm << "\t";
		
		semi_file << s.str();
		semi_file.flush();
		s.str("");
		s.clear();

		s << c3.Usr()/iterations_bm << "\t";

		compl_file << s.str();
		compl_file.flush();
		s.str("");
		s.clear();
	}


	void chain_benchmarkG1( int elements, 
							int dim, 
							int iterations_bm, 
							std::fstream & def_file,
							std::fstream & semi_file,
							std::fstream & compl_file )
//-------------------------------------- benchmarking, comparison to stl map ----------------------------------
	{
		
		
		std::ostringstream s;
		Clock c1;
		Clock c2;
		Clock c3;
		
		for( int times = 0; times < iterations_bm; ++ times )
		{
			//struct tm * act_time_format;
			//time_t act_time_simple;
			//time( &act_time_simple );
			//act_time_format = localtime( &act_time_simple );

			//char * time_string = new char[10];
			//sprintf( time_string, "%2i%2i_%2i%2i%2i", act_time_format->tm_mday,
			//										act_time_format->tm_mon,
			//										act_time_format->tm_hour,
			//										act_time_format->tm_min,
			//										act_time_format->tm_sec );
			//std::cout << time_string << std::endl;
			//std::ostringstream file_name;
			//file_name << "G1" << time_string;

			//std::fstream frag_file;
			//std::ostringstream frag_file_name;
			//frag_file_name << "fragments.cha";
			//frag_file.open( frag_file_name.str().c_str(), std::ios::out );

			//std::fstream data_file;
			//std::ostringstream data_file_name;
			//data_file_name << "g0.data";
			//data_file.open( data_file_name.str().c_str(), std::ios::out );

			//std::fstream def_chain_file;
			//std::ostringstream def_file_name;
			//def_file_name << "defg0.data.cha.dbf";
			//def_chain_file.open( def_file_name.str().c_str(), std::ios::out );

			//std::fstream semi_chain_file;
			//std::ostringstream semi_file_name;
			//semi_file_name << "semig0.data.cha.dbf";
			//semi_chain_file.open( semi_file_name.str().c_str(), std::ios::out );

			//std::fstream compl_chain_file;
			//std::ostringstream compl_file_name;
			//compl_file_name << "complg0.data.cha.dbf";
			//compl_chain_file.open( compl_file_name.str().c_str(), std::ios::out );


			String< Fragment<int > > data;
			reserve( data, elements );
			
			_generateRandomFrags( data, elements, 1, 3 * elements, 10, 20, dim );
		
			//writeFragments( frag_file, data, FragFile() );
			//writeFragments( data_file, data, FragGnuplot() );

				
			String< Fragment<int> > chain1;
			reserve( chain1, elements );
			benchChainerG1( data, chain1, Deferred(), c1 );

			//writeFragments( def_chain_file, chain1, ChainStatistics() );

			String< Fragment<int> > chain2;
			reserve( chain2, elements );
			benchChainerG1( data, chain2, SemiDeferred(), c2 );

			//writeFragments( semi_chain_file, chain2, ChainStatistics() );

			String< Fragment<int> > chain3;
			reserve( chain3, elements );
			benchChainerG1( data, chain3, Complete(), c3 );

			//writeFragments( compl_chain_file, chain3, ChainStatistics() );*/

			//data_file.flush();
			//data_file.close();
			//frag_file.flush();
			//frag_file.close();
			//def_chain_file.flush();
			//def_chain_file.close();
			//semi_chain_file.flush();
			//semi_chain_file.close();
			//compl_chain_file.flush();
			//compl_chain_file.close();
		}
		
		s << c1.Usr()/iterations_bm << "\t";

		def_file << s.str();
		def_file.flush();
		s.str("");
		s.clear();

		s << c2.Usr()/iterations_bm << "\t";
		
		semi_file << s.str();
		semi_file.flush();
		s.str("");
		s.clear();

		s << c3.Usr()/iterations_bm << "\t";

		compl_file << s.str();
		compl_file.flush();
		s.str("");
		s.clear();
		
	}


	void chain_benchmarkGSoP( int elements, 
							int dim, 
							int iterations_bm, 
							std::fstream & def_file,
							std::fstream & semi_file,
							std::fstream & compl_file )
//-------------------------------------- benchmarking, comparison to stl map ----------------------------------
	{
		
		
		std::ostringstream s;
		Clock c1;
		Clock c2;
		Clock c3;
		
		for( int times = 0; times < iterations_bm; ++ times )
		{
			//struct tm * act_time_format;
			//time_t act_time_simple;
			//time( &act_time_simple );
			//act_time_format = localtime( &act_time_simple );

			//char * time_string = new char[10];
			//sprintf( time_string, "%2i%2i_%2i%2i%2i", act_time_format->tm_mday,
			//										act_time_format->tm_mon,
			//										act_time_format->tm_hour,
			//										act_time_format->tm_min,
			//										act_time_format->tm_sec );
			//std::cout << time_string << std::endl;
			//std::ostringstream file_name;
			//file_name <<  "GSoP" << time_string;

			//std::fstream frag_file;
			//std::ostringstream frag_file_name;
			//frag_file_name << elements << "_fragments.cha";
			//frag_file.open( frag_file_name.str().c_str(), std::ios::out );

			//std::fstream data_file;
			//std::ostringstream data_file_name;
			//data_file_name << elements << ".data";
			//data_file.open( data_file_name.str().c_str(), std::ios::out );

			//std::fstream def_chain_file;
			//std::ostringstream def_file_name;
			//def_file_name << "defg0.data.cha.dbf";
			//def_chain_file.open( def_file_name.str().c_str(), std::ios::out );

			//std::fstream semi_chain_file;
			//std::ostringstream semi_file_name;
			//semi_file_name << "semig0.data.cha.dbf";
			//semi_chain_file.open( semi_file_name.str().c_str(), std::ios::out );

			//std::fstream compl_chain_file;
			//std::ostringstream compl_file_name;
			//compl_file_name << elements << "_compl.data.cha.dbf";
			//compl_chain_file.open( compl_file_name.str().c_str(), std::ios::out );


			String< Fragment<int > > data;
			reserve( data, elements );
			
			_generateRandomFrags( data, elements, 1, 3 * elements, 10, 20, dim );
		
			//writeFragments( frag_file, data, FragFile() );
			//writeFragments( data_file, data, FragGnuplot() );

			String< Fragment<int> > chain1;
			reserve( chain1, elements );
			benchChainerGSoP( data, chain1, Deferred(), c1 );

			//writeFragments( def_chain_file, chain1, ChainStatistics() );

			String< Fragment<int> > chain2;
			reserve( chain2, elements );
			benchChainerGSoP( data, chain2, SemiDeferred(), c2 );

			//writeFragments( semi_chain_file, chain2, ChainStatistics() );

			String< Fragment<int> > chain3;
			reserve( chain3, elements );
			benchChainerGSoP( data, chain3, Complete(), c3 );

			//writeFragments( compl_chain_file, chain3, ChainStatistics() );

			//data_file.flush();
			//data_file.close();
			//frag_file.flush();
			//frag_file.close();
			//def_chain_file.flush();
			//def_chain_file.close();
			//semi_chain_file.flush();
			//semi_chain_file.close();
			//compl_chain_file.flush();
			//compl_chain_file.close();
		}
		
		s << c1.Usr()/iterations_bm << "\t";

		def_file << s.str();
		def_file.flush();
		s.str("");
		s.clear();

		s << c2.Usr()/iterations_bm << "\t";
		
		semi_file << s.str();
		semi_file.flush();
		s.str("");
		s.clear();

		s << c3.Usr()/iterations_bm << "\t";

		compl_file << s.str();
		compl_file.flush();
		s.str("");
		s.clear();
	}


	

	void timetest_chaining()
	{
		srand(static_cast<unsigned>(time(0)));
		typedef int TKey;
						
		TKey iterations_bm = 5;

		struct tm * act_time_format;
		time_t act_time_simple;
		time( &act_time_simple );
		act_time_format = localtime( &act_time_simple );

		char * time_string = new char[10];
		sprintf( time_string, "%2i%2i_%2i%2i%2i", act_time_format->tm_mday,
													act_time_format->tm_mon,
													act_time_format->tm_hour,
													act_time_format->tm_min,
													act_time_format->tm_sec );
		std::cout << time_string << std::endl;
		std::ostringstream file_name;
		file_name << time_string;

		std::fstream def_file;
		std::ostringstream def_file_name;
		def_file_name << "def_searches_h" << _getMaximalSLTowerHeight( 1 ) << "_t" << _rt_thresh << "_d" << time_string << ".txt";
		def_file.open( def_file_name.str().c_str(), std::ios::out );

		std::fstream semi_file;
		std::ostringstream semi_file_name;
		semi_file_name << "semi_searches_h" << _getMaximalSLTowerHeight( 1 ) << "_t" << _rt_thresh << "_d" << time_string << ".txt";
		semi_file.open( semi_file_name.str().c_str(), std::ios::out );

		std::fstream compl_file;
		std::ostringstream compl_file_name;
		compl_file_name << "compl_searches_h" << _getMaximalSLTowerHeight( 1 ) << "_t" << _rt_thresh << "_d" << time_string << ".txt";
		compl_file.open( compl_file_name.str().c_str(), std::ios::out );

		//TKey keyval_pairs_bm[] = { 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000};
		TKey keyval_pairs_bm[] = { 10000, 50000, 100000, 200000, 500000 };
		for( int i = 0; i < 5; ++i ){
			def_file << keyval_pairs_bm[ i ] << "\t";
			semi_file << keyval_pairs_bm[ i ] << "\t";
			compl_file << keyval_pairs_bm[ i ] << "\t";
		}
		def_file << std::endl;
		semi_file << std::endl;
		compl_file << std::endl;

		def_file.flush();
		semi_file.flush();
		compl_file.flush();

		for( int i = 0; i < 4; ++i )
		{
			def_file << dimensions[ i ] << "\t";
			semi_file << dimensions[ i ] << "\t";
			compl_file << dimensions[ i ] << "\t";
		
			for( int j = 0; j < 5; ++j ){
			#ifdef CHAINGG0
				chain_benchmarkG0( keyval_pairs_bm[ j ], dimensions[ i ], iterations_bm, def_file, semi_file, compl_file );
			#endif
			#ifdef CHAING1
				chain_benchmarkG1( keyval_pairs_bm[ j ], dimensions[ i ], iterations_bm, def_file, semi_file, compl_file );
			#endif
			#ifdef CHAINSOP
				chain_benchmarkGSoP( keyval_pairs_bm[ j ], dimensions[ i ], iterations_bm, def_file, semi_file, compl_file );
			#endif	
			}	

			def_file << std::endl;
			semi_file << std::endl;
			compl_file << std::endl;
				
		}
			
		def_file.close();
		semi_file.close();
		compl_file.close();
	}




}

#endif

