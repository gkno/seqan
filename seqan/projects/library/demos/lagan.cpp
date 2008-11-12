#include <seqan/seeds.h>
#include <seqan/file.h>
#include <fstream>

using namespace seqan;
using namespace std;

template<typename TText, typename TScore, typename TValue, typename TChain>
void
chaos(int ktup,
	  String<TText> &query_ori, 
	  String<TText> &database_ori, 
	  TValue s1,
	  TValue s2,
	  TValue e1,
	  TValue e2,
	  Score<TScore, Simple> &scoreMatrix, 
	  int gap_threshold,
	  TChain &chain)
{

	TScore score_threshold=4*ktup;
	if ((e2-s2 > 500 ) && (e1-s1 > 500)){ 
		typedef Index< String<TText>, Index_QGram<SimpleShape > > TQGramIndex;
		typedef SeedSet<int, SimpleSeed, DefaultScore> TSeedSet;
		Segment<String<TText>, InfixSegment> query(query_ori, s1, e1+1);
		Segment<String<TText>, InfixSegment> database(database_ori, s2, e2+1);

		TQGramIndex index_qgram(database);
		TSeedSet seedCollector(15, score_threshold, scoreMatrix);

		while ((length(seedCollector) == 0) && (ktup > 5)){
			resize(indexShape(index_qgram), ktup);
			Finder<TQGramIndex> finder(index_qgram);
			
			String<TText> qgram;
			int qLength =  length(query)-ktup+1;
			for (int i = 0; i < qLength; ++i){

				qgram= infix(query, i, i+ktup);
				while (find(finder, qgram)){
					int pos = position(finder);
					if (!addSeed(seedCollector, s1+i, s2+pos, ktup, 0, Merge()))
						if (!addSeed(seedCollector,s1+i, s2+pos, ktup, query_ori, database_ori, 5, Chaos()))
							addSeed(seedCollector, s1+i, s2+pos, ktup, Single());
				}
				clear(finder);
			}
			ktup-=1;
		}

		if (length(seedCollector) != 0)
		{
			globalChaining(seedCollector, chain);
			clear(seedCollector);

			if (ktup >=3){
				TChain chain2;
				typename TChain::iterator it = chain.begin();
				typename TChain::iterator it2 = it;
				++it2;

			
				typename TChain::iterator it_end = chain.end();

				if (((e1 - rightDim0(*it)) >= (int)gap_threshold) && ((e2 - rightDim1(*it)) >= (int)gap_threshold))
				{	
					chaos(ktup, query_ori, database_ori, rightDim0(*it)+1, rightDim1(*it)+1, e1, e2, scoreMatrix, gap_threshold, chain2);
					chain.splice(it, chain2);
				}

				for (;it2 != it_end; ++it2)
				{
					it = it2;
					--it;
					if (((leftDim0(*it) - rightDim0(*it2)) >= gap_threshold) && ((leftDim1(*it) - rightDim1(*it2)) >= gap_threshold)){
						list<Seed<int, SimpleSeed> > chain2;
						chaos(ktup, query_ori, database_ori, rightDim0(*it2)+1, rightDim1(*it2)+1, leftDim0(*it)-1, leftDim1(*it)-1, scoreMatrix, gap_threshold, chain2);
						chain.splice(it2, chain2);
					}
				}
					it = it2;
					--it;
				if (((leftDim0(*it) - s1) >= gap_threshold) && ((leftDim1(*it) - s2) >= gap_threshold)){
					chaos(ktup, query_ori, database_ori, s1, s2, leftDim0(*it)-1, leftDim1(*it)-1, scoreMatrix, gap_threshold, chain2);
					chain.splice(it2, chain2);
				}

			}
		}
	}
}




template<typename TValue>
void
appendValue(list<TValue> vec, TValue value)
{
	vec.push_back(value);
}

template<typename TText, typename TScore>
void
lagan(int ktup, 
	  String<TText> query,
	  String<TText> database,
	  Score<TScore, Simple> scoreMatrix, 
	  int gap_threshold,
	  int score_threshold)
{

	typedef Index< String<TText>, Index_QGram<SimpleShape > > TQGramIndex;
	typedef SeedSet<int, SimpleSeed, DefaultScore> TSeedSet;

	TQGramIndex index_qgram(database);
	Finder<TQGramIndex> finder(index_qgram);
	TSeedSet seedCollector(6, score_threshold, scoreMatrix);
	String<TText> qgram;

	while ((length(seedCollector) == 0)&&(ktup > 5)){
		resize(indexShape(index_qgram), ktup);	
		for (int i = 0; i < (int)length(query)-ktup+1; ++i)
		{
			qgram= infix(query, i, i+ktup);
			
			while (find(finder, qgram)){
				
				int pos = position(finder);
				if (!addSeed(seedCollector, i, pos, ktup, 0, Merge()))
					if (!addSeed(seedCollector,i, pos, ktup, query, database, 5, Chaos()))
						addSeed(seedCollector, i, pos, ktup, Single());
			}
			clear(finder);
		}
		ktup -= 3;
	}

	if (length(seedCollector) != 0)
	{
		typedef list<Seed<int, SimpleSeed> > TChain;
		TChain chain;
		TChain chain2;
		globalChaining(seedCollector, chain);
	
		clear(seedCollector);
		if (ktup >=3){
			TChain::iterator it2 = chain.begin();
			TChain::iterator it = it2;
			++it2;

			if (((static_cast<int>(length(query))-1 - rightDim0(*it)) >= gap_threshold) && ((static_cast<int>(length(database))-1 - rightDim1(*it)) >= gap_threshold)){	
				chaos(ktup, query, database, rightDim0(*it)+1, rightDim1(*it)+1, static_cast<int>(length(query)), static_cast<int>(length(database)), scoreMatrix, gap_threshold, chain2);
				chain.splice(it, chain2);
			}

			TChain::iterator it_end = chain.end();

			for (; it2 != it_end; ++it2)
			{
				it = it2;
				--it;
				if (((leftDim0(*it) - rightDim0(*it2)) >= gap_threshold) && ((leftDim1(*it) - rightDim1(*it2)) >= gap_threshold)){
					chaos(ktup, query, database, rightDim0(*it2)+1, rightDim1(*it2)+1, leftDim0(*it)-1, leftDim1(*it)-1, scoreMatrix, gap_threshold, chain2);
					chain.splice(it2, chain2);
				}
				
			}
			it = it2;
			--it;
			if ((leftDim0(*it) >= gap_threshold) && (leftDim1(*it) >= gap_threshold)){
				chaos(ktup, query, database, 0, 0, leftDim0(*it)-1, leftDim1(*it)-1, scoreMatrix, gap_threshold, chain2);
				chain.splice(it2, chain2);
			 }
		}
		
		TChain::iterator chainend = chain.end();

		Align<String<char>,ArrayGaps> alignment;
		resize(rows(alignment), 2);
		assignSource(row(alignment, 0), query);
		assignSource(row(alignment, 1), database);
		cout << "Score: " << bandedChainAlignment(chain, 7, alignment, scoreMatrix) << endl;
		cout << alignment << endl;
	}
}



int
main()
{
	fstream fstrm;
	fstrm.open("./lagan1.fasta", ios_base::in | ios_base::binary);
    String<char> fasta_tag;
    String<Dna> fasta_seq_query;
    readMeta(fstrm, fasta_tag, Fasta());

    read(fstrm, fasta_seq_query, Fasta());
    fstrm.close(); 
	fstream fstrm2;
	fstrm2.open("./lagan2.fasta", ios_base::in | ios_base::binary);
    String<Dna> fasta_seq_database;
    readMeta(fstrm2, fasta_tag, Fasta());
    read(fstrm2, fasta_seq_database, Fasta());
    fstrm2.close(); 

	//Fasta file
	Score<int, Simple> scoreMatrix(3,-2,-1,-3);
	if ((length(fasta_seq_query) > 0) && (length(fasta_seq_database) > 0))
		lagan(13, fasta_seq_query, fasta_seq_database, scoreMatrix, 500, 30);
	else
		cout << "Error - File problem" << endl;
}
