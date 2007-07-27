#ifndef SEQAN_HEADER_TEST_BLAST_PARSING_H
#define SEQAN_HEADER_TEST_BLAST_PARSING_H

using namespace std;
//using namespace seqan;

namespace SEQAN_NAMESPACE_MAIN
{

template<typename T>
void Test_BlastHsp() {
	// Windows
#ifdef PLATFORM_WINDOWS
	//String<char> blast_path("C:\\emde\\blast\\bin\\");
	//String<char> db_path("C:\\emde\\blast\\bin\\data\\");
	//String<char> query_path("C:\\emde\\blast\\bin\\data\\");
	//String<char> out_path("C:\\emde\\blast\\bin\\data\\");
	String<char> blast_path("D:\\emde\\blast\\bin\\");
	String<char> db_path("D:\\emde\\blast\\bin\\data\\");
	String<char> query_path("D:\\emde\\blast\\bin\\data\\");
	String<char> out_path("D:\\emde\\blast\\bin\\data\\");
#else
	// Linux
	String<char> blast_path("/home/bude2/emde/blast/blast-2.2.16/");
	String<char> db_path("/export/local-1/public/emde/");
	String<char> query_path("/export/local-1/public/emde/");
	String<char> out_path("/export/local-1/public/emde/");
#endif


	std::stringstream s;
	s << out_path << "drosophp.out"; 


	typedef BlastHsp<BlastP, FullInfo> TBlastHsp;
	typedef ::std::fstream TFile;
	typedef BlastReport<TBlastHsp, StreamReport<TFile> > TBlastReport;
	typedef typename Hit<TBlastReport>::Type TBlastHit;

	typedef typename Iterator<TBlastReport,HitIterator>::Type THitIterator;
	typedef typename Iterator<TBlastHit,HspIterator>::Type THspIterator;

	typedef StringSet<String<AminoAcid> > TStringSet;
 	typedef Graph<Alignment<TStringSet> > TAliGraph;
	
	std::fstream strm;
	strm.open((s.str()).c_str(),ios_base::in | ios_base::binary);

	int ali_count = 0;
	TBlastReport blast;
	while(!atEnd(strm,blast)) 
	{
		int hitcount = 0;
		int hspcount = 0;
		read(strm,blast,Blast());
		THitIterator hit_it(blast); 
		for(; !atEnd(strm,hit_it); goNext(strm,hit_it)) 
		{
			++hitcount;
			TBlastHit hit = getValue(strm,hit_it);
			THspIterator hsp_it(hit);
			for(; !atEnd(strm,hsp_it); goNext(strm,hsp_it)) 
			{
				++hspcount;
				++ali_count;	
				TBlastHsp hsp = getValue(strm,hsp_it);
				TAliGraph ali_g;
				getAlignment(hsp,ali_g);
				//cout << "QueryBegin: " << queryBegin(hsp)<< "=" << getQueryBegin(hsp)<<"\n";
				//cout << "DataBBegin: " << databaseBegin(hsp)<< "=" <<getDatabaseBegin(hsp)<<"\n";
				//cout << "QueryEnd  : " << queryEnd(hsp)<<"="<<getQueryEnd(hsp)<<"\n";
				//cout << "DataBEnd  : " << databaseEnd(hsp)<< "=" <<getDatabaseEnd(hsp)<<"\n";
				//queryAlignmentString(hsp);
				//getQueryAlignmentString(hsp);
				//databaseAlignmentString(hsp);
				//databaseQueryAlignmentString(hsp);
				//cout << "EValue    : "<<eValue(hsp)<<"="<<getEValue(hsp)<<"\n";
				//cout << "Score     : "<<score(hsp)<<"="<<getScore(hsp)<<"\n";
				//cout << "Bits      : "<<bitScore(hsp)<<"="<<getBitScore(hsp)<<"\n";
				//cout << "Identity  : "<<percentIdentity(hsp)<<"="<<getPercentIdentity(hsp)<<"\n";
				//cout << "Gaps      : "<<percentGaps(hsp)<<"="<<getPercentGaps(hsp)<<"\n";
				//cout << "QueryStran: "<<queryOrientationPlus(hsp)<<"\n";
				//cout << "DBStrand  : "<<databaseOrientationPlus(hsp)<<"\n";
//				cout << "QueryFrame: "<<queryFrame(hsp)<<"="<<getQueryFrame(hsp)<<"\n";
//				cout << "DBFrame   : "<<databaseFrame(hsp)<<"="<<getDatabaseFrame(hsp)<<"\n";
//				cout << "Positives : "<<percentPositives(hsp)<<"="<<getPercentPositives(hsp)<<"\n\n";
			}
	//		cout <<"\n";
		}
		cout << "\nhits: "<<hitcount <<"\n"<<std::flush;
		cout << "\nhsps: "<<hspcount <<"\n"<<std::flush;
	}


	//cout << "NumHsps: "<< ali_count<<"\n";

//37 42
//33 39

}



//////////////////////////////////////////////////////////////////////////////
template<typename T>
void Test_BlastStoreReport() {

	// Windows
#ifdef PLATFORM_WINDOWS
	//String<char> blast_path("C:\\emde\\blast\\bin\\");
	//String<char> db_path("C:\\emde\\blast\\bin\\data\\");
	//String<char> query_path("C:\\emde\\blast\\bin\\data\\");
	//String<char> out_path("C:\\emde\\blast\\bin\\data\\");
	String<char> blast_path("D:\\emde\\blast\\bin\\");
	String<char> db_path("D:\\emde\\blast\\bin\\data\\");
	String<char> query_path("D:\\emde\\blast\\bin\\data\\");
	String<char> out_path("D:\\emde\\blast\\bin\\data\\");
#else
	// Linux
	String<char> blast_path("/home/bude2/emde/blast/blast-2.2.16/");
	String<char> db_path("/export/local-1/public/emde/");
	String<char> query_path("/export/local-1/public/emde/");
	String<char> out_path("/export/local-1/public/emde/");
#endif

	typedef BlastHsp<BlastN,BasicInfo> TBlastHsp;
	typedef BlastReport<TBlastHsp,StoreReport<FullInfo> > TBlastReport;
	typedef typename Hit<TBlastReport>::Type TBlastHit;
	typedef typename Hsp<TBlastReport>::Type TBlastHsp;
	typedef typename Iterator<TBlastReport,HitIterator>::Type THitIterator;
	typedef typename Iterator<TBlastHit,HspIterator>::Type THspIterator;

	// the blast report that is supposed to be parsed
	//std::stringstream s;
	//s << out_path << "test2mac.out"; 

	// default: reads the whole report and stores only the basic information (default assumes blastn report)
	std::fstream strm;
	//strm.open((s.str()).c_str(),ios_base::in | ios_base::binary);
	strm.open(TEST_PATH "test2_crlf.out",ios_base::in | ios_base::binary);
	TBlastReport blast;
	
	int numhsps = 0;
	while(!atEnd(strm,blast)) 
	{
		read(strm,blast,Blast());   // complete report is now parsed and all the hits (and hsps) are stored
		cout << "Query  : " << queryName(blast)<<" = "<<	getQueryName(blast) <<"\n";
		cout << "DB     : " << databaseName(blast)<<" = "<< getDatabaseName(blast)<<"\n";
		cout << "NumHits: "<< numHits(blast)<<"  \nNumHsps: "<< numHsps(blast) <<"\n";
		numhsps+=numHsps(blast);

		THitIterator hit_it(blast); 
		for(; !atEnd(hit_it); goNext(hit_it)) 
		{
			TBlastHit hit = getValue(hit_it);
			THspIterator hsp_it(value(hit_it));
			for(; !atEnd(hsp_it); goNext(hsp_it)) 
			{
				TBlastHsp hsp = getValue(hsp_it);
				//
			}
		}
		if(atEnd(strm,blast))
		{
			cout << "NumHsps insgesamt: "<<numhsps<<"\n";
			cout << "ECutoff: "<< eValueCutoff(blast)<< " = "<< getEValueCutoff(blast)<< "\n";
			cout << "Matrix : "<< matrixName(blast)<<" = "<< getMatrixName(blast)<< "\n";
			if(gapsAllowed(blast))
			{
				cout << "GapOpen: "<< gapOpen(blast) << " = " << getGapOpen(blast) <<"\n";
				cout << "GapExt : "<< gapExtension(blast)<<" = "<<getGapExtension(blast)<<"\n";
			}
			//cout << "Lambda : " <<lambda(blast)<<" = " <<getLambda(blast) << "\n";
			//cout << "K : "<< kappa(blast)<< " = "<<getKappa(blast)<<"  H : " << entropy(blast)<< " = "<<getEntropy(blast)<<"\n";
			//cout << "GLambda: " <<gappedLambda(blast)<<" = " <<getGappedLambda(blast) << "\n";
			//cout << "GK: "<< gappedKappa(blast)<< " = "<<getGappedKappa(blast)<<"  GH: " << gappedEntropy(blast)<< " = "<<getGappedEntropy(blast)<<"\n";
		}
	}


	

}

//////////////////////////////////////////////////////////////////////////////
template<typename T>
void Test_BlastParsing() {

	// Windows
#ifdef PLATFORM_WINDOWS
	//String<char> blast_path("C:\\emde\\blast\\bin\\");
	//String<char> db_path("C:\\emde\\blast\\bin\\data\\");
	//String<char> query_path("C:\\emde\\blast\\bin\\data\\");
	//String<char> out_path("C:\\emde\\blast\\bin\\data\\");
	String<char> blast_path("D:\\emde\\blast\\bin\\");
	String<char> db_path("D:\\emde\\blast\\bin\\data\\");
	String<char> query_path("D:\\emde\\blast\\bin\\data\\");
	String<char> out_path("D:\\emde\\blast\\bin\\data\\");
#else
	// Linux
	String<char> blast_path("/home/bude2/emde/blast/blast-2.2.16/");
	String<char> db_path("/export/local-1/public/emde/");
	String<char> query_path("/export/local-1/public/emde/");
	String<char> out_path("/export/local-1/public/emde/");
#endif


	typedef BlastHsp<BlastN, FullInfo> TBlastHsp;
	typedef ::std::fstream TFile;
	typedef BlastReport<TBlastHsp, StreamReport<TFile> > TBlastReport;
	typedef typename Hit<TBlastReport>::Type TBlastHit;
	typedef typename Iterator<TBlastReport,HitIterator>::Type THitIterator;
	typedef typename Iterator<TBlastHit,HspIterator>::Type THspIterator;
	typedef StringSet<String<Dna> > TStringSet;
 	typedef Graph<Alignment<TStringSet> > TAliGraph;



	std::stringstream s;
	s << out_path << "testecoli.out"; 

	std::fstream strm2;
	strm2.open((s.str()).c_str(),ios_base::in | ios_base::binary);

	int ali_count = 0;
	TBlastReport blast2;
	while(!atEnd(strm2,blast2)) 
	{
		int hitcount = 0;
		int hspcount = 0;
		read(strm2,blast2,Blast());
		THitIterator hit_it(blast2); 
		for(; !atEnd(strm2,hit_it); goNext(strm2,hit_it)) 
		{
			++hitcount;
			TBlastHit hit = getValue(strm2,hit_it);
			THspIterator hsp_it(hit);
			for(; !atEnd(strm2,hsp_it); goNext(strm2,hsp_it)) 
			{
				++hspcount;
				TBlastHsp hsp = getValue(strm2,hsp_it);
	//			if(eValue(hsp) < 0.02)
	//			{
					Align< String<Dna>, ArrayGaps> ali;
					getAlignment(hsp,ali,UnknownSource());
			//		cout << ali <<"\n";
					//TAliGraph ali_g(str);
					//getAlignment(hsp,ali_g,0,1); //hit ID
					TAliGraph ali_g;
					getAlignment(hsp,ali_g); //hit ID
					++ali_count;
					//cout << percentGaps(hsp) << "\n";
			//		if(percentGaps(hsp) > 0)
			//		{	
						//cout << queryAlignmentString(hsp)<< "\n";
						//cout << databaseAlignmentString(hsp) << "\n\n";
						//cout << ali_g <<"\n";
						//cout << ali <<"\n";
			//		}
	//			}
			}
		}
//		cout << "\nhits: "<<hitcount <<"\n"<<std::flush;
//		cout << "\nhsps: "<<hspcount <<"\n"<<std::flush;
	}

	cout << ali_count<< " alignments parsed\n";
}


}

#endif

