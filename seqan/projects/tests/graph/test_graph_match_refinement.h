#ifndef SEQAN_HEADER_TEST_GRAPH_MATCH_REFINEMENT_H
#define SEQAN_HEADER_TEST_GRAPH_MATCH_REFINEMENT_H

using namespace std;
using namespace seqan;

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

int Test_ConvertSequences(String<char> const in_path, String<char> const in_file, String<char> const path, String<char> const file_prefix) {
	typedef String<Dna5, External<ExternalConfig<File<>, 64*1024> > > TString;
	
	// count sequences
	unsigned seqCount = 0;

	ifstream file;
	std::stringstream input;
	input << in_path << in_file;
	file.open(input.str().c_str(), ios_base::in | ios_base::binary);
	if (!file.is_open()) return 0;
	while (!_streamEOF(file)) {
		String<char> id;
		readID(file, id, Fasta());
		std::cout << id << std::endl;
		goNext(file, Fasta());
		++seqCount;
	}
	std::cout << "Number of sequences: " << seqCount << std::endl;

	// import sequences
	file.clear();
	file.seekg(0, ios_base::beg);
	for(unsigned i = 0; (i < seqCount) && !_streamEOF(file); ++i) 
	{
		TString str;
		//String<TraceBack, External<> > trace;
		//open(trace, "D:\\seqan.dat");
		std::stringstream s;
		s << path << file_prefix << i;
		open(str, s.str().c_str());
		read(file, str, Fasta());
	}
    file.close();

	return seqCount;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TVal1, typename TVal2>
inline bool 
Test_ReadSequences(String<char> const path, String<char> const file_prefix, TStringSet& str, TVal1 const start, TVal2 const nseq) {
	for(unsigned i = start; i < start + nseq; ++i) {
		std::stringstream s;
		s << path << file_prefix << i - start;
		bool f = open(str[i], s.str().c_str());
		if (!f) return false;
	}
	return true;
}

//////////////////////////////////////////////////////////////////////////////

void Test_MatchRefinement() {
	// Sequences
	typedef String<Dna5, External<ExternalConfig<File<>, 64*1024> > > TString;
	typedef StringSet<TString, Owner<> > TStringSet;
	typedef Id<TStringSet>::Type TId;
	typedef Size<TStringSet>::Type TSize;

	// Matches
	typedef String<Fragment<>, External<> > TFragmentString;
	
	// Windows
#ifdef PLATFORM_WINDOWS
	String<char> in_path("Z:\\matches\\");
	String<char> out_path("Z:\\matches\\out\\");
#else
	// Linux
	String<char> in_path("/home/takifugu/rausch/matches/");
	String<char> out_path("/home/takifugu/rausch/matches/out/");
#endif
	

	TSize hSeq = 24;
	TSize wSeq = 24;
	TSize bSeq = 24;
	

	// Convert all sequences only once
	//TSize tmp = 0;
	//tmp = Test_ConvertSequences(in_path, "HUREF6CHROM.fasta",out_path,"H.chr.");
	//if (tmp != hSeq) {
	//	std::cerr << "Did not read all HUREF sequences." << std::endl;
	//	exit(1);
	//}
	//tmp = Test_ConvertSequences(in_path, "WGSACHROM.fasta",out_path,"W.chr.");
	//if (tmp != wSeq) {
	//	std::cerr << "Did not read all WGSA sequences." << std::endl;
	//	exit(1);
	//}
	//tmp = Test_ConvertSequences(in_path, "B36LCCHROM.fasta",out_path,"B.chr.");
	//if (tmp != bSeq) {
	//	std::cerr << "Did not read all B36LC sequences." << std::endl;
	//	exit(1);
	//}

	// Read all the sequences
	TStringSet str;
	resize(str, hSeq + wSeq + bSeq);
	bool f = Test_ReadSequences(out_path,"H.chr.", str, 0, hSeq);
	if (!f) {
		std::cerr << "Error importing HUREF sequences." << std::endl;
		exit(1);
	}
	f = Test_ReadSequences(out_path,"W.chr.", str, hSeq, wSeq);
	if (!f) {
		std::cerr << "Error importing WGSA sequences." << std::endl;
		exit(1);
	}
	f = Test_ReadSequences(out_path,"B.chr.", str, hSeq + wSeq, bSeq);
	if (!f) {
		std::cerr << "Error importing B36LC sequences." << std::endl;
		exit(1);
	}

	// Build a map:
	// SeqId -> Identifier
	typedef std::map<TId, String<char> > TIdToNameMap;
	TIdToNameMap idToName;
	for(TId i = 0; i < length(str); ++i) {
		TSize index = 0;
		std::stringstream s;
		if (i < 24) {
			s << "H";
			index = i;
		}
		else if (i < 48) {
			s << "W";
			index = i - hSeq;
		}
		else {
			s << "B";
			index = i - (hSeq + wSeq);
		}
		s << ":" << index;
		idToName.insert(std::make_pair(i, s.str().c_str()));
	}

	// Just to check that everything worked
	std::cout << "Number of sequences: " << length(str) << std::endl;
	for(TIdToNameMap::const_iterator pos =  idToName.begin(); pos != idToName.end(); ++pos) {
		std::cout << pos->second << ") ";
		for(TSize i=0; i<10;++i) {
			std::cout << str[pos->first][i];
		}
		std::cout << std::endl;
	}

	// Access the matches
	TFragmentString matches;
	std::stringstream strstream;
	strstream << out_path << "matchesTest.dat"; // 10 Matches
	//strstream << out_path << "matches10000.dat"; // 10 Matches
	open(matches, strstream.str().c_str());

	// Convert the matches to an external string
	//for(TSize i = 1; i<4; ++i) {
	//	fstream strm; 
	//	std::stringstream s;
	//	if (i==0 ) s << in_path << "TvsT.atac";
	//	else if (i==1 ) s << in_path << "BvsH.atac";
	//	else if (i==2 ) s << in_path << "BvsW.atac";
	//	else if (i==3 ) s << in_path << "WvsH.atac";
	//	strm.open(s.str().c_str(), ios_base::in);
	//	read(strm, matches, 10000, AtacMatches());
	//	strm.close();
	//}

	// Print all matches
	std::cout << "Number of matches: " << length(matches) << std::endl;
	typedef Infix<TString>::Type TInfix;
	for(TSize i = 0; i < length(matches); ++i) {
		TId seqId1 = sequenceId(matches[i],0);
		TId seqId2 = sequenceId(matches[i],1);
		TSize seqBegin1 = fragmentBegin(matches[i], seqId1);
		TSize seqBegin2 = fragmentBegin(matches[i], seqId2);
		TSize len = fragmentLength(matches[i], seqId1);
		TIdToNameMap::const_iterator pos1 =  idToName.find(seqId1);
		TIdToNameMap::const_iterator pos2 =  idToName.find(seqId2);
		std::cout << pos1->second << ") ";
		TInfix infix1 = infix(str[seqId1], seqBegin1, seqBegin1+len);
		std::cout << infix1 << std::endl;
		std::cout << pos2->second << ") ";
		TInfix infix2 = infix(str[seqId2], seqBegin2, seqBegin2+len);
		std::cout << infix2 << std::endl;
		std::cout << std::endl;
	}

	/*
	Score<int> score_type = Score<int>(1,-1,-2,0) ;
	typedef Graph<Alignment<TStringSet> > TAliGraph;
	
	TAliGraph ali_graph(str);

	matchRefinement(matches,str,score_type,ali_graph);//,StoreEdges());
	std::cout << "\nnumEdges: "<<numEdges(ali_graph)<<"\n";
	std::cout << "\nnumVertices: "<<numVertices(ali_graph)<<"\n";
	*/
}

}

#endif

