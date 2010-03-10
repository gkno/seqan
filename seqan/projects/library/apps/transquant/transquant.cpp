#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/store.h>
#include <seqan/misc/misc_cmdparser.h>

// TODO(holtgrew): This raises a warning with Boost 1.42. Deactivate warnings, activate again afterwards. The correct #pragma has to be used for each supported compiler.
#include <boost/math/distributions/normal.hpp>

using namespace seqan;
using namespace std;

///////////////////////////////////////////////////////////////////////////////

struct Count
{
	String<int>	ids;
	double		count;
};

struct Match {
	unsigned	contigId;
	int			posBegin;
	int			posEnd;
	bool		forward;	// true..forward strand, false..reverse complement strand
	std::string line;
};

struct Stats
{
	String<Pair<double, int> >	contigStats;			// MContig id -> (count,length)
	String<double>	histSingleContigs;	// histogram of the number of exons covered by per read
	String<double>	histMPContigs;		// histogram of the number of exons covered by per read
	String<double>	histSingleAmbig;	// histogram of the ambiguity per read
	String<double>	histMPAmbig;		// histogram of the ambiguity per read
	String<int>		histMatchErrors;	// histogram of the number of errors per match
	double			fragMean;			// empirical mean of fragment size
	double			fragStd;			// empirical std.dev of fragment size using *command line mean*
	unsigned long	fragments;			// to divide fragMean and fragStd by
	unsigned long	brokenMPs;
// stats directly taken from the result file
	unsigned long	matches;	
	unsigned long	reads;
	unsigned long	matepairs;
};

Stats stats;

// DContig - Dave's contig is a contig (e.g. chromosome) in the fragment store
// MContig - Marcel's contig is a node in the de Bruijn graph

String<String<int> >		transToDContig;	// (locus,transcript) -> DContig id
String<CharString>			dContigToTrans;	// DContig id -> original name, e.g. Locus_0_Transcript_0_Confidence_0.800

String<int>					dContigToAnno;	// DContig id -> its first node annotation
String<String<Count> >		counts;			// node[0] -> ((node[1..],count), (node[1..],count), ...)

String<int>					nodeToMContig;	// node number -> MContig id (compact renumeration)
String<int>					mContigToNode;	// MContig id -> node number (reverse map)
int							lastMContig = 0;
	
typedef std::map<string,int>	TTransName2Num;
TTransName2Num					transName2Num;

bool disablePairedEnds;
boost::math::normal fragmentDistr;
double width = 0.5;
unsigned histSize = 18;

void initStats()
{
	stats.fragMean = 0;
	stats.fragStd = 0;
	stats.fragments = 0;
	stats.brokenMPs = 0;

	stats.matches = 0;
	stats.reads = 0;
	stats.matepairs = 0;
}

template <typename TString>
void sumStats(TString &s)
{
	double sum = 0;
	double prod = 0;
	unsigned i;
	for (i = 0; i < length(s); ++i)
	{
		sum += s[i];
		prod += i*s[i];
	}
	if (i < histSize + 2)
		fill(s, histSize + 2, 0);
	s[histSize] = sum;
	s[histSize + 1] = prod;
}

void printStats()
{
	if (stats.fragments != 0)
	{
		stats.fragMean /= stats.fragments;
		stats.fragStd = sqrt(stats.fragStd / stats.fragments);
	}
	sumStats(stats.histSingleContigs);
	sumStats(stats.histMPContigs);
	sumStats(stats.histSingleAmbig);
	sumStats(stats.histMPAmbig);
	sumStats(stats.histMatchErrors);
	
	cout << "Matches:               \t" << stats.matches << endl;
	cout << "Reads:                 \t" << stats.reads << endl;
	cout << "Mate-Pairs:            \t" << stats.matepairs << endl;

	cout << "Valid Mate-Pairs:      \t" << stats.fragments << endl;
	cout << "Broken Mate-Pairs:     \t" << stats.brokenMPs << endl;
	cout << "Fragment size mean:    \t" << stats.fragMean << endl;
	cout << "Fragment size std:     \t" << stats.fragStd << endl;
	cout << "                       ";
	for (unsigned i = 0; i < histSize; ++i)   cout << '\t' << i;
	cout << "\tSum\tWghtd Sum";
	cout << "\nContigs per read:    ";
	for (unsigned i = 0; i < histSize+2; ++i) cout << '\t' << std::fixed << stats.histSingleContigs[i];
	cout << "\nContigs per MP:      ";
	for (unsigned i = 0; i < histSize+2; ++i) cout << '\t' << std::fixed << stats.histMPContigs[i];
	cout << "\nMatches per read:    ";
	for (unsigned i = 0; i < histSize+2; ++i) cout << '\t' << std::fixed << stats.histSingleAmbig[i];
	cout << "\nMP-Matches per MP:   ";
	for (unsigned i = 0; i < histSize+2; ++i) cout << '\t' << std::fixed << stats.histMPAmbig[i];
	cout << "\nErrors per match:   ";
	for (unsigned i = 0; i < histSize+2; ++i) cout << '\t' << stats.histMatchErrors[i];
	cout << '\n';
}

template <typename TCounts, typename TCount, typename TValue>
inline void
addCount(TCounts &cnts, TCount &cnt, TValue val)
{
	typedef typename Value<TCounts>::Type TCountLine;
	typename Suffix<TCount>::Type suf = suffix(cnt, 1);
	if (cnt[0] >= (int)length(cnts))
		resize(cnts, cnt[0] + 1, Generous());
	
	TCountLine &c = cnts[cnt[0]];
	unsigned i;
	for (i = 0; i < length(c); ++i)
	{
		if (c[i].ids > suf) break;
		if (c[i].ids == suf)
		{
			c[i].count += val;
			return;
		}
	}
	Count a;
	a.ids = suf;
	a.count = val;
	insertValue(c, i, a);
}

template <typename TSpec, typename TConfig>
inline void
loadTranscriptAnnotation(FragmentStore<TSpec, TConfig> & store, CharString const &fileName)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type TAnnotation;
	typedef std::fstream TFile;
	
	TFile file;
	file.open(toCString(fileName), ios_base::in | ios_base::binary);

	typename Value<TFile>::Type c = _streamGet(file);
	CharString header_;
	std::string header;
	std::string transName;

	unsigned contigId = 0;
	if (!empty(transToDContig))
		contigId = back(back(transToDContig)) + 1;
	if (empty(dContigToAnno))
		appendValue(dContigToAnno, 0);

	while (!_streamEOF(file))
	{
		while (c != '>' && !_streamEOF(file))
			_parse_skipLine(file, c);
		if (_streamEOF(file)) break;
		
		c = _streamGet(file);
		header_ = _parse_readIdentifier(file, c);
		assign(header, suffix(header_, 6));
		size_t confPos = header.find("_Confidence_");
		if (confPos != header.npos)
			header.erase(confPos);		

		std::istringstream iss(header);
		_parse_skipLine(file, c);

		char d;
		int locusNum = -1;
		int transNum = -1;
//		double confidence = 0;

		iss >> locusNum >> d >> d;
		if (d != 'T') continue;
		for (int i = 0; i < 10; ++i) iss >> d;
		iss >> transName;
//		for (int i = 0; i < 12; ++i) iss >> d;
//		iss >> confidence;

		if (locusNum >= (int)length(transToDContig))
			resize(transToDContig, locusNum + 1, Generous());

		// next 2 lines patched
		transNum = length(transToDContig[locusNum]);
		transName2Num[transName] = transNum;
		if (transNum >= (int)length(transToDContig[locusNum]))
			fill(transToDContig[locusNum], transNum + 1, 0, Generous());

		transToDContig[locusNum][transNum] = contigId;
		appendValue(dContigToTrans, header_, Generous());

		int pos = 0;
		while (!_streamEOF(file))
		{
			int nodeId = _parse_readNumber(file, c);			
			if (c != ':') 
				std::cerr<<"HUH1?"<<std::endl;

			c = _streamGet(file);	// :
			// quick fix for Hugues' bug
			//int endPos = _parse_readNumber(file, c);
			int endPos = pos + _parse_readNumber(file, c);
			if (nodeId < 0) nodeId = -nodeId;

			// rename nodeIds
			if (nodeId >= (int)length(nodeToMContig))
				fill(nodeToMContig, nodeId + 1, -1, Generous());
			if (nodeToMContig[nodeId] == -1)
			{
				appendValue(stats.contigStats, Pair<double, int>(0, endPos-pos));
				nodeToMContig[nodeId] = lastMContig++;
				appendValue(mContigToNode, nodeId, Generous());
			}
			nodeId = nodeToMContig[nodeId];

			TAnnotation a;
			a.beginPos = pos;
			a.endPos = endPos;
			a.parentId = 0;
			a.contigId = contigId;
			a.countId = nodeId;
			appendValue(store.annotationStore, a, Generous());

			if (_streamEOF(file) || c != '-')
				break;

			c = _streamGet(file);						
			if (c != '(') 
				std::cerr<<"HUH3?"<<std::endl;

			c = _streamGet(file);
			int gapLen = _parse_readNumber(file, c);	
			if (c != ')') 
				std::cerr<<"HUH4?"<<std::endl;

			c = _streamGet(file);	// )
			c = _streamGet(file);	// -
			c = _streamGet(file);	// >

			pos = endPos + gapLen;
		}
		++contigId;
		appendValue(dContigToAnno, length(store.annotationStore), Generous());
	}
	file.close();
}

template<typename TSpec, typename TConfig, typename TMatches>
inline void
countSingleMatches(FragmentStore<TSpec, TConfig> &store, TMatches const &matches)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type						TAnnotation;
	typedef typename Iterator<TMatches>::Type							TMIter;
	typedef typename Iterator<TAnnotationStore>::Type					TIter;
	
	String<int> ids;
	unsigned ambig = length(matches);
	TMIter m = begin(matches, Standard());
	TMIter mEnd = end(matches, Standard());
	for (; m != mEnd; ++m)
	{
		TIter it = begin(store.annotationStore, Standard()) + dContigToAnno[(*m).contigId];
		TIter itEnd = begin(store.annotationStore, Standard()) + dContigToAnno[(*m).contigId + 1];
		clear(ids);
		
		for (; it != itEnd; ++it)
		{
			if ((*m).posEnd <= (int)(*it).beginPos) break;
			if ((*m).posBegin >= (int)(*it).endPos) continue;
			appendValue(ids, (*it).countId, Generous());
			stats.contigStats[(*it).countId].i1 += 1.0 / ambig;
		}
		if (empty(ids))
		{
			std::cerr << "ERROR: Read completely mapped in a gap (";
			std::cerr << (*m).line;
			std::cerr << ')' << std::endl;
		} else
			addCount(counts, ids, 1.0 / ambig);
		
		unsigned i = length(ids);
		if (i >= length(stats.histSingleContigs))
			fill(stats.histSingleContigs, i + 1, 0.0, Generous());
		stats.histSingleContigs[i] += 1.0 / ambig;
	}
	if (ambig >= length(stats.histSingleAmbig))
		fill(stats.histSingleAmbig, ambig + 1, 0, Generous());
	++stats.histSingleAmbig[ambig];
}


template<typename TSpec, typename TConfig, typename TMatches>
inline void
countPairedMatches(FragmentStore<TSpec, TConfig> &store, TMatches const matches[2])
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type						TAnnotation;
	typedef typename Iterator<TMatches>::Type							TMIter;
	typedef typename Iterator<TAnnotationStore>::Type					TIter;
	
	if (disablePairedEnds)
	{
		++stats.brokenMPs;
		countSingleMatches(store, matches[0]);
		countSingleMatches(store, matches[1]);
		return;
	}
	
	String<int> ids;
	TMIter m0;
	TMIter m1;
	TMIter m0End = end(matches[0], Standard());
	TMIter m1End = end(matches[1], Standard());
	
	double fragSize = 0;
	double sum = 0;
	unsigned ambig = 0;
	for (m0 = begin(matches[0], Standard()); m0 != m0End; ++m0)
		for (m1 = begin(matches[1], Standard()); m1 != m1End; ++m1)
			// reads must be mapped on different strands of the same contig
			if ((*m0).contigId == (*m1).contigId && ((*m0).forward != (*m1).forward))
			{
				int beg0 = (*m0).posBegin;
				int end0 = (*m0).posEnd;
				int beg1 = (*m1).posBegin;
				int end1 = (*m1).posEnd;

				// the left read must be on the forward and the right on the backward strand
				if ((*m0).forward == (beg0 > beg1))
					continue;
				
				if (beg1 < beg0)
				{
					int tmp = beg0;
					beg0 = beg1;
					beg1 = tmp;
					tmp = end0;
					end0 = end1;
					end1 = tmp;
				}
				double x = end1 - beg0;
				sum += cdf(fragmentDistr, x + width) - cdf(fragmentDistr, x - width);
				fragSize += x;
				++ambig;
			}

	if (sum == 0)
	{
		std::cerr<<matches[0][0].line<<" broken"<<std::endl;
		++stats.brokenMPs;
		countSingleMatches(store, matches[0]);
		countSingleMatches(store, matches[1]);
		return;
	}

	fragSize /= ambig;
	stats.fragMean += fragSize;
	stats.fragStd += (fragSize - fragmentDistr.mean())*(fragSize - fragmentDistr.mean());
	++stats.fragments;
	
	m0 = begin(matches[0], Standard());
	m1 = begin(matches[1], Standard());
	for (m0 = begin(matches[0], Standard()); m0 != m0End; ++m0)
		for (m1 = begin(matches[1], Standard()); m1 != m1End; ++m1)
			if ((*m0).contigId == (*m1).contigId && ((*m0).forward != (*m1).forward))
			{
				int beg0 = (*m0).posBegin;
				int end0 = (*m0).posEnd;
				int beg1 = (*m1).posBegin;
				int end1 = (*m1).posEnd;

				// the left read must be on the forward and the right on the backward strand
				if ((*m0).forward == (beg0 > beg1))
					continue;
				
				if (beg1 < beg0)
				{
					int tmp = beg0;
					beg0 = beg1;
					beg1 = tmp;
					tmp = end0;
					end0 = end1;
					end1 = tmp;					
				}
				
				double x = end1 - beg0;
				double val = (cdf(fragmentDistr, x + width) - cdf(fragmentDistr, x - width)) / sum;

				TIter it = begin(store.annotationStore, Standard()) + dContigToAnno[(*m0).contigId];
				TIter itEnd = begin(store.annotationStore, Standard()) + dContigToAnno[(*m0).contigId + 1];
				clear(ids);
				
				for (; it != itEnd; ++it)
				{
					if ((*m0).posEnd <= (int)(*it).beginPos) break;
					if ((*m0).posBegin >= (int)(*it).endPos) continue;
					appendValue(ids, (*it).countId, Generous());
					stats.contigStats[(*it).countId].i1 += val;
				}

				if (empty(ids))
				{
					std::cerr << "ERROR: Read completely mapped in a gap (";
					std::cerr << (*m0).line;
					std::cerr << ')' << std::endl;
					continue;
				}
				appendValue(ids, -1, Generous());
				unsigned old = length(ids);

				it = begin(store.annotationStore, Standard()) + dContigToAnno[(*m1).contigId];
				itEnd = begin(store.annotationStore, Standard()) + dContigToAnno[(*m1).contigId + 1];

				for (; it != itEnd; ++it)
				{
					if ((*m1).posEnd <= (int)(*it).beginPos) break;
					if ((*m1).posBegin >= (int)(*it).endPos) continue;
					appendValue(ids, (*it).countId, Generous());
					stats.contigStats[(*it).countId].i1 += val;
				}

				if (old == length(ids))
				{
					std::cerr << "ERROR: Read completely mapped in a gap (";
					std::cerr << (*m1).line;
					std::cerr << ')' << std::endl;
				}
				else
					addCount(counts, ids, val);
				
				unsigned i = length(ids);
				if (i >= length(stats.histMPContigs))
					fill(stats.histMPContigs, i + 1, 0.0, Generous());
				stats.histMPContigs[i] += val;
			}
	if (ambig >= length(stats.histMPAmbig))
		fill(stats.histMPAmbig, ambig + 1, 0, Generous());
	++stats.histMPAmbig[ambig];
}



template<typename TSpec, typename TConfig>
inline void
loadAlignments(FragmentStore<TSpec, TConfig> &store, CharString const &fileName)
{
	typedef std::fstream TFile;
	
	TFile file;
	file.open(toCString(fileName), ios_base::in | ios_base::binary);

	CharString header_;
	std::string line, readName, lastReadName, transName;
	
	String<Match> matches[2];
	Match m;

	while (!_streamEOF(file))
	{
		// Parse Match
	
		getline(file, line);
		if (empty(line) || line[0] != '>')
			continue;

		std::istringstream iss1(line);
		char d;
		int locusNum = -1;
		int transNum = -1;
		int ambig = -1;
		unsigned errors = 0;
//		double confidence = 0;
		m.posBegin = -1;
		m.posEnd = -1;

		iss1 >> d >> m.posBegin >> d >> m.posEnd;

		size_t posId = line.find("id=");
		size_t posFragId = line.find(",fragId=");
		size_t posContigId = line.find("contigId=");
		size_t posAmbig = line.find("ambiguity=");
		size_t posErrors = line.find("errors=");
		size_t posConf = line.find("_Confidence_");
		if (posId == line.npos || posFragId == line.npos || posContigId == line.npos || posAmbig == line.npos || posErrors == line.npos || posConf == line.npos)
		{
			std::cerr << "Error parsing: " << line << std::endl;
			break;
		}
		
		readName = line.substr(posId + 3, posFragId - (posId + 3));	// skip "id=" and stop before ",contigId="
		if (readName.length() > 2 && readName[readName.length() - 2] == '_')
			readName.resize(readName.length() - 2);
		
		int mateNum = 0;
		if (readName.length() > 2 && readName[readName.length() - 2] == '/')
		{
			if (readName[readName.length() - 1] == '2')
				mateNum = 1;
			readName.resize(readName.length() - 2);
		}
		
		std::istringstream iss2(line.substr(posContigId + 9 + 6, posConf - (posContigId + 9 + 6)));	// skip "contigId=" and "Locus_"
		iss2 >> locusNum >> d;
		for (int i = 0; i < 11; ++i) iss2 >> d;
		iss2 >> transName;
		for (int i = 0; i < 12; ++i) iss2 >> d;
//		iss2 >> confidence;

		std::istringstream iss3(line.substr(posAmbig + 10));		// skip "ambiguity="
		iss3 >> ambig;

		std::istringstream iss4(line.substr(posErrors + 7));		// skip "errors="
		iss4 >> errors;
		if (errors >= length(stats.histMatchErrors))
		fill(stats.histMatchErrors, errors + 1, 0, Generous());
		++stats.histMatchErrors[errors];


		if (m.posBegin > m.posEnd)
		{
			int temp = m.posBegin;
			m.posBegin = m.posEnd;
			m.posEnd = temp;
			m.forward = false;
		} else
			m.forward = true;
		
		TTransName2Num::iterator it = transName2Num.find(transName);
		if (it == transName2Num.end())
		{
			std::cerr << "ERROR: No annotation for transcript: " << transName << std::endl;
			continue;
		}
		
		transNum = it->second;
		
		if (locusNum >= (int)length(transToDContig) || transNum >= (int)length(transToDContig[locusNum]))
		{
			std::cerr << "ERROR: No annotation for read: " << line << std::endl;
			continue;
		}
		m.contigId = transToDContig[locusNum][transNum];
		m.line = line;
		
		// Count Matches
		
		if (lastReadName != readName)
		{
			if (!empty(matches[0]))
			{
				++stats.reads;
				if (!empty(matches[1]))
				{
					++stats.reads;
					++stats.matepairs;
					countPairedMatches(store, matches);
				} else
					countSingleMatches(store, matches[0]);
			}
			else
			{
				if (!empty(matches[1]))
				{
					++stats.reads;
					countSingleMatches(store, matches[1]);
				}
			}
			clear(matches[0]);
			clear(matches[1]);
			lastReadName = readName;
		}
		appendValue(matches[mateNum], m, Generous());
		++stats.matches;
	}
	if (!empty(matches[0]))
	{
		++stats.reads;
		if (!empty(matches[1]))
		{
			++stats.reads;
			++stats.matepairs;
			countPairedMatches(store, matches);
		} else
			countSingleMatches(store, matches[0]);
	}
	else
	{
		if (!empty(matches[1]))
		{
			++stats.reads;
			countSingleMatches(store, matches[1]);
		}
	}
	file.close();
}

template<typename TSpec, typename TConfig>
inline void
dumpResults(FragmentStore<TSpec, TConfig> &store, CharString const &prefix)
{
	std::fstream file;
	file << fixed;
	
	// dump node names
	CharString fileName = prefix;
	append(fileName, ".node");
	file.open(toCString(fileName), ios_base::out | ios_base::binary);
	for (unsigned i = 0; i < length(mContigToNode); ++i)
		file << "Node_" << mContigToNode[i] << std::endl;
	file.close();

	// dump transcript names
	fileName = prefix;
	append(fileName, ".trans");
	file.open(toCString(fileName), ios_base::out | ios_base::binary);
	for (unsigned i = 0; i < length(dContigToTrans); ++i)
		file << dContigToTrans[i] << std::endl;
	file.close();	
	
	// dump indicator matrices
	String<int> headerNodes, nodes;
	for (unsigned i = 0; i < length(transToDContig); ++i)
	{
		std::ostringstream ss;
		ss.fill('0');
		ss << prefix << ".mat/" << setw(3) << (i / 1000) << "xxx/" << setw(6) << i << ".mat";
		
		file.open(ss.str().c_str(), ios_base::out | ios_base::binary);
		
		clear(headerNodes);
		for (unsigned j = 0; j < length(transToDContig[i]); ++j)
		{
			unsigned l = transToDContig[i][j];
			for (int k = dContigToAnno[l]; k < dContigToAnno[l + 1]; ++k)
				appendValue(headerNodes, store.annotationStore[k].countId, Generous());
			std::sort(begin(headerNodes, Standard()), end(headerNodes, Standard()));
		}

		file << 'T';	
		for (unsigned j = 0; j < length(headerNodes); ++j)
		{
			if (j > 0 && headerNodes[j] == headerNodes[j-1]) continue;
			file << '\t' << headerNodes[j];
		}
		file << std::endl;
		
		for (unsigned j = 0; j < length(transToDContig[i]); ++j)
		{
			unsigned l = transToDContig[i][j];
			file << l;
			clear(nodes);
			for (int k = dContigToAnno[l]; k < dContigToAnno[l + 1]; ++k)
				appendValue(nodes, store.annotationStore[k].countId, Generous());
			std::sort(begin(nodes, Standard()), end(nodes, Standard()));

			unsigned k = 0;
			for (unsigned m = 0; m < length(headerNodes); ++m)
			{
				if (m > 0 && headerNodes[m] == headerNodes[m-1]) continue;
				if (k < length(nodes) && headerNodes[m] == nodes[k])
				{
					++k;
					file << "\t1";
				} else
					file << "\t0";
			}
			file << std::endl;
		}
		file.close();
	}

	// dump statistics
	fileName = prefix;
	append(fileName, ".stat");
	file.open(toCString(fileName), ios_base::out | ios_base::binary);
	file << "contig\tcount\tlength" << std::endl;
	for (unsigned i = 0; i < length(stats.contigStats); ++i)
		file << i << '\t' << setprecision(2) << stats.contigStats[i].i1 << '\t' << stats.contigStats[i].i2 << std::endl;
	file.close();

	// dump all connections with counts
	fileName = prefix;
	append(fileName, ".cnt");
	file.open(toCString(fileName), ios_base::out | ios_base::binary);
	file << "connection\tposterior-count" << std::endl;
	for (unsigned i = 0; i < length(counts); ++i)
	{
		for (unsigned j = 0; j < length(counts[i]); ++j)
		{
			file << i;
			for (unsigned k = 0; k < length(counts[i][j].ids); ++k)
			{
				if (counts[i][j].ids[k] == -1)	// ^ gap in mate pair
				{
					file << '^';
					++k;
				} else
					file << '-';
				file << counts[i][j].ids[k];
			}
			file << '\t' << setprecision(2) << counts[i][j].count << std::endl;
		}
	}
	file.close();
}

/*
template<typename TSpec, typename TConfig>
inline void
loadAlignments(FragmentStore<TSpec, TConfig> & store, CharString const &fileName)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type TAnnotation;
	typedef std::fstream TFile;
	
	TFile file;
	file.open(toCString(fileName), ios_base::in | ios_base::binary);

	typename Value<TFile>::Type c = _streamGet(file);
	CharString header_;
	std::string header;

	while (!_streamEOF(file))
	{
		_parse_readIdentifier(file, c);		// read name
		_parse_skipWhitespace(file, c);
		_parse_readNumber(file, c);			// read begin
		_parse_skipWhitespace(file, c);
		_parse_readNumber(file, c);			// read end
		_parse_skipWhitespace(file, c);
		_streamGet(file);					// F/R
		_parse_skipWhitespace(file, c);

		header_ = _parse_readIdentifier(file, c);
		assign(header, suffix(header_, 7));
		std::istringstream iss(header);

		char d;
		int locusNum = -1;
		int transNum = -1;
		double confidence = 0;

		iss >> locusNum >> d;
		for (int i = 0; i < 11; ++i) iss >> d;
		iss >> transNum >> d;
		for (int i = 0; i < 11; ++i) iss >> d;
		iss >> confidence;

		_parse_skipWhitespace(file, c);
		int posBeg = _parse_readNumber(file, c);	// begin position in transcript
		_parse_skipWhitespace(file, c);
		int posEnd = _parse_readNumber(file, c);	// end position in transcript

		_parse_skipLine(file, c);
	}
	std::string header;

	file.close();
}
*/

///////////////////////////////////////////////////////////////////////////////
////// main
///////////////////////////////////////////////////////////////////////////////

int main( int argc, const char *argv[] ) 
{	
	CommandLineParser	parser;
	FragmentStore<>		store;			// stores all of the tables

	//////////////////////////////////////////////////////////////////////////////
	// Define options
	addTitleLine(parser, "*****************************************");
	addTitleLine(parser, "***       Transcript Quantifier       ***");
	addTitleLine(parser, "*** (c) Copyright 2009 by David Weese ***");
	addTitleLine(parser, "*****************************************");
	addUsageLine(parser, "[OPTION]... <mapping result file>");
	
	addOption(parser, CommandLineOption("d", "definition",    "read contig ordering file", OptionType::String | OptionType::Label | OptionType::List));
	addOption(parser, CommandLineOption("o", "output-prefix", "prefix for all output files (.mat, .stat, .cnt, .node, .trans)", OptionType::String | OptionType::Label, "result"));
	addOption(parser, CommandLineOption("m", "mean",          "fragment size mean value", OptionType::Double | OptionType::Label, 220));
	addOption(parser, CommandLineOption("s", "std",           "fragment size standard deviation", OptionType::Double | OptionType::Label, 30));
	addOption(parser, CommandLineOption("w", "width",         "pdf integral radius", OptionType::Double | OptionType::Label, width));
	addOption(parser, CommandLineOption("",  "single",        "interpret paired-end as single matches", OptionType::Bool));
	addHelpLine(parser, "");
	
	if (argc == 1)
	{
		shortHelp(parser, cerr);	// print short help and exit
		return 0;
	}
	
	!parse(parser, argc, argv, cerr);
	initStats();

	double mean = 0;
	double sd = 0;
	disablePairedEnds = isSetLong(parser, "single");
	getOptionValueLong(parser, "mean", mean);
	getOptionValueLong(parser, "std", sd);
	getOptionValueLong(parser, "width", width);
	fragmentDistr = boost::math::normal(mean, sd);

	unsigned defCount = length(getOptionValuesShort(parser, "d"));
	for (unsigned i = 0; i < defCount; ++i)
		loadTranscriptAnnotation(store, getOptionValuesShort(parser, "d")[i]);

//	printAnnotation(store);
	unsigned argCount = argumentCount(parser);
	for (unsigned i = 0; i < argCount; ++i)
		loadAlignments(store, getArgumentValues(parser)[i]);

	if (argCount > 0)
		printStats();

	CharString prefix;
	getOptionValueShort(parser, "o", prefix);
	dumpResults(store, prefix);

	return 0;
}
