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

//#define RENAME_NODES
#define MANY_BINS
#define SINGLE_MAT_FILE
#define OMIT_UNCOVERED_MATS

using namespace seqan;
using namespace std;

///////////////////////////////////////////////////////////////////////////////

struct Count
{
	String<int>	ids;
	double		count;
};

struct Match {
//	unsigned		contigId;
	unsigned		locusId;
	unsigned		transId;
	int				posBegin;
	int				posEnd;
	unsigned char	errors;
	bool			forward;	// true..forward strand, false..reverse complement strand
	std::string		line;
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

struct CoverageBin
{
	unsigned			binLength;		// size in bps all transcripts with length >= binLength are counted
	unsigned			transcripts;	// coverages of how many transcripts were counted
#ifdef MANY_BINS
	String<double, Array<64> >	coverage;		// stores the sum of coverages for each position i=0,...,binLength-1
#else
	String<double>		coverage;		// stores the sum of coverages for each position i=0,...,binLength-1
#endif
	
	CoverageBin():
		binLength(0),
		transcripts(0) {}

	CoverageBin(unsigned _binLength):
		binLength(_binLength),
		transcripts(0) 
	{
#ifdef MANY_BINS
		resize(coverage, 64, 0);
#else
		resize(coverage, binLength, 0);
#endif
	}
};

String<CoverageBin> bins;

typedef StringSet<CharString>	TNames;
typedef NameStoreCache<TNames>	TNameCache;
typedef std::map<string,Pair<int,int> >	TTransNameCache;

// DContig - Dave's contig is a contig (e.g. chromosome) in the fragment store
// MContig - Marcel's contig is a node in the de Bruijn graph

TNames						locusNames;
TNames						transNames;		// DContig id -> original name, e.g. Locus_0_Transcript_0_Confidence_0.800
String<String<int> >		transToDContig;	// (locus,transcript) -> DContig id

TNameCache					locusNameCache(locusNames);
TTransNameCache				transNameCache;

String<int>					dContigToAnno;	// DContig id -> its first node annotation
String<int>					dContigLength;	// transcript length
String<String<Count> >		counts;			// node[0] -> ((node[1..],count), (node[1..],count), ...)

String<int>					nodeToMContig;	// node number -> MContig id (compact renumeration)
String<int>					mContigToNode;	// MContig id -> node number (reverse map)
int							lastMContig = 0;

bool disablePairedEnds;
boost::math::normal fragmentDistr;
double width = 0.5;
unsigned histSize = 18;


template <typename TSequence>
void removeEqualElements(TSequence &seq)
{
	typedef typename Iterator<TSequence, Standard>::Type TIter;

	TIter itBegin = begin(seq, Standard());
	TIter itEnd = end(seq, Standard());
	TIter src = itBegin;
	TIter dst = itBegin;
	for (; src != itEnd; ++src)
	{
		if (*dst != *src)
		{
			++dst;
			*dst = *src;
		}
	}
	if (itBegin != itEnd)
		resize(seq, dst - itBegin + 1);
}

template <typename TId, typename TName>
inline void 
appendLocus(TId & locusId, TName const & locusName)
{
	if (!getIdByName(locusNames, locusName, locusId, locusNameCache))
	{
		locusId = length(locusNames);
		appendName(locusNames, locusName, locusNameCache);
	}
}

template <typename TId, typename TName>
inline int 
appendTrans(TId & locusId, TId & transId, TName const & transName)
{
	if (locusId >= (int)length(transToDContig))
		resize(transToDContig, locusId + 1, Generous());		
	transId = length(transToDContig[locusId]);
	if (transId >= (int)length(transToDContig[locusId]))
		resize(transToDContig[locusId], transId + 1, 0, Generous());
	
	int contigId = length(transNames);
	transToDContig[locusId][transId] = contigId;
	transNameCache[transName] = Pair<int>(locusId, transId);
	appendValue(transNames, transName, Generous());
	return contigId;
}

template <typename TId>
inline bool 
getTransId(TId & locusId, TId & transId, std::string const & transName)
{
	TTransNameCache::iterator it = transNameCache.find(transName);
	if (it != transNameCache.end())
	{
		locusId = it->second.i1;
		transId = it->second.i2;
		return true;
	}
	return false;
}


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
		resize(s, histSize + 2, 0);
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
	std::string transName;

	if (empty(dContigToAnno))
		appendValue(dContigToAnno, 1);

	while (!_streamEOF(file))
	{
		while (c != '>' && !_streamEOF(file))
			_parseSkipLine(file, c);
		if (_streamEOF(file)) break;
		
		c = _streamGet(file);
		assign(transName, _parseReadIdentifier(file, c));
		size_t locusPos = transName.find("Locus_");
		size_t transPos = transName.find("Transcript_");
		if (locusPos == transName.npos || transPos == transName.npos)
		{
			std::cerr << "Error parsing: " << transName << std::endl;
			break;
		}

		_parseSkipLine(file, c);

		int locusNum = -1;
		int transNum = -1;
		
		appendLocus(locusNum, transName.substr(locusPos + 6, transPos - 1 - (locusPos + 6)));
		unsigned contigId = appendTrans(locusNum, transNum, transName);

		int beginPos = 0;
		int endPos = 0;
		while (!_streamEOF(file))
		{
			if (!_parseIsDigit(c) && c !='-')
			{
				std::cerr << "Ignoring entry " << transName << std::endl;
				break;
			}
			int nodeId = _parseReadNumber(file, c);			
			if (c != ':')
				std::cerr << "HUH1? " << transName << std::endl;

			c = _streamGet(file);	// :
			// quick fix for Hugues' bug
			endPos = _parseReadNumber(file, c);
			//int endPos = beginPos + _parseReadNumber(file, c);
			if (nodeId < 0) nodeId = -nodeId;

			// rename nodeIds
			if (nodeId >= (int)length(nodeToMContig))
				resize(nodeToMContig, nodeId + 1, -1, Generous());
#ifdef RENAME_NODES
			if (nodeToMContig[nodeId] == -1)
			{
				nodeToMContig[nodeId] = lastMContig++;
				appendValue(mContigToNode, nodeId, Generous());
			}
			nodeId = nodeToMContig[nodeId];
#endif
			if ((int)length(stats.contigStats) <= nodeId)
				resize(stats.contigStats, nodeId + 1, Pair<double, int>(0, 0));
			stats.contigStats[nodeId] = Pair<double, int>(0, endPos - beginPos);

			TAnnotation a;
			a.beginPos = beginPos;
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
			int gapLen = _parseReadNumber(file, c);	
			if (c != ')') 
				std::cerr<<"HUH4?"<<std::endl;

			c = _streamGet(file);	// )
			c = _streamGet(file);	// -
			c = _streamGet(file);	// >

			beginPos = endPos + gapLen;
		}
		++contigId;
		appendValue(dContigLength, endPos, Generous());
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
		unsigned contigId = transToDContig[(*m).locusId][(*m).transId];
		TIter it = begin(store.annotationStore, Standard()) + dContigToAnno[contigId];
		TIter itEnd = begin(store.annotationStore, Standard()) + dContigToAnno[contigId + 1];
		clear(ids);
		
		int transPosBegin = 0;
		int transPosEnd = 0;
		double overlapRelBegin = 0.0;
		double overlapRelEnd = 0.0;
		for (; it != itEnd; ++it)
		{
			if ((*m).posEnd <= (int)(*it).beginPos) 
			{
				// node is right of our read
				break;
			}
			if ((*m).posBegin >= (int)(*it).endPos) 
			{
				// node is left of our read
				transPosBegin += (int)(*it).endPos - (int)(*it).beginPos;
				transPosEnd += (int)(*it).endPos - (int)(*it).beginPos;
				continue;
			}
			
			// here we have an overlap
		
			if ((*m).posBegin >= (int)(*it).beginPos)
			{
				// read begin is in the node
				transPosBegin += (*m).posBegin - (int)(*it).beginPos;
			}
			
			if ((*m).posEnd <= (int)(*it).endPos)
			{
				// read end is in the node
				transPosEnd += (*m).posEnd - (int)(*it).beginPos;
				overlapRelEnd = 1.0;
			} else
			{
				// read end is right of the node
				transPosEnd += (int)(*it).endPos - (int)(*it).beginPos;
				overlapRelEnd = ((int)(*it).endPos - (*m).posBegin) / (double)((*m).posEnd - (*m).posBegin);
			}
			
			appendValue(ids, (*it).countId, Generous());
			stats.contigStats[(*it).countId].i1 += (overlapRelEnd - overlapRelBegin) / ambig;
			overlapRelBegin = overlapRelEnd;
		}
		if (empty(ids))
		{
			std::cerr << "ERROR: Read completely mapped in a gap (";
			std::cerr << (*m).line;
			std::cerr << ')' << std::endl;
		} else
		{
			addCount(counts, ids, 1.0 / ambig);
			if (length(transToDContig[(*m).locusId]) == 1 && (*m).errors < 2 && ambig == 1)
			{
				unsigned transLen = dContigLength[contigId];
//				std::cout << (*m).line << '\t'<<transPosBegin<<'\t'<<transPosEnd<<'\t'<<transLen<<std::endl;

				unsigned binNo;
#ifdef MANY_BINS
				binNo = transLen;
				if (length(bins) <= binNo)
					resize(bins, binNo + 1);
				bins[binNo].binLength = binNo;
				resize(bins[binNo].coverage, 64, 0);
#else
				for (binNo = 0; binNo < (int)length(bins); ++binNo)
					if (bins[binNo].binLength > transLen)
						break;
				if (binNo != 0) --binNo;
#endif				
				CoverageBin &bin = bins[binNo];
				transPosBegin = (transPosBegin * length(bin.coverage)) / transLen;
				transPosEnd = (transPosEnd * length(bin.coverage)) / transLen;
				for (; transPosBegin < transPosEnd; ++transPosBegin)
					bin.coverage[transPosBegin] += 1;
			}
		}
		
		unsigned i = length(ids);
		if (i >= length(stats.histSingleContigs))
			resize(stats.histSingleContigs, i + 1, 0.0, Generous());
		stats.histSingleContigs[i] += 1.0 / ambig;
	}
	if (ambig >= length(stats.histSingleAmbig))
		resize(stats.histSingleAmbig, ambig + 1, 0, Generous());
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
			if ((*m0).locusId == (*m1).locusId && (*m0).transId == (*m1).transId && (*m0).forward != (*m1).forward)
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
			if ((*m0).locusId == (*m1).locusId && (*m0).transId == (*m1).transId && (*m0).forward != (*m1).forward)
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
				unsigned contigId0 = transToDContig[(*m0).locusId][(*m0).transId];

				TIter it = begin(store.annotationStore, Standard()) + dContigToAnno[contigId0];
				TIter itEnd = begin(store.annotationStore, Standard()) + dContigToAnno[contigId0 + 1];
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
				unsigned contigId1 = transToDContig[(*m1).locusId][(*m1).transId];

				it = begin(store.annotationStore, Standard()) + dContigToAnno[contigId1];
				itEnd = begin(store.annotationStore, Standard()) + dContigToAnno[contigId1 + 1];

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
					resize(stats.histMPContigs, i + 1, 0.0, Generous());
				stats.histMPContigs[i] += val;
			}
	if (ambig >= length(stats.histMPAmbig))
		resize(stats.histMPAmbig, ambig + 1, 0, Generous());
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
		size_t posErrors = line.find(",errors=");
//		size_t posConf = line.find("Confidence_");
		if (posId == line.npos || posFragId == line.npos || posContigId == line.npos || posAmbig == line.npos || posErrors == line.npos)
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
		
		transName = line.substr(posContigId + 9, posErrors - (posContigId + 9));	// skip "contigId="
/*		
		std::istringstream iss2(line.substr(posContigId + 9 + 6, posConf - (posContigId + 9 + 6)));	// skip "contigId=" and "Locus_"
		iss2 >> locusNum >> d;
		for (int i = 0; i < 11; ++i) iss2 >> d;
		iss2 >> transName;
		for (int i = 0; i < 12; ++i) iss2 >> d;
//		iss2 >> confidence;
*/
		std::istringstream iss3(line.substr(posAmbig + 10));		// skip "ambiguity="
		iss3 >> ambig;

		std::istringstream iss4(line.substr(posErrors + 8));		// skip ",errors="
		iss4 >> errors;
		if (errors >= length(stats.histMatchErrors))
		resize(stats.histMatchErrors, errors + 1, 0, Generous());
		++stats.histMatchErrors[errors];

		m.errors = errors;
		if (m.posBegin > m.posEnd)
		{
			int temp = m.posBegin;
			m.posBegin = m.posEnd;
			m.posEnd = temp;
			m.forward = false;
		} else
			m.forward = true;
		
		if (!getTransId(locusNum, transNum, transName))
		{
			std::cerr << "ERROR: No annotation for transcript: " << transName << std::endl;
			continue;
		}
				
		if (locusNum >= (int)length(transToDContig) || transNum >= (int)length(transToDContig[locusNum]))
		{
			std::cerr << "ERROR: No annotation for read: " << line << std::endl;
			continue;
		}
		m.locusId = locusNum;
		m.transId = transNum;
//		m.contigId = transToDContig[locusNum][transNum];
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
#ifdef RENAME_NODES
	append(fileName, ".node");
	file << "id\tname\tnodes\tlength" << std::endl;
	file.open(toCString(fileName), ios_base::out | ios_base::binary);
	for (unsigned i = 0; i < length(mContigToNode); ++i)
		file << "Node_" << mContigToNode[i] << std::endl;
	file.close();
#endif

	// dump transcript names
	fileName = prefix;
	append(fileName, ".locus");
	file.open(toCString(fileName), ios_base::out | ios_base::binary);
	file << "id\tname\tisoforms\tnodes\tcoverage" << std::endl;

	String<int> nodes;
	for (unsigned i = 0; i < length(locusNames); ++i)
	{
		file << i << '\t' << locusNames[i] << '\t' << length(transToDContig[i]);

		// get and sort all nodes of the locus i
		clear(nodes);
		for (unsigned j = 0; j < length(transToDContig[i]); ++j)
		{
			unsigned l = transToDContig[i][j];
			for (int k = dContigToAnno[l]; k < dContigToAnno[l + 1]; ++k)
				appendValue(nodes, store.annotationStore[k].countId, Generous());
		}

		// sort and remove nodes that are contained in multiple isoforms
		std::sort(begin(nodes, Standard()), end(nodes, Standard()));
		removeEqualElements(nodes);

		// count coverage over all segments of the locus
		double locusCoverage = 0.0;
		for (unsigned j = 0; j < length(nodes); ++j)
			locusCoverage += stats.contigStats[nodes[j]].i1;

		file << '\t' << length(nodes) << '\t' << locusCoverage << std::endl;
	}
	file.close();	
	
	// dump transcript names
	fileName = prefix;
	append(fileName, ".trans");
	file.open(toCString(fileName), ios_base::out | ios_base::binary);
	file << "id\tname\tnodes\tlength\tcoverage" << std::endl;
	for (unsigned i = 0; i < length(transNames); ++i)
	{
		file << i << '\t' << transNames[i] << '\t' << (dContigToAnno[i+1]-dContigToAnno[i]) << '\t' << dContigLength[i];

		// count coverage over all segments of the isoform
		double transCoverage = 0.0;
		for (int k = dContigToAnno[i]; k < dContigToAnno[i + 1]; ++k)
			transCoverage += stats.contigStats[store.annotationStore[k].countId].i1;
		file << '\t' << transCoverage << std::endl;
	}
	file.close();	
	
	// dump indicator matrices
	String<int> headerNodes;
#ifdef SINGLE_MAT_FILE
	fileName = prefix;
	append(fileName, ".mat");
	file.open(toCString(fileName), ios_base::out | ios_base::binary);
#endif
	for (unsigned i = 0; i < length(transToDContig); ++i)
	{
		// get and sort all nodes of the locus i
		clear(headerNodes);
		for (unsigned j = 0; j < length(transToDContig[i]); ++j)
		{
			unsigned l = transToDContig[i][j];
			for (int k = dContigToAnno[l]; k < dContigToAnno[l + 1]; ++k)
				appendValue(headerNodes, store.annotationStore[k].countId, Generous());
		}

		// sort and remove nodes that are contained in multiple isoforms
		std::sort(begin(headerNodes, Standard()), end(headerNodes, Standard()));
		removeEqualElements(headerNodes);

		// count coverage over all segments of the locus
		double locusCoverage = 0.0;
		for (unsigned j = 0; j < length(headerNodes); ++j)
			locusCoverage += stats.contigStats[headerNodes[j]].i1;

#ifdef OMIT_UNCOVERED_MATS
		if (locusCoverage == 0.0) continue;
#endif		

		
#ifdef SINGLE_MAT_FILE
		if (empty(transToDContig[i])) continue;
//		file << ">" << i << "\tcvrg=" << locusCoverage << std::endl;
		file << ">" << i << std::endl;
#else
		std::ostringstream ss;
		ss.resize('0');
		ss << prefix << ".mat/" << setw(3) << (i / 1000) << "xxx/" << setw(6) << i << ".mat";
		
		file.open(ss.str().c_str(), ios_base::out | ios_base::binary);
#endif		

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
#ifndef SINGLE_MAT_FILE
		file.close();
#endif
	}
#ifdef SINGLE_MAT_FILE
	file.close();
#endif

	// dump statistics
	fileName = prefix;
	append(fileName, ".stat");
	file.open(toCString(fileName), ios_base::out | ios_base::binary);
	file << "contig\tcount\tlength" << std::endl;
	for (unsigned i = 0; i < length(stats.contigStats); ++i)
		file << i << '\t' << setprecision(2) << stats.contigStats[i].i1 << '\t' << stats.contigStats[i].i2 << std::endl;
	file.close();
	
	// dump sequencing bias
	fileName = prefix;
	append(fileName, ".bias");
	file << "bin-length\tcoverage" << std::endl;
	file.open(toCString(fileName), ios_base::out | ios_base::binary);
	// begin with bin no.2 and always output the same number of coverage points
	for (unsigned i = 2, s = 1; i < length(bins); ++i)
		if (bins[i].binLength > 0)
		{
			file << bins[i].binLength;
			for (unsigned j = 0; j < length(bins[i].coverage); j+=s)
				file << '\t' << bins[i].coverage[j];
			file << std::endl;
#ifndef MANY_BINS
			s<<=1;
#endif
		}
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
		_parseReadIdentifier(file, c);		// read name
		_parseSkipWhitespace(file, c);
		_parseReadNumber(file, c);			// read begin
		_parseSkipWhitespace(file, c);
		_parseReadNumber(file, c);			// read end
		_parseSkipWhitespace(file, c);
		_streamGet(file);					// F/R
		_parseSkipWhitespace(file, c);

		header_ = _parseReadIdentifier(file, c);
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

		_parseSkipWhitespace(file, c);
		int posBeg = _parseReadNumber(file, c);	// begin position in transcript
		_parseSkipWhitespace(file, c);
		int posEnd = _parseReadNumber(file, c);	// end position in transcript

		_parseSkipLine(file, c);
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

	std::string rev = "$Revision$";
	addVersionLine(parser, "TransQuant version 1.0 20100901 [" + rev.substr(11, 4) + "]");

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
	addOption(parser, CommandLineOption("w", "width",         "probability distribution function integral radius", OptionType::Double | OptionType::Label, width));
	addOption(parser, CommandLineOption("",  "single",        "interpret paired-end as single matches", OptionType::Bool));
	addHelpLine(parser, "");
	requiredArguments(parser, 1);
	
	if (!parse(parser, argc, argv, cerr)) return 0;
	initStats();

	double mean = 0;
	double sd = 0;

#ifndef MANY_BINS
	for (unsigned i = 4, l = 16; i < 16; ++i, l<<=1)
		append(bins, CoverageBin(l));
#endif
	
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
