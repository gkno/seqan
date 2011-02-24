//#define SEQAN_ENABLE_DEBUG 1
#include <fstream>
#include <iostream>
#include <string>
#include <seqan/store.h>
#include <seqan/misc/misc_cmdparser.h>

using namespace seqan;

struct Match {
	unsigned		contigId;
	unsigned		transId;
	int				posBegin;
	int				posEnd;
	int				genomeBegin;
	int				genomeEnd;
    int             mateDelta;
    int             matePos;
	unsigned char	errors;
	bool			transForward;	// true..forward strand, false..reverse complement strand
    bool            genomeForward;
	std::string		line;
    std::string     readSeq;
    std::string     readId;
    String<int>     cigar;
};

bool disablePairedEnds;

// Get Exon Boundaries

template <typename TExonBoundString, typename TFragmentStore, typename TId>
void getExonBoundaries(TExonBoundString &exonBounds, TFragmentStore &store, TId transId)
{
	typedef typename TFragmentStore::TAnnotationStore TAnnotationStore;
	typedef typename Iterator<TFragmentStore, AnnotationTree<> >::Type TAnnoTreeIter;
	typedef typename Value<TAnnotationStore>::Type TAnnotation;
	
	clear(exonBounds);
	TAnnoTreeIter exonIter = begin(store, AnnotationTree<>());
    goTo(exonIter, transId);
    if (!goDown(exonIter))
        return;
    
	while (!atEnd(exonIter))
	{
		TAnnotation &anno = getAnnotation(exonIter);
		if (anno.typeId == TFragmentStore::ANNO_EXON)
		{
			appendValue(exonBounds, anno.beginPos);
			appendValue(exonBounds, anno.endPos);
		} 
		goNextRight(exonIter);
	}
    std::sort(begin(exonBounds, Standard()), end(exonBounds, Standard()));
}

template<typename TSpec, typename TConfig, typename TOutFile, typename TMatch>
inline void
outputMatch(FragmentStore<TSpec, TConfig> &store, TOutFile &outFile, TMatch const &match, int mateNo)
{
    int bitfield = match.genomeForward? 0x00: 0x10;
    if (mateNo != -1)
    {
        if (mateNo == 0)
            bitfield = 0x40 + 0x02 + 0x01;
        else
            bitfield = 0x80 + 0x02 + 0x01;
        bitfield |= match.genomeForward? 0x20: 0x00;
    }
    
    // read name
    outFile << match.readId << '\t';

    // bitfield
    outFile << bitfield << '\t';

    // reference name
    outFile << store.contigNameStore[match.contigId] << '\t';

    // forward begin pos
    outFile << match.genomeBegin + 1 << '\t';

    // mapping qual (not avail.)
    outFile << 255 << '\t';

    // cigar string
    for (unsigned i = 0; i < length(match.cigar); ++i)
    {
        outFile << match.cigar[i];
        if ((i & 1) == 0)
            outFile << 'M';
        else
            outFile << 'N';
    }
    outFile << '\t';
    
    if (mateNo == -1)
    {
        outFile << "*\t0\t0\t";
    }
    else
    {
        outFile << "=\t";
        outFile << match.matePos + 1 << '\t';
        outFile << match.mateDelta << '\t';
    }

    outFile << match.readSeq << "\t*\n";
}

template<typename TSpec, typename TConfig, typename TOutFile, typename TMatches>
inline void
outputMatches(FragmentStore<TSpec, TConfig> &store, TOutFile &outFile, TMatches const &matches)
{
	typedef typename Iterator<TMatches>::Type TMIter;	
	TMIter m = begin(matches, Standard());
	TMIter mEnd = end(matches, Standard());

	for (; m != mEnd; ++m)
        outputMatch(store, outFile, *m, -1);
}

template<typename TSpec, typename TConfig, typename TOutFile, typename TMatches>
inline void
outputPairedMatches(FragmentStore<TSpec, TConfig> &store, TOutFile &outFile, TMatches const matches[2])
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type						TAnnotation;
	typedef typename Iterator<TMatches>::Type							TMIter;
	typedef typename Iterator<TAnnotationStore>::Type					TIter;
	
	if (disablePairedEnds)
	{
		outputMatches(store, outFile, matches[0]);
		outputMatches(store, outFile, matches[1]);
		return;
	}

	TMIter m0;
	TMIter m1;
	TMIter m0End = end(matches[0], Standard());
	TMIter m1End = end(matches[1], Standard());

	for (m0 = begin(matches[0], Standard()); m0 != m0End; ++m0)
		for (m1 = begin(matches[1], Standard()); m1 != m1End; ++m1)
			// reads must be mapped on different strands of the same contig
			if ((*m0).transId == (*m1).transId && (*m0).genomeForward != (*m1).genomeForward)
			{
				// the left read must be on the forward and the right on the backward strand
				if ((*m0).genomeForward == ((*m0).genomeBegin > (*m1).genomeBegin))
					continue;
                
                (*m0).matePos = (*m1).genomeBegin;
                (*m1).matePos = (*m0).genomeBegin;

                if ((*m0).genomeBegin < (*m1).genomeBegin)
                {
                    (*m0).mateDelta = (*m1).genomeEnd - (*m0).genomeBegin;
                    (*m1).mateDelta = (*m0).genomeBegin - (*m1).genomeEnd;
                } else 
                {
                    (*m0).mateDelta = (*m1).genomeBegin - (*m0).genomeEnd;
                    (*m1).mateDelta = (*m0).genomeEnd - (*m1).genomeBegin;
                }
                outputMatch(store, outFile, (*m0), 0);
                outputMatch(store, outFile, (*m1), 1);
			}
}

template<typename TSpec, typename TConfig, typename TMatches>
inline void
transformMatches(FragmentStore<TSpec, TConfig> &store, TMatches const &matches)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type						TAnnotation;
	typedef typename Iterator<TMatches>::Type							TMIter;
	typedef typename Iterator<TAnnotationStore>::Type					TIter;
	
	TMIter m = begin(matches, Standard());
	TMIter mEnd = end(matches, Standard());
    String<int> exonBounds;
	for (; m != mEnd; ++m)
	{
        clear(exonBounds);
        getExonBoundaries(exonBounds, store, (*m).transId);
		
		int transBegin = 0;
		int transEnd = 0;
        (*m).genomeBegin = -1;
        (*m).genomeEnd = -1;
        clear((*m).cigar);
        bool transOnForward = store.annotationStore[(*m).transId].beginPos < store.annotationStore[(*m).transId].endPos;
        (*m).genomeForward = ((*m).transForward == transOnForward);
		for (unsigned e=0; e<length(exonBounds); e+=2)
		{
            int transPosBegin = 0;
            int transPosEnd = exonBounds[e+1] - exonBounds[e];

            transBegin = transEnd;
            transEnd += transPosEnd;

			if ((*m).posEnd <= transBegin) 
			{
				// exon is right of our read
				break;
			}
			if (transEnd <= (*m).posBegin)
			{
				// node is left of our read
				continue;
			}
            
            // tb     tetb     tetb       te
            // 11111111122222222233333333333
            //       rrrrrrr
            //       mb   me
            //
            // 111111111           222222222        33333333333
            //       rrr           rrrr
			
            // here holds:  transPosBegin < (*m).posEnd  and  (*m).posBegin < transPosEnd
			// here we have an overlap
		
			if (transBegin <= (*m).posBegin)
            {
				// read begin is in the node
                transPosBegin = (*m).posBegin - transBegin;
				(*m).genomeBegin = exonBounds[e] + transPosBegin;
			}
			
			if ((*m).posEnd <= transEnd)
			{
				// read end is in the node
                transPosEnd = (*m).posEnd - transBegin;
				(*m).genomeEnd = exonBounds[e] + transPosEnd;
			}
            if (!empty((*m).cigar))
            {
                // append intron length
                appendValue((*m).cigar, exonBounds[e] - exonBounds[e-1]);
            }
                        
            // append exon overlap length
            appendValue((*m).cigar, transPosEnd - transPosBegin);
            if (empty((*m).cigar))
            {
				std::cerr << "ERROR: No exons overlap with read: " << (*m).line << std::endl;
            }
            
            (*m).mateDelta = 0;
        }
	}
}


template <typename TId, typename TFragmentStore>
inline bool 
getTransId(TId & contigId, TId & transId, std::string const & transName, TFragmentStore &fragStore)
{
    if (getIdByName(fragStore.annotationNameStore, transName, transId, fragStore.annotationNameStoreCache))
    {
        contigId = fragStore.annotationStore[transId].contigId;
        return true;
    }
	return false;
}


template <typename TRNAStream, typename TDNAStream, typename TFragmentStore>
void transformAlignments(TRNAStream &rnaAlignFile, TDNAStream &dnaAlignFile, TFragmentStore &store)
{
	CharString header_;
	std::string line, readName, lastReadName, refName;
	std::string transName;
    
	String<Match> matches[2];
	Match m;
    Match *lastMatch = NULL;
	
	while (!_streamEOF(rnaAlignFile))
	{
		// Parse Match
		
		getline(rnaAlignFile, line);
		if (empty(line) || line[0] != '>')
        {
            if (lastMatch != NULL)
            {
                lastMatch->readSeq = line;
                lastMatch = NULL;
            }
			continue;
        }
		
		std::istringstream iss1(line);
		char d;
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
        
        m.readId = readName;

        int mateNum = 0;
        if (readName.length() > 2 && (readName[readName.length() - 2] == '/' || readName[readName.length() - 2] == '_'))
        {
            if (readName[readName.length() - 1] == '2')
                mateNum = 1;
            readName.resize(readName.length() - 2);
        }
        
        refName = line.substr(posContigId + 9, posErrors - (posContigId + 9));	// skip "contigId="
        
        size_t posLocusEnd = refName.find("_Transcript_");
        size_t posTransEnd = refName.find("_Confidence_");
        if (posLocusEnd == refName.npos || posTransEnd == refName.npos)
        {
            std::cerr << "Error parsing: " << refName << std::endl;
            break;
        }
        
        transName = refName.substr(posLocusEnd + 12, posTransEnd - (posLocusEnd + 12));
        
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
        
        m.errors = errors;
        if (m.posBegin > m.posEnd)
        {
            int temp = m.posBegin;
            m.posBegin = m.posEnd;
            m.posEnd = temp;
            m.transForward = false;
        } else
            m.transForward = true;
        
        if (!getTransId(m.contigId, m.transId, transName, store))
        {
            std::cerr << "ERROR: No annotation for transcript: " << transName << std::endl;
            continue;
        }
        
        m.line = line;
        
        // Count Matches
        
        if (lastReadName != readName)
        {
            transformMatches(store, matches[0]);
            transformMatches(store, matches[1]);
            
            if (empty(matches[0]) || empty(matches[1]))
            {
                outputMatches(store, dnaAlignFile, matches[0]);
                outputMatches(store, dnaAlignFile, matches[1]);
            }
            else
                outputPairedMatches(store, dnaAlignFile, matches);

            clear(matches[0]);
            clear(matches[1]);
            lastReadName = readName;
        }
        appendValue(matches[mateNum], m, Generous());
        lastMatch = &back(matches[mateNum]);
    }
    transformMatches(store, matches[0]);
    transformMatches(store, matches[1]);
    if (empty(matches[0]) || empty(matches[1]))
    {
        outputMatches(store, dnaAlignFile, matches[0]);
        outputMatches(store, dnaAlignFile, matches[1]);
    }
    else
        outputPairedMatches(store, dnaAlignFile, matches);
}

int main(int argc, char const * argv[])
{
	typedef	FragmentStore<> TFragmentStore;
	typedef String<int> TContigExonBounds;
	typedef String<TContigExonBounds> TExonBounds;
	typedef Id<TFragmentStore>::Type TId;

	CommandLineParser parser;
	std::string rev = "$Revision: 8423 $";
	addVersionLine(parser, "TransSplice version 1.0 20100901 [" + rev.substr(11, 4) + "]");

	//////////////////////////////////////////////////////////////////////////////
	// Define options
	addTitleLine(parser, "*****************************************");
	addTitleLine(parser, "***  Transform RNA to DNA Alignments  ***");
	addTitleLine(parser, "*** (c) Copyright 2010 by David Weese ***");
	addTitleLine(parser, "*****************************************");
	addUsageLine(parser, "[OPTION]... <alignment.results> <annotation.gff> <output.sam>");
	addUsageLine(parser, "[OPTION]... <alignment.results> <knownGene.txt> <knownIsoforms.txt> <output.sam>");

	addOption(parser, CommandLineOption("",  "single",        "interpret paired-end as single matches", OptionType::Bool));
	addLine(parser, "");
//	addLine(parser, "Splice transcripts of a given annotation according to a given genome.");
//	addLine(parser, "The resulting transcripts are written to <output-prefix>/transcripts.fa.");
//	addLine(parser, "Overlapping exons on the same strand are divided and annotated as subexons (<output-prefix>/refined.gff).");
//	addLine(parser, "The sequences of subexon ids for every transcript are written to <output-prefix>/ordering.fa.");
	requiredArguments(parser, 3);
	
	if (!parse(parser, argc, argv, std::cerr)) return 0;
	
	double t0 = sysTime();
	TFragmentStore store;

	CharString rnaAlignFileName;
	CharString dnaAlignFileName;
		
	disablePairedEnds = isSetLong(parser, "single");
	
	//////////////////////////////////////////////////////////////////////////////
	// Step 1: Reading the annotation
	if (argumentCount(parser) == 3) // Gtf/Gff?
	{
		std::ifstream annotationFile(toCString(getArgumentValue(parser, 1)), ::std::ios_base::in | ::std::ios_base::binary);
		if (!annotationFile.is_open())
		{
			std::cerr << "Could not open " << getArgumentValue(parser, 1) << std::endl;
			return 1;
		}
		read(annotationFile, store, Gff());
        rnaAlignFileName = getArgumentValue(parser, 0);
        dnaAlignFileName = getArgumentValue(parser, 2);
	} 
	else 
	{
		std::ifstream knownGenes(toCString(getArgumentValue(parser, 1)), ::std::ios_base::in | ::std::ios_base::binary);
		std::ifstream knownIsoforms(toCString(getArgumentValue(parser, 2)), ::std::ios_base::in | ::std::ios_base::binary);
		if (!knownGenes.is_open())
		{
			std::cerr << "Could not open " << getArgumentValue(parser, 1) << std::endl;
			return 1;
		}
		if (!knownIsoforms.is_open())
		{
			std::cerr << "Could not open " << getArgumentValue(parser, 2) << std::endl;
			return 1;
		}
		read(knownGenes, store, Ucsc());
		read(knownIsoforms, store, Ucsc());
        rnaAlignFileName = getArgumentValue(parser, 0);
        dnaAlignFileName = getArgumentValue(parser, 3);
	}
	double t1 = sysTime();		std::cout << "Reading the annotation took "<< t1-t0 << " seconds." << std::endl;

	//////////////////////////////////////////////////////////////////////////////
	// Step 2: Transform alignments
	std::ifstream rnaAlignFile(toCString(rnaAlignFileName), ::std::ios_base::in | ::std::ios_base::binary);
	std::ofstream dnaAlignFile(toCString(dnaAlignFileName), ::std::ios_base::out | ::std::ios_base::binary);

    if (!rnaAlignFile.is_open())
    {
        std::cerr << "Could not open input file " << rnaAlignFileName << std::endl;
        return 1;
    }

    if (!dnaAlignFile.is_open())
    {
        std::cerr << "Could not open output file " << dnaAlignFileName << std::endl;
        return 1;
    }

	transformAlignments(rnaAlignFile, dnaAlignFile, store);
	rnaAlignFile.close();
	dnaAlignFile.close();
	double t2 = sysTime();		std::cout << "Writing transformed alignments took "<< t2-t1 << " seconds." << std::endl;

	return 0;
}
