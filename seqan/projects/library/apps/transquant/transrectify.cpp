#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/store.h>
#include <seqan/misc/misc_cmdparser.h>

using namespace seqan;
using namespace std;

struct GTFLine
{
	CharString	contig;
	CharString	app;
	CharString	type;
	__int64		beginPos, endPos;
	CharString	extra1;
	char		orientation; 
	CharString	extra2;
	CharString	ids;
};

template <typename TStream>
void printLine(TStream &stream, GTFLine const &line)
{
	stream << line.contig << '\t' << line.app << '\t' << line.type << '\t' << line.beginPos << '\t' << line.endPos << '\t';
	stream << line.extra1 << '\t' << line.orientation << '\t' << line.extra2 << '\t' << line.ids << std::endl;
}

template <typename TExons, typename TFragStore>
void merge(TExons &lines, int beginExon, TFragStore &fragStore, int minDistance)
{
	int merged = 1;
	int contigId = -1;
	for (int j = beginExon + 1; j < (int)length(lines); ++j)
	{
		GTFLine &prevLine = lines[j - 1];
		GTFLine &line = lines[j];
		if (prevLine.beginPos < line.beginPos && prevLine.endPos + minDistance > line.beginPos)
		{
			if (getIdByName(fragStore.contigNameStore, line.contig, contigId))
			{
				std::cerr << "remove gap " << infix(fragStore.contigStore[contigId].seq, prevLine.endPos, prevLine.endPos + 2);
				std::cerr << "..." << infix(fragStore.contigStore[contigId].seq, line.beginPos - 3, line.beginPos - 1);
				std::cerr << "  length:" << (line.beginPos - prevLine.endPos - 1) << std::endl;
			}
			prevLine.endPos = line.endPos;
			erase(lines, j);
			--j;
			++merged;
		}
		else if (prevLine.beginPos > line.beginPos && line.endPos + minDistance > prevLine.beginPos)
		{
			if (getIdByName(fragStore.contigNameStore, line.contig, contigId))
			{
				std::cerr << "remove gap " << infix(fragStore.contigStore[contigId].seq, line.endPos, line.endPos + 2);
				std::cerr << "..." << infix(fragStore.contigStore[contigId].seq, prevLine.beginPos - 3, prevLine.beginPos - 1);
				std::cerr << "  length:" << (prevLine.beginPos - line.endPos - 1) << std::endl;
			}
			prevLine.beginPos = line.beginPos;
			erase(lines, j);
			--j;
			++merged;
		}
	}
	if (merged > 1)
		std::cerr << "merged " << merged << " exons of " << lines[beginExon].ids << std::endl;
}

template <typename TExons, typename TFragStore>
int rectify(TExons &lines, int beginExon, int endExon, TFragStore &fragStore, int minDifference)
{
	if (beginExon == endExon) return 0;
	
	int forwardIntrons = 0;
	int reverseIntrons = 0;
	int contigId = -1;
	for (int j = beginExon; j < endExon; ++j)
	{
		GTFLine &line = lines[j];
		if (getIdByName(fragStore.contigNameStore, line.contig, contigId))
		{
			if (line.beginPos > 3)
			{
				if (infix(fragStore.contigStore[contigId].seq, line.beginPos - 3, line.beginPos - 1) == "AG")
					forwardIntrons += (j == beginExon)? 0: 3;
				if (infix(fragStore.contigStore[contigId].seq, line.beginPos - 3, line.beginPos - 1) == "AC")
					reverseIntrons += (j == beginExon)? 0: 3;
			}
			if (line.endPos + 2 < (int)length(fragStore.contigStore[contigId].seq))
			{
				if (infix(fragStore.contigStore[contigId].seq, line.endPos, line.endPos + 2) == "GT")
					forwardIntrons += (j + 1 == endExon)? 0: 3;
				if (infix(fragStore.contigStore[contigId].seq, line.endPos, line.endPos + 2) == "CT")
					reverseIntrons += (j + 1 == endExon)? 0: 3;
			}
		}
	}
	if (-minDifference < forwardIntrons - reverseIntrons && forwardIntrons - reverseIntrons < minDifference)
	{
		std::cerr << "no orientation (fwd:" << forwardIntrons << " bwd:" << reverseIntrons << ") for " << lines[beginExon].ids << std::endl;
		return 0;
	}
	else 
	{
		int changed = 0;
		char orientation = (forwardIntrons > reverseIntrons)? '+': '-';
		for (int j = beginExon; j < endExon; ++j)
		{
			if (lines[j].orientation != orientation) 
			{
				++changed;
				lines[j].orientation = orientation;
			}
		}
		if (changed == 0)
			std::cerr << "good (" << orientation << ')';
		else
			std::cerr << "modified (" << ((orientation=='+')? '-':'+') << " -> " << orientation << ')';
		std::cerr << " (fwd:" << forwardIntrons << " bwd:" << reverseIntrons << ") for transcript:" << lines[beginExon].ids << std::endl;
		if (orientation == '+')
			return 1;
		else
			return -1;
	}
}

///////////////////////////////////////////////////////////////////////////////
template <typename TFragStore>
int rectify(CharString const &src, CharString const &dst, TFragStore &fragStore, int minDistance)
{
	std::ifstream srcFile;
	std::ofstream dstFile;
	
	srcFile.open(toCString(src), ios_base::in | ios_base::binary);
	dstFile.open(toCString(dst), ios_base::out | ios_base::binary);


	string rest;
	char c = _streamGet(srcFile);
	
	
	String<GTFLine> transcriptLines;
	GTFLine line;
	
	CharString lastLocus, lastTranscript;
	CharString locus, transcript;
	string ids;
	int firstExon = 0;
	int forwardTrans = 0;
	int reverseTrans = 0;
	bool hasIntrons = false;
	
	while (!_streamEOF(srcFile))
	{
		line.contig = _parseReadIdentifier(srcFile, c);
		_parseSkipBlanks(srcFile, c);
		line.app = _parseReadIdentifier(srcFile, c);
		_parseSkipBlanks(srcFile, c);
		line.type = _parseReadIdentifier(srcFile, c);
		_parseSkipBlanks(srcFile, c);
		line.beginPos = _parseReadNumber(srcFile, c);
		_parseSkipBlanks(srcFile, c);
		line.endPos = _parseReadNumber(srcFile, c);
		_parseSkipBlanks(srcFile, c);
		line.extra1 = _parseReadIdentifier(srcFile, c);
		_parseSkipBlanks(srcFile, c);
		line.orientation = c;
		c = _streamGet(srcFile);
		_parseSkipBlanks(srcFile, c);
		line.extra2 = _parseReadIdentifier(srcFile, c);
		_parseSkipBlanks(srcFile, c);
		line.ids = _parseReadFilepath(srcFile, c);
		assign(ids, line.ids);
		_parseSkipLine(srcFile, c);
		
		if (line.type != "Exon")
		{
			printLine(dstFile, line);
			continue;
		}

		size_t locusP = ids.find("gene_id \"");
		size_t transP = ids.find("\"; transcript_id \"");
		if (locusP == ids.npos || transP == ids.npos)
		{
			std::cerr << "Error parsing: " << ids << std::endl;
			break;
		}
		locus = ids.substr(locusP + 9, transP - (locusP + 9));
		transcript = ids.substr(transP + 18, ids.size() - (transP + 18 + 2));
//		std::cerr<<"##"<<line.contig<<"##"<<line.type <<"##"<<locus<<"##"<<transcript<<"##"<<std::endl;

		if (transcript != lastTranscript)
		{
			merge(transcriptLines, firstExon, fragStore, minDistance);
			if (length(transcriptLines) - firstExon > 1)
			{
				hasIntrons = true;
				int r = rectify(transcriptLines, firstExon, length(transcriptLines), fragStore, 2);
				if (r > 0) ++forwardTrans;
				if (r < 0) ++reverseTrans;
			}
			firstExon = length(transcriptLines);
			lastTranscript = transcript;
		}
		if (locus != lastLocus)
		{
			if (!empty(transcriptLines) && hasIntrons)
			{
				if (forwardTrans != 0 && reverseTrans != 0)
					std::cerr << "ambig. locus orientation (fwdTrans:" << forwardTrans << " bwdTrans:" << reverseTrans << ") for locus:" << locus << std::endl;
				else if (forwardTrans == 0 && reverseTrans == 0)
					std::cerr << "no locus orientation (fwdTrans:" << forwardTrans << " bwdTrans:" << reverseTrans << ") for locus:" << locus << std::endl;
			}
			for (unsigned i = 0; i < length(transcriptLines); ++i)
				printLine(dstFile, transcriptLines[i]);

			clear(transcriptLines);
			lastLocus = locus;
			firstExon = 0;
			forwardTrans = 0;
			reverseTrans = 0;
			hasIntrons = false;
		}
		appendValue(transcriptLines, line);
		
		int contigId = -1;
		if (getIdByName(fragStore.contigNameStore, line.contig, contigId))
		{
			std::cerr << "\tLeft:"<<infix(fragStore.contigStore[contigId].seq, line.beginPos-5, line.beginPos + 1)<<" Right:"<<infix(fragStore.contigStore[contigId].seq, line.endPos-2, line.endPos+4)<<" Locus:"<<locus<<" Trans:"<<transcript<<" Begin:"<<line.beginPos<<" End:"<<line.endPos<<std::endl;
		} else {
			std::cerr << "Ref. contig not found:" << line.contig << std::endl;
		}

	}

	if (!empty(transcriptLines))
	{
		if (length(transcriptLines) - firstExon > 1)
		{
			hasIntrons = true;
			int r = rectify(transcriptLines, firstExon, length(transcriptLines), fragStore, 2);
			if (r > 0) ++forwardTrans;
			if (r < 0) ++reverseTrans;
		}

		if (hasIntrons)
		{
			if (forwardTrans != 0 && reverseTrans != 0)
				std::cerr << "ambig. locus orientation (fwdTrans:" << forwardTrans << " bwdTrans:" << reverseTrans << ") for locus:" << locus << std::endl;
			if (forwardTrans == 0 && reverseTrans == 0)
				std::cerr << "no locus orientation (fwdTrans:" << forwardTrans << " bwdTrans:" << reverseTrans << ") for locus:" << locus << std::endl;
		}
		for (unsigned i = 0; i < length(transcriptLines); ++i)
			printLine(dstFile, transcriptLines[i]);
	}
	
	srcFile.close();
	dstFile.close();
	return 0;
}

///////////////////////////////////////////////////////////////////////////////
////// main
///////////////////////////////////////////////////////////////////////////////

int main( int argc, const char *argv[] ) 
{	
	CommandLineParser	parser;
	FragmentStore<>		store;			// stores all of the tables

	std::string rev = "$Revision$";
	addVersionLine(parser, "TransRectify version 1.0 20100901 [" + rev.substr(11, 4) + "]");

	//////////////////////////////////////////////////////////////////////////////
	// Define options
	addTitleLine(parser, "*****************************************");
	addTitleLine(parser, "***       Transcript Rectifier        ***");
	addTitleLine(parser, "*** (c) Copyright 2010 by David Weese ***");
	addTitleLine(parser, "*****************************************");
	addUsageLine(parser, "[OPTION]... <transcript annotation Gff file>");
	
	addOption(parser, CommandLineOption("r",  "reference",        "reference sequence file", OptionType::String | OptionType::Label | OptionType::List));
	addOption(parser, CommandLineOption("o",  "output",           "rectified annotation output file", OptionType::String | OptionType::Label));
	addOption(parser, CommandLineOption("md", "minimal-distance", "minimal distance between adjacent exons", OptionType::Int | OptionType::Label, 20));
	addLine(parser, "");
	requiredArguments(parser, 1);
	
	if (!parse(parser, argc, argv, cerr)) return 0;
	
	FragmentStore<> fragStore;
	for (unsigned i = 0; i < length(getOptionValuesShort(parser, "r")); ++i)
		loadContigs(fragStore, getOptionValuesShort(parser, "r")[i]);
		
	CharString dstFileName;
	int minDistance = 20;
	getOptionValueShort(parser, "md", minDistance);
	getOptionValueShort(parser, "o", dstFileName);
	rectify(getArgumentValues(parser)[0], dstFileName, fragStore, minDistance);

	return 0;
}
