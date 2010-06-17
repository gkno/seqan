///A tutorial about global alignments.
#include <fstream>
#include <iostream>
#include <seqan/align.h>
#include <seqan/graph_align.h>
#include <seqan/store.h>
#include <seqan/store/store_io_gff.h>
#include <seqan/store/store_io_ucsc.h>

using namespace seqan;

namespace seqan
{

struct SVGFile
{
	std::fstream file;

	Pair<int> cursor;
	String<CharString> style;
	
	int dpSequence;
	int dpMatrix;
	int dpTrace;
	
	int text;
	
	int readForward;
	int readReverse;
	int readText;
	int rulerTextTicks;
	int rulerTextLabel;

	friend inline void svgWriteHeader(SVGFile &svg);
	friend inline void svgWriteFooter(SVGFile &svg);
	
	SVGFile(char const *fileName):
		cursor(0,0) 
	{
		readText = dpSequence = length(style);
		appendValue(style, "style=\"text-anchor:middle; font-family:Verdana,sans-serif;\" font-size=\"18px\"");
		text = length(style);
		appendValue(style, "style=\"text-anchor:middle; font-family:Courier New,Verdana,sans-serif; font-weight:bold;\" font-size=\"20px\"");
		rulerTextTicks = length(style);
		appendValue(style, "style=\"text-anchor:middle; font-family:Verdana,sans-serif;\" font-size=\"6px\"");
		rulerTextLabel = length(style);
		appendValue(style, "style=\"text-anchor:left; font-family:Verdana,sans-serif; font-weight:bold;\" font-size=\"8px\"");

		dpMatrix = length(style);
		appendValue(style, "style=\"stroke:lightgray;stroke-width:1;\" marker-end=\"url(#startMarkerNormal)\"");
		dpTrace = length(style);
		appendValue(style, "style=\"stroke:black;stroke-width:2;\" marker-end=\"url(#startMarkerBold)\"");

		readForward = length(style);
		appendValue(style, "style=\"stroke:darkblue;stroke-width:3;\"");
		appendValue(style, "style=\"stroke:lightblue;stroke-width:1;\"");
		appendValue(style, "style=\"stroke:darkblue;stroke-width:3;\" marker-end=\"url(#startMarkerForward)\"");
		appendValue(style, "style=\"stroke:lightblue;stroke-width:1;\" marker-end=\"url(#startMarkerForward)\"");
		readReverse = length(style);
		appendValue(style, "style=\"stroke:darkred;stroke-width:3;\"");
		appendValue(style, "style=\"stroke:salmon;stroke-width:1;\"");
		appendValue(style, "style=\"stroke:darkred;stroke-width:3;\" marker-end=\"url(#startMarkerReverse)\"");
		appendValue(style, "style=\"stroke:salmon;stroke-width:1;\" marker-end=\"url(#startMarkerReverse)\"");

		file.open(fileName, std::ios_base::out);
		svgWriteHeader(*this);
	}

	~SVGFile()
	{
		svgWriteFooter(*this);
		file.close();
	}
};

inline void svgWriteHeader(SVGFile &svg)
{
	svg.file << "<?xml version=\"1.0\"?>" << std::endl;
	svg.file << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"" << std::endl;
	svg.file << "  \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << std::endl;
	svg.file << std::endl;
	svg.file << "<svg version=\"1.1\"" << std::endl;
	svg.file << "    xmlns=\"http://www.w3.org/2000/svg\"" << std::endl;
	svg.file << "    xmlns:xlink=\"http://www.w3.org/1999/xlink\">" << std::endl;
    svg.file << "<defs>" << std::endl;

	// trace back arrow markers
    svg.file << "	<g id=\"arrowMarker\" stroke-linecap=\"round\">" << std::endl;
    svg.file << "		<line x1=\"-6\" y1=\"1.5\" x2=\"0\" y2=\"0\" />" << std::endl;
    svg.file << "		<line x1=\"-6\" y1=\"-1.5\" x2=\"0\" y2=\"0\" />" << std::endl;
    svg.file << "	</g>" << std::endl;
    svg.file << "	<marker id=\"startMarkerNormal\" markerWidth=\"6\" markerHeight=\"3\" style=\"overflow:visible\" orient=\"auto\" markerUnits=\"userSpaceOnUse\">" << std::endl;
    svg.file << "		<use xlink:href=\"#arrowMarker\" stroke=\"lightgray\" />" << std::endl;
    svg.file << "	</marker>" << std::endl;
    svg.file << "	<marker id=\"startMarkerBold\" markerWidth=\"6\" markerHeight=\"4\" style=\"overflow:visible\" orient=\"auto\" markerUnits=\"userSpaceOnUse\">" << std::endl;
    svg.file << "		<use xlink:href=\"#arrowMarker\" stroke-width=\"2\" stroke=\"black\" />" << std::endl;
    svg.file << "	</marker>" << std::endl;

	// read mapping arrow markers
    svg.file << "	<marker id=\"startMarkerForward\" markerWidth=\"10\" markerHeight=\"4\" style=\"overflow:visible\" orient=\"auto\" markerUnits=\"userSpaceOnUse\">" << std::endl;
    svg.file << "		<polyline points=\"-8,0 -9,-6 3,0 -9,6 -8,0\" fill=\"darkblue\" />" << std::endl;
    svg.file << "	</marker>" << std::endl;
    svg.file << "	<marker id=\"startMarkerReverse\" markerWidth=\"10\" markerHeight=\"4\" style=\"overflow:visible\" orient=\"auto\" markerUnits=\"userSpaceOnUse\">" << std::endl;
    svg.file << "		<polyline points=\"-8,0 -9,-6 3,0 -9,6 -8,0\" fill=\"darkred\" />" << std::endl;
    svg.file << "	</marker>" << std::endl;
	
    svg.file << " </defs>" << std::endl;
	svg.file << std::endl;
}

inline void svgWriteFooter(SVGFile &svg)
{
	svg.file << std::endl;
	svg.file << "</svg>" << std::endl;
}


template <typename TChar>
inline void
_streamPut(SVGFile & svg, TChar character)
{
	if (convert<char>(character) == '\n')
	{
		++svg.cursor.i2;
		svg.cursor.i1 = 0;
	} else if (convert<char>(character) == '\t')
	{
		svg.cursor.i1 = (svg.cursor.i1 & ~7) + 8;
	} else
	{
		if (convert<char>(character) != ' ')
		{
			svg.file << "<g transform=\"translate(" << svg.cursor.i1*20+10 << ',' << svg.cursor.i2*20+10 << ")\"><text y=\"0.3em\" " << svg.style[svg.text] << '>';
			_streamPut(svg.file, convert<char>(character));
			svg.file << "</text></g>" << std::endl;
		}
		++svg.cursor.i1;
	}
}

template <typename TFormatTag, typename TContigGaps, typename TContigName>
inline void _printContig(
	SVGFile &svg, 
	Tag<TFormatTag> const &format,
	AlignedReadLayout &, 
	TContigGaps &contigGaps,
	TContigName const &contigName)
{
	typedef typename Iterator<TContigGaps, Standard>::Type TContigIterator;

	TContigIterator cit = begin(contigGaps, Standard());
	TContigIterator citEnd = end(contigGaps, Standard());	
	for (__int64 ofs = 0; cit != citEnd; ++cit, ++ofs)
	{
		if (!isGap(cit))
		{
			if (ofs == 0)
			{
				svg.file << "<g transform=\"translate(" << ofs*20+2 << ',' << svg.cursor.i2*20+10 << ")\"><text y=\"0.3em\" " << svg.style[svg.rulerTextLabel] << '>';
				svg.file << contigName << "</text></g>" << std::endl;
			}
			
			__int64 seqPos = cit.current.seqPos + 1;
			if (seqPos % 5 == 0)
			{
				if (seqPos % 10 == 0)
				{
					if (ofs >= 5)
					{
						svg.file << "<g transform=\"translate(" << ofs*20+10 << ',' << svg.cursor.i2*20+10 << ")\"><text y=\"0.3em\" " << svg.style[svg.rulerTextTicks] << '>';
						svg.file << seqPos << "</text></g>" << std::endl;
					}
					
					svg.file << "<line x1=\"" << ofs*20+10 << "\" x2=\"" << ofs*20+10 << "\" ";
					svg.file << "y1=\"" << svg.cursor.i2*20+12 << "\" y2=\"" << svg.cursor.i2*20+15 << "\" ";
				} else {
					svg.file << "<line x1=\"" << ofs*20+10 << "\" x2=\"" << ofs*20+10 << "\" ";
					svg.file << "y1=\"" << svg.cursor.i2*20+12 << "\" y2=\"" << svg.cursor.i2*20+15 << "\" ";
				}
				svg.file << "stroke-width=\"1\" stroke=\"gray\" />" << std::endl;
			}
		}
	}
	_streamPut(svg, '\n');
	
	int savedStyle = svg.text;
	svg.text = svg.readText;
	write(svg, contigGaps, "", format);
	svg.text = savedStyle;
}

template <typename TFormatTag, typename TContigGaps, typename TReadGaps, typename TAlignedRead, typename TLine>
inline void _printRead(
	SVGFile &svg, 
	Tag<TFormatTag> const &,
	AlignedReadLayout &layout, 
	TContigGaps &contigGaps,
	TReadGaps &readGaps,
	TAlignedRead &alignedRead,
	TLine line)
{
	typedef typename Iterator<TContigGaps, Standard>::Type TContigIterator;
	typedef typename Iterator<TReadGaps, Standard>::Type TIterator;

	__int64 x;
	__int64 xEnd;
	int style, arrow;
	const char *first;
	const char *second;
	if (alignedRead.beginPos < alignedRead.endPos)
	{
		xEnd = alignedRead.beginPos * 20;
		x = xEnd + 5;
		first = "<line x1=\"";
		second = "\" x2=\"";
		style = svg.readForward;
		arrow = 0;
	} else {
		xEnd = alignedRead.endPos * 20;
		x = xEnd + 10;
		first = "<line x2=\"";
		second = "\" x1=\"";
		style = svg.readReverse;
		arrow = 2;
	}
	line = svg.cursor.i2 * 20 + 10;
	
	if (length(layout.mateCoords) <= alignedRead.pairMatchId)
		fill(layout.mateCoords, alignedRead.pairMatchId + 1, Pair<int>(-1,-1));
	else
	{
		if (layout.mateCoords[alignedRead.pairMatchId].i2 != -1)
		{
			Pair<__int64> a(alignedRead.beginPos * 20, line);
			Pair<__int64> b(layout.mateCoords[alignedRead.pairMatchId]);
			if (a.i1 < b.i1)
			{
				Pair<__int64> tmp = a;
				a = b;
				b = tmp;
			}
			__int64 dx = (b.i1 - a.i1);
			__int64 dy = (b.i2 - a.i2);
			
			svg.file << "<path d=\"M " << a.i1 << ',' << a.i2;
			svg.file << " C " << a.i1+dy/10 << ',' << a.i2-dx/10;
			svg.file << ' ' << b.i1+dy/10 << ',' << b.i2-dx/10;
			svg.file << ' ' << b.i1 << ',' << b.i2;
			svg.file << "\" stroke-width=\"2\" stroke=\"black\" stroke-opacity=\"0.2\" fill=\"none\"/>";
		}
		else
			layout.mateCoords[alignedRead.pairMatchId] = Pair<int>(alignedRead.beginPos * 20, line);
	}
	

	TContigIterator cit = begin(contigGaps, Standard()) + _min(alignedRead.beginPos, alignedRead.endPos);
	TIterator it = begin(readGaps, Standard());
	TIterator itEnd = end(readGaps, Standard());
	int lastWasGap = -1;
	int inGap;
	
	for (; it != itEnd; ++it, ++cit, xEnd += 20)
	{
		inGap = isGap(it);
		if (lastWasGap != inGap || inGap != isGap(cit) || (!inGap && convert<Dna5>(*cit) != convert<Dna5>(*it)))
		{
			if (x < xEnd && lastWasGap != -1)
			{
				svg.file << first << x << "\" y1=\"" << line << second << xEnd; 
				svg.file << "\" y2=\"" << line << "\" " << svg.style[style + arrow + lastWasGap] << " />" << std::endl;
				arrow = 0;
			}
			lastWasGap = inGap;
			x = xEnd;
			if (!inGap && convert<Dna5>(*cit) != convert<Dna5>(*it))
			{
				svg.file << "<g transform=\"translate(" << xEnd + 10 << ',' << line << ")\"><text y=\"0.3em\" " << svg.style[svg.readText] << '>';
				_streamPut(svg.file, convert<char>(*it));
				svg.file << "</text></g>" << std::endl;
				x += 20;
				arrow = 0;
			}
		}
	}
	if (x < xEnd && lastWasGap != -1)
	{
		if (alignedRead.beginPos < alignedRead.endPos)
		{
			arrow = 2;
			xEnd -= 10;
		} else
			xEnd -= 5;
		svg.file << first << x << "\" y1=\"" << line << second << xEnd;
		svg.file << "\" y2=\"" << line << "\" " << svg.style[style + arrow + lastWasGap] << " />" << std::endl;
	}
}


//////////////////////////////////////////////////////////////////////////////
// stream operators
//////////////////////////////////////////////////////////////////////////////
/*
template <typename TSource>
inline SVGFile &
operator << (SVGFile & target, 
			 TSource  source)
{
SEQAN_CHECKPOINT
	write(target, source);
	return target;
}
*/

inline SVGFile &
operator << (SVGFile & target, char source)
{
SEQAN_CHECKPOINT
	write(target, source);
	return target;
}

inline SVGFile &
operator << (SVGFile & target, char const *source)
{
SEQAN_CHECKPOINT
	write(target, source);
	return target;
}


template <typename TStringSet, typename TTrace, typename TIndexPair>
void
_align_needleman_wunsch_matrix(SVGFile& svg,
							  TStringSet const& str,
							  TTrace const& trace,
							  TIndexPair const&)
{
SEQAN_CHECKPOINT
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename Value<TTrace>::Type TTraceValue;

	// TraceBack values
	const int Diagonal = 0; 
	const int Horizontal = 1; 
	const int Vertical = 2;
	
	// Initialization
	TSize numCols = length(str[0]);
	TSize numRows = length(str[1]);

	// Print trace matrix	
	for (TSize pos0 = 0; pos0 < numCols; ++pos0)
	{
		for (TSize pos1 = 0; pos1 < numRows; ++pos1)
		{
			int tv =(int)trace[pos0*numRows + pos1];
				if (tv & (1 << Diagonal))
					svg.file << "<line x2=\"" << pos0*20+10 << "\" y2=\"" << pos1*20+10 << "\"" << " x1=\"" << (pos0+1)*20+10 << "\" y1=\"" << (pos1+1)*20+10 << "\" " << svg.style[svg.dpMatrix] << " />" << std::endl;
				
				if (tv & (1 << Horizontal))
					svg.file << "<line x2=\"" << pos0*20+10 << "\" y2=\"" << (pos1+1)*20+10 << "\"" << " x1=\"" << (pos0+1)*20+10 << "\" y1=\"" << (pos1+1)*20+10 << "\" " << svg.style[svg.dpMatrix] << " />" << std::endl;

				if (tv & (1 << Vertical))
					svg.file << "<line x2=\"" << (pos0+1)*20+10 << "\" y2=\"" << pos1*20+10 << "\"" << " x1=\"" << (pos0+1)*20+10 << "\" y1=\"" << (pos1+1)*20+10 << "\" " << svg.style[svg.dpMatrix] << " />" << std::endl;
		}
	}

	// Print sequences
	for (TSize pos0 = 0; pos0 < numCols; ++pos0)
		svg.file << "<g transform=\"translate(" << pos0*20+30 << ",10)\"><text y=\"0.3em\" " << svg.style[svg.dpSequence] << '>' << str[0][pos0] << "</text></g>" << std::endl;

	for (TSize pos1 = 0; pos1 < numRows; ++pos1)
		svg.file << "<g transform=\"translate(10," << pos1*20+30 << ")\"><text y=\"0.3em\" " << svg.style[svg.dpSequence] << '>' << str[1][pos1] << "</text></g>" << std::endl;
}

template <typename TStringSet, typename TId, typename TPos, typename TTraceValue>
inline void
_align_trace_print(SVGFile& svg,
				   TStringSet const&,
				   TId const,
				   TPos pos1,
				   TId const,
				   TPos pos2,
				   TPos const segLen,
				   TTraceValue const tv)
{
	// TraceBack values
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;

	svg.file << "<line x2=\"" << pos1*20+10 << "\" y2=\"" << pos2*20+10 << "\"";
	if (tv == Diagonal) {
		pos1 += segLen; pos2 += segLen;
	} else if (tv == Horizontal) {
		pos1 += segLen;
	} else if (tv == Vertical) {
		pos2 += segLen;
	}
	svg.file << " x1=\"" << pos1*20+10 << "\" y1=\"" << pos2*20+10;
	svg.file << "\" " << svg.style[svg.dpTrace] << " />" << std::endl;
}


}

int main()
{
///Two DNA sequences that shall be aligned.
	typedef String<char> TSequence;
	TSequence seq1 = "atcgaatgcgimmiga";
	TSequence seq2 = "actcgttgca";
///Scoring objects are used to define a scoring scheme.
///In this case, affine gap costs with match = 0, mismatch = -1, gapextend = -1 and gapopen = -2.
	Score<int> score(0, -1, -1, -1);
///Example 1: We use @Class.Align@ to align the two sequences. 
///Since we do not specify an @Tag.Global Alignment Algorithms|algorithm tag@ when we call @Function.globalAlignment@, 
///a suitable algorithm (@Tag.Global Alignment Algorithms|Gotoh@) is automatically choosen.
	Align<TSequence, ArrayGaps> align;
	resize(rows(align), 2);
	assignSource(row(align, 0), seq1);
	assignSource(row(align, 1), seq2);

	SVGFile svg("/Volumes/weese/small_case.svg");

	typedef	FragmentStore<> TFragmentStore;
	typedef  TFragmentStore::TAnnotationStore TAnnotationStore;
	typedef   Iterator<TFragmentStore, AnnotationTree<> >::Type TAnnoTreeIter;
	typedef   Value<TAnnotationStore>::Type TAnnotation;
	
	TFragmentStore store;
//	loadContigs(store, "/Users/weese/Documents/Development/samtools/examples/ex1.fa", true);
//	loadContigs(store, "/Volumes/weese/seqan/projects/library/cmake/make/tmp/drosophila-part.fasta", true);
	loadContigs(store, "/Volumes/weese/small_case.fa", true);

	double t1=sysTime();
	std::ifstream file1("/tmp/knownIsoforms.txt", ::std::ios_base::in | ::std::ios_base::binary);
	std::ifstream file1b("/tmp/knownGene.txt", ::std::ios_base::in | ::std::ios_base::binary);
	read(file1, store, UCSC());
	double t1b=sysTime();
	read(file1b, store, UCSC());
	double t2=sysTime();
	std::cout << "time for reading knownGene: "<< t1b-t1 << std::endl;
	std::cout << "time for reading knownToLocusLink: "<< t2-t1b << std::endl;

	std::ofstream file2c("/tmp/output.gff", ::std::ios_base::out | ::std::ios_base::binary);
	write(file2c, store, GFF());
	file2c.close();

	std::ofstream file2("/tmp/output.ucsc", ::std::ios_base::out | ::std::ios_base::binary);
	write(file2, store, UCSC());
	file2.close();
	std::cout << "time for writing GFF: "<< sysTime()-t2 << std::endl;
	
	return 0;
	
//	std::ifstream file("/Users/weese/Documents/Development/samtools/examples/ex1.sam", ::std::ios_base::in | ::std::ios_base::binary);
//	std::ifstream file("/Volumes/weese/seqan/projects/library/cmake/make/tmp/SRR027007.828.1-part.sam", ::std::ios_base::in | ::std::ios_base::binary);	
	std::ifstream file("/Volumes/weese/small_case.sam", ::std::ios_base::in | ::std::ios_base::binary);	
	read(file, store, SAM());
	
	AlignedReadLayout layout;
	layoutAlignment(layout, store);
	for (unsigned i=0;i<length(store.contigStore);++i)
	{
		printAlignment(std::cout, Raw(), layout, store, i, 280, 320, 0, 100);
		printAlignment(svg, Raw(), layout, store, i, 0, 2000, 0, 100);
	}
/*
	TAnnoTreeIter it = begin(store, AnnotationTree<>());
	
	int line=svg.cursor.i2;
	while (!atEnd(it))
	{
		TAnnotation &anno = getAnnotation(it);
		if (anno.contigId == length(store.contigStore)-1 && anno.typeId > 1)
		{			
			if (anno.typeId == TFragmentStore::ANNO_GENE)
				line=svg.cursor.i2;
			
			if (anno.typeId == TFragmentStore::ANNO_MRNA)
				++line;
				
			if (anno.typeId == TFragmentStore::ANNO_EXON)
			{
				anno.beginPos -=41260;
				anno.endPos -=41260;
				svg.file << "<polygon style=\"stroke:black;stroke-width:1;fill:lightgreen\" points=\"";
				svg.file << anno.beginPos*20 << ' ' << line*20+3  << ", ";
				svg.file << anno.beginPos*20 << ' ' << line*20+17 << ", ";
				svg.file << anno.endPos*20   << ' ' << line*20+17 << ", ";
				svg.file << anno.endPos*20   << ' ' << line*20+3  << ", ";
				svg.file << anno.beginPos*20 << ' ' << line*20+3  << "\" />" << std::endl;
			}
		}
		
		goNext(it);
	}

*/
	SVGFile svgTrace("trace1.svg");
	
	_globalAlignment(svgTrace, stringSet(align), score, AlignConfig<>(), NeedlemanWunsch());
//	globalAlignment(align, score);
	::std::cout << align << ::std::endl;

	return 0;

	::std::cout << "Score = " << globalAlignment(align, score) << ::std::endl;
	::std::cout << align << ::std::endl;
///Example 2: We now choose explicitely the algorithm @Tag.Global Alignment Algorithms|MyersHirschberg@.
///Since this algorithm always works on Levenshtein distance, $score$ is ignored here.
///Therefore, this algorithm computes a different alignment and returns a different score.
	::std::cout << "Score = " << globalAlignment(align, score, MyersHirschberg()) << ::std::endl;
	::std::cout << align << ::std::endl;
///Example 3: We now do the same as in case 1, but now we use an @Spec.Alignment Graph@ for storing the alignment.
///Here we use @Tag.Global Alignment Algorithms|Gotoh's algorithm@.
	typedef StringSet<TSequence, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;

	TStringSet string_set;
	appendValue(string_set, seq1);
	appendValue(string_set, seq2);
	TAlignmentGraph alignment_graph(string_set);

	::std::cout << "Score = " << globalAlignment(alignment_graph, score, Gotoh()) << ::std::endl;
	::std::cout << alignment_graph << ::std::endl;
	return 0;
}
