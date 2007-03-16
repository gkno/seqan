#ifndef SEQAN_HEADER_FIND_AHOCORASICK_H
#define SEQAN_HEADER_FIND_AHOCORASICK_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// AhoCorasick
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.AhoCorasick:
..summary: Multiple exact string matching using Aho-Corasick.
..general:Class.Pattern
..cat:Pattern Matching
..signature:Pattern<TNeedle, AhoCorasick>
..param.TNeedle:The needle type, a string of keywords.
...type:Class.String
..remarks.text:The types of the keywords in the needle container and the haystack have to match.
..remarks.text:Matching positions do not come in order because we report beginning positions of matches.
..remarks.text:Likewise, if multiple keywords match at a given position no pre-specified order is guaranteed.
*/

///.Class.Pattern.param.TSpec.type:Spec.AhoCorasick

struct _AhoCorasick;
typedef Tag<_AhoCorasick> AhoCorasick;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, AhoCorasick> {
//____________________________________________________________________________
private:
	Pattern(Pattern const& other);
	Pattern const& operator=(Pattern const & other);

//____________________________________________________________________________
public:
	typedef typename Size<TNeedle>::Type TSize;
	typedef typename Value<TNeedle>::Type TKeyword;
	typedef typename Value<TKeyword>::Type TAlphabet;
	typedef Graph<Automaton<TAlphabet> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	
	Holder<TNeedle> data_needle;
	String<TVertexDescriptor> data_supplyMap;
	String<String<TSize> > data_terminalStateMap;
	TGraph data_graph;

	// To restore the automaton after a hit
	String<TSize> data_endPositions;	// All remaining keyword indices
	TSize data_keywordIndex;			// Current keyword that produced a hit
	TSize data_needleLength;			// Last length of needle to reposition finder
	TVertexDescriptor data_lastState;   // Last state in the trie

//____________________________________________________________________________

	Pattern() {
	}

	template <typename TNeedle2>
	Pattern(TNeedle2 const & ndl)
	{
		SEQAN_CHECKPOINT
		setHost(*this, ndl);
	}

	~Pattern() {
		SEQAN_CHECKPOINT
	}
//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Host Metafunctions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
struct Host< Pattern<TNeedle, AhoCorasick> >
{
	typedef TNeedle Type;
};

template <typename TNeedle>
struct Host< Pattern<TNeedle, AhoCorasick> const>
{
	typedef TNeedle const Type;
};


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
inline void
_createAcTrie(Pattern<TNeedle, AhoCorasick> & me)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TNeedle>::Type TSize;
	typedef typename Value<TNeedle>::Type TKeyword;
	typedef typename Value<TKeyword>::Type TAlphabet;
	typedef Graph<Automaton<TAlphabet> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();

	// Create regular trie
	createTrie(me.data_graph,me.data_terminalStateMap, host(me));

	// Create parent map
	String<TVertexDescriptor> parentMap;
	String<TAlphabet> parentCharMap;
	initVertexMap(me.data_graph,parentMap);
	initVertexMap(me.data_graph,parentCharMap);
	for(unsigned int i = 0;i<length(parentMap);++i) {
		assignProperty(parentMap, i, nilVal);
	}
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	TEdgeIterator itEd(me.data_graph);
	for(;!atEnd(itEd);goNext(itEd)) {
		assignProperty(parentMap, targetVertex(itEd), sourceVertex(itEd));
		assignProperty(parentCharMap, targetVertex(itEd), label(itEd));
	}

	// Build AC
	TVertexDescriptor root = getRoot(me.data_graph);
	initVertexMap(me.data_graph,me.data_supplyMap);
	assignProperty(me.data_supplyMap, root, nilVal);

	// Bfs Traversal
	typedef typename Iterator<TGraph, BfsIterator>::Type TBfsIterator;
	TBfsIterator it(me.data_graph,root);
	for(;!atEnd(it);goNext(it)) {
		if (atBegin(it)) continue;
		TVertexDescriptor parent = getProperty(parentMap, *it);
		TAlphabet sigma = getProperty(parentCharMap, *it);
		TVertexDescriptor down = getProperty(me.data_supplyMap, parent);
		while ((down != nilVal) &&
			(getSuccessor(me.data_graph, down, sigma) == nilVal)) 
		{
			down = getProperty(me.data_supplyMap, down);
		}
		if (down != nilVal) {
			assignProperty(me.data_supplyMap, *it, getSuccessor(me.data_graph, down, sigma));
			String<unsigned int> endPositions = getProperty(me.data_terminalStateMap, getProperty(me.data_supplyMap, *it));
			if (!empty(endPositions)) {
				String<unsigned int> endPositionsCurrent = getProperty(me.data_terminalStateMap, *it);
				typedef typename Iterator<String<unsigned int> >::Type TStringIterator;
				TStringIterator sit = begin(endPositions);
				for(;!atEnd(sit);goNext(sit)) {
					appendValue(endPositionsCurrent, *sit);
				}
				assignProperty(me.data_terminalStateMap, *it, endPositionsCurrent);
			}
		} else {
			assignProperty(me.data_supplyMap, *it, root);
		}

	}
}


template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, AhoCorasick> & me, TNeedle2 const & needle) {
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(!empty(needle));
	me.data_needle = needle;
	clear(me.data_graph);
	clear(me.data_supplyMap);
	clear(me.data_terminalStateMap);
	clear(me.data_endPositions);
	me.data_keywordIndex = 0;
	_createAcTrie(me);
	me.data_needleLength = 0;
	me.data_lastState = getRoot(me.data_graph);

	/*
	fstream strm;
	strm.open(TEST_PATH "my_trie.dot", ios_base::out | ios_base::trunc);
	String<String<char> > nodeMap;
	_createTrieNodeNames(me.data_graph, me.data_terminalStateMap, nodeMap);
	String<String<char> > edgeMap;
	_createEdgeNames(me.data_graph,edgeMap);
	write(strm,me.data_graph,nodeMap,edgeMap,DotDrawing());
	strm.close();
	// Supply links
	for(unsigned int i=0;i<length(me.data_supplyMap);++i) {
		std::cout << i << "->" << getProperty(me.data_supplyMap,i) << ::std::endl;
	}
	*/
}

template <typename TNeedle, typename TNeedle2>
inline void 
setHost (Pattern<TNeedle, AhoCorasick> & me, TNeedle2 & needle)
{
	setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}

//____________________________________________________________________________

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, AhoCorasick>const>::Type & 
host(Pattern<TNeedle, AhoCorasick> & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, AhoCorasick>const>::Type & 
host(Pattern<TNeedle, AhoCorasick> const & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

//____________________________________________________________________________


template <typename TNeedle>
inline typename Size<TNeedle>::Type
position(Pattern<TNeedle, AhoCorasick> & me)
{
	return me.data_keywordIndex;
}


template <typename TFinder, typename TNeedle>
inline bool find(TFinder & finder, Pattern<TNeedle, AhoCorasick> & me) {
	SEQAN_CHECKPOINT
	typedef typename Size<TNeedle>::Type TSize;
	typedef typename Value<TNeedle>::Type TKeyword;
	typedef typename Value<TKeyword>::Type TAlphabet;
	typedef Graph<Automaton<TAlphabet> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	
	if (empty(finder))
		goBegin(finder);
	else {
		finder += me.data_needleLength;
		++finder; // Set forward the finder
	}

	// Process left-over hits
	if (!empty(me.data_endPositions)) {
		--finder; // Set back the finder
		me.data_keywordIndex = me.data_endPositions[length(me.data_endPositions)-1];
		me.data_needleLength = length(getValue(host(me), me.data_keywordIndex))-1;
		if (length(me.data_endPositions) > 1) resize(me.data_endPositions, (length(me.data_endPositions)-1));
		else clear(me.data_endPositions);
		finder -= me.data_needleLength;
		return true;
	}



	TVertexDescriptor current = me.data_lastState;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	while (!atEnd(finder)) {
		while ((getSuccessor(me.data_graph, current, *finder) == nilVal) &&
			(getProperty(me.data_supplyMap, current) != nilVal))
		{
			current = getProperty(me.data_supplyMap,current);
		}
		if (getSuccessor(me.data_graph, current, *finder) != nilVal) {
			current = getSuccessor(me.data_graph, current, *finder);
		}
		else {
			current = getRoot(me.data_graph);
		}
		me.data_endPositions = getProperty(me.data_terminalStateMap,current);
		if (!empty(me.data_endPositions)) {
			me.data_keywordIndex = me.data_endPositions[length(me.data_endPositions)-1];
			me.data_needleLength = length(getValue(host(me), me.data_keywordIndex))-1;
			if (length(me.data_endPositions) > 1) resize(me.data_endPositions, length(me.data_endPositions)-1);
			else clear(me.data_endPositions);
			me.data_lastState = current;
			finder -= me.data_needleLength;
			return true;
		}
		++finder;
	}
	return false;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_AHOCORASICK_H
