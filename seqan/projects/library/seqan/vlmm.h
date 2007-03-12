#ifndef SEQAN_HEADER_VLMM_H
#define SEQAN_HEADER_VLMM_H

namespace SEQAN_NAMESPACE_MAIN
{

//this file is the initial version of the construction of
// a variable length markov model (ConstrainedTraversal)


// Defining specialized iterators

	// specializations for the ConstrainedTraversal iterator
	struct Absolute;
	struct Relative;

template <typename TSpec = Absolute>
struct ConstrainedTraversal;



template < typename TIndex,typename TSpec>
class Iter< TIndex, VSTree< TopDown< ParentLinks<ConstrainedTraversal<TSpec> > > > >:
	public Iter< TIndex, VSTree< TopDown< ParentLinks<> > > >
{
	typedef Iter< TIndex, VSTree< TopDown< ParentLinks<> > > > Base;

public:
	// although both variables are hold by the iterator
	// actually either one of them is used during traversal
	// indicated by the specialization Absolute or Relative
	unsigned AbsoluteThreshold;
	float    RelativeThreshold;

	bool	 Down;
	unsigned Up;

	Iter(TIndex &__index):
		Base(__index){}

	// use constructor for an absolute threshold
	Iter(TIndex &__index,int threshold):
		Base(__index),AbsoluteThreshold(threshold),RelativeThreshold(0.0),Down(false),Up(0){}

	// use constructor for a relative threshold
	Iter(TIndex &__index,float threshold):
		Base(__index),RelativeThreshold(threshold),AbsoluteThreshold(0),Down(false),Up(0){}

	Iter(Iter const &_origin):
		Base(_origin) {}

	Iter(Iter const &_origin, typename Base::TPair const &_childRange):
		Base(_origin, _childRange) {}

};



// go down the leftmost edge
template < typename TIndex, class TSpec >
inline bool goDown(Iter< TIndex, VSTree< TopDown<ParentLinks<ConstrainedTraversal<TSpec > > > > > &it) {
	if (isLeaf(it)) return false;
	_historyPush(it, value(it));
		typename Size<TIndex>::Type i = _getUp(value(it).i2, container(it));
	if (value(it).i1 < i && i < value(it).i2)
		value(it).i2 = i;
	else
		value(it).i2 = _getDown(value(it).i1, container(it));

	// set the the flag for going down
	it.Down = true;
	return true;
}

// go up one edge (returns false if in root node) and changes it.Down
template < typename TIndex, class TSpec >
inline bool goUp(Iter< TIndex, VSTree< TopDown< ParentLinks<ConstrainedTraversal<TSpec> > > > > &it) {
	if (!isRoot(it)) {
		value(it) = top(it.history);
		pop(it.history);

		// whenever we go up it is necessary to set it.Down to false
		// and if we did not go down before we increase it.Up by one that we know
		// how far we went up(this can then be used to go up in the new automata)
		if(it.Down)
				it.Down = false;
		else
				++it.Up;
				
		return true;
	}
	return false;
}


template < typename TText, typename TSpec >
inline void goNext(Iter< Index<TText, TSpec>, VSTree< TopDown< ParentLinks<ConstrainedTraversal<Absolute> > > > > &it) {
	unsigned walk_down=0,not_finished = 1;
	//cout << "do the iterator with AbsoluteeThreshol:"<<it.AbsoluteThreshold<<endl;
	// resembles the construction of the VLMM-core of the context tree
	// only takes nodes that occur at least AbsoluteThreshold times
	// note: that we do not enforce that the tree is balanced

	if((countOccurences(it)) >= it.AbsoluteThreshold) 
		walk_down = 1;
	

	it.Down = false;
	it.Up = 0;
	while(not_finished)
	{
		// if threshold is reached walk down else go always right or up
		if(walk_down){
			if(!goDown(it) && !goRight(it))
				while (goUp(it) && !goRight(it));
		}
		else
			if (!goRight(it))
				while (goUp(it) && !goRight(it));
		
		if( countOccurences(it) >= it.AbsoluteThreshold)
			not_finished = 0;
		else
			walk_down = 0;

		// stop if we reach the root
		if (isRoot(it)) {
			clear(it);
			break;
		}

	}

}

template < typename TText, typename TSpec >
inline void goNext(Iter< Index<TText, TSpec>, VSTree< TopDown< ParentLinks<ConstrainedTraversal<Relative> > > > > &it) {
	unsigned walk_down=0,not_finished = 1;
	
	// resembles the construction of the VLMM-core of the PST algorithm of Ron et.al
	// only take nodes where their counts or their parents occur more 
	// often than RelativeThreshold times by considering the actual context length
	// and the text length. This is known as relative frequency or empirical probability.
	// note: that we do not enforce that the tree is balanced
	
	if(countOccurences(it) >= floor(it.RelativeThreshold*(float)(length(container(it))-repLength(it)+1))) 
		walk_down = 1;

	it.Down = false;
	it.Up = 0;
	while(not_finished)
	{
		// if threshold is reached walk down else go always right or up
		if(walk_down){
			// this is temporary for avoiding the implicit $-edges
			if(!isLeaf(it) && isRightTerminal(it)){
				// this must be possible
				goDown(it);
				goRight(it);
			}
			else
				if (!goDown(it)&&!goRight(it))
					while (goUp(it) && !goRight(it));
		}
		else
			if (!goRight(it))
				while (goUp(it) && !goRight(it));
		
		/*while(length(parentEdgeLabel(it)) == 0){
			walk_down = 0;
						goRight(it);	
		}*/
	/*	if(length(parentEdgeLabel(it)) == 0){ 
			walk_down = 0;
			continue;
		}*/
		if(!isLeaf(it) && (countOccurences(it) >= floor(it.RelativeThreshold*(float)(length(container(it))-repLength(it)+1))))
			not_finished = 0;
		else
			walk_down = 0;

		// stop if we reach the root
		if (isRoot(it)) {
			clear(it);
			break;
		}

	}


}



// Definition of the Graph Specialization for the VLMM

struct ContextTree{
	float threshold ;	
	
};
	struct PST{
		float threshold ;				   // r
		float minEmpiricalProbability ;    // P min
		float minConditionalProbability ;  // y min
		float alpha ;					   // alpha :-)
		
	};

void setParameters(PST & params,float t,float minE,float minC,float alpha){
	params.threshold = t;
	params.minEmpiricalProbability = minE;
	params.minConditionalProbability = minC;
	params.alpha = alpha;
}
void setParameters(ContextTree & params,float t){
	params.threshold = t;
}

template<typename TSpec = ContextTree>
struct VLMM;


////////////////
// VLMM
////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec>
class Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec> 

{
	public:
		typedef typename Id<Graph>::Type TIdType;
		typedef typename EdgeType<Graph>::Type TEdge;
		typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec > > >, TGraphSpec> TGraph;
		typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;


		String<AutomatonEdgeArray<TEdge, TAlphabet> > data_vertex;		// List of tables
		IdManager<TIdType> data_id_managerV;
		TVertexDescriptor root;
		String<String<TAlphabet> > data_edge_label;
		String<TVertexDescriptor>  data_suffix_link;
		// oder besser automaton edgearray
		String<AutomatonEdgeArray<TEdge, TAlphabet> > data_reverse_suffix_link;
		//String<TVertexDescriptor>  data_auxiliary_link;
		String<TVertexDescriptor>  data_father;
		String<bool>	   data_marked;
		String<String<float> >	   data_probability_vector;

//____________________________________________________________________________


		Graph() {
			SEQAN_CHECKPOINT
		}


		~Graph() {
			SEQAN_CHECKPOINT
			clear(*this);
		}

		Graph(Graph const & _other) 
		{
			SEQAN_CHECKPOINT
			_copyGraph(_other, *this);
		}
	
		Graph const& operator = (Graph const & _other) {
			SEQAN_CHECKPOINT
			if (this == &_other) return *this;
			_copyGraph(_other, *this);
			return *this;
		}
};

// new addVertex Specialization which enlarges data_edge_label,data_count,data_father at once
template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM< TSpec> > >, TGraphSpec> >::Type 
addVertex(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec>& g) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec> TGraph;
	typedef typename Size<TGraph>::Type TSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	unsigned num = numVertices(g)+1;
	resize(g.data_edge_label, (num) * table_length,Generous());
	appendValue(g.data_father, nilVal);
	appendValue(g.data_marked, false);
	resize(g.data_probability_vector,num,Generous());
	//append Alphabetsize many 0s
	for (int i = 0;i<table_length;++i)	
		append(g.data_probability_vector[num-1],0);
		

	return _addAutomatonVertex(g);
}

// new addVertex Specialization which enlarges data_edge_label,data_count,data_father at once
template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM< TSpec> > >, TGraphSpec> >::Type 
addAdditionalVertex(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec>& g) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec> TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph >::Type TEdge;
	typedef typename Size<TAlphabet>::Type TAlphabetSize;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	TSize table_length = ValueSize<TAlphabet>::VALUE;

	TVertexDescriptor vd = _addAutomatonVertex(g);
	// if new vertex vd is at the end,i.e. points to the last cell
	// a new entry has been created/appended in addAutomatonVertex

if (vd == length(g.data_vertex)-1) {
	String<String<TAlphabet> > edgeLabel;
	resize(g.data_edge_label, (length(g.data_edge_label)+1) * table_length, Generous());
	//appendValue(g.data_edge_label, edgeLabel);
	//resize(g.data_father,length(g.data_vertex));
	//g.data_father[vd] = nilVal;
	appendValue(g.data_father, nilVal);
	//resize(g.data_count,length(g.data_vertex));
	//g.data_count[vd] = 0;
	appendValue(g.data_marked,false);
	appendValue(g.data_suffix_link,nilVal);
	appendValue(g.data_reverse_suffix_link, AutomatonEdgeArray<TEdge, TAlphabet>());
	unsigned int Length = length(g.data_vertex);
	resize(g.data_probability_vector,Length,Generous());
	//append Alphabetsize many 0s
	for (int i = 0;i<table_length;++i)	
		append(g.data_probability_vector[Length-1],0);
}

// else they filled a gap in between -> no new space to allocate

return vd;
}

// init tables suffix_link  and reverse_suffix_link
template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec>
inline void
initMaps(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec>& g) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec> TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	typedef typename EdgeType<TGraph >::Type TEdge;
	unsigned num = numVertices(g);
	resize(g.data_reverse_suffix_link,num);
	resize(g.data_suffix_link, num);
	//resize(g.data_auxiliary_link, num);

	// init table data_suffix_link with NIL value
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	for(unsigned int i = 0;i<num;++i){
			g.data_suffix_link[i] = nilVal;
			value(g.data_reverse_suffix_link, i) =  AutomatonEdgeArray<TEdge, TAlphabet>();
	}
	

	return;
}

template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec,typename TVertexDescriptor>
inline  void
setFather(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec>& wg,
		  TVertexDescriptor& father,
		  TVertexDescriptor& child) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec> TGraph;

	wg.data_father[child] = father;
	return;
}

template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec,typename TVertexDescriptor>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec> >::Type 
getFather(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec>& wg,
		  TVertexDescriptor & child) 
{
	return wg.data_father[child] ;
}


template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec,typename TVertexDescriptor>
inline  void
setMarked(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec>& wg,
		  TVertexDescriptor& vertex,
		  bool state) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec> TGraph;

	wg.data_marked[vertex] = state;
	return;
}

template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec,typename TVertexDescriptor>
inline bool
isMarked(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec>& wg,
		  TVertexDescriptor & vertex) 
{

	return wg.data_marked[vertex] ;
}

template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec,typename TVertexDescriptor>
inline bool
isMarked(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec>const& wg,
		  TVertexDescriptor & vertex) 
{

	return wg.data_marked[vertex] ;
}


template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec,typename TVertexDescriptor>
inline  void
setSuffixLink(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec>& wg,
		  TVertexDescriptor& source,
		  TVertexDescriptor& target) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec> TGraph;

	wg.data_suffix_link[source] = target;
	return;
}

template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec,typename TVertexDescriptor>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec> >::Type 
getSuffixLink(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec>& wg,
		  TVertexDescriptor & source) 
{

	return wg.data_suffix_link[source] ;
}


template<typename TAlphabet, typename TCargo,typename TSpec,typename TGraphSpec, typename TVertexDescriptor, typename TChar>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec> >::Type 
getReverseSuffixLink(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec>& g,
			 TVertexDescriptor & vertex,
			 TChar const c) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec> TGraph;
	typedef typename Size<TAlphabet>::Type TSize;
	TAlphabet letter(c);
	return g.data_reverse_suffix_link[vertex].data_edge[(TSize) letter].data_target;
}


template<typename TAlphabet, typename TCargo, typename TSpec,typename TGraphSpec,  typename TVertexDescriptor, typename TChar>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec> >::Type 
setReverseSuffixLink(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec>& g,
			 TVertexDescriptor source,
			 TVertexDescriptor target,
			 TChar const c) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, source) == true)
	//typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec> TGraph;
	typedef typename Size<TAlphabet>::Type TSize;
	TAlphabet letter(c);
	return g.data_reverse_suffix_link[source].data_edge[(TSize) letter].data_target = target;
}


//Graph<Automaton<TText, TCargo, WordGraph< VLMM < TSpec > > >, TGraphSpec> 
// currently this function can only be called for index of
// specialization: Index_ESA
template<typename TIndex,typename TIterSpec,typename TSpec,typename TCargo,typename TAlphabet,typename TGraphSpec,typename TVertexDescriptor>
inline void
initProbabilityVector(Iter<TIndex, VSTree< TopDown< ParentLinks<ConstrainedTraversal<TIterSpec> > > > > &it,
						 Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > >,TGraphSpec> &target,
						 TVertexDescriptor & node)
{
	//typedef Graph<Automaton<TAlphabet,TCargo,WordGraph< VLMM < TSpec > > >,TGraphSpec> TVlmm;
	
	typedef Iter<TIndex, VSTree< TopDown< > > >  TIter;

	TIter childs(it);
	goDown(childs);
	TAlphabet startCharacter;
	unsigned fatherLength = repLength(it);
	//cout <<"Cilds:" << representative(childs) << endl;
	if( repLength(childs) > fatherLength){

		startCharacter = value(representative(childs),fatherLength);
		setProbability(target,node,startCharacter,(float)countOccurences(childs));
	}
		while(goRight(childs)){

			// cout << representative(childs) << endl;
			if(repLength(childs) > fatherLength){
			startCharacter = value(representative(childs),fatherLength);
			setProbability(target,node,startCharacter,(float)countOccurences(childs));
			}
		}

}


//Graph<Automaton<TText, TCargo, WordGraph< VLMM < TSpec > > >, TGraphSpec> 
// currently this function can only be called for index of
// specialization: Index_ESA
template<typename TIndex,typename TIterSpec,typename TSpec,typename TCargo,typename TAlphabet,typename TGraphSpec>
inline void 
buildSuffixTreeFromIndex(Iter<TIndex, VSTree< TopDown< ParentLinks<ConstrainedTraversal<TIterSpec> > > > > &it,
						 Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > >,TGraphSpec> &target)
{
	typedef Graph<Automaton<TAlphabet,TCargo,WordGraph< VLMM < TSpec > > >,TGraphSpec> TVlmm;
	typedef typename VertexDescriptor<TVlmm>::Type TVertexDescriptor;
	typedef Iter<TIndex, VSTree< TopDown< > > >  TIter;
	typedef Iter<TIndex, VSTree< TopDown< ParentLinks<ConstrainedTraversal<TIterSpec> > > > >  TIter2;
	
	TVertexDescriptor root = addVertex(target);
	TVertexDescriptor child,father;
	child = root;

	initProbabilityVector(it,target,root);
	assignRoot(target,root);
	cout << "first chars are done.."<<endl;
	//go to the first node
	goNext(it);
	while(!atEnd(it)){
		
		
		// it.Down == true: the previous child becomes the new father
		// neither it.Down == true nor it.Up > 0: the father remains the same
		// it.Up == x: the new father is x levels up in the target

		if(it.Down){
			father = child;
		}
		while(it.Up>0){
			father = getFather(target,father);
			--it.Up;
		}

		child = addVertex(target);
		cout << "Node:"<<child<<"  "<<value(it) << " = " << repLength(it)<< " " << representative(it) << "  toFather:"<<parentEdgeLabel(it)<<"  hits: "<<length(getOccurences(it))<<endl;
		setFather(target,father,child);
		cout <<"set Father\t";
		// set child relation with edgelabel
		//cout <<"ParentEdgeLAbel:"<<parentEdgeLabel(it)<<endl;
		addEdge(target,father,child,parentEdgeLabel(it));
		cout <<"set Edge\t";
		initProbabilityVector(it,target,child);
		//// provisional set prob-vector of nodes
		//TIter childs(it);
		//goDown(childs);
		//unsigned fatherLength = 0;
		//TAlphabet startCharacter;
		//if(isRightTerminal(childs)){
		//fatherLength = repLength(it);
		//startCharacter = value(representative(childs),fatherLength);
		//setProbability(target,child,startCharacter,(float)countOccurences(childs));
		//}
		//while(goRight(childs)){
		//	if(isRightTerminal(childs)){
		//	startCharacter = value(representative(childs),fatherLength);
		//	setProbability(target,child,startCharacter,(float)countOccurences(childs));
		//	}
		//}
		cout << "added node:"<<child<<endl;
		// go to next valid node and add it to the graph
		goNext(it);

	}

	// All valid nodes are now part of the vlmm
	return;
}

template<typename TSpec,typename TCargo,typename TAlphabet,typename TGraphSpec,typename TVertexDescriptor>
inline TAlphabet 
getChildCharacter(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > >,TGraphSpec> &vlmm,
			   TVertexDescriptor &father,
			   TVertexDescriptor &child)
{
	typedef Graph<Automaton<TAlphabet,TCargo,WordGraph< VLMM < TSpec > > >,TGraphSpec> TVlmm;
	typedef typename Iterator<TVlmm, OutEdgeIterator<> >::Type TOutEdgeIterator;

	TOutEdgeIterator itout(vlmm,father);
	while(!atEnd(itout)){
		if(targetVertex(vlmm, getValue(itout)) == child)
			break;

	goNext(itout);
	}
	return(itout.data_pos);
// look in the EdgeAutomaton of father where the son/daughter is

}


template<typename TSpec,typename TCargo,typename TAlphabet,typename TGraphSpec,typename TVertexDescriptor>
inline void 
getChildLabel(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > >,TGraphSpec> &vlmm,
			   TVertexDescriptor &father,
			   TVertexDescriptor &child,
			   String<TAlphabet> &Label)
{
	typedef Graph<Automaton<TAlphabet,TCargo,WordGraph< VLMM < TSpec > > >,TGraphSpec> TVlmm;
	typedef typename EdgeDescriptor<TVlmm>::Type TEdgeDescriptor;

	TAlphabet childPosition = getChildCharacter(vlmm,father,child);

	append(Label,childPosition);
	append(Label,getProperty(vlmm.data_edge_label, TEdgeDescriptor(father, childPosition)));

	return;
}

template<typename TSpec,typename TCargo,typename TAlphabet,typename TChar,typename TGraphSpec,typename TVertexDescriptor>
inline float
getProbability(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > >,TGraphSpec> &vlmm,
			   TVertexDescriptor &father,
			   TChar pos)
{
	TAlphabet letter(pos);
	return value(vlmm.data_probability_vector[father],(int)letter);
}

template<typename TSpec,typename TCargo,typename TAlphabet,typename TChar,typename TGraphSpec,typename TVertexDescriptor>
inline void
setProbability(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > >,TGraphSpec> &vlmm,
			   TVertexDescriptor &father,
			   TChar pos,
			   float prob)
{
	TAlphabet letter(pos);
	value(vlmm.data_probability_vector[father],(int)letter) = prob;
}


// splitPosition is the position in the array string edgelabel
// e.g. edgelabel = ACGT , the new edge should be GT than splitPosition is 2
template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec, typename TVertexDescriptor, typename TChar,typename TPos>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec> >::Type 
splitEdge(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec> & vlmm,
			TVertexDescriptor & father,
			TChar  childCharacter,
			TPos  splitPosition)
{
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM <TSpec> > >, TGraphSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TGraph>::Type TSize;
	TAlphabet letter(childCharacter);
	TVertexDescriptor child = vlmm.data_vertex[father].data_edge[(TSize) letter].data_target;
	String<TAlphabet> edgeString = getProperty(vlmm.data_edge_label, TEdgeDescriptor(father, letter));
	TVertexDescriptor newNode = addAdditionalVertex(vlmm);
	String<TAlphabet> newEdgeString = childCharacter;
	append(newEdgeString,prefix(edgeString,splitPosition),Exact());
	addEdge(vlmm,father,newNode,newEdgeString );

	setFather(vlmm,father,newNode);
	setProbability(vlmm,newNode,value(edgeString,splitPosition),1);
	addEdge(vlmm,newNode,child,suffix(edgeString,splitPosition));
	setFather(vlmm,newNode,child);

	return newNode;
}

template<typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec, typename TVertexDescriptor, typename TLabel>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec> >::Type 
parseString2(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec> & g,
			TVertexDescriptor & vertex,
			TLabel & label)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM <TSpec> > >, TGraphSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TGraph>::Type TSize;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	TVertexDescriptor succ = vertex;
	TSize pos = 0;
	//cout <<"start parseString at node:"<< succ<<endl;
	TAlphabet letter = value(label,pos);
	TVertexDescriptor tmp = g.data_vertex[succ].data_edge[(TSize) letter].data_target;
	++pos;
	while (tmp != nilVal) {
		String<TAlphabet> edgeString = getProperty(g.data_edge_label, TEdgeDescriptor(succ, letter));
		if( (pos == length(label)) && (length(edgeString) == 0))
				return tmp;

		if(length(edgeString) > 0){
			//cout <<"at  node:" << succ<<endl;
			//cout << "infix: "<<infix(label,pos,pos+length(edgeString))<<" EdgeString:"<<edgeString<<endl;
			if ( (pos+length(edgeString)-1 < length(label) ) && edgeString == infix(label,pos,pos+length(edgeString))){
				pos += length(edgeString);
				if(pos == length(label))
						return tmp;
			}
			else
				return nilVal;
		}

		succ = tmp;
		letter = value(label,pos);
		tmp = g.data_vertex[succ].data_edge[(TSize) letter].data_target;
		++pos;
	}
	return nilVal;
}



template<typename TSpec,typename TCargo,typename TAlphabet,typename TGraphSpec>
inline void 
addSuffixLinks(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > >,TGraphSpec> &vlmm)
{
	typedef Graph<Automaton<TAlphabet,TCargo,WordGraph< VLMM < TSpec > > >,TGraphSpec> TVlmm;
	typedef typename VertexDescriptor<TVlmm>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TVlmm>::Type TEdgeDescriptor;
	typedef typename Iterator<TVlmm, OutEdgeIterator<> >::Type TOutEdgeIterator;
	typedef typename Size<TVlmm>::Type TSize;
	
	TVertexDescriptor root = getRoot(vlmm);
	TVertexDescriptor father,target;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	TOutEdgeIterator itout(vlmm,root);
	setSuffixLink(vlmm,root,root);
	TAlphabet letter;
	for(;!atEnd(itout);goNext(itout)) {
		TVertexDescriptor v = targetVertex(vlmm, getValue(itout));
		letter = itout.data_pos;
		// means the suffix link must point to the root
		if(length(getProperty(vlmm.data_edge_label, TEdgeDescriptor(root, letter))) == 0){
				setSuffixLink(vlmm,v,root);
				setReverseSuffixLink(vlmm,root,v,letter);
				
			}
		else{
				// getProperty(vlmm.data_edge_label, TEdgeDescriptor(root, letter) == edgelabel without first char
				target = parseString2(vlmm,root,getProperty(vlmm.data_edge_label, TEdgeDescriptor(root, letter)));
				SEQAN_ASSERT(target!=nilVal)
				setSuffixLink(vlmm,v,target);
				setReverseSuffixLink(vlmm,target,v,letter);
			}
	}
// Having set these trivial cases we can treat every node the same way


deque<TVertexDescriptor> queue;
TOutEdgeIterator children(vlmm,root);
while(!atEnd(children)){
	TAlphabet childCharacter = children.data_pos;
	queue.push_back(targetVertex(vlmm, getValue(children)));
	
	// Bfs
	while (!queue.empty()) {
		TVertexDescriptor n = queue.front();
		queue.pop_front();
		
		if(getSuffixLink(vlmm,n) == nilVal){

			father = getFather(vlmm,n);
			String<TAlphabet> walkDown;
			getChildLabel(vlmm,father,n,walkDown);
			//cout << "go for SL of father:" << father<< " String:"<< walkDown<<endl;
			TVertexDescriptor fatherSuffixLink = getSuffixLink(vlmm,father);
			target = parseString2(vlmm,fatherSuffixLink,walkDown);
			SEQAN_ASSERT(target!=nilVal)
			setSuffixLink(vlmm,n,target);
			setReverseSuffixLink(vlmm,target,n,childCharacter);
		}

		TOutEdgeIterator itout(vlmm,n);
		for(;!atEnd(itout);goNext(itout)) {
			TVertexDescriptor v = targetVertex(vlmm, getValue(itout));
			queue.push_back(v);
			
		}
	}
	goNext(children);
}
}


/**********************
* Pruning
**********************/

// spezialization of DifferenceOperator are PST and ContextTree like for VLMM 
// really needed ?????
template <typename TSpec = Default>
struct DifferenceOperator
;

template<typename TCargo,typename TAlphabet,typename TSpec,typename TGraphSpec,typename TVertexDescriptor>
inline void
turnNodeCountsIntoProbability(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > >,TGraphSpec> &vlmm,
			   TVertexDescriptor & node) 
{
	unsigned size = ValueSize<TAlphabet>::VALUE;
	float sum = 0;;
	for(int pos = 0;pos< size ;++pos)
		sum = sum + getProbability(vlmm,node,pos);
	SEQAN_ASSERT(sum != 0);
	for(int pos = 0;pos< size;++pos)
		setProbability(vlmm,node,pos, (getProbability(vlmm,node,pos)/sum) );
}



template<typename TCargo,typename TAlphabet,typename TGraphSpec,typename TChar,typename TVertexDescriptor>
inline bool
extendNode(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < PST > > >,TGraphSpec> &vlmm,
			   TVertexDescriptor & node,
			   TChar childCharacter,
			   PST & parameters) 
{
	
		TAlphabet letter(childCharacter);

		float gammaMin = (1 + parameters.alpha)*parameters.minConditionalProbability;
		
		if( 1 >= gammaMin){
			SEQAN_ASSERT(getProbability(vlmm,node,childCharacter) != 0);
			
			float ratio = 1/getProbability(vlmm,node,childCharacter);

			if(ratio >= parameters.threshold  || ratio <= 1/parameters.threshold )
				return true;
		}
	return false;

}

template<typename TCargo,typename TAlphabet,typename TGraphSpec,typename TVertexDescriptor>
inline bool
pruneNode(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < PST > > >,TGraphSpec> &vlmm,
			   TVertexDescriptor &son,
			   PST & parameters) 
{
	// the original pruning criteria defined by Ron et.al. 1996
	float gammaMin = (1 + parameters.alpha)*parameters.minConditionalProbability;
	float ratio = 0;
	TVertexDescriptor father = getSuffixLink(vlmm,son);
	for(int pos = 0;pos< ValueSize<TAlphabet>::VALUE;++pos){
		
		if(getProbability(vlmm,son,pos) >= gammaMin){
			if(getProbability(vlmm,father,pos) == 0)
				ratio = 0;
			else
				ratio = getProbability(vlmm,son,pos)/getProbability(vlmm,father,pos);

			if(ratio >= parameters.threshold  || ratio <= 1/parameters.threshold )
				return false;
		}
	}
 return true;
}

template<typename TCargo,typename TAlphabet,typename TGraphSpec,typename TVertexDescriptor>
inline void
smoothNode(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < PST > > >,TGraphSpec> &vlmm,
			   TVertexDescriptor & node,
			   PST & parameters) 
{
	typedef typename Size<TAlphabet>::Type TSize;
	TSize alphaSize = ValueSize<TAlphabet>::VALUE;
	// the original pruning criteria defined by Ron et.al. 1996
	// (1-|Sigma|*ymin)
	float gammaFactor = (1- ( alphaSize * parameters.minConditionalProbability) );

	for(int pos = 0;pos< ValueSize<TAlphabet>::VALUE;++pos)
		setProbability(vlmm,node,pos, ( getProbability(vlmm,node,pos)*(gammaFactor)  + parameters.minConditionalProbability ) );
}

template<typename TAlphabet,typename TCargo,typename TSpec,typename TGraphSpec,typename TVLMMSPec>
inline void
pruneTree(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > >,TGraphSpec> &vlmm,
			   TVLMMSPec & parameters) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator<> >::Type TVertexIterator;

	TVertexDescriptor root = getRoot(vlmm);
	
	String<bool> original;
	// remember which nodes have been there
	initVertexMap(vlmm, original);
	TVertexIterator it(vlmm);
	for(;!atEnd(it);goNext(it)) {
		TVertexDescriptor node =  getValue(it);
		assignProperty(original, node, true);
		// change counts of all nodes into probabilities
		turnNodeCountsIntoProbability(vlmm,node);
	}
	
	// root cannot be pruned
	setMarked(vlmm,root,true);
	//assignProperty(marked,root,true);

	pruneTreeRecursively(vlmm,root,original,parameters);

		for(unsigned i =0;i<numVertices(vlmm);++i){
				cout << isMarked(vlmm,i)<< "  ";
		}
		cout <<endl;
		for(unsigned i =0;i<length(original);++i){
				cout << original[i]<< "  ";
		}
		cout <<endl;

cout << "ready with the function\n";

 return;
}

// recursive walk over reverse suffix links
template<typename TAlphabet,typename TCargo,typename TSpec,typename TGraphSpec,typename TVertexDescriptor,typename TMarked,typename TVLMMSpec>
inline void
pruneTreeRecursively(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > >,TGraphSpec> &vlmm,
				TVertexDescriptor & node,
				TMarked & original,
			    TVLMMSpec & parameters) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TAlphabet>::Type TSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	

	TVertexDescriptor next;
	// check all outgoing reverse suffix links
	for(int pos = 0;pos< table_length;++pos){
		next = getReverseSuffixLink(vlmm,node,pos);
		if(next != nilVal){
			cout <<" start Recursion from node:"<<node<<" char:"<<pos<<endl;
			pruneTreeRecursively(vlmm,next,original,parameters);
		}
		
	}
	// the root is excluded from further processing
	if(isRoot(vlmm,node)){ 
	
		return ;
			
	}

	TVertexDescriptor father = getFather(vlmm,node);
	String<TAlphabet> childLabel;
	getChildLabel(vlmm,father,node,childLabel);
	// check for possible extended 
	if(length(childLabel) > 1){
		// seek the node in the tree which could be the suffix link father of the new node
		
		if(! isRoot(vlmm,father) ){
				TVertexDescriptor target = getSuffixLink(vlmm,father);

			cout <<"..check potential nodes above node:" <<node<<endl;
			TVertexDescriptor potVertex;
			unsigned lastVertex = 0;
			for(unsigned pos = 1;pos < length(childLabel);++pos ){
				cout << "check string:" << prefix(childLabel,pos)<<" from node:"<< target<<endl;
				// potVertex = potentialVertex is the potential SuffixLinkFather of the node to be created
				potVertex = parseString2(vlmm,target,prefix(childLabel,pos));
				if( (potVertex != nilVal) && (potVertex < length(original)) && original[potVertex] && extendNode(vlmm,potVertex,value(childLabel,pos),parameters) ){
					// find character for Reverse Suffix Link
					TVertexDescriptor faader = getFather(vlmm,father);
					TVertexDescriptor vorfaader = father;
					while(! isRoot(vlmm,faader) ){
						vorfaader = faader;
						faader = getFather(vlmm,faader);
						
					}
					TAlphabet letter = getChildCharacter(vlmm,faader,vorfaader);

					father = splitEdge(vlmm,father,value(childLabel,lastVertex),pos-1-lastVertex);
					smoothNode(vlmm,father,parameters);
					setSuffixLink(vlmm,father,potVertex);
					setReverseSuffixLink(vlmm,potVertex,father,letter);
					setMarked(vlmm,father,true);
					lastVertex = pos;
					
				}
						
			}
		} // if father == root check only for the first possible extension
		  // extension longer than one are not checked, because it is already known that there is no 
		  // node one context shorter elsewhere in the tree
		else{
				TVertexDescriptor root = getRoot(vlmm);
				TAlphabet startChar = value(childLabel,0);
				if(extendNode(vlmm,root,startChar,parameters))
				{
					cout <<"startChar" << startChar<<endl;
					father = splitEdge(vlmm,father,startChar,0);
					smoothNode(vlmm,father,parameters);
					setSuffixLink(vlmm,father,root);
					setReverseSuffixLink(vlmm,root,father,startChar);
					setMarked(vlmm,father,true);
				}
			
		} // else case
	} // Label > 1
	cout << "check keeping of node:" << node;
	// all potential nodes are build
	if( isMarked(vlmm,node) || (! pruneNode(vlmm,node,parameters)) ){
		// node should be kept
		cout <<" keep it"<<endl;
		setMarked(vlmm,node,true);
		// bubble up and set all nodes on the way to the root true
		TVertexDescriptor fatherSuffixLink = getSuffixLink(vlmm,node);
		while( !isRoot(vlmm,fatherSuffixLink) ){
			if(! isMarked(vlmm, fatherSuffixLink) ){
				setMarked(vlmm,node,true);
			}
			else
				break;
			fatherSuffixLink = getSuffixLink(vlmm,fatherSuffixLink);
		}
		// smooth node
		smoothNode(vlmm,node,parameters);
	}
	else
	{
	// delete node
		cout <<"  delete it"<<endl;
		setMarked(vlmm,node,false);
	//assignProperty(marked,node,false);
	}

	//cout <<"finished node:"<<node<<endl;
 return;
}

				// recursive walk over reverse suffix links
template<typename TIndexType,typename TAlphabet,typename TCargo,typename TGraphSpec>
inline void
buildPST(Index<TIndexType, Index_ESA<> > & index,
		 Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < PST > > >,TGraphSpec> &vlmm,
			    float threshold,
				float minEmpiricalProbability,
				float minConditionalProbability,
				float alpha) 
{
	typedef Index<TIndexType, Index_ESA<> > TIndex;
	PST parameters;
	setParameters(parameters,threshold,minEmpiricalProbability,minConditionalProbability,alpha);

	Iter< TIndex, VSTree< TopDown< ParentLinks<ConstrainedTraversal<Relative> > > > > it(index,parameters.minEmpiricalProbability);
	//Iter< TIndex, VSTree< TopDown< ParentLinks<ConstrainedTraversal<Absolute> > > > > it(index,2);

	buildSuffixTreeFromIndex(it,vlmm);
	cout << "constructed suffix tree core" <<endl;
	initMaps(vlmm);
	cout << "init all maps" <<endl;
	addSuffixLinks(vlmm);
	cout << "added suffix links and reverse suffix links" <<endl;
	cout <<vlmm;
	pruneTree(vlmm,parameters);
	cout << "pruned the PST" <<endl;
	cout << vlmm;
	cout << "READY!" <<endl;
}

//write the vlmm to a file in Easy-readable format,also used by the  << Operator
template<typename TFile, typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec, typename TIDString>
inline void
write(TFile & target,
	  Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec>const & g,
	  TIDString const &,
	  Raw)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Size<TAlphabet>::Type TSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();

	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const>::Type TIterConst;
	_streamWrite(target,"VLMM - Probability Vector of each node:\n");
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		if (!idInUse(g.data_id_managerV, position(it))) continue;
		TVertexDescriptor sourceVertex = position(it);
			_streamPutInt(target, sourceVertex);
			_streamPut(target, ' ');
			_streamPut(target, ' ');
			_streamWrite(target, "SLink: ");
			cout << g.data_suffix_link[sourceVertex];
			_streamWrite(target,"  (");

			// wollte Wahrscheinlichleiten ausgeben
		for(int i = 0;i<ValueSize<TAlphabet>::VALUE;++i){
			cout <<value(g.data_probability_vector[sourceVertex],i);
				//_streamPut(target,(float)value(g.data_probability_vector[sourceVertex],i));
				if(i < (ValueSize<TAlphabet>::VALUE -1))
						_streamPut(target, ' ');	
		}
		_streamPut(target, ')');
		if(g.data_marked[getRoot(g)])
			if(! g.data_marked[sourceVertex])
				_streamWrite(target, "  deleted");
		_streamWrite(target, "\n");
	}
		_streamWrite(target,"VLMM - EdgeList:\n");

	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		if (!idInUse(g.data_id_managerV, position(it))) continue;
		TVertexDescriptor sourceVertex = position(it);
		for(TSize i=0;i<table_length;++i) {
			if (g.data_vertex[sourceVertex].data_edge[i].data_target ==  nilVal) continue;
			_streamPutInt(target, sourceVertex);
			_streamWrite(target,"->");
			_streamPutInt(target, g.data_vertex[sourceVertex].data_edge[i].data_target);
			_streamPut(target, ' ');
			_streamPut(target, ' ');
			_streamWrite(target, "Label: ");
			_streamPut(target, TAlphabet(i));
			_streamWrite(target, getProperty(g.data_edge_label, TEdgeDescriptor(sourceVertex, TAlphabet(i))));
		
			_streamPut(target, '\n');
		}
	}
}

//write the vlmm to a file in Dot-Format
template <typename TFile,typename TAlphabet,typename TCargo,typename TGraphSpec, typename TSpec, typename TNodeMap, typename TEdgeMap>
void write(TFile & file, 
	   Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > >,TGraphSpec> const& g,
	   TNodeMap const& nodeMap,
	   TEdgeMap const& edgeMap,
	   DotDrawing) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > >,TGraphSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;

	_streamWrite(file, "digraph G {\n");
	_streamPut(file, '\n');
	_streamWrite(file, "/* Graph Attributes */\n");
	_streamWrite(file, "graph [rankdir = LR];\n");
	_streamPut(file, '\n');
	_streamWrite(file, "/* Node Attributes */\n");
	_streamWrite(file, "node [shape = ellipse, fillcolor = lightgrey, style = filled, fontname = \"Times-Italic\"];\n");
	_streamPut(file, '\n');
	_streamWrite(file, "/* Edge Attributes */\n");
	_streamWrite(file, "edge [fontname = \"Times-Italic\", arrowsize = 0.75, fontsize = 16];\n");
	_streamPut(file, '\n');

	_streamWrite(file, "/* Nodes */\n");
	typedef typename Iterator<TGraph, VertexIterator<> >::Type TConstIter;
	TConstIter it(g);
	for(;!atEnd(it);goNext(it)) {
		_streamWrite(file, getProperty(nodeMap, *it));
		/* wollte Wahrscheinlichleiten ausgeben
		_streamWrite(file, "|");
		for(int i = 0;i<ValueSize<TAlphabet>::VALUE;++i){
			//cout <<"write:"<<value(g.data_probability_vector[*it],i)<<endl;
				_streamWrite(file,(int)value(g.data_probability_vector[*it],i));
				if(i < (ValueSize<TAlphabet>::VALUE -1))
						_streamPut(file, ':');	
		}
		*/
		_streamPut(file, ';');
		_streamPut(file, '\n');
	}
	_streamPut(file, '\n');

	_streamWrite(file, "/* Edges */\n");
	typedef typename Iterator<TGraph, EdgeIterator<> >::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		TVertexDescriptor sc = sourceVertex(itEd);
		TVertexDescriptor tr = targetVertex(itEd);
		_streamWrite(file, getProperty(nodeMap, sc));
		_streamWrite(file, " -> ");
		_streamWrite(file, getProperty(nodeMap, tr));
		_streamWrite(file, " [label = \"");
		_streamWrite(file, getProperty(edgeMap, *itEd));
		_streamWrite(file, "\"];\n");
	}
	//NEW Reverse Suffix Links
	/*
	goBegin(it);
	for(;!atEnd(it);++it) {

		TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
		typedef typename Size<TAlphabet>::Type TSize;
		TSize table_length = ValueSize<TAlphabet>::VALUE;
		TSize pos = 0;
		for(pos = 0;pos<table_length;++pos){
			TVertexDescriptor sc = g.data_reverse_suffix_link[*it].data_edge[pos].data_target;
			if(sc != nilVal){
	
				TVertexDescriptor tr = *it;
				_streamWrite(file, getProperty(nodeMap, tr));
				_streamWrite(file, " -> ");
				_streamWrite(file, getProperty(nodeMap, sc));
				_streamWrite(file, " [label = \"");
				_streamWrite(file, (TAlphabet)pos);
				_streamWrite(file, "\"];\n");
			}
		}
	}
	*/
	_streamPut(file, '\n');

	_streamWrite(file, "}\n");
}

template <typename TFile,typename TAlphabet,typename TCargo,typename TGraphSpec, typename TSpec>
void write(TFile & file, 
	   Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > >,TGraphSpec> const& g, 
	   DotDrawing) 
{
	SEQAN_CHECKPOINT
	String<String<char> > nodeMap;
	_createNodeNames(g,nodeMap);
	String<String<char> > edgeMap;
	_createEdgeNames(g,edgeMap);
	write(file,g,nodeMap,edgeMap,DotDrawing());
}

//write the vlmm not as a suffix tree but reveersed edge labels to a file in Dot-Format
template <typename TFile,typename TAlphabet,typename TCargo,typename TGraphSpec, typename TSpec, typename TNodeMap, typename TEdgeMap>
void writeAsVLMM(TFile & file, 
	   Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > >,TGraphSpec> const& g,
	   TNodeMap const& nodeMap,
	   TEdgeMap const& edgeMap,
	   DotDrawing) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > >,TGraphSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;

	_streamWrite(file, "digraph G {\n");
	_streamPut(file, '\n');
	_streamWrite(file, "/* Graph Attributes */\n");
	_streamWrite(file, "graph [rankdir = LR];\n");
	_streamPut(file, '\n');
	_streamWrite(file, "/* Node Attributes */\n");
	_streamWrite(file, "node [shape = ellipse, fillcolor = lightgrey, style = filled, fontname = \"Times-Italic\"];\n");
	_streamPut(file, '\n');
	_streamWrite(file, "/* Edge Attributes */\n");
	_streamWrite(file, "edge [fontname = \"Times-Italic\", arrowsize = 0.75, fontsize = 16];\n");
	_streamPut(file, '\n');

	_streamWrite(file, "/* Nodes */\n");
	typedef typename Iterator<TGraph, VertexIterator<> >::Type TConstIter;
	TConstIter it(g);
	for(;!atEnd(it);++it) {
		_streamWrite(file, getProperty(nodeMap, *it));
		_streamPut(file, ';');
		_streamPut(file, '\n');
	}
	_streamPut(file, '\n');

	_streamWrite(file, "/* Edges */\n");
	typedef typename Iterator<TGraph, EdgeIterator<> >::Type TConstEdIter;
	TConstEdIter itEd(g);
		//NEW Reverse Suffix Links
	goBegin(it);
	for(;!atEnd(it);++it) {

		TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
		typedef typename Size<TAlphabet>::Type TSize;
		TSize table_length = ValueSize<TAlphabet>::VALUE;
		TSize pos = 0;
		for(pos = 0;pos<table_length;++pos){
			TVertexDescriptor sc = g.data_reverse_suffix_link[*it].data_edge[pos].data_target;
			if(sc != nilVal){
	
				TVertexDescriptor tr = *it;
				_streamWrite(file, getProperty(nodeMap, tr));
				_streamWrite(file, " -> ");
				_streamWrite(file, getProperty(nodeMap,sc));
				_streamWrite(file, " [label = \"");
				_streamWrite(file, (TAlphabet)pos);
				_streamWrite(file, "\"];\n");
			}
		}
	}
	_streamPut(file, '\n');

	_streamWrite(file, "}\n");
}

template <typename TFile,typename TAlphabet,typename TCargo,typename TGraphSpec, typename TSpec>
void writeAsVLMM(TFile & file, 
	   Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > >,TGraphSpec> const& g, 
	   DotDrawing) 
{
	SEQAN_CHECKPOINT
	String<String<char> > nodeMap;
	_createNodeNames(g,nodeMap);
	String<String<char> > edgeMap;
	_createEdgeNames(g,edgeMap);
	writeAsVLMM(file,g,nodeMap,edgeMap,DotDrawing());
}

} // End Namespace

#endif

