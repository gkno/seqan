#ifndef SEQAN_HEADER_VLMM_H
#define SEQAN_HEADER_VLMM_H
#define SEQAN_PROFILE 
// for time stopping

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void inittheclock(void);
double gettheruntime(void);
int gettheclockticks(void); 

/*
  This module implements function to measures the running time
  of programs or program parts.
*/

/*
  The following values store the the clockticks at start time 
  and stop time of the clock.
*/

static clock_t startclock, 
               stopclock;

/*EE
  The following function initializes the clock.
*/

void inittheclock(void)
{ 
  startclock = clock(); 
}

/*EE
  The following function delivers the time since the 
  clock was initialized. The time is reported in seconds
  as a floating point value.
*/

double gettheruntime(void)
{
   stopclock = clock();
   return (double) (stopclock-startclock) / (double) CLOCKS_PER_SEC;
}

/*EE
  The following function delivers the clock ticks betwenn 
  \texttt{startclock} to \texttt{stopclock}.
*/

int gettheclockticks(void)
{
   stopclock = clock();
   return (int) (stopclock-startclock);
}


using namespace std;
namespace SEQAN_NAMESPACE_MAIN
{

//this file is the initial version of the construction of
// a variable length markov model (ConstrainedTraversal)


// Defining specialized iterators

	// specializations for the ConstrainedTraversal iterator
	struct Absolute;
	struct Relative;
	struct Support;

template <typename TSpec = Absolute>
struct ConstrainedTraversal;


  template < typename TSTree, typename T >
    struct GetVSTreeIteratorTraits< Iter< TSTree, VSTree< TopDown<
ParentLinks< ConstrainedTraversal<T> > > > > >{
        typedef Preorder Type;
    };



template < typename TIndex,typename TSpec>
class Iter< TIndex, VSTree< TopDown< ParentLinks<ConstrainedTraversal<TSpec> > > > >:
	public Iter< TIndex, VSTree< TopDown< ParentLinks<> > > >
{
	typedef Iter< TIndex, VSTree< TopDown< ParentLinks<> > > >TBase;
	typedef	Pair<typename Size<TIndex>::Type>	TStackEntry, TRange;
		typedef String<TStackEntry, Block<> >		TStack;
		typedef Iter								iterator;
public:
	TStack			history;	// contains all previously visited intervals (allows to go up)
	// although both variables are hold by the iterator
	// actually either one of them is used during traversal
	// indicated by the specialization Absolute or Relative
	unsigned AbsoluteThreshold;
	float    RelativeThreshold;
	unsigned MaxDepth;
	bool	 Down;
	unsigned Up;

	
	Iter(TIndex &__index):
		TBase(__index){}

	// use constructor for an absolute threshold
	Iter(TIndex &__index,unsigned threshold,unsigned maxdepth):
		TBase(__index),AbsoluteThreshold(threshold),RelativeThreshold(0.0),MaxDepth(maxdepth),Down(false),Up(0){}

	// use constructor for a relative threshold
	Iter(TIndex &__index,float threshold,unsigned maxdepth):
		TBase(__index),AbsoluteThreshold(0),RelativeThreshold(threshold),MaxDepth(maxdepth),Down(false),Up(0){}

		// brauch man die überhaupt
		Iter(Iter const &_origin):
			TBase((TBase const &)_origin),
			history(_origin.history) {}

 		Iter(Iter const &_origin, TRange const &_childRange):
			TBase((TBase const &)_origin, _childRange),
			history(_origin.history)
		{
			push(history, value(_origin.i1));
		}


};


//
//// go down the leftmost edge
//template < typename TIndex, class TSpec >
//inline bool goDown(Iter< TIndex, VSTree< TopDown<ParentLinks<ConstrainedTraversal<TSpec > > > > > &it) {
//	if (isLeaf(it) || (length(representative(it)) > it.MaxDepth) ) return false;
//	_historyPush(it, value(it).i1);
//
//		typename Size<TIndex>::Type i = _getUp(value(it).i1.i2, container(it));
//	if (value(it).i1.i1 < i && i < value(it).i1.i2)
//		value(it).i1.i2 = i;
//	else
//		value(it).i1.i2 = _getDown(value(it).i1.i1, container(it));
//
//	// set the the flag for going down
//	it.Down = true;
//	return true;
//}

// go down the leftmost edge (skip empty $-edges)
	template < typename TIndex, class TSpec >
	inline bool goDown(Iter< TIndex, VSTree< TopDown<ParentLinks<ConstrainedTraversal<TSpec > > > > > &it)
	{
		if (isLeaf(it) || (length(representative(it)) > it.MaxDepth) ) return false;
		//set to to true, if you go up again it will be changed in goUp(it)
		it.Down = true;
		_historyPush(it, value(it).i1);

		TIndex const &index = container(it);

		typename Size<TIndex>::Type lval = _getUp(value(it).i1.i2, index);
		if (!(value(it).i1.i1 < lval && lval < value(it).i1.i2))
			lval = _getDown(value(it).i1.i1, index);
		value(it).i1.i2 = lval;

		typename Size<TIndex>::Type lcp = lcpAt(lval - 1, index);
		//typename typename StringSetLimits<TIndex const>::Type &limits = stringSetLimits(index);
		while (isLeaf(it)) {
			typename SAValue<TIndex>::Type pos = getOccurence(it);
			if (getSeqOffset(pos, stringSetLimits(index)) + lcp == sequenceLength(getSeqNo(pos, stringSetLimits(index)), index)) {
				if (!goRight(it)) {
					goUp(it);
					return false;
				}
			} else
				break;
		}
		
		return true;
	}

// go up one edge (returns false if in root node) and changes it.Down
template < typename TIndex, class TSpec >
inline bool goUp(Iter< TIndex, VSTree< TopDown< ParentLinks<ConstrainedTraversal<TSpec> > > > > &it) {

	if (!empty(it.history)) {
			value(it).i1 = top(it.history);
			pop(it.history);
			// whenever we go up it is necessary to set it.Down to false
			// and if we did not go down before we increase it.Up by one that we know
			// how far we went up(this can then be used to go up in the new automata)
			if(it.Down)
					it.Down = false;
			else
					++it.Up;

			if (!empty(it.history))
				value(it).i2 = top(it.history).i2;	// copy right boundary of parent's range
			return true;
		}
		return false;
}


template < typename TText, typename TSpec >
inline void goNext(Iter< Index<TText, TSpec>, VSTree< TopDown< ParentLinks<ConstrainedTraversal<Absolute> > > > > &it) {
	unsigned walk_down=0,not_finished = 1;
	//std::cout << "do the iterator with AbsoluteeThreshol:"<<it.AbsoluteThreshold<<std::endl;
	// resembles the construction of the VLMM-core of the context tree
	// only takes nodes that occur at least AbsoluteThreshold times
	// note: that we do not enforce that the tree is balanced

	if(repLength(it) < it.MaxDepth){
		if((countOccurences(it)) >= it.AbsoluteThreshold) 
			walk_down = 1;
	}

	it.Down = false;
	it.Up = 0;
	while(not_finished)
	{
		// if threshold is reached walk down else go always right or up
		if(walk_down){
			// for avoiding the implicit $-edges
			
			//if(!isLeaf(it) && isRightTerminal(it)){
			//	// this must be possible
			//	
			//	goDown(it);
			//	goRight(it);
			//	while( (repLength(container(it), nodeUp(it)) == repLength(it)) && goRight(it))
			//		
			//	// if the current node remains to be right terminal
			//	// than we have found a node where multiple strings end
			//	// goDown was wrong and we have to go up again
			//	if(isLeaf(it) && (repLength(container(it), nodeUp(it)) == repLength(it))){
			//		goUp(it);
			//		//cout << value(it) << " = " << length(representative(it)) << " " << representative(it) << "  toFather:"<<parentEdgeLabel(it)<<"  hits: "<<length(getOccurences(it))<<endl;	
			//		if (!goRight(it))
			//			while (goUp(it) && !goRight(it));
			//		}
			//}
			//else
				if(!goDown(it) && !goRight(it))
					while (goUp(it) && !goRight(it));
		}
		else
			if (!goRight(it))
				while (goUp(it) && !goRight(it));
		
		if( countOccurences(it) >= it.AbsoluteThreshold){
			not_finished = 0;
			// we have to check if the current node can be extended with at least one letter if the following
			// condition is true:
			if(isRightTerminal(it) && length(parentEdgeLabel(it)) == 1){
				//cout << "In Abs Clausel  "<<value(it) << " = " << length(representative(it)) << " " << representative(it) << "  toFather:"<<parentEdgeLabel(it)<<"  hits: "<<length(getOccurences(it))<<endl;
				if(isLeaf(it)){
						walk_down=0;
						not_finished =1;
				}
				else
				{
						// test with the topdown iterator (without considering $-Edges) if non-$-Leaf
						Iter<Index<TText, TSpec>, VSTree< TopDown< > > >  copy(container(it),value(it));
						if( !goDown(copy)){
							not_finished=1;
							walk_down=0;
						}
				}
					
			}
		}
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
inline void goNext(Iter< Index<TText, TSpec>, VSTree< TopDown< ParentLinks<ConstrainedTraversal<Support> > > > > &it) {
	unsigned walk_down=0,not_finished = 1;
	// this iterator pecializations is needed for the assessment of nodes that are
	// only takes nodes that occur at least AbsoluteThreshold times
	// note: that we do not enforce that the tree is balanced

	if(repLength(it) < it.MaxDepth){
		if((unsigned)(getFrequency(it)) >= it.AbsoluteThreshold) 
			walk_down = 1;
	}

	it.Down = false;
	it.Up = 0;
	while(not_finished)
	{
		// if threshold is reached walk down else go always right or up
		if(walk_down){
			if(!goDown(it) && !goRight(it))
				while (goUp(it) && !goRight(it));
		}
		else{
			if (!goRight(it))
				while (goUp(it) && !goRight(it));
		}
		if( (unsigned)getFrequency(it) >= it.AbsoluteThreshold){
			not_finished = 0;
			// we have to check if the current node can be extended with at least one letter if the following
			// condition is true:
			if(isRightTerminal(it) && length(parentEdgeLabel(it)) == 1){
				//cout << "In Abs Clausel  "<<value(it) << " = " << length(representative(it)) << " " << representative(it) << "  toFather:"<<parentEdgeLabel(it)<<"  hits: "<<length(getOccurences(it))<<endl;
				if(isLeaf(it)){
						walk_down=0;
						not_finished =1;
				}
				else
				{
						// test with the topdown iterator (without considering $-Edges) if non-$-Leaf
						Iter<Index<TText, TSpec>, VSTree< TopDown< > > >  copy(container(it),value(it));
						if( !goDown(copy)){
							not_finished=1;
							walk_down=0;
						}
				}
					
			}
		}
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
	unsigned d = 0;
	// resembles the construction of the VLMM-core of the PST algorithm of Ron et.al
	// only take nodes where their counts or their parents occur more 
	// often than RelativeThreshold times by considering the actual context length
	// and the text length. This is known as relative frequency or empirical probability.
	// note: that we do not enforce that the tree is balanced
	// countSequences(container(it));
	//sequenceLength(SeqNr,index);
	if(repLength(it) < it.MaxDepth){

		d = repLength(it)-1;
		std::cout<<"AT TOP goNext.representative::"<<representative(it)<<"  relative Tresh:"<<it.RelativeThreshold<<"  length containeer:"<<length(container(it))<<"  d:"<<d<<"  countSeqs:" <<countSequences(container(it))<<endl;
		std::cout<<" Threshold: "<<floor(it.RelativeThreshold*(float)(length(container(it))-
			(d)*countSequences(container(it)) ))<<"  countOcc:" <<countOccurences(it)<<endl;
		if(countOccurences(it) >= floor(it.RelativeThreshold*(float)(length(container(it))-
		(d)*countSequences(container(it)) ))) 
			walk_down = 1;
	}

	it.Down = false;
	it.Up = 0;
	while(not_finished)
	{
		// if threshold is reached walk down else go always right or up
		if(walk_down){
			// this is temporary for avoiding the implicit $-edges
			/*
			if(!isLeaf(it) && isRightTerminal(it)){
				goDown(it);
				goRight(it);
				while(isRightTerminal(it) && goRight(it));
					
			// if the current node remains to be right terminal
			// than we have found a node where multiple strings end
			// goDown is not possible here
			if(isRightTerminal(it))
				goUp(it);
				if (!goRight(it))
				while (goUp(it) && !goRight(it));

			}
			skippe alle kinder mit leerem parentedgelabel
			repLength(container(it), nodeUp(it)) == repLength(it))
			*/
			if(!isLeaf(it) && isRightTerminal(it)){
				// this must be possible
				
				goDown(it);
				std::cout << value(it) << " = " << length(representative(it)) << " " << representative(it) << "  toFather:"<<parentEdgeLabel(it)<<"  hits: "<<length(getOccurences(it))<<endl;
				goRight(it);
				std::cout << value(it) << " = " << length(representative(it)) << " " << representative(it) << "  toFather:"<<parentEdgeLabel(it)<<"  hits: "<<length(getOccurences(it))<<endl;
				while( (repLength(container(it), nodeUp(it)) == repLength(it)) && goRight(it)){
					
				std::cout << value(it) << " = " << length(representative(it)) << " " << representative(it) << "  toFather:"<<parentEdgeLabel(it)<<"  hits: "<<length(getOccurences(it))<<endl;	
				}
					
				// if the current node remains to be right terminal
				// than we have found a node where multiple strings end
				// goDown was wrong and we have to go up again
				if(isLeaf(it) && (repLength(container(it), nodeUp(it)) == repLength(it))){
					goUp(it);
					std::cout << value(it) << " = " << length(representative(it)) << " " << representative(it) << "  toFather:"<<parentEdgeLabel(it)<<"  hits: "<<length(getOccurences(it))<<endl;	
					if (!goRight(it))
						while (goUp(it) && !goRight(it));
					}
			}
			else
				if (!goDown(it)&&!goRight(it))
					while (goUp(it) && !goRight(it));
		}
		else
			if (!goRight(it))
				while (goUp(it) && !goRight(it));
		
		if(!isLeaf(it)){
			d = min(repLength(it),it.MaxDepth)-1;}
		else{
			d = min(repLength(it)-1,it.MaxDepth)-1;}
		std::cout<<"AT BOTTOM goNext.representative::"<<representative(it)<<"  relative Tresh:"<<it.RelativeThreshold<<"  length containeer:"<<length(container(it))<<"  d:"<<d<<"  countSeqs:" <<countSequences(container(it))<<endl;
		std::cout<<" Threshold: "<<floor(it.RelativeThreshold*(float)(length(container(it))-
			(d)*countSequences(container(it)) ))<<"  countOcc:" <<countOccurences(it)<<endl;
	if( (countOccurences(it) >= floor(it.RelativeThreshold*(float)(length(container(it))-
	(d)*countSequences(container(it)) )) )  &&    ( !isRightTerminal(it) || length(parentEdgeLabel(it))>1 ))
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
	unsigned threshold;	
	float K;    //the name originates from the paper of Bühlmann (1999) 
	unsigned MaxDepth;
	float alpha;	// pseudocount for smoothing
	
};

struct BioPST{
	float threshold;		// r
	unsigned minSupport;	// N min
	float minConditionalProbability ;  // y
	unsigned MaxDepth;		// L
};

struct PST{
	float threshold ;				   // r
	float minEmpiricalProbability ;    // P min
	float minConditionalProbability ;  // y min
	float alpha ;					   // alpha :-)
	unsigned MaxDepth;				   // L
		
};

void setParameters(PST & params,float t,float minE,float minC,float alpha,unsigned d ){
	params.threshold = t;
	params.minEmpiricalProbability = minE;
	params.minConditionalProbability = minC;
	params.alpha = alpha;
	params.MaxDepth = d;
}
void setParameters(ContextTree & params,unsigned t,float K, unsigned d, float alpha){
	params.threshold = t;
	params.K = K;
	params.MaxDepth = d;
	params.alpha = alpha;
}


void setParameters(BioPST & params,float t, unsigned minSupport,float y,unsigned d){
	
	params.threshold = t;		
	params.minSupport = minSupport;	
	params.minConditionalProbability = y; 
	params.MaxDepth = d;		
}
template<typename TSpec = ContextTree>
struct VLMM;


////////////////
// VLMM
////////////////

template<typename TAlphabet,  typename TSpec>
class Graph<Automaton<TAlphabet, String<TAlphabet>, WordGraph<VLMM<TSpec> > > > 

{
	public:
		typedef typename VertexIdHandler<Graph>::Type TVertexIdManager;
		typedef typename EdgeIdHandler<Graph>::Type TEdgeIdManager;
		typedef typename VertexDescriptor<Graph>::Type TVertexDescriptor;
		typedef typename EdgeType<Graph>::Type TEdge;


		String<AutomatonEdgeArray<TEdge, TAlphabet> > data_vertex;		// List of tables
		TVertexIdManager data_id_managerV;
		TEdgeIdManager data_id_managerE;
		TVertexDescriptor data_root;
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
template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM< TSpec> > > > >::Type 
addIncompleteVertex(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > >& g) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	unsigned num = numVertices(g)+1;
	//resize(g.data_edge_label, (num) * table_length,Generous());
	appendValue(g.data_father, nilVal);
	appendValue(g.data_marked, false);
	resize(g.data_probability_vector,num,Generous());
	//append Alphabetsize many 0s
	for (unsigned int i = 0;i<table_length;++i)	
		append(g.data_probability_vector[num-1],0);
		

	return addVertex(g);
}

// new addVertex Specialization which should be called if the vlmm is set up
// it enlarges all tables at once
template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM< TSpec> > > > >::Type 
addAdditionalVertex(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > >& g) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph >::Type TEdge;
	typedef typename Size<TAlphabet>::Type TAlphabetSize;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	TSize table_length = ValueSize<TAlphabet>::VALUE;

	TVertexDescriptor vd = addVertex(g);
	// if new vertex vd is at the end,i.e. points to the last cell
	// a new entry has been created/appended in addAutomatonVertex

if (vd == length(g.data_vertex)-1) {
	unsigned int Length = length(g.data_vertex);
	//resize(g.data_edge_label, Length * table_length, Generous());
	//resize(g.data_father,length(g.data_vertex));
	//g.data_father[vd] = nilVal;
	appendValue(g.data_father, nilVal);
	//resize(g.data_count,length(g.data_vertex));
	//g.data_count[vd] = 0;
	appendValue(g.data_marked,false);
	appendValue(g.data_suffix_link,nilVal);
	appendValue(g.data_reverse_suffix_link, AutomatonEdgeArray<TEdge, TAlphabet>());
	
	resize(g.data_probability_vector,Length,Generous());
	//append Alphabetsize many 0s
	for (unsigned int i = 0;i<table_length;++i)	
		append(g.data_probability_vector[Length-1],0);
}

// else they filled a gap in between -> no new space to allocate

return vd;
}

// init tables suffix_link  and reverse_suffix_link
template<typename TAlphabet, typename TCargo, typename TSpec>
inline void
initMaps(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > >& g) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph >::Type TEdge;
	unsigned num = numVertices(g);
	resize(g.data_reverse_suffix_link,num);
	resize(g.data_suffix_link, num);
	//resize(g.data_auxiliary_link, num);

	// init table data_suffix_link with NIL value
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	for(unsigned int i = 0;i<num;++i){
			g.data_suffix_link[i] = nilVal;
			value(g.data_reverse_suffix_link, i) =  AutomatonEdgeArray<TEdge, TAlphabet>();
	}
	

	return;
}

// init tables suffix_link  and reverse_suffix_link
template<typename TAlphabet, typename TCargo, typename TSpec,typename TSize>
inline void
initGraph(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > >& vlmm,
		 TSize &size) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph >::Type TEdge;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	vlmm.data_root = 0;
	resize(vlmm.data_father,size,Exact());
	resize(vlmm.data_suffix_link,size,Exact());
	resize(vlmm.data_marked,size,Exact());
	resize(vlmm.data_id_managerV.data_in_use,size,Exact());
	resize(vlmm.data_probability_vector,size,Exact());
	for(TSize i = 0;i<size;++i){
		TVertexDescriptor dummy = i;
		assignValue(vlmm.data_id_managerV.data_in_use, dummy, false);
		setFather(vlmm,nilVal,dummy);
		assignValue(vlmm.data_suffix_link,dummy, nilVal);
		setMarked(vlmm,dummy,false);
		appendValue(vlmm.data_vertex, AutomatonEdgeArray<TEdge, TAlphabet>());
		appendValue(vlmm.data_reverse_suffix_link, AutomatonEdgeArray<TEdge, TAlphabet>());
		
	
		//append Alphabetsize many 0s
		for (unsigned int j = 0;j<table_length;++j)	
			append(vlmm.data_probability_vector[dummy],0);

	}
	return;
}

template<typename TAlphabet, typename TCargo, typename TSpec,typename TVertexDescriptor>
inline  void
setFather(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > >& wg,
		  TVertexDescriptor& father,
		  TVertexDescriptor& child) 
{
	SEQAN_CHECKPOINT

	wg.data_father[child] = father;
	return;
}

template<typename TAlphabet, typename TCargo, typename TSpec,typename TVertexDescriptor>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > > >::Type 
getFather(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > >& wg,
		  TVertexDescriptor & child) 
{
	return wg.data_father[child] ;
}


template<typename TAlphabet, typename TCargo, typename TSpec,typename TVertexDescriptor>
inline  void
setMarked(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > >& wg,
		  TVertexDescriptor& vertex,
		  bool state) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > > TGraph;
	wg.data_marked[vertex] = state;
	return;
}

template<typename TAlphabet, typename TCargo, typename TSpec ,typename TVertexDescriptor>
inline bool
isMarked(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > >& wg,
		  TVertexDescriptor & vertex) 
{

	return wg.data_marked[vertex] ;
}

template<typename TAlphabet, typename TCargo, typename TSpec ,typename TVertexDescriptor>
inline bool
isMarked(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > >const& wg,
		  TVertexDescriptor & vertex) 
{

	return wg.data_marked[vertex] ;
}


template<typename TAlphabet, typename TCargo, typename TSpec ,typename TVertexDescriptor>
inline  void
setSuffixLink(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > >& wg,
		  TVertexDescriptor& source,
		  TVertexDescriptor& target) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > > TGraph;
	SEQAN_ASSERT(idInUse(wg.data_id_managerV, source) == true)
	wg.data_suffix_link[source] = target;
	return;
}

template<typename TAlphabet, typename TCargo, typename TSpec ,typename TVertexDescriptor>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > > >::Type 
getSuffixLink(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > >& wg,
		  TVertexDescriptor & source) 
{
	SEQAN_ASSERT(idInUse(wg.data_id_managerV, source) == true)
	return wg.data_suffix_link[source] ;
}


template<typename TAlphabet, typename TCargo,typename TSpec , typename TVertexDescriptor, typename TChar>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > > >::Type 
getReverseSuffixLink(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > >& g,
			 TVertexDescriptor & vertex,
			 TChar const c) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > > TGraph;
	typedef typename Size<TAlphabet>::Type TSize;
	TAlphabet letter(c);
	return g.data_reverse_suffix_link[vertex].data_edge[(TSize) letter].data_target;
}


template<typename TAlphabet, typename TCargo, typename TSpec ,  typename TVertexDescriptor, typename TChar>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > > >::Type 
setReverseSuffixLink(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > >& g,
			 TVertexDescriptor &source,
			 TVertexDescriptor &target,
			 TChar const c) 
{
	SEQAN_CHECKPOINT

	//typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > >, TGraphSpec> TGraph;
	typedef typename Size<TAlphabet>::Type TSize;
	TAlphabet letter(c);
	return g.data_reverse_suffix_link[source].data_edge[(TSize) letter].data_target = target;
}

template<typename TSpec,typename TCargo,typename TAlphabet ,typename TVertexDescriptor>
inline TAlphabet 
getReverseSuffixLinkCharacter(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > > > &vlmm,
			   TVertexDescriptor &SLfather,
			   TVertexDescriptor &SLchild)
{
	typedef typename Size<TAlphabet>::Type TSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TSize i;
	for(i=0;i<table_length;++i){
		if(getReverseSuffixLink(vlmm,SLfather,i) == SLchild)
			return((TAlphabet)i);
	
	}
	SEQAN_ASSERT(getReverseSuffixLink(vlmm,SLfather,i) == SLchild);
	return ((TAlphabet)0);
}


template<typename TAlphabet, typename TCargo, typename TSpec ,  typename TVertexDescriptor>
inline void 
removeAllReverseSuffixLinks(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > >& g,
			 TVertexDescriptor &source) 
{
	SEQAN_CHECKPOINT

	typedef typename Size<TAlphabet>::Type TSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	TVertexDescriptor dummy;
	for(TSize i = 0;i<table_length;++i){
		if(getReverseSuffixLink(g,source,i) != nilVal){
				dummy = getReverseSuffixLink(g,source,i);
				setReverseSuffixLink(g,source,nilVal,i);
				setSuffixLink(g,dummy,nilVal);
		}
	}
}


//Graph<Automaton<TText, TCargo, WordGraph< VLMM < TSpec > > >, TGraphSpec> 
// currently this function can only be called for index of
// specialization: Index_ESA
template<typename TIndex,typename TIterSpec,typename TSpec,typename TCargo,typename TAlphabet ,typename TVertexDescriptor>
inline bool
initProbabilityVector(Iter<TIndex, VSTree< TopDown< ParentLinks<ConstrainedTraversal<TIterSpec> > > > > &it,
						 Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > > > &target,
						 TVertexDescriptor & node)
{
	//typedef Graph<Automaton<TAlphabet,TCargo,WordGraph< VLMM < TSpec > > >,TGraphSpec> TVlmm;
	
	typedef Iter<TIndex, VSTree< TopDown< > > >  TIter;

	TIter childs(it);
	if(!goDown(childs))
		return false;
	TAlphabet startCharacter;
	unsigned fatherLength = repLength(it);
	unsigned count = 0;
	//std::cout <<"Cilds:" << representative(childs) << std::endl;
	if( repLength(childs) > fatherLength){

		startCharacter = value(representative(childs),fatherLength);
		setProbability(target,node,startCharacter,(float)countOccurences(childs));
		count += countOccurences(childs);
	}
		while(goRight(childs)){

			// std::cout << representative(childs) << std::endl;
			if(repLength(childs) > fatherLength){
			startCharacter = value(representative(childs),fatherLength);
			setProbability(target,node,startCharacter,(float)countOccurences(childs));
			count += countOccurences(childs);
			}
		}
	if(count >0)
		return true;
	return false;

}

//Graph<Automaton<TText, TCargo, WordGraph< VLMM < TSpec > > >, TGraphSpec> 
// currently this function can only be called for index of
// specialization: Index_ESA
template<typename TIndex,typename TIterSpec,typename TSpec,typename TCargo,typename TAlphabet ,typename TVertexDescriptor,typename TChar>
inline void
initProbabilityVectorForLeaf(Iter<TIndex, VSTree< TopDown< ParentLinks<ConstrainedTraversal<TIterSpec> > > > > &it,
						 Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > > > &target,
						 TVertexDescriptor & node,
						 TChar & letter)
{
	//only one letter can be there, so set it to the number of counts found in the node
	setProbability(target,node,letter,(float)countOccurences(it));
	
}

// this function cuts the EdgeLabel of the node pointed to by the iterator.
// It is used in 2 ways, for leafs the parameter diff == 1, and for inner nodes diff is length(EdgeLabel-it.MaxDepth)
template<typename TIndex,typename TIterSpec,typename TSpec,typename TCargo,typename TAlphabet ,typename TVertexDescriptor,typename TDiff>
inline void
cutEdgeForLeaf(Iter<TIndex, VSTree< TopDown< ParentLinks<ConstrainedTraversal<TIterSpec> > > > > &it,
						 Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > > > &target,
						 TVertexDescriptor & father,
						 TVertexDescriptor & child,
						 TDiff & diff)
{
					String<TAlphabet> EdgeLabel = parentEdgeLabel(it);
					TAlphabet letter = value(EdgeLabel,length(EdgeLabel)-diff );
					String<TAlphabet> pref = prefix( EdgeLabel, length(EdgeLabel)-diff );
					//cout << "prefix of edgelabel "<<pref<<endl;
					addEdge(target,father,child,  pref);
					initProbabilityVectorForLeaf(it,target,child,letter);
					SEQAN_ASSERT(parseString2(target,father,pref)==child)
}

//Graph<Automaton<TText, TCargo, WordGraph< VLMM < TSpec > > >, TGraphSpec> 
// currently this function can only be called for index of
// specialization: Index_ESA
template<typename TIndex,typename TIterSpec,typename TSpec,typename TCargo,typename TAlphabet>
inline void 
buildSuffixTreeFromIndex(Iter<TIndex, VSTree< TopDown< ParentLinks<ConstrainedTraversal<TIterSpec> > > > > &it,
						 Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > > > &target)
{
	typedef Graph<Automaton<TAlphabet,TCargo,WordGraph< VLMM < TSpec > > > > TVlmm;
	typedef typename VertexDescriptor<TVlmm>::Type TVertexDescriptor;
	typedef Iter<TIndex, VSTree< TopDown< > > >  TIter;
	typedef Iter<TIndex, VSTree< TopDown< ParentLinks<ConstrainedTraversal<TIterSpec> > > > >  TIter2;
	
	TVertexDescriptor root = addIncompleteVertex(target);
	TVertexDescriptor child,father;
	child = root;

	initProbabilityVector(it,target,root);
	assignRoot(target,root);
	//go to the first node
	goNext(it);
	while(!atEnd(it)){
		
		
		// it.Down == true: the previous child becomes the new father
		// neither it.Down == true nor it.Up > 0: the father remains the same
		// it.Up == x: the new father is x levels up in the target
		//cout << "Before father: "<< father<< "  child:" <<child;
		if(it.Down){
			father = child;
		}
		while(it.Up>0){
			father = getFather(target,father);
			--it.Up;
		}
		//cout << "After father: "<< father<< "  child:" <<child<<endl;;

		
		child = addIncompleteVertex(target);
		// check if iterator
		//cout << "Node:"<<child<<"  "<<value(it) << " = " << repLength(it)<< " " << representative(it) << "  toFather:"<<parentEdgeLabel(it)<<"  hits: "<<length(getOccurences(it))<<std::endl;
		//DAVIDDEBUGcout <<  " " << representative(it) <<endl;  
		
		setFather(target,father,child);
		unsigned diff = 1;
		//std::cout <<"set Father\t";
		// set child relation with edgelabel
		//std::cout <<"ParentEdgeLAbel:"<<parentEdgeLabel(it)<<std::endl;
		// only if not a leaf
		if(!isLeaf(it)){
			if(length(representative(it)) <= it.MaxDepth){
				//std::cout<<"case1\t"<<representative(it)<<" ";
				if(initProbabilityVector(it,target,child))
				{
				String<TAlphabet> EdgeLabel = parentEdgeLabel(it);
				addEdge(target,father,child,EdgeLabel);
				}
				else{// this is actually a leaf in the sense that no further extension
					 // of the edgelabel is possible, still this sequence occurs more than once
					/*String<TAlphabet> EdgeLabel = parentEdgeLabel(it);
					TAlphabet letter = value(EdgeLabel,length(EdgeLabel)-1 );
					String<TAlphabet> pref = prefix( EdgeLabel, length(EdgeLabel)-1 );
					addEdge(target,father,child,  pref);
					initProbabilityVectorForLeaf(it,target,child,letter);*/
					cutEdgeForLeaf(it,target,father,child,diff); // diff is 1 by default
				}
			}
			else{
				diff = length(representative(it)) - it.MaxDepth;
				//std::cout<<"case2\t"<<representative(it)<<" ";
				//addEdge with one char less , da muss man auch nicht zählen
				//String<TAlphabet> EdgeLabel = parentEdgeLabel(it);
				////TAlphabet letter = value(EdgeLabel,length(EdgeLabel)-1 );
				//// for the non leaf case 
				//TAlphabet letter = value(EdgeLabel,length(EdgeLabel)-diff );
				////if(child == 13 || child ==4)std::cout<<"EdgeLabel: "<<EdgeLabel<<"RepLength: "<<repLength(it)<<"  prefix"<<prefix( EdgeLabel, length(EdgeLabel)-1 ) <<endl;
				//String<TAlphabet> pref = prefix( EdgeLabel, length(EdgeLabel)-diff );
				//addEdge(target,father,child,  pref);
				////if(child == 13 || child ==4)std::cout<<"addedEdge"<<"letter:"<<letter<<endl;
				//initProbabilityVectorForLeaf(it,target,child,letter);
				////if(child == 13 || child ==4) std::cout<<"initVector"<<endl;
				cutEdgeForLeaf(it,target,father,child,diff); 
			}
		}
		else{
			////std::cout<<"case3\t"<<representative(it)<<" ";
			////addEdge with one char less , da muss man auch nicht zählen
			//String<TAlphabet> EdgeLabel = parentEdgeLabel(it);
			//TAlphabet letter = value(EdgeLabel,length(EdgeLabel)-1 );
			//// for the non leaf case 
			////if(child == 13 || child ==4)std::cout<<"EdgeLabel: "<<EdgeLabel<<"RepLength: "<<repLength(it)<<"  prefix"<<prefix( EdgeLabel, length(EdgeLabel)-1 ) <<endl;
			//String<TAlphabet> pref = prefix( EdgeLabel, length(EdgeLabel)-1 );
			//addEdge(target,father,child,  pref);
			////if(child == 13 || child ==4)std::cout<<"addedEdge"<<"letter:"<<letter<<endl;
			//initProbabilityVectorForLeaf(it,target,child,letter);
			////if(child == 13 || child ==4) std::cout<<"initVector"<<endl;
			cutEdgeForLeaf(it,target,father,child,diff); // diff is 1 by default
		}
		//std::cout<<"created node:"<<child<<std::endl;
		//std::cout <<"set Edge\t";
		// go to next valid node and add it to the graph
		goNext(it);
	}
	// All valid nodes are now part of the vlmm
	return;
}

template<typename TSpec,typename TCargo,typename TAlphabet ,typename TVertexDescriptor>
inline TAlphabet 
getChildCharacter(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > > > &vlmm,
			   TVertexDescriptor &father,
			   TVertexDescriptor &child)
{
	typedef Graph<Automaton<TAlphabet,TCargo,WordGraph< VLMM < TSpec > > > > TVlmm;
	typedef typename Iterator<TVlmm, OutEdgeIterator>::Type TOutEdgeIterator;

	TOutEdgeIterator itout(vlmm,father);
	while(!atEnd(itout)){
		if(targetVertex(vlmm, getValue(itout)) == child)
			break;
	// *(TOutEdgeIterator) 

	goNext(itout);
	}
	return(itout.data_pos);
// look in the EdgeAutomaton of father where the son/daughter is

}

template<typename TSpec,typename TCargo,typename TAlphabet ,typename TVertexDescriptor,typename TChar>
inline void 
getSuffixChildLabel(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > > > &vlmm,
			   TVertexDescriptor &father,
			   TChar  Letter,
			   String<TAlphabet> &Label)
{
	typedef Graph<Automaton<TAlphabet,TCargo,WordGraph< VLMM < TSpec > > > > TVlmm;
	typedef typename EdgeType<TVlmm>::Type TEdgeStump;
	
//typedef typename Iterator<TVlmm, OutEdgeIterator >::Type TOutEdgeIterator;
	TEdgeStump* ed = findEdge(vlmm,father, Letter);
	append(Label,getCargo(ed),Exact());
/*	TOutEdgeIterator itout(vlmm,father);
	while(!atEnd(itout)){
		if(targetVertex(vlmm, getValue(itout)) == child)
			break;

	goNext(itout);
	}
    SEQAN_ASSERT(targetVertex(vlmm, getValue(itout)) == child)
	append(Label,getCargo(*(itout)));
	//append(Label,getProperty(vlmm.data_edge_label,*(itout)));
	//append(Label,getProperty(vlmm.data_edge_label, TEdgeDescriptor(father, childPosition)));
*/
	return;
}

template<typename TSpec,typename TCargo,typename TAlphabet ,typename TVertexDescriptor>
inline void 
getSuffixChildLabel(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > > > &vlmm,
			   TVertexDescriptor &father,
			   TVertexDescriptor &child,
			   String<TAlphabet> &Label)
{
	typedef Graph<Automaton<TAlphabet,TCargo,WordGraph< VLMM < TSpec > > > > TVlmm;
	//typedef typename EdgeDescriptor<TVlmm>::Type TEdgeDescriptor;
typedef typename Iterator<TVlmm, OutEdgeIterator >::Type TOutEdgeIterator;

	TOutEdgeIterator itout(vlmm,father);
	while(!atEnd(itout)){
		if(targetVertex(vlmm, getValue(itout)) == child)
			break;

	goNext(itout);
	}
    SEQAN_ASSERT(targetVertex(vlmm, getValue(itout)) == child)
	append(Label,getCargo(*(itout)));
	//append(Label,getProperty(vlmm.data_edge_label,*(itout)));
	//append(Label,getProperty(vlmm.data_edge_label, TEdgeDescriptor(father, childPosition)));

	return;
}

template<typename TSpec,typename TCargo,typename TAlphabet ,typename TVertexDescriptor>
inline void 
getSuffixChildLabel(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > > > const &vlmm,
			   TVertexDescriptor &father,
			   TVertexDescriptor &child,
			   String<TAlphabet> &Label)
{
	typedef Graph<Automaton<TAlphabet,TCargo,WordGraph< VLMM < TSpec > > > > TVlmm;
	//typedef typename EdgeDescriptor<TVlmm>::Type TEdgeDescriptor;
typedef typename Iterator<TVlmm, OutEdgeIterator >::Type TOutEdgeIterator;

	TOutEdgeIterator itout(vlmm,father);
	while(!atEnd(itout)){
		if(targetVertex(vlmm, getValue(itout)) == child)
			break;

	goNext(itout);
	}
	SEQAN_ASSERT(targetVertex(vlmm, getValue(itout)) == child)
	append(Label,getCargo(*(itout)));
	
	//append(Label,getProperty(vlmm.data_edge_label, TEdgeDescriptor(father, childPosition)));

	return;
}

template<typename TSpec,typename TCargo,typename TAlphabet ,typename TVertexDescriptor>
inline void 
getChildLabel(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > > > &vlmm,
			   TVertexDescriptor &father,
			   TVertexDescriptor &child,
			   String<TAlphabet> &Label)
{
	typedef Graph<Automaton<TAlphabet,TCargo,WordGraph< VLMM < TSpec > > > > TVlmm;
	//typedef typename EdgeDescriptor<TVlmm>::Type TEdgeDescriptor;
typedef typename Iterator<TVlmm, OutEdgeIterator >::Type TOutEdgeIterator;

	TOutEdgeIterator itout(vlmm,father);
	while(!atEnd(itout)){
		if(targetVertex(vlmm, getValue(itout)) == child)
			break;

	goNext(itout);
	}
    SEQAN_ASSERT(targetVertex(vlmm, getValue(itout)) == child)
	//TAlphabet letter = itout.data_pos;
	//position(it) = = itout.data_pos??
	append(Label,itout.data_pos);
	append(Label,getCargo(*(itout)));
	//append(Label,getProperty(vlmm.data_edge_label, TEdgeDescriptor(father, childPosition)));

	return;
}

template<typename TSpec,typename TCargo,typename TAlphabet,typename TChar ,typename TVertexDescriptor>
inline float
getProbability(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > > > &vlmm,
			   TVertexDescriptor &father,
			   TChar pos)
{
	TAlphabet letter(pos);
	return value(vlmm.data_probability_vector[father],(int)letter);
}

template<typename TSpec,typename TCargo,typename TAlphabet,typename TChar, typename TVertexDescriptor>
inline void
setProbability(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > > > &vlmm,
			   TVertexDescriptor &father,
			   TChar pos,
			   float prob)
{
	TAlphabet letter(pos);
	value(vlmm.data_probability_vector[father],(int)letter) = prob;
}


// splitPosition is the position in the array string edgelabel
// e.g. edgelabel = ACGT , the new edge should be GT than splitPosition is 2
template<typename TAlphabet, typename TCargo, typename TSpec , typename TVertexDescriptor, typename TChar,typename TPos>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > > >::Type 
splitEdge(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > > & vlmm,
			TVertexDescriptor & father,
			TChar  childCharacter,
			TPos  splitPosition)
{
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM <TSpec> > > > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TGraph>::Type TSize;
	TAlphabet letter(childCharacter);
	TVertexDescriptor child = vlmm.data_vertex[father].data_edge[(TSize) letter].data_target;
	//std::cout << "split edge at pos:"<< splitPosition<<" character:" <<(int)childCharacter<<" father" <<father<<endl;
	String<TAlphabet> edgeString;
	float number = getProbability(vlmm,father,childCharacter);
	getSuffixChildLabel(vlmm,father,letter,edgeString);
	//String<TAlphabet> edgeString = getProperty(vlmm.data_edge_label, TEdgeDescriptor(father, letter));
	//std::cout<<"before new node created, NumVertc:"<<numVertices(vlmm)<<std::endl;
	TVertexDescriptor newNode = addAdditionalVertex(vlmm);
	//std::cout<<"after node created, NumVertc:"<<numVertices(vlmm)<<std::endl;
	String<TAlphabet> newEdgeString = childCharacter;
		//std::cout << "eedestring:"<<edgeString<<" length: "<<length(edgeString)<<endl;
	SEQAN_ASSERT(splitPosition < (TPos)length(edgeString))
	if(splitPosition >= (TPos)length(edgeString)){
		std::cout<<"splitPosition >= edgeString"<<std::endl;
		exit(1);
	
	}
	append(newEdgeString,prefix(edgeString,splitPosition),Exact());
	addEdge(vlmm,father,newNode,newEdgeString );
	setFather(vlmm,father,newNode);
	setProbability(vlmm,newNode,value(edgeString,splitPosition),number);
	String<TAlphabet> EdgeLabel = suffix(edgeString,splitPosition);
	addEdge(vlmm,newNode,child,EdgeLabel);
	setFather(vlmm,newNode,child);

	return newNode;
}

template<typename TAlphabet, typename TCargo, typename TSpec , typename TVertexDescriptor, typename TLabel>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > > >::Type 
parseString2(Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > > & g,
			TVertexDescriptor & vertex,
			TLabel & label)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM <TSpec> > > > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TGraph>::Type TSize;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	TVertexDescriptor succ = vertex;
	TSize pos = 0;
	//std::cout <<"start parseString at node:"<< succ<<std::endl;
	TAlphabet letter = value(label,pos);
	TVertexDescriptor tmp = g.data_vertex[succ].data_edge[(TSize) letter].data_target;
	++pos;
	while (tmp != nilVal) {
		// can be substituted by using getChildLabel
		//edgeString = getProperty(g.data_edge_label, Edge);
		String<TAlphabet> edgeString;
		getSuffixChildLabel(g,succ,letter,edgeString);
		
		
		if( (pos == length(label)) && (length(edgeString) == 0))
				return tmp;

		if(length(edgeString) > 0){
			//std::cout <<"at  node:" << succ<<std::endl;
			//std::cout << "infix: "<<infix(label,pos,pos+length(edgeString))<<" EdgeString:"<<edgeString<<std::endl;
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



template<typename TSpec,typename TCargo,typename TAlphabet>
inline void 
addSuffixLinks(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > > > &vlmm)
{
	typedef Graph<Automaton<TAlphabet,TCargo,WordGraph< VLMM < TSpec > > > > TVlmm;
	typedef typename VertexDescriptor<TVlmm>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TVlmm>::Type TEdgeDescriptor;
	typedef typename Iterator<TVlmm, OutEdgeIterator >::Type TOutEdgeIterator;
	typedef typename Size<TVlmm>::Type TSize;
	
	TVertexDescriptor root = getRoot(vlmm);
	TVertexDescriptor father,target;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	TOutEdgeIterator itout(vlmm,root);
	setSuffixLink(vlmm,root,root);
	TAlphabet letter;
	for(;!atEnd(itout);goNext(itout)) {
		TVertexDescriptor v = targetVertex(vlmm, getValue(itout));
		letter = itout.data_pos;
		// means the suffix link must point to the root
		String<TAlphabet> edgeString;
		getSuffixChildLabel(vlmm,root,letter,edgeString);
		if(length(edgeString) == 0){
				setSuffixLink(vlmm,v,root);
				setReverseSuffixLink(vlmm,root,v,letter);
				
			}
		else{
				// getProperty(vlmm.data_edge_label, TEdgeDescriptor(root, letter) == edgelabel without first char
				target = parseString2(vlmm,root,edgeString);
				// it might occur that a suffix link target for n does not exist
				// because n is leaf which edge has be pruned
				if(target != nilVal){
					setSuffixLink(vlmm,v,target);
					setReverseSuffixLink(vlmm,target,v,letter);
				}
			}
	}
// Having set these trivial cases we can treat every node the same way


std::deque<TVertexDescriptor> queue;
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
			//std::cout << "go for SL of father:" << father<< " String:"<< walkDown<<std::endl;
			TVertexDescriptor fatherSuffixLink = getSuffixLink(vlmm,father);
			target = parseString2(vlmm,fatherSuffixLink,walkDown);
			// it might occur that a suffix link target for n does not exist
			// because n is leaf which edge has be pruned
			if(target != nilVal){
				setSuffixLink(vlmm,n,target);
				setReverseSuffixLink(vlmm,target,n,childCharacter);
			}
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

template<typename TCargo,typename TAlphabet,typename TSpec ,typename TVertexDescriptor>
inline void
turnNodeCountsIntoProbability(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > > > &vlmm,
			   TVertexDescriptor & node) 
{
	unsigned size = ValueSize<TAlphabet>::VALUE;
	float sum = 0;
	for(unsigned int pos = 0;pos< size ;++pos)
		sum = sum + getProbability(vlmm,node,pos);
	SEQAN_ASSERT(sum != 0);
	for(unsigned int pos = 0;pos< size;++pos)
		setProbability(vlmm,node,pos, (getProbability(vlmm,node,pos)/sum) );
}

 

template<typename TCargo,typename TAlphabet,typename TChar,typename TVertexDescriptor>
inline bool
extendNode(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < PST > > > > &vlmm,
			   TVertexDescriptor & node,
			   TChar childCharacter,
			   float countChar,
			   PST & parameters) 
{
	
		TAlphabet letter(childCharacter);
		//std::cout<<"extendNode on node "<<node<<std::endl;
		float gammaMin = (1 + parameters.alpha)*parameters.minConditionalProbability;
		
		if( 1 >= gammaMin){
			SEQAN_ASSERT(getProbability(vlmm,node,childCharacter) != 0);
			
			float ratio = 1/getProbability(vlmm,node,childCharacter);

			if(ratio >= parameters.threshold  || ratio <= 1/parameters.threshold )
				return true;
		}
	return false;

}

// check if we should create the new node
template<typename TCargo,typename TAlphabet,typename TChar,typename TVertexDescriptor>
inline bool
extendNode(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ContextTree> > > > &vlmm,
			   TVertexDescriptor & father,
			   TChar childCharacter,
			   float countChar,
			   ContextTree & param) 
{
	
		TAlphabet letter(childCharacter);
		unsigned size = ValueSize<TAlphabet>::VALUE;
	// the pruning criterion mentioned by Bühlmann & Mächler, taken from Rissanen 1983
		float fathersum = 0,sonsum=0,difference=0;
    
	for(unsigned int pos = 0;pos< size ;++pos){
		fathersum = fathersum + getProbability(vlmm,father,pos);
	}
	SEQAN_ASSERT(fathersum != 0)
	float pseudocounts = param.alpha * (float)size;
	fathersum += pseudocounts; 
	sonsum += pseudocounts + countChar;
	//std::cout<<"extendNode on node "<<node<<std::endl;
	//std::cout << " diff= "<<countChar * log(1/(getProbability(vlmm,father,childCharacter)/fathersum))<<endl;
	for(unsigned pos=0;pos<size;++pos){
		if(pos == (unsigned)letter)
				difference += countChar*((countChar+param.alpha)/sonsum) * log( ((countChar+param.alpha)/sonsum)/((getProbability(vlmm,father,letter)+param.alpha)/fathersum));
		else
				difference += param.alpha * (param.alpha/sonsum)*
						  log( (param.alpha/sonsum)/( (getProbability(vlmm,father,pos)+param.alpha)/fathersum));

	}
	if ( difference  <= param.K )
			return false;

	return true;
}


template<typename TCargo,typename TChar,typename TVertexDescriptor>
inline bool
extendNode(Graph<Automaton<AminoAcid, TCargo , WordGraph < VLMM < BioPST > > > > &vlmm,
			   TVertexDescriptor & father,
			   TChar childCharacter,
			   float countChar,
			   BioPST & parameters) 
{
	
	AminoAcid letter(childCharacter);
	//std::cout<<"extendNode(BioPST) on node "<<node<<std::endl;
	unsigned alphaSize = ValueSize<AminoAcid>::VALUE;
	float fathersum = 0;
    
	for(unsigned int pos = 0;pos< alphaSize ;++pos){
		fathersum = fathersum + getProbability(vlmm,father,pos);
	}
	SEQAN_ASSERT(fathersum != 0)
		
	SEQAN_ASSERT(getProbability(vlmm,father,letter) != 0);
			
	float ratio = 1/getProbability(vlmm,father,letter)/fathersum;

	if(ratio >= parameters.threshold  || ratio <= 1/parameters.threshold )
		return true;
		
	return false;

}


/* 
	In general this function can only be called for nodes which still contain the pattern counts and
	have not be converted into probabilities yet.
*/
template<typename TCargo,typename TAlphabet ,typename TVertexDescriptor>
inline bool
pruneNode(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ContextTree > > > > &vlmm,
			   TVertexDescriptor &son,
			   ContextTree & param) 
{
    unsigned size = ValueSize<TAlphabet>::VALUE;
	// the pruning criterion mentioned by Bühlmann, taken from Rissanen 1983
	float fathersum = 0,sonsum = 0, difference = 0;
    
	SEQAN_ASSERT(getSuffixLink(vlmm,son) != getNil<TVertexDescriptor>())
	TVertexDescriptor father = getSuffixLink(vlmm,son);

	for(unsigned int pos = 0;pos< size ;++pos){
		sonsum = sonsum + getProbability(vlmm,son,pos);
		fathersum = fathersum + getProbability(vlmm,father,pos);
	}
	SEQAN_ASSERT(sonsum != 0)
	SEQAN_ASSERT(fathersum != 0)
	//add the pseudocounts to avoid zero counts
	fathersum += param.alpha * (float)size; 
	sonsum += param.alpha * (float)size; 
	for(int pos = 0;pos< ValueSize<TAlphabet>::VALUE;++pos){
		//if(getProbability(vlmm,father,pos) > 0 && getProbability(vlmm,son,pos) > 0) it is checked in trainvlmm that alpha >0
			difference += (getProbability(vlmm,son,pos)+param.alpha) * ((getProbability(vlmm,son,pos)+param.alpha)/sonsum)*
						  log(( (getProbability(vlmm,son,pos)+param.alpha)/sonsum)/( (getProbability(vlmm,father,pos)+param.alpha)/fathersum));
		
	}
	if(difference <= param.K)
		return true;

 return false;
}


template<typename TCargo,typename TVertexDescriptor>
inline bool
pruneNode(Graph<Automaton<AminoAcid, TCargo , WordGraph < VLMM < BioPST > > > > &vlmm,
			   TVertexDescriptor &son,
			   BioPST & parameters) 
{
	// the pruning criteria defined by Bejeran & Yona 2001 for BioPST
	unsigned alphaSize = ValueSize<AminoAcid>::VALUE;
	float fathersum = 0,sonsum = 0,ratio = 0;
    
	SEQAN_ASSERT(getSuffixLink(vlmm,son) != getNil<TVertexDescriptor>())
	TVertexDescriptor father = getSuffixLink(vlmm,son);

	for(unsigned int pos = 0;pos< alphaSize ;++pos){
		sonsum = sonsum + getProbability(vlmm,son,pos);
		fathersum = fathersum + getProbability(vlmm,father,pos);
	}
	SEQAN_ASSERT(sonsum != 0)
	SEQAN_ASSERT(fathersum != 0)

	
	for(unsigned pos = 0;pos< alphaSize;++pos){
		
		if(getProbability(vlmm,son,pos) >= parameters.minConditionalProbability){
			if(getProbability(vlmm,father,pos) == 0)
				ratio = 0;
			else
				ratio = (getProbability(vlmm,son,pos)/sonsum)/(getProbability(vlmm,father,pos)/fathersum);

			if(ratio >= parameters.threshold  || ratio <= 1/parameters.threshold )
				return false;
		}
	}
 return true;
}

template<typename TCargo,typename TAlphabet ,typename TVertexDescriptor>
inline bool
pruneNode(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < PST > > > > &vlmm,
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

template<typename TCargo,typename TAlphabet ,typename TVertexDescriptor>
inline void
smoothNode(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < PST > > > > &vlmm,
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


template<typename TCargo,typename TAlphabet ,typename TVertexDescriptor, typename TVLMMSpec>
inline void
smoothNode(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec > > > > &vlmm,
			   TVertexDescriptor & node,
			   TVLMMSpec & parameters) 
{
	// simple pseudocounts scalable with alpha
	unsigned alphaSize = ValueSize<TAlphabet>::VALUE;
	float sum = 0;
	for(unsigned int pos = 0;pos< alphaSize ;++pos)
		sum = sum + getProbability(vlmm,node,pos);
	SEQAN_ASSERT(sum != 0);
	// add the pseudocounts
	sum += parameters.alpha * (float)alphaSize; 
	for(unsigned int pos = 0;pos< alphaSize;++pos)
		setProbability(vlmm,node,pos, ( (getProbability(vlmm,node,pos)+parameters.alpha)/sum) );
}

float Q[24][24]=  // Dayhoff substitution matrix
{
{0.986693	, 0.000109553	, 0.000398041	, 0.00056237	, 0.000120508	, 0.000339613	, 0.000971367	, 0.00211437	, 7.66868e-005	, 0.000241016	, 0.000346917	, 0.00020815	, 0.000105901	, 7.30351e-005	, 0.00125986	, 0.00281915	, 0.00215454	, 0.			, 7.30351e-005	, 0.00133289	,0., 0., 0., 0.},
{0.000233351, 0.99135		, 0.000132232	, 0.			, 7.77837e-005	, 0.000933404	, 0				,7.77837e-005	, 0.000801172	, 0.000233351	, 0.000132232	, 0.00371028	, 0.000132232	, 5.44486e-005	, 0.000521151	, 0.00106564	, 0.000155567	, 0.000210016	, 2.33351e-005	, 0.000155567	,0., 0., 0., 0.},
{0.000857731, 0.000133775	, 0.982169		, 0.00418636	, 0.			, 0.000393455	, 0.000739695	, 0.00122758	, 0.00177842	, 0.000283287	, 0.000291157	, 0.00253385	, 0.			, 5.50837e-005	, 0.000212466	, 0.00339945	, 0.00132988	, 2.36073e-005	, 0.000283287	, 0.000102298	,0., 0., 0., 0.},
{0.00104535	, 0.0			, 0.0036112		, 0.985895		, 0.			, 0.000515886	, 0.00564081	, 0.00109965	, 0.000291883	, 8.82437e-005	, 0.			, 0.000576978	, 0.			, 0				, 6.78798e-005	, 0.000665222	, 0.000386915	, 0.			, 0.			, 0.000115396	,0., 0., 0., 0.},
{0.000313665, 9.50499e-005	, 0.			, 0.			, 0.997339		,		0.		, 0.			, 9.50499e-005	, 9.50499e-005	, 0.000161585	, 0.			, 0.			, 0.			, 0.			, 9.50499e-005	, 0.00111208	, 9.50499e-005	, 0.			, 0.00028515	, 0.000313665	,0., 0., 0., 0.},
{0.000773469, 0.000998024	, 0.000415844	, 0.000632082	, 0.			, 0.987624		, 0.00350972	, 0.000249506	, 0.002021		, 6.6535e-005	, 0.000623765	, 0.00122258	, 0.000166337	, 0.			, 0.000773469	, 0.000390893	, 0.000307724	, 0.			, 0.			, 0.000224555	,0., 0., 0., 0.},
{0.00170869	, 0.0			, 0.000603821	, 0.00533804	, 0.			, 0.00271077	, 0.986427		, 0.000719447	, 0.000147744	, 0.000224827	, 9.63545e-005	, 0.000668058	, 4.49654e-005	, 0.			, 0.000256945	, 0.000552432	, 0.000199133	, 0.			, 6.42363e-005	, 0.000237674	,0., 0., 0., 0.},
{0.00207892	, 3.59054e-005	, 0.000560125	, 0.000581668	, 3.59054e-005	, 0.000107716	, 0.000402141	, 0.99348		, 3.59054e-005	, 0.			, 6.10392e-005	, 0.000215433	, 2.51338e-005	, 6.10392e-005	, 0.000175937	, 0.00161574	, 0.000179527	, 0.			, 0.			, 0.000348283	,0., 0., 0., 0.},
{0.000198745, 0.000974795	, 0.00213887	, 0.000406953	, 9.46402e-005	, 0.00229976	, 0.000217673	, 9.46402e-005	, 0.991217		, 2.83921e-005	, 0.000378561	, 0.000217673	, 0.			, 0.00018928	, 0.000473201	, 0.000246065	, 0.000132496	, 2.83921e-005	, 0.000378561	, 0.000283921	,0., 0., 0., 0.},
{0.000569298, 0.000258772	, 0.000310526	, 0.000112135	, 0.000146637	, 6.90059e-005	, 0.000301901	, 0.			, 2.58772e-005	, 0.987225		, 0.00218231	, 0.000370907	, 0.000491667	, 0.000776316	, 6.03801e-005	, 0.000172515	, 0.00111272	, 0.			, 0.000112135	, 0.00570161	,0., 0., 0., 0.},
{0.000354108, 6.33666e-005	, 0.000137916	, 0.			, 0.			, 0.000279559	, 5.59117e-005	, 6.33666e-005	, 0.000149098	, 0.000943045	, 0.994677		, 0.000145371	, 0.000771582	, 0.000622484	, 0.00016028	, 0.000119278	, 0.000193827	, 4.84568e-005	, 8.57313e-005	, 0.00112942	,0., 0., 0., 0.},
{0.000225336, 0.00188571	, 0.00127295	, 0.000336028	, 0.			, 0.000581131	, 0.00041114	, 0.000237196	, 9.09252e-005	, 0.000169991	, 0.000154178	, 0.992548		, 0.000355794	, 0.			, 0.000169991	, 0.00066415	, 0.000790654	, 0.			, 3.95327e-005	, 7.72056e-005	,0., 0., 0., 0.},
{0.000625429, 0.000366631	, 0.			, 0.			, 0.			, 0.00043133	, 0.000150966	, 0.000150966	, 0.			, 0.00122929	, 0.00446427	, 0.00194099	, 0.987491		, 0.000366631	, 8.6266e-005	, 0.00043133	, 0.000603862	, 0.			, 0.			, 0.00166062	,0., 0., 0., 0.},
{0.000159996, 5.59986e-005	, 5.59986e-005	, 0.			, 0.			, 0.			, 0.			, 0.000135997	, 0.000159996	, 0.000719982	, 0.00133597	, 0.			, 0.000135997	, 0.994544		, 5.59986e-005	, 0.000319992	, 7.9998e-005	, 7.9998e-005	, 0.00207995	, 7.9998e-005	,0., 0., 0., 0.},
{0.00216589	, 0.000420622	, 0.000169504	, 6.27794e-005	, 6.27794e-005	, 0.000583849	, 0.000251118	, 0.000307619	, 0.000313897	, 4.39456e-005	, 0.000269951	, 0.000269951	, 2.51118e-005	, 4.39456e-005	, 0.992548		, 0.00168877	, 0.00045829	, 0.			, 0.			, 0.000313897	,0., 0., 0., 0.},
{0.00353024	, 0.00062648	, 0.00197547	, 0.000448139	, 0.000535023	, 0.000214924	, 0.000393265	, 0.00205778	, 0.000118894	, 9.1457e-005	, 0.000146331	, 0.000768239	, 9.1457e-005	, 0.000182914	, 0.0012301		, 0.984032		, 0.0031827		, 7.77384e-005	, 0.000100603	, 0.000196633	,0., 0., 0., 0.},
{0.00320656	, 0.000108697	, 0.000918491	, 0.000309787	, 5.43486e-005	, 0.00020109	, 0.000168481	, 0.000271743	, 7.6088e-005	, 0.000701096	, 0.000282613	, 0.00108697	, 0.000152176	, 5.43486e-005	, 0.000396744	, 0.00378266	, 0.987092		, 0.			, 0.000125002	, 0.00101088	,0., 0., 0., 0.},
{0.			, 0.000818633	, 9.09592e-005	, 0.			, 0.			, 0.			, 0.			, 0.			, 9.09592e-005	, 0.			, 0.000394156	, 0.			, 0.			, 0.000303197	, 0.			, 0.000515435	, 0.			, 0.997605		, 0.000181918	, 0				,0., 0., 0., 0.},
{0.000212704, 3.19057e-005	, 0.000382868	, 0.			, 0.000319057	, 0.			, 0.000106352	, 0.			, 0.000425409	, 0.000138258	, 0.00024461	, 0.000106352	, 0.			, 0.00276516	, 0.			, 0.000233975	, 0.00024461	, 6.38113e-005	, 0.994544		, 0.000180799	,0., 0., 0., 0.},
{0.00179442	, 9.83243e-005	, 6.39108e-005	, 8.35756e-005	, 0.000162235	, 0.000132738	, 0.0001819		, 0.000476873	, 0.000147486	, 0.00324962	, 0.00148961	, 8.35756e-005	, 0.000378548	, 4.91621e-005	, 0.000245811	, 0.000211397	, 0.000914416	, 0.			, 8.35756e-005	, 0.990153		,0., 0., 0., 0.},
{0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
{0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
{0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
{0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}
};

// this is the position/node specific acdlculation of pseudocounts depending on amino acid
// substitution probabilities
template<typename TCargo,typename TVertexDescriptor>
inline void
smoothNode(Graph<Automaton<AminoAcid, TCargo , WordGraph < VLMM <BioPST > > > > &vlmm,
			   TVertexDescriptor & node,
			   BioPST & parameters) 
{

	unsigned alphaSize = 20,pos=0;
	float Bs =5; // == m taken from Henikoff & Henikoff 1996 or == mü in Bejerano and Yona 2001
	float Rs =0;//will count how many aminoacids have a non-zero entry at node
	float ba[20];  // these are the relative pseudocounts for every amino acid
	float sum = 0;
	for(;pos< alphaSize ;++pos){
		sum = sum + getProbability(vlmm,node,pos);
		if(getProbability(vlmm,node,pos) >0)
			++Rs;
	}
	SEQAN_ASSERT(sum != 0);
	Bs = Bs * Rs;     //Bs = mü * Rs
	//calculate the ba`s,i.e. the proportion of pseudocounts/Bs for every amino acid
	for(unsigned aminoacid = 0;aminoacid< alphaSize;++aminoacid){
		ba[aminoacid] *= Bs;
		for(pos = 0;pos< alphaSize;++pos){
			ba[aminoacid] = getProbability(vlmm,node,pos)/sum*Q[aminoacid][pos];
		}
		
	}
	
	sum += Bs;  // add the overall pseudocounts to the node
	for(pos = 0;pos< alphaSize;++pos)
		setProbability(vlmm,node,pos, ( (getProbability(vlmm,node,pos)+ba[pos])/sum) );
}

template<typename TAlphabet,typename TCargo,typename TSpec ,typename TVLMMSPec>
inline void
pruneTree(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > > > &vlmm,
			   TVLMMSPec & parameters) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator >::Type TVertexIterator;

	TVertexDescriptor root = getRoot(vlmm);
	
	String<bool> original;
	// remember which nodes have been there
	resizeVertexMap(vlmm, original);
	TVertexIterator it(vlmm);
	for(;!atEnd(it);goNext(it)) {
		TVertexDescriptor node =  getValue(it);
		assignProperty(original, node, true);
		// change counts of all nodes into probabilities
		//std::cout<<" at node:"<<node<<endl;
		//turnNodeCountsIntoProbability(vlmm,node);
	}
	
	// root cannot be pruned
	setMarked(vlmm,root,true);
	//assignProperty(marked,root,true);
	TVertexDescriptor dummy = 0;
	//std::cout<<" start pruning from the root"<<std::endl;
	pruneTreeRecursivelyFast(vlmm,root,original,parameters,dummy);

		/*for(unsigned i =0;i<numVertices(vlmm);++i){
				std::cout << isMarked(vlmm,i)<< "  ";
		}
		std::cout <<std::endl;
		for(unsigned i =0;i<length(original);++i){
				std::cout << original[i]<< "  ";
		}
		std::cout <<std::endl;

std::cout << "ready with the function\n";*/

 return;
}

// recursive walk over reverse suffix links
template<typename TAlphabet,typename TCargo,typename TSpec ,typename TVertexDescriptor,typename TMarked,typename TVLMMSpec>
inline void
pruneTreeRecursively(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > > > &vlmm,
				TVertexDescriptor & node,
				TMarked & original,
			    TVLMMSpec & parameters) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TAlphabet>::Type TSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	

	TVertexDescriptor next;
	// check all outgoing reverse suffix links
	for(unsigned int pos = 0;pos< table_length;++pos){
		next = getReverseSuffixLink(vlmm,node,pos);
		if(next != nilVal){
			//std::cout <<" start Recursion from node:"<<node<<" char:"<<pos<<std::endl;
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

			//std::cout <<"..check potential nodes above node:" <<node<<std::endl;
			TVertexDescriptor potVertex;
			unsigned lastVertex = 0;
			for(unsigned pos = 1;pos < length(childLabel);++pos ){
				//std::cout << "check string:" << prefix(childLabel,pos)<<" from node:"<< target<<std::endl;
				// potVertex = potentialVertex is the potential SuffixLinkFather of the node to be created
				String<TAlphabet> prefixChildLabel = prefix(childLabel,pos);
				potVertex = parseString2(vlmm,target,prefixChildLabel);
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
					std::cout <<"startChar" << startChar<<std::endl;
					father = splitEdge(vlmm,father,startChar,0);
					smoothNode(vlmm,father,parameters);
					setSuffixLink(vlmm,father,root);
					setReverseSuffixLink(vlmm,root,father,startChar);
					setMarked(vlmm,father,true);
				}
			
		} // else case
	} // Label > 1
	
	//std::cout << "check keeping of node:" << node;
	// all potential nodes are build
	if( isMarked(vlmm,node) || (! pruneNode(vlmm,node,parameters)) ){
		// node should be kept
		//std::cout <<" keep it"<<std::endl;
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
		//std::cout <<"  delete it"<<std::endl;
		setMarked(vlmm,node,false);
	//assignProperty(marked,node,false);
	}

	//std::cout <<"finished node:"<<node<<std::endl;
 return;
}


// recursive walk over reverse suffix links
template<typename TAlphabet,typename TCargo,typename TSpec ,typename TVertexDescriptor,typename TMarked,typename TVLMMSpec,typename TChar>
inline void
pruneTreeRecursivelyFast(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > > > &vlmm,
				TVertexDescriptor & node,
				TMarked & original,
			    TVLMMSpec & parameters,
				TChar & letter) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TAlphabet>::Type TSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	

	TVertexDescriptor next;
	// check all outgoing reverse suffix links
	for(unsigned int pos = 0;pos< table_length;++pos){
		next = getReverseSuffixLink(vlmm,node,pos);
		if(next != nilVal){
			//std::cout <<" start Recursion from node:"<<node<<" char:"<<pos<<std::endl;
			pruneTreeRecursivelyFast(vlmm,next,original,parameters,pos);
		}
		
	}

	// the root is excluded from further processing
	if(isRoot(vlmm,node)){ 
	
		return ;
			
	}

	TVertexDescriptor father = getFather(vlmm,node);
	String<TAlphabet> childLabel;
	TAlphabet childCharacter;
	TVertexDescriptor target;
	getChildLabel(vlmm,father,node,childLabel);
	childCharacter = value(childLabel,0);
	float countChar = getProbability(vlmm,father,childCharacter);
	if(! isRoot(vlmm,father) ){
	
	//NOTE we don´t check if the suffixlink is valid, because we assume that up to here 
	// all inner nodes are consistent, i.e. all nodes which do not have a valid suffixlink taget
	//have been deleted before
	target = getSuffixLink(vlmm,father);
	}
	else{
	  // if father == root check only for the first possible extension
		  // extension longer than one are not checked, because it is already known that there is no 
		  // node one context shorter elsewhere in the tree

		  if(length(childLabel)>1){
					TVertexDescriptor root = getRoot(vlmm);
					childCharacter = value(childLabel,1);
					if(extendNode(vlmm,father,childCharacter,countChar,parameters))
					{
						
						father = splitEdge(vlmm,father,value(childLabel,0),0);
						//smoothNode(vlmm,father,parameters);
						setSuffixLink(vlmm,father,root);
						setReverseSuffixLink(vlmm,root,father,value(childLabel,0));
						setMarked(vlmm,father,true);
						if(father < length(original))
									original[father] = false;

					}
					//supply with childLabel
					
					childLabel = suffix(childLabel,1);
					childCharacter = value(childLabel,0);
					target = root;
		  }
	} // else case for root
	// check for possible extended 
	if(length(childLabel) > 1){
	

		
		// seek the node in the tree which could be the suffix link father of the potential new node
		
		//if(! isRoot(vlmm,father) ){
			
		
			// we first check if there exist a node on the edge of SLfather->SLnode
			// if not, nothing needs to be done
			unsigned lastVertex = 0,lastSplit=0;
			unsigned pos;
			unsigned labelLength = length(childLabel);
			//std::cout <<"..check potential nodes above node:" <<node<<" chidLAbel"<<childLabel<<std::endl;
			TVertexDescriptor potVertex = vlmm.data_vertex[target].data_edge[(TSize) childCharacter].data_target;
			SEQAN_ASSERT(potVertex != nilVal)
			// we have a node which may lead to a new node on the edge father->node
			// how far down is this node ?
			String<TAlphabet> edgeString;
			getSuffixChildLabel(vlmm,target,potVertex,edgeString);
			pos = 1 + length(edgeString);
			//pos = 1 + length(getProperty(vlmm.data_edge_label,TEdgeDescriptor(target,value(childLabel,lastVertex))));

			//while(! potVertex == getSuffixLink(vlmm,node) ){
			SEQAN_ASSERT(pos<=labelLength)
			while(pos != labelLength){
			
					//std::cout <<"another run on the edge"<<std::endl;
					if( (potVertex < length(original)) && original[potVertex] && extendNode(vlmm,potVertex,value(childLabel,pos),countChar,parameters) ){
						//std::cout <<"split edge at node: "<<node<<std::endl;
						father = splitEdge(vlmm,father,childCharacter,pos-1-lastSplit);
						//remember where the last node has been split
						lastSplit=pos;
						//smoothNode(vlmm,father,parameters);
						setSuffixLink(vlmm,father,potVertex);
						setReverseSuffixLink(vlmm,potVertex,father,letter);
						setMarked(vlmm,father,true);
						//std::cout <<"created node: "<<father<<" by checking above node:"<<node<<std::endl;
						if(father < length(original))
								original[father] = false;

						// we created a new node so have to retain the starting character for node
						// from the new father
						childCharacter = value(childLabel,pos);
					
				 }
					lastVertex = pos;
					// update the next vertex on edge SLfather->SLnode
					//std::cout << " childLabel:" << childLabel<<std::endl;
					target = potVertex;
					potVertex = vlmm.data_vertex[target].data_edge[(TSize) value(childLabel,lastVertex)].data_target;
					// we have a node which may lead to a new node on the edge father->node
					// how far down is this node ?
					String<TAlphabet> edgeString2;
					getSuffixChildLabel(vlmm,target,potVertex,edgeString2);
					pos += 1 + length(edgeString2);
					//pos += 1 + length(getProperty(vlmm.data_edge_label,TEdgeDescriptor(target,value(childLabel,lastVertex))));
			}

			
		//}
	} // Label > 1

	//std::cout << "check keeping of node:" << node;
	// all potential nodes are build
	if( (getSuffixLink(vlmm,node) != nilVal) && (isMarked(vlmm,node) || (! pruneNode(vlmm,node,parameters))) ){
		// node should be kept
		//std::cout <<" keep it"<<std::endl;
		setMarked(vlmm,node,true);
		// bubble up and set all nodes on the way to the root true
		TVertexDescriptor fatherSuffixLink = getSuffixLink(vlmm,node);
		while( !isRoot(vlmm,fatherSuffixLink) ){
			if(! isMarked(vlmm, fatherSuffixLink) ){
				setMarked(vlmm,fatherSuffixLink,true);
			}
			else
				break;
			fatherSuffixLink = getSuffixLink(vlmm,fatherSuffixLink);
		}
		// smooth node during deletion of nodes
		//smoothNode(vlmm,node,parameters);
	}
	else
	{	// we are sure thatr we have to remove the whole subtree below
		// node, because if node has no valid suffixLink means that either
		// node only exits because the edge has been pruned or the potential
		// suffixLink target did not reach the threshold (not in the monotonous decreasing case)
		if(getSuffixLink(vlmm,node) == nilVal)
				removeSubtree(vlmm,node);
	// delete node
		//std::cout <<"  delete it"<<std::endl;
		//setMarked(vlmm,node,false);    default is false
	//assignProperty(marked,node,false);
	}

	//std::cout <<"finished node:"<<node<<std::endl;
 return;
}

template<typename TAlphabet,typename TCargo,typename TVLMMSpec,typename TVertexDescriptor >
inline void
deleteProbabilityVector( Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec > > > > &vlmm,
					 TVertexDescriptor const trashNode )
{
	typedef typename Size<TAlphabet>::Type TSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	for(TSize i=0;i<table_length;i++)
	{
		setProbability(vlmm,trashNode,i,0);
	}
}

template<typename TAlphabet,typename TCargo,typename TVLMMSpec,typename TVertexDescriptor >
inline void
removeVertex( Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec > > > > &vlmm,
					 TVertexDescriptor &trashNode )
{
	SEQAN_CHECKPOINT
	typedef typename Size<TAlphabet>::Type TSize;
	
	typedef Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec > > > > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	SEQAN_ASSERT(idInUse(vlmm.data_id_managerV, trashNode) == true)
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	if(trashNode < length(vlmm.data_vertex)-1)
	{
		removeOutEdges(vlmm,trashNode); // Remove all outgoing edges
		TVertexDescriptor dummy = getFather(vlmm,trashNode);
		TAlphabet letter = getChildCharacter(vlmm,dummy,trashNode);
		TEdgeDescriptor ed = &vlmm.data_vertex[dummy].data_edge[(TSize)letter];
		assignTarget(ed, nilVal);
		releaseId(vlmm.data_id_managerE, _getId(ed));
		setFather(vlmm,nilVal,trashNode);
		dummy = getSuffixLink(vlmm,trashNode);
		if(dummy != nilVal){
			letter = getReverseSuffixLinkCharacter(vlmm,dummy,trashNode);
			setReverseSuffixLink(vlmm,dummy,nilVal,letter);
			setSuffixLink(vlmm,trashNode,nilVal);
		}
		removeAllReverseSuffixLinks(vlmm,trashNode);
		setMarked(vlmm,trashNode,false);
		deleteProbabilityVector(vlmm,trashNode);
		}
	else{//delete entries in the tables
		unsigned int Length = length(vlmm.data_vertex)-1;
		TVertexDescriptor dummy = getFather(vlmm,trashNode);
		TAlphabet letter = getChildCharacter(vlmm,dummy,trashNode);
		vlmm.data_vertex[dummy].data_edge[(TSize) letter].data_target = nilVal;
		dummy = getSuffixLink(vlmm,trashNode);
		if(dummy != nilVal){
			letter = getReverseSuffixLinkCharacter(vlmm,dummy,trashNode);
			setReverseSuffixLink(vlmm,dummy,nilVal,letter);
		}
		resize(vlmm.data_vertex,Length,Generous());
		resize(vlmm.data_father,Length,Generous());
		resize(vlmm.data_marked,Length,Generous());
		resize(vlmm.data_suffix_link,Length,Generous());
		resize(vlmm.data_reverse_suffix_link,Length,Generous());
		resize(vlmm.data_probability_vector,Length,Generous());
	}

	releaseId(vlmm.data_id_managerV, trashNode); // Release id

}

template<typename TAlphabet,typename TCargo,typename TVLMMSpec,typename TVertexDescriptor >
inline void
removeSubtree( Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec > > > > &vlmm,
			  TVertexDescriptor &head)
{
	// recursive removal of all nodes in the subtree starting at the head node
	typedef Graph<Automaton<TAlphabet,TCargo,WordGraph< VLMM < TVLMMSpec > > > > TVlmm;
	typedef typename Iterator<TVlmm, OutEdgeIterator>::Type TOutEdgeIterator;
	SEQAN_ASSERT(idInUse(vlmm.data_id_managerV, head) == true)
	TOutEdgeIterator itout(vlmm,head);
	TVertexDescriptor dummy;
	while(!atEnd(itout)){
		if(targetVertex(vlmm, getValue(itout)) != getNil<TVertexDescriptor>())
		{
			dummy = targetVertex(vlmm, getValue(itout));
			removeSubtree(vlmm,dummy);
		}
	// *(TOutEdgeIterator) 

	goNext(itout);
	}
	removeVertex(vlmm,head);

}

// returns 1 if the node will be kept

template<typename TAlphabet,typename TCargo,typename TVLMMSpec, typename TVertexDescriptor >
inline int
removeRedundantNodes( Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec > > > > &vlmm,
					 TVertexDescriptor &start,
					 TVLMMSpec &parameters)
{
//top down traversal to figure out which nodes can be deleted
typedef Graph<Automaton<TAlphabet,TCargo,WordGraph< VLMM < TVLMMSpec > > > > TVlmm;
	typedef typename Iterator<TVlmm, OutEdgeIterator>::Type TOutEdgeIterator;
	int sum = 0;

	TOutEdgeIterator itout(vlmm,start);
	TVertexDescriptor dummy;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	while(!atEnd(itout)){
		dummy = targetVertex(vlmm, getValue(itout));
		if(dummy != nilVal){
			sum += removeRedundantNodes(vlmm,dummy,parameters);

		}
		goNext(itout);
	}
	if(!isMarked(vlmm,start) && sum == 0)
	{
		removeVertex(vlmm,start);
		return 0;
	}
	if(isMarked(vlmm,start))
		smoothNode(vlmm,start,parameters);

	return 1;
}




template<typename TIndexType,typename TAlphabet,typename TCargo >
inline void
buildContextTree(Index<TIndexType, Index_ESA<> > & index,
		 Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ContextTree > > > > &vlmm,
			    unsigned threshold,
				float K,
				unsigned d,
				float alpha) 
{
	typedef Index<TIndexType, Index_ESA<> > TIndex;
	typedef Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ContextTree > > > > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	ContextTree parameters;
	setParameters(parameters,threshold,K,d,alpha);
	Iter< TIndex, VSTree< TopDown< ParentLinks<ConstrainedTraversal<Absolute> > > > > it(index,threshold,d);
	std::cout << "create the core Suffix Tree from the suffix array"<<endl;
	SEQAN_PROTIMESTART(create_sa_time);
	buildSuffixTreeFromIndex(it,vlmm);
    std::cout << "suffix array creation + sufix tree building took " << SEQAN_PROTIMEDIFF(create_sa_time) << " seconds" << endl;

	std::cout << "in initMaps:";
	initMaps(vlmm);
	std::cout << "Size of vlmm after suffix core:"<<numVertices(vlmm)<<std::endl;
	SEQAN_PROTIMESTART(addSuffixLinks);
    addSuffixLinks(vlmm);
	std::cout << "added suffix links and reverse suffix links: " <<SEQAN_PROTIMEDIFF(addSuffixLinks)<<" seconds"<<std::endl;
	//std::cout <<vlmm;
	SEQAN_PROTIMESTART(pruneTree);
	pruneTree(vlmm,parameters);
	std::cout << "Size of vlmm after prune Tree:"<<numVertices(vlmm)<< " Time: "<<SEQAN_PROTIMEDIFF(pruneTree)<<std::endl;
	std::cout << "pruned the ContextTree" <<std::endl;
	TVertexDescriptor root = getRoot(vlmm);
	SEQAN_PROTIMESTART(removeNodes);
	removeRedundantNodes(vlmm,root,parameters);
	std::cout << "Time removal nodes: "<<SEQAN_PROTIMEDIFF(removeNodes)<< " Number of nodes left:" <<numVertices(vlmm)<<endl;
	std::cout << "READY!" <<std::endl;
}

template<typename TIndexType,typename TAlphabet,typename TCargo >
inline void
buildBioPST(Index<TIndexType, Index_ESA<> > & index,
		 Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < BioPST > > > > &vlmm,
			    float threshold,
				unsigned minSupport,
				float minCondProb,
				unsigned d
				) 
{
	typedef Index<TIndexType, Index_ESA<> > TIndex;
	typedef Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < BioPST > > > > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	BioPST parameters;
	setParameters(parameters,threshold,minSupport,minCondProb,d);
	Iter< TIndex, VSTree< TopDown< ParentLinks<ConstrainedTraversal<Support> > > > > it(index,minSupport,d);
	std::cout << "create the core Suffix Tree from the suffix array"<<endl;
	buildSuffixTreeFromIndex(it,vlmm);
	std::cout << "in initMaps:";
	initMaps(vlmm);
	std::cout << "Size of vlmm after suffix core:"<<numVertices(vlmm)<<std::endl;
    addSuffixLinks(vlmm);
	std::cout << "added suffix links and reverse suffix links" <<std::endl;
	//std::cout <<vlmm;
	pruneTree(vlmm,parameters);
	std::cout << "Size of vlmm after prune Tree:"<<numVertices(vlmm)<<std::endl;
	std::cout << "pruned the BioPST" <<std::endl;
	TVertexDescriptor root = getRoot(vlmm);
	removeRedundantNodes(vlmm,root,parameters);
	std::cout <<"after remove redundant nodes: "<<numVertices(vlmm)<<" remain in the tree"<<endl;
	std::cout << "READY!" <<std::endl;
}

template<typename TIndexType,typename TAlphabet,typename TCargo >
inline void
buildPST(Index<TIndexType, Index_ESA<> > & index,
		 Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < PST > > > > &vlmm,
			    float threshold,
				float minEmpiricalProbability,
				float minConditionalProbability,
				float alpha,
				unsigned d) 
{
	typedef Index<TIndexType, Index_ESA<> > TIndex;
	typedef Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < PST > > > > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	if(minConditionalProbability > 1/(float)ValueSize<TAlphabet>::VALUE){
		cerr<<" minimum conditional probability is too large. It is not possible\n to distribute "<<minConditionalProbability<<
			" over "<<ValueSize<TAlphabet>::VALUE<< " entries "<<endl;
		std::cout<<" Can be at most:" <<1/ValueSize<TAlphabet>::VALUE<<endl;
		exit(1);
	}
	PST parameters;
	setParameters(parameters,threshold,minEmpiricalProbability,minConditionalProbability,alpha,d);

	Iter< TIndex, VSTree< TopDown< ParentLinks<ConstrainedTraversal<Relative> > > > > it(index,parameters.minEmpiricalProbability,d);
	inittheclock();
	
	indexRequire(index, ESA_SA());
	double buildESA = gettheruntime();
	
	indexRequire(index, ESA_LCP());
	double buildLCP = gettheruntime()-buildESA;		
	
	indexRequire(index, ESA_ChildTab());
	double buildChildTab = gettheruntime()-buildESA-buildLCP;
	
	std::cout << "ESA: " <<buildESA<< "  LCP: " << buildLCP << "  ChildTab: "<< buildChildTab<<std::endl;
	double build = gettheruntime();
	std::cout << "in buildSuffixTreeFromIndex:";
	buildSuffixTreeFromIndex(it,vlmm);
	std::cout << "constructed suffix tree core:" <<gettheruntime()-build<<std::endl;
	build = gettheruntime();
	std::cout << "in initMaps:";
	initMaps(vlmm);
	std::cout << " Size of vlmm after suffix core:"<<numVertices(vlmm)<<std::endl;
	std::cout << "init all maps" <<gettheruntime()-build<<std::endl;
	build = gettheruntime();
	addSuffixLinks(vlmm);
	std::cout << "added suffix links and reverse suffix links" <<gettheruntime()-build<<std::endl;
	std::cout <<vlmm;
	build = gettheruntime();
	std::cout <<"in Prune Tree";
	pruneTree(vlmm,parameters);
	std::cout << " Size of vlmm after prune Tree:"<<numVertices(vlmm)<<std::endl;
	std::cout << "pruned the PST" <<gettheruntime()-build<<std::endl;
	TVertexDescriptor root = getRoot(vlmm);
	removeRedundantNodes(vlmm,root,parameters);
	std::cout <<"after remove redundant nodes: "<<numVertices(vlmm)<<" remain in the tree"<<endl;
	std::cout << "READY!" <<std::endl;
}

/*******

Training

********/

/**
*  Likelihood Estimation : works such that the reverse suffix links are walked starting from the root
*  whenever a node is not marked(or a leaf) is reached the walking down is finished and the deepest possible 
*  context for the estimation is identified
*/



template<typename TAlphabet,typename TCargo,typename TVLMMSpec>
inline float
estimateLikelihood( Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec > > > > &vlmm,
					String<TAlphabet> &text )
{
	float result = 0;
	typedef typename Iterator<String<TAlphabet> >::Type TIter;
	TIter it = begin(text);

	for(goBegin(it);!atEnd(it);goNext(it))
	{
		
		result += log(getProbabilityForLongestContext(vlmm,it));
		//std::cout <<" prob for letter: "<<value(it)<< " is: "<<getProbabilityForLongestContext(vlmm,it)<<endl;
	}
	return result;
}

// estimates the highest scoring window of size windowSize
template<typename TAlphabet,typename TCargo,typename TVLMMSpec>
inline float
estimateLikelihoodWindow( Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec > > > > &vlmm,
					String<TAlphabet> &text,
					unsigned windowSize)
{
	if(length(text)<windowSize){
		cerr << "Text smaller than window size. In estimateLikelihoodWindow\n";
		return (float)0;
	}
	float best,result = 0;
	typedef typename Iterator<String<TAlphabet> >::Type TIter;
	TIter windowEnd = begin(text),windowStart = begin(text);
	goBegin(windowEnd);
	goBegin(windowStart);
	
	for(unsigned i = 0;i<windowSize;++i,goNext(windowEnd))
			result += log(getProbabilityForLongestContext(vlmm,windowEnd));
			
	best = result;
	std::cout <<"window score:" <<result<<endl;
	for(;!atEnd(windowEnd);goNext(windowEnd),goNext(windowStart))
	{
		result += -log(getProbabilityForLongestContext(vlmm,windowStart));
		result += log(getProbabilityForLongestContext(vlmm,windowEnd));
		//std::cout <<"window score:" <<result<<endl;
		if(result > best)
			best = result;
	}
	return best;
}


// estimates the highest scoring window of size windowSize and the likelihood on the whole sequence
//NOTE if the window size is greater than the sequence length than the bestWindow scorwe is the same as the score
//of the whole sequence
template<typename TAlphabet,typename TCargo,typename TVLMMSpec>
inline float
estimateLikelihoodWindowAndWhole( Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec > > > > &vlmm,
					String<TAlphabet> &text,
					unsigned windowSize,
					float &bestWindow)
{

	float wholeSequenceScore=0,result = 0;
	typedef typename Iterator<String<TAlphabet> >::Type TIter;
	TIter windowEnd = begin(text),windowStart = begin(text);
	goBegin(windowEnd);
	goBegin(windowStart);
	float dummy;
	for(unsigned i = 0;i<windowSize && !atEnd(windowEnd);++i,goNext(windowEnd)){
			wholeSequenceScore += log(getProbabilityForLongestContext(vlmm,windowEnd));
		}

	result = wholeSequenceScore;
	bestWindow = result;
	
	for(;!atEnd(windowEnd);goNext(windowEnd),goNext(windowStart))
	{
		dummy = log(getProbabilityForLongestContext(vlmm,windowEnd));
		result += -log(getProbabilityForLongestContext(vlmm,windowStart));
		result += dummy;
		wholeSequenceScore += dummy;
		if(result > bestWindow)
			bestWindow = result;
	}
	std::cout <<" about to return, wholeSequenceScore: "<<wholeSequenceScore<<" and bestWindow: "<<bestWindow<<endl;
	return wholeSequenceScore;
}

// estimates the likelihood for every sequences in the input file and save it in the outout file
template<typename TAlphabet,typename TCargo,typename TVLMMSpec, typename TFile>
inline void estimateLikelihoodOnFile(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec > > > > &vlmm,
									 String<char> &sequenceFile,
									 TFile &outFile,
									 unsigned windowSize)
{
		String<String<char> > ids;
		String<String<TAlphabet> > sequences;
		createInputString(sequenceFile,sequences,ids);
		//for(unsigned i=0;i<length(sequences);++i)
			//std::cout <<" first sequence:"<<sequences[0]<< " and the id: "<<ids[0]<<endl;
			//std::cout <<" second sequence:"<<sequences[1]<< " and the id: "<<ids[1]<<endl;
		typedef typename Iterator<String<String<TAlphabet> > >::Type TIter;
		typedef typename Iterator<String<String<char> > >::Type TIterChar;
		TIter it = begin(sequences);
		TIterChar id = begin(ids);
		//put first line for file
		_streamWrite(outFile,"Id\tSequenceLength\tLikelihood(Sequence)\tNormalizedLikelihood\tLikelihood of BestWindow of size ");
		_streamPutInt(outFile,windowSize);
		_streamWrite(outFile,")\n");
		//std::cout <<" value(it): "<<value(it)<<endl;
		for(;!atEnd(it);goNext(it),goNext(id))
		{
			float bestWindow=0,wholeSequenceScore=0;
			wholeSequenceScore =estimateLikelihoodWindowAndWhole(vlmm,value(it),windowSize,bestWindow);
			//std::cout <<"wSS" << wholeSequenceScore<<" bestWindow "<<bestWindow<<endl;
			saveSequenceLikelihood(outFile,windowSize,length(value(it)),value(id),wholeSequenceScore,bestWindow);
		}
}

template<typename TAlphabet,typename TCargo,typename TVLMMSpec,typename TIter>
inline float
getProbabilityForLongestContext( Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec > > > > &vlmm,
					TIter &it )
{
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TVLMMSpec> > > > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TIter copy = it;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	TVertexDescriptor node = getRoot(vlmm);


	while(!atBegin(copy) ){
		goPrevious(copy);
		if(getReverseSuffixLink(vlmm,node,value(copy)) != nilVal)
		{
				node = getReverseSuffixLink(vlmm,node,value(copy));
				if(!isMarked(vlmm,node)){
						node = getSuffixLink(vlmm,node);
						break;
				}
		}
		else
			break;
		
	}
	return getProbability(vlmm,node,value(it));
}



/*************

save the graph
**************/

template<typename TFile, typename TCargo >
inline void
writeHead(Graph<Automaton<AminoAcid, TCargo , WordGraph < VLMM < ContextTree > > > > &vlmm,
	   TFile & target)
{	
		_streamWrite(target,"VLMM\tContextTree\tAminoAcid\t");
		_streamPutInt(target,length(vlmm.data_marked));
}

template<typename TFile, typename TCargo >
inline void
writeHead(Graph<Automaton<Dna, TCargo , WordGraph < VLMM < ContextTree > > > > &vlmm,
	   TFile & target)
{	
		_streamWrite(target,"VLMM\tContextTree\tDna\t");
		_streamPutInt(target,length(vlmm.data_marked));
}

template<typename TFile,  typename TCargo >
inline void
writeHead(Graph<Automaton<AminoAcid, TCargo , WordGraph < VLMM < BioPST > > > > &vlmm,
	   TFile & target)
{	
	_streamWrite(target,"VLMM\tBio-PST\tAminoAcid\t");
	_streamPutInt(target,length(vlmm.data_marked));
}

template<typename TFile,  typename TCargo >
inline void
writeHead(Graph<Automaton<Dna, TCargo , WordGraph < VLMM < BioPST > > > > &vlmm,
	   TFile & target)
{	
	_streamWrite(target,"VLMM\tBio-PST\tDna\t");
	_streamPutInt(target,length(vlmm.data_marked));
		
}

template <typename TStream>
inline void
_streamPutFloat(TStream & target,
			  float number, 
			  char const * format_string)
{
SEQAN_CHECKPOINT
	char str[32];
	sprintf(str, format_string, number);
	_streamWrite(target, str);
}
template <typename TStream>
inline void
_streamPutFloat(TStream & target,
			  float number)
{
SEQAN_CHECKPOINT
	_streamPutFloat(target, number, "%f");
}


//save/export the vlmm in a file for reading it again using the import funcion
template<typename TFile, typename TAlphabet, typename TCargo, typename TVLMMSpec ,typename TVertexDescriptor>
inline void
saveNode(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec > > > > &vlmm,
		   TVertexDescriptor &node,
			TFile & target)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TVLMMSpec> > > > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Size<TAlphabet>::Type TSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	// get to the next line in the file
	_streamPut(target,'\r');
	_streamPutInt(target, (unsigned)node);
	_streamPut(target, '\t');
	_streamPutInt(target, (unsigned)getFather(vlmm,node));
	_streamPut(target, '\t');
	_streamPutInt(target, (unsigned)getSuffixLink(vlmm,node));
	_streamPut(target, '\t');
	_streamPutInt(target, (unsigned)isMarked(vlmm,node));
	_streamPut(target, '\t');
//Childrens with label
	
	for(TSize i = 0;i<table_length;++i) {
			TVertexDescriptor child = vlmm.data_vertex[node].data_edge[i].data_target;
			_streamPutInt(target,(unsigned )child);
			_streamPut(target, '\t');
			if(child != nilVal )
			{
				String<TAlphabet> edgeString;
				getChildLabel(vlmm,node,child,edgeString);
				_streamWrite(target, edgeString);
				_streamPut(target, '\t');
			}
	
	}
	// Probability vector
		for(TSize i = 0;i<table_length;++i){
				_streamPutFloat(target,getProbability(vlmm,node,i));
				_streamPut(target,'\t');
		}

//ReverseSuffixLinks
	
		for(TSize i = 0;i< (table_length-1) ;++i){
				_streamPutInt(target,getReverseSuffixLink(vlmm,node,i));
				_streamPut(target,'\t');
		}
		// last entry should not end with a '\t' instead with a  '\n'
		_streamPutInt(target,getReverseSuffixLink(vlmm,node,table_length-1));
		
}


/* Export plain file format for
VLMM\t	Type[ContextTree,PST]\t	Alphabet[Dna,Dna5,AminoAcids]\t	Params[gamma:threshold:pMin ... etc]
Node\t	Father\t	SuffixLink\t	Marked\t	Children[1 .. N]\t	ChildLabel[1 .. N]\t		ProbabilityVector[1 .. N]\t	ReverseSuffixLink[1 .. N]	

	



//save/export the vlmm in a file for reading it again using the import funcion

*/


// save function where a filehandle has to be created before 
template<typename TFile, typename TAlphabet, typename TCargo, typename TVLMMSpec>
inline void
save(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec > > > > &vlmm,
	   TFile & target)
{

SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TVLMMSpec> > > > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > >::Type TIterConst;
	std::cout<<"Save final vlmm"<<endl;
	writeHead(vlmm,target);
	for(TIterConst it = begin(vlmm.data_vertex);!atEnd(it);goNext(it)) {
		if (!idInUse(vlmm.data_id_managerV, position(it))) continue;
		TVertexDescriptor dummy = position(it); 
		saveNode(vlmm,dummy,target);
	}



}
// save function where a filehandle is created for the filename
template<typename TFile, typename TAlphabet, typename TCargo, typename TVLMMSpec>
inline void
save(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec > > > > &vlmm,
	   String<char> & filename)
{

SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TVLMMSpec> > > > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > >::Type TIterConst;
	fstream target;
	target.open(toCString(filename), ios_base::out | ios_base::trunc);
	if (!target.is_open()) {
			cerr << "Import of sequence " << filename << " failed." << endl;
			exit(1);
	}
	writeHead(vlmm,target);
	for(TIterConst it = begin(vlmm.data_vertex);!atEnd(it);goNext(it)) {
		if (!idInUse(vlmm.data_id_managerV, position(it))) continue;
		TVertexDescriptor dummy = position(it); 
		saveNode(vlmm,dummy,target);
	}
	target.close();

}


// save function where a filehandle is created for the filename
template<typename TFile>
inline void saveSequenceLikelihood(TFile &outFile,
					   unsigned windowSize,
					   unsigned SequenceLength,
					   String<char> & id,
					   float wholeSequenceScore,
					   float bestWindow)
{
SEQAN_CHECKPOINT
	_streamWrite(outFile,id);
	_streamPut(outFile,'\t');
	_streamPutInt(outFile,SequenceLength);
	_streamPut(outFile,'\t');
	_streamPutFloat(outFile,wholeSequenceScore);
	_streamPut(outFile,'\t');
	_streamPutFloat(outFile,(wholeSequenceScore/(float)SequenceLength));
	_streamPut(outFile,'\t');
	_streamPutFloat(outFile,bestWindow);
	_streamPut(outFile,'\n');



}


/*************

read the graph
**************/

/* this function can be called to scan every entry of
	the vlmm output file
*/
template <typename TFile,typename TAlphabet>
inline bool
_scanNextEntry(TFile & file,
			   String<TAlphabet> &entry)
{
SEQAN_CHECKPOINT
	resize(entry,0);
	//SEQAN_ASSERT(!_streamEOF(file))
	if(_streamEOF(file))
		return false;
	
	while (true)
	{
		typename Value<TFile>::Type c = _streamGet(file);

		if (_streamEOF(file))
		{
			if(length(entry) == 0)
				return false;
			return true;
		}

		if ((c == '\t') ||  (c == '\n') || (c == '\r'))
		{
			if(length(entry) == 0)
				return false;
			return true;
		}

		
			appendValue(entry,c,Generous());
	}
}


template <typename TFile>
inline bool
_removeTailingNewlines(TFile & file)
{
SEQAN_CHECKPOINT
	
	if(_streamEOF(file))
		return false;
	bool anotherRound = true;
	while (anotherRound)
	{
		typename Value<TFile>::Type c = _streamGet(file);

		if (_streamEOF(file))
		{
			return false;
		}
		bool anotherRound = false;
		if (c == '\n')
		{
			anotherRound = true;
		}
		

	}
}

template <typename TFile>
inline int
_scanNextIntEntry(TFile & file)
{
	String<char> entry;
	reserve(entry,40);
	_scanNextEntry(file,entry);

	return atoi((char*)toCString(entry));

}

template <typename TFile>
inline float
_scanNextFloatEntry(TFile & file)
{
	String<char> entry;
	reserve(entry,40);
	_scanNextEntry(file,entry);
	return atof((char*)toCString(entry));

}


template<typename TFile, typename TAlphabet, typename TCargo>
inline void
readGraph(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ContextTree > > > >&vlmm,
	   TFile & file)
{
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<ContextTree> > > > TGraph;
	typedef typename Size<TAlphabet>::Type TSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdge;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > >::Type TIterConst;

	unsigned size = (unsigned)_scanNextIntEntry(file);
	std::cout <<"size= "<<size<<endl;
	//init the vlmm
	initGraph(vlmm,size);
	//size = _scanNextIntEntry(file);  // remove the tailing newline character
	//_removeTailingNewlines(file);
	// go for every node
	
	do{
		//std::cout <<"go for next node:"<<endl;
		TVertexDescriptor node = (unsigned)_scanNextIntEntry(file);
		assignValue(vlmm.data_id_managerV.data_in_use, node, true);
		TVertexDescriptor father = (unsigned)_scanNextIntEntry(file);
		if(node != getRoot(vlmm))
				setFather(vlmm,father,node);
		TVertexDescriptor suffixLink = (unsigned)_scanNextIntEntry(file);
		if(suffixLink != nilVal)
				setSuffixLink(vlmm,node,suffixLink);
		int yes = _scanNextIntEntry(file);
		if(yes)
			setMarked(vlmm,node,true);
		//std::cout <<"node:" << node<< " father:  "<< father<< " suffixLinkTarget: "<<suffixLink<< " Marked?: "<<yes<<endl;

		for(TSize i =0;i<table_length;++i)
		{
			TVertexDescriptor child = (unsigned)_scanNextIntEntry(file);
			if( child  != nilVal)
			{	
				//set child true for function adEdge
				assignValue(vlmm.data_id_managerV.data_in_use, child, true);
				String<TAlphabet> edgeString;
				_scanNextEntry(file,edgeString);
				//std::cout << "edge string to node:"<<child<< " is:" <<edgeString<<endl;
				addEdge(vlmm,node,child,edgeString);
			}


		}
		
		for(TSize i =0;i<table_length;++i)
		{
			float set = _scanNextFloatEntry(file);
			//std::cout << "float: "<<set;
			setProbability(vlmm,node,i,set);
		}
		//cout <<endl;
		for(TSize i =0;i<table_length;++i)
		{
			TVertexDescriptor target = (unsigned)_scanNextIntEntry(file);
			if( target  != nilVal)
			{	
				//set child true for function addEdge
				setReverseSuffixLink(vlmm,node,target,i);
				//cout << "reverseSL to node:"<<target<< " is:" <<(TAlphabet)i<<endl;
			}
		}
		//_scanNextIntEntry(file);
	} while(!_streamEOF(file));
    
	for(TIterConst it = begin(vlmm.data_vertex);!atEnd(it);goNext(it)) {
		if (!idInUse(vlmm.data_id_managerV, position(it))){
			//std::cout << "id:" <<position(it)<<" is not used\n";
			appendValue(vlmm.data_id_managerV.data_freeIds, position(it));
		}
		else{
			//std::cout << "id:" <<position(it)<<" is used\n";
		}

	}
	//std::cout<<"laenge: "<<length(vlmm.data_id_managerV.data_freeIds);

}

template<typename TAlphabet>
void createInputString(String<char>  &filename,
					   String<String<TAlphabet> > &sequences,
					   String<String<char> > &ids){

	//Read in sequences

	ifstream file2,file;
	file2.open(toCString(filename), ios_base::in | ios_base::binary );
	if (!file2.is_open()) {
				cerr << "Import of sequence " << filename << " failed." << endl;
				exit(1);
		}
	int count = 0;
	while(!_streamEOF(file2))
	{
		goNext(file2,Fasta());
		++count;
	}
	file2.close();
	std::cout<<"There are "<<count<<" seqs in the file "<<filename<<endl;
	file.open(toCString(filename), ios_base::in | ios_base::binary );
	resize(ids,count);
	resize(sequences,count);
	for(int i = 0;i<count;++i){
		readID(file, ids[i], Fasta());
		read(file, sequences[i],  Fasta());


	}
	file.close();
}

template<typename TFile>
inline void
read( TFile & file)
{
	
	
	String<char> entry;
	_scanNextEntry(file,entry);
	_scanNextEntry(file,entry);
	_scanNextEntry(file,entry);
	if(entry == "Dna"){
		Graph<Automaton<Dna, String<Dna> , WordGraph < VLMM < ContextTree > > > > vlmmDna;
		readGraph(vlmmDna,file);
		std::cout<< vlmmDna;
	}
	if(entry == "AminoAcid"){
		Graph<Automaton<AminoAcid, String<AminoAcid> , WordGraph < VLMM < ContextTree > > > > vlmmProtein;
		readGraph(vlmmProtein,file);

	}
	else{
		cerr<<"Alphabet in file "<<file<<" not supported for vlmm.Maybe incorrect input file format.\n";
		exit(1);
	}
	
	//typedef Graph<Automaton<TAlphabet, String<TAlphabet> , WordGraph < VLMM < ContextTree > > > TVLMM;
	//TVLMM vlmm;
	//readGraph(vlmm,file);

}

template<typename TInputFile,typename TFile>
inline void
readForLikelihoodEstimate(TInputFile			& file,
						  String<char>	& sequenceFile,
						  TFile			& outFile,
						  unsigned		windowSize)
{
	
	
	String<char> entry;
	_scanNextEntry(file,entry);
	_scanNextEntry(file,entry);
	_scanNextEntry(file,entry);
	if(entry == "Dna"){
		Graph<Automaton<Dna, String<Dna> , WordGraph < VLMM < ContextTree > > > > vlmmDna;
		readGraph(vlmmDna,file);
		estimateLikelihoodOnFile(vlmmDna,sequenceFile,outFile,windowSize);

	}
	if(entry == "AminoAcid"){
		Graph<Automaton<AminoAcid, String<AminoAcid> , WordGraph < VLMM < ContextTree > > > > vlmmProtein;
		readGraph(vlmmProtein,file);
		estimateLikelihoodOnFile(vlmmProtein,sequenceFile,outFile,windowSize);
	}
	else{
		cerr<<"Alphabet in file "<<file<<" not supported for vlmm.Maybe incorrect input file format.\n";
		exit(1);
	}


}




//write the vlmm to a file in Easy-readable format,also used by the  << Operator
template<typename TFile, typename TAlphabet, typename TCargo, typename TSpec , typename TIDString>
inline void
write(TFile & target,
	  Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > >const & g,
	  TIDString const &,
	  Raw)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TSpec> > > > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Size<TAlphabet>::Type TSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();

	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const>::Type TIterConst;
	_streamWrite(target,"VLMM - Probability Vector of each node:\n");
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		if (!idInUse(g.data_id_managerV, position(it))) continue;
		TVertexDescriptor sourceVertex = position(it);
			_streamPutInt(target, sourceVertex);
			_streamPut(target, ' ');
			_streamPut(target, ' ');
			_streamWrite(target, "SLink: ");
			std::cout << g.data_suffix_link[sourceVertex];
			_streamWrite(target, " Father: ");
			std::cout << g.data_father[sourceVertex];
			_streamWrite(target,"  (");

			// wollte Wahrscheinlichleiten ausgeben
		for(int i = 0;i<ValueSize<TAlphabet>::VALUE;++i){
			std::cout <<value(g.data_probability_vector[sourceVertex],i);
				//_streamPut(target,(float)value(g.data_probability_vector[sourceVertex],i));
				if(i < (ValueSize<TAlphabet>::VALUE -1))
						_streamPut(target, ' ');	
		}
		_streamPut(target, ')');
		if(g.data_marked[getRoot(g)])
			if(! g.data_marked[sourceVertex])
				_streamWrite(target, " unmarked");
		_streamWrite(target, "\n");
	}
		_streamWrite(target,"VLMM - Directed:\n");

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
			
			String<TAlphabet> edgeString;
			TVertexDescriptor child = g.data_vertex[sourceVertex].data_edge[i].data_target;
			getSuffixChildLabel(g,sourceVertex,child,edgeString);
			_streamWrite(target, edgeString);
			//_streamWrite(target, getProperty(g.data_edge_label, TEdgeDescriptor(sourceVertex, TAlphabet(i))));
		
			_streamPut(target, '\n');
		}
	}
}

//write the vlmm to a file in Dot-Format
template <typename TFile,typename TAlphabet,typename TCargo , typename TSpec, typename TNodeMap, typename TEdgeMap>
void write(TFile & file, 
	   Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > > > const& g,
	   TNodeMap const& nodeMap,
	   TEdgeMap const& edgeMap,
	   DotDrawing) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > > > TGraph;
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
	typedef typename Iterator<TGraph, VertexIterator >::Type TConstIter;
	TConstIter it(g);
	for(;!atEnd(it);goNext(it)) {
		_streamWrite(file, getProperty(nodeMap, *it));
		/* wollte Wahrscheinlichleiten ausgeben
		_streamWrite(file, "|");
		for(int i = 0;i<ValueSize<TAlphabet>::VALUE;++i){
			//std::cout <<"write:"<<value(g.data_probability_vector[*it],i)<<std::endl;
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
	typedef typename Iterator<TGraph, EdgeIterator >::Type TConstEdIter;
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

		TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
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

template <typename TFile,typename TAlphabet,typename TCargo, typename TSpec>
void write(TFile & file, 
	   Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > > > const& g, 
	   DotDrawing) 
{
	SEQAN_CHECKPOINT
	String<String<char> > nodeMap;
	_createNodeAttributes(g,nodeMap);
	String<String<char> > edgeMap;
	_createEdgeAttributes(g,edgeMap);
	write(file,g,nodeMap,edgeMap,DotDrawing());
}

//write the vlmm not as a suffix tree but reveersed edge labels to a file in Dot-Format
template <typename TFile,typename TAlphabet,typename TCargo, typename TSpec, typename TNodeMap, typename TEdgeMap>
void writeAsVLMM(TFile & file, 
	   Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > > > const& g,
	   TNodeMap const& nodeMap,
	   TEdgeMap const& edgeMap,
	   DotDrawing) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > > > TGraph;
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
	typedef typename Iterator<TGraph, VertexIterator >::Type TConstIter;
	TConstIter it(g);
	for(;!atEnd(it);++it) {
		_streamWrite(file, getProperty(nodeMap, *it));
		_streamPut(file, ';');
		_streamPut(file, '\n');
	}
	_streamPut(file, '\n');

	_streamWrite(file, "/* Edges */\n");
	typedef typename Iterator<TGraph, EdgeIterator >::Type TConstEdIter;
	TConstEdIter itEd(g);
		//NEW Reverse Suffix Links
	goBegin(it);
	for(;!atEnd(it);++it) {

		TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
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

template <typename TFile,typename TAlphabet,typename TCargo, typename TSpec>
void writeAsVLMM(TFile & file, 
	   Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TSpec > > > > const& g, 
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


