#ifndef SEQAN_HEADER_VLMM_H
#define SEQAN_HEADER_VLMM_H

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
		TBase(__index),AbsoluteThreshold(threshold),RelativeThreshold(0.0),Down(false),Up(0),MaxDepth(maxdepth){}

	// use constructor for a relative threshold
	Iter(TIndex &__index,float threshold,unsigned maxdepth):
		TBase(__index),AbsoluteThreshold(0),RelativeThreshold(threshold),Down(false),Up(0),MaxDepth(maxdepth){}

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



// go down the leftmost edge
template < typename TIndex, class TSpec >
inline bool goDown(Iter< TIndex, VSTree< TopDown<ParentLinks<ConstrainedTraversal<TSpec > > > > > &it) {
	if (isLeaf(it) || (length(representative(it)) > it.MaxDepth) ) return false;
	_historyPush(it, value(it).i1);

		typename Size<TIndex>::Type i = _getUp(value(it).i1.i2, container(it));
	if (value(it).i1.i1 < i && i < value(it).i1.i2)
		value(it).i1.i2 = i;
	else
		value(it).i1.i2 = _getDown(value(it).i1.i1, container(it));

	// set the the flag for going down
	it.Down = true;
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
			
			if(!isLeaf(it) && isRightTerminal(it)){
				// this must be possible
				
				goDown(it);
				goRight(it);
				while( (repLength(container(it), nodeUp(it)) == repLength(it)) && goRight(it))
					
				// if the current node remains to be right terminal
				// than we have found a node where multiple strings end
				// goDown was wrong and we have to go up again
				if(isLeaf(it) && (repLength(container(it), nodeUp(it)) == repLength(it))){
					goUp(it);
					it.Up += -1;
					//cout << value(it) << " = " << length(representative(it)) << " " << representative(it) << "  toFather:"<<parentEdgeLabel(it)<<"  hits: "<<length(getOccurences(it))<<endl;	
					if (!goRight(it))
						while (goUp(it) && !goRight(it));
					}
			}
			else
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
			if(repLength(it) == it.MaxDepth && length(parentEdgeLabel(it)) == 1){
				//cout << "In Abs Clausel  "<<value(it) << " = " << length(representative(it)) << " " << representative(it) << "  toFather:"<<parentEdgeLabel(it)<<"  hits: "<<length(getOccurences(it))<<endl;
				if(isLeaf(it)){
						walk_down=0;
						not_finished =1;
				}
				else
					if(isRightTerminal(it)){ // check if node can be extended
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
		if((getFrequency(it)) >= it.AbsoluteThreshold) 
			walk_down = 1;
	}

	it.Down = false;
	it.Up = 0;
	while(not_finished)
	{
		// if threshold is reached walk down else go always right or up
		if(walk_down){
			// for avoiding the implicit $-edges
			
			if(!isLeaf(it) && isRightTerminal(it)){
				// this must be possible
				
				goDown(it);
				goRight(it);
				while( (repLength(container(it), nodeUp(it)) == repLength(it)) && goRight(it))
					
				// if the current node remains to be right terminal
				// than we have found a node where multiple strings end
				// goDown was wrong and we have to go up again
				if(isLeaf(it) && (repLength(container(it), nodeUp(it)) == repLength(it))){
					goUp(it);
					it.Up += -1;
					//cout << value(it) << " = " << length(representative(it)) << " " << representative(it) << "  toFather:"<<parentEdgeLabel(it)<<"  hits: "<<length(getOccurences(it))<<endl;	
					if (!goRight(it))
						while (goUp(it) && !goRight(it));
					}
			}
			else
				if(!goDown(it) && !goRight(it))
					while (goUp(it) && !goRight(it));
		}
		else
			if (!goRight(it))
				while (goUp(it) && !goRight(it));
		
		if( getFrequency(it) >= it.AbsoluteThreshold){
			not_finished = 0;
			// we have to check if the current node can be extended with at least one letter if the following
			// condition is true:
			if(repLength(it) == it.MaxDepth && length(parentEdgeLabel(it)) == 1){
				//cout << "In Abs Clausel  "<<value(it) << " = " << length(representative(it)) << " " << representative(it) << "  toFather:"<<parentEdgeLabel(it)<<"  hits: "<<length(getOccurences(it))<<endl;
				if(isLeaf(it)){
						walk_down=0;
						not_finished =1;
				}
				else
					if(isRightTerminal(it)){ // check if node can be extended
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
				cout << value(it) << " = " << length(representative(it)) << " " << representative(it) << "  toFather:"<<parentEdgeLabel(it)<<"  hits: "<<length(getOccurences(it))<<endl;
				goRight(it);
				cout << value(it) << " = " << length(representative(it)) << " " << representative(it) << "  toFather:"<<parentEdgeLabel(it)<<"  hits: "<<length(getOccurences(it))<<endl;
				while( (repLength(container(it), nodeUp(it)) == repLength(it)) && goRight(it)){
					
				cout << value(it) << " = " << length(representative(it)) << " " << representative(it) << "  toFather:"<<parentEdgeLabel(it)<<"  hits: "<<length(getOccurences(it))<<endl;	
				}
					
				// if the current node remains to be right terminal
				// than we have found a node where multiple strings end
				// goDown was wrong and we have to go up again
				if(isLeaf(it) && (repLength(container(it), nodeUp(it)) == repLength(it))){
					goUp(it);
					cout << value(it) << " = " << length(representative(it)) << " " << representative(it) << "  toFather:"<<parentEdgeLabel(it)<<"  hits: "<<length(getOccurences(it))<<endl;	
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
	unsigned threshold;	 // minimum support
	float K;    //the name originates from the paper of Bühlmann (1999) 
	unsigned MaxDepth;
	float alpha;	// pseudocount for smoothing
	
};

struct PST{
	float threshold ;				   // r
	float minEmpiricalProbability ;    // P min
	float minConditionalProbability ;  // y min
	float alpha ;					   // alpha :-)
	unsigned MaxDepth;				   // d
		
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

void setParameters(BioPST & params,unsigned t,float K, unsigned d, float alpha){
	params.threshold = t;
	params.K = K;
	params.MaxDepth = d;
	params.alpha = alpha;
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
	int i;
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
	goDown(childs);
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
					cout << "prefix of edgelabel "<<pref<<endl;
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

		if(it.Down){
			father = child;
		}
		while(it.Up>0){
			father = getFather(target,father);
			--it.Up;
		}
		

		
		child = addIncompleteVertex(target);
		// check if iterator
		std::cout << "Node:"<<child<<"  "<<value(it) << " = " << repLength(it)<< " " << representative(it) << "  toFather:"<<parentEdgeLabel(it)<<"  hits: "<<length(getOccurences(it))<<std::endl;
		
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
	cout << "split edge at pos:"<< splitPosition<<" character:" <<(int)childCharacter<<" father" <<father<<endl;
	String<TAlphabet> edgeString;
	float number = getProbability(vlmm,father,childCharacter);
	getSuffixChildLabel(vlmm,father,letter,edgeString);
	//String<TAlphabet> edgeString = getProperty(vlmm.data_edge_label, TEdgeDescriptor(father, letter));
	//std::cout<<"before new node created, NumVertc:"<<numVertices(vlmm)<<std::endl;
	TVertexDescriptor newNode = addAdditionalVertex(vlmm);
	//std::cout<<"after node created, NumVertc:"<<numVertices(vlmm)<<std::endl;
	String<TAlphabet> newEdgeString = childCharacter;
		//cout << "eedestring:"<<edgeString<<" length: "<<length(edgeString)<<endl;
	SEQAN_ASSERT(splitPosition < length(edgeString))
	if(splitPosition >= length(edgeString)){
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
			std::cout << "go for SL of father:" << father<< " String:"<< walkDown<<std::endl;
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
			   ContextTree & parameters) 
{
	
		TAlphabet letter(childCharacter);
		unsigned size = ValueSize<TAlphabet>::VALUE;
	// the pruning criterion mentioned by Bühlmann & Mächler, taken from Rissanen 1983
		float fathersum = 0;
    
	for(unsigned int pos = 0;pos< size ;++pos){
		fathersum = fathersum + getProbability(vlmm,father,pos);
	}
	SEQAN_ASSERT(fathersum != 0)
	//std::cout<<"extendNode on node "<<node<<std::endl;
	//cout << " diff= "<<countChar * log(1/(getProbability(vlmm,father,childCharacter)/fathersum))<<endl;
		if ( (countChar * log(1/(getProbability(vlmm,father,childCharacter)/fathersum)))   <= parameters.K )
			return false;

	return true;
}

/* 
	In general this function can only be called for nodes which still contain the pattern counts and
	have not be converted into probabilities yet.
*/
template<typename TCargo,typename TAlphabet ,typename TVertexDescriptor>
inline bool
pruneNode(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ContextTree > > > > &vlmm,
			   TVertexDescriptor &son,
			   ContextTree & parameters) 
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

	for(int pos = 0;pos< ValueSize<TAlphabet>::VALUE;++pos){
		if(getProbability(vlmm,father,pos) > 0 && getProbability(vlmm,son,pos) > 0)
			difference += getProbability(vlmm,son,pos) * (getProbability(vlmm,son,pos)/sonsum)*
						  log((getProbability(vlmm,son,pos)/sonsum)/(getProbability(vlmm,father,pos)/fathersum));
		
	}
	if(difference <= parameters.K)
		return true;

 return false;
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
		//cout<<" at node:"<<node<<endl;
		//turnNodeCountsIntoProbability(vlmm,node);
	}
	
	// root cannot be pruned
	setMarked(vlmm,root,true);
	//assignProperty(marked,root,true);
	TVertexDescriptor dummy = 0;
	std::cout<<" start pruning from the root"<<std::endl;
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
	for(int i=0;i<table_length;i++)
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
	cout << "create the core Suffix Tree from the suffix array"<<endl;
	buildSuffixTreeFromIndex(it,vlmm);

	std::cout << "in initMaps:";
	initMaps(vlmm);
	std::cout << " Size of vlmm after suffix core:"<<numVertices(vlmm)<<std::endl;
    addSuffixLinks(vlmm);
	std::cout << "added suffix links and reverse suffix links" <<std::endl;
	std::cout <<vlmm;
	pruneTree(vlmm,parameters);
	std::cout << " Size of vlmm after prune Tree:"<<numVertices(vlmm)<<std::endl;
	std::cout << "pruned the ContextTree" <<std::endl;
	TVertexDescriptor root = getRoot(vlmm);
	removeRedundantNodes(vlmm,root,parameters);
	std::cout <<" remove redundant nodes: "<<endl<< vlmm;
	std::cout << "READY!" <<std::endl;
}

template<typename TIndexType,typename TAlphabet,typename TCargo >
inline void
buildBioPST(Index<TIndexType, Index_ESA<> > & index,
		 Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < BioPST > > > > &vlmm,
			    unsigned threshold,
				float K,
				unsigned d,
				float alpha) 
{
	typedef Index<TIndexType, Index_ESA<> > TIndex;
	typedef Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < BioPST > > > > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	BioPST parameters;
	setParameters(parameters,threshold,K,d,alpha);
	Iter< TIndex, VSTree< TopDown< ParentLinks<ConstrainedTraversal<Support> > > > > it(index,threshold,d);
	cout << "create the core Suffix Tree from the suffix array"<<endl;
	buildSuffixTreeFromIndex(it,vlmm);
	std::cout << "in initMaps:";
	initMaps(vlmm);
	std::cout << " Size of vlmm after suffix core:"<<numVertices(vlmm)<<std::endl;
    addSuffixLinks(vlmm);
	std::cout << "added suffix links and reverse suffix links" <<std::endl;
	std::cout <<vlmm;
	pruneTree(vlmm,parameters);
	std::cout << " Size of vlmm after prune Tree:"<<numVertices(vlmm)<<std::endl;
	std::cout << "pruned the BioPST" <<std::endl;
	TVertexDescriptor root = getRoot(vlmm);
	removeRedundantNodes(vlmm,root,parameters);
	std::cout <<" remove redundant nodes: "<<endl<< vlmm;
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
	PST parameters;
	setParameters(parameters,threshold,minEmpiricalProbability,minConditionalProbability,alpha,d);

	Iter< TIndex, VSTree< TopDown< ParentLinks<ConstrainedTraversal<Relative> > > > > it(index,parameters.minEmpiricalProbability,d);
	//Iter< TIndex, VSTree< TopDown< ParentLinks<ConstrainedTraversal<Absolute> > > > > it(index,2,d);
	
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
	//std::cout << vlmm;
	std::cout << "READY!" <<std::endl;
}

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
	Iterator<String<TAlphabet> >::Type it = begin(text);

	for(goBegin(it);!atEnd(it);goNext(it))
	{
		
		result += log(getProbabilityForLongestContext(vlmm,it));
		//cout <<" prob for letter: "<<value(it)<< " is: "<<getProbabilityForLongestContext(vlmm,it)<<endl;
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
	if(length(text)<windowSize)
		return (float)0;
	float best,result = 0;
	Iterator<String<TAlphabet> >::Type windowEnd = begin(text),windowStart = begin(text);
	goBegin(windowEnd);
	goBegin(windowStart);
	
	for(int i = 0;i<windowSize;++i,goNext(windowEnd))
			result += log(getProbabilityForLongestContext(vlmm,windowEnd));
	goNext(windowEnd);
	best = result;
	for(;!atEnd(windowEnd);goNext(windowEnd),goNext(windowStart))
	{
		result += -log(getProbabilityForLongestContext(vlmm,windowStart));
		result += log(getProbabilityForLongestContext(vlmm,windowEnd));
		//cout <<" prob for letter: "<<value(it)<< " is: "<<getProbabilityForLongestContext(vlmm,it)<<endl;
		if(result > best)
			best = result;
	}
	return best;
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





template<typename TFile, typename TAlphabet, typename TCargo >
inline void
writeHead(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ContextTree > > > > &vlmm,
	   TFile & target)
{	
		_streamWrite(target,"VLMM\tContextTree\tAlphabet\t\n");
		
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
exportNode(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec > > > > &vlmm,
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

	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const>::Type TIterConst;
	_streamPutInt(target, (unsigned)node);
	_streamPut(target, '\t');
	_streamPutInt(target, (unsigned)getFather(vlmm,node));
	_streamPut(target, '\t');
	_streamPutInt(target, (unsigned)getSuffixLink(vlmm,node));
	_streamPut(target, '\t');
	_streamPutInt(target, (unsigned)isMarked(vlmm,node));
	_streamPut(target, '\t');
//Childrens with label
	
	for(int i = 0;i<ValueSize<TAlphabet>::VALUE;++i) {
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
		for(int i = 0;i<ValueSize<TAlphabet>::VALUE;++i){
				_streamPutFloat(target,getProbability(vlmm,node,i));
				_streamPut(target,'\t');
		}

//ReverseSuffixLinks
	
		for(int i = 0;i<ValueSize<TAlphabet>::VALUE;++i){
				_streamPutInt(target,getReverseSuffixLink(vlmm,node,i));
				_streamPut(target,'\t');
		}
		_streamPut(target,'\n');
}


/* Export plain file format for
VLMM\t	Type[ContextTree,PST]\t	Alphabet[Dna,Dna5,AminoAcids]\t	Params[gamma:threshold:pMin ... etc]
Node\t	Father\t	SuffixLink\t	Marked\t	Children[1 .. N]\t	ChildLabel[1 .. N]\t		ProbabilityVector[1 .. N]\t	ReverseSuffixLink[1 .. N]	

	



//save/export the vlmm in a file for reading it again using the import funcion

*/


template<typename TFile, typename TAlphabet, typename TCargo, typename TVLMMSpec>
inline void
exportVLMM(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec > > > > &vlmm,
	   TFile & target)
{

SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TVLMMSpec> > > > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Size<TAlphabet>::Type TSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();

	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > >::Type TIterConst;

	writeHead(vlmm,target);
	for(TIterConst it = begin(vlmm.data_vertex);!atEnd(it);goNext(it)) {
		if (!idInUse(vlmm.data_id_managerV, position(it))) continue;
		TVertexDescriptor dummy = position(it); 
		exportNode(vlmm,dummy,target);
		cout << "eported node: " <<dummy<<endl;
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


