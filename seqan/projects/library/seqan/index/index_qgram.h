/*
 *  Created by Anne-Katrin Emde on 3.1.2007.
 */

#ifndef SEQAN_HEADER_INDEX_QGRAM_H
#define SEQAN_HEADER_INDEX_QGRAM_H

namespace SEQAN_NAMESPACE_MAIN
{



/**
.Function.createQGramIndex:
..summary:Builds a q-gram index on a sequence. 
..cat:Index
..signature:createQGramIndex(text,shape,pos,index)
..param.text:The sequence.
..param.shape:The shape to be used.
...type:Class.Shape
..param.pos:The resulting array that indicates at which position in index the corresponding q-grams can be found.
..param.index:The resulting list in which all q-grams are sorted alphabetically.
..returns:Index contains the sorted list of qgrams. For each possible q-gram pos contains the first position in index that corresponds to this q-gram. 
*/
template <typename TSequence, typename TArray, typename TList, typename TString, typename TShapeType>
void
createQGramIndex(TSequence & text, 
			 Shape<TString, TShapeType> & shape,
			 TArray & pos, 
			 TList & index)
{
SEQAN_CHECKPOINT
	typedef typename Position<TSequence >::Type TPosition;
	typedef typename Value<Shape<TString, TShapeType> >::Type TValue;
	typedef typename Iterator<TSequence,Standard >::Type TIterator;
	typedef typename Size<TSequence>::Type TSize;
	
	typename Iterator<TArray,Standard>::Type pos_it = begin(pos,Standard());
	typename Iterator<TArray,Standard>::Type pos_end = end(pos,Standard());
	while(pos_it != pos_end)
	{
		*pos_it = 0;
		++pos_it;
	}

	//count q-grams
	TIterator qgram_1_it = begin(text,Standard());
	TIterator qgram_2_it = qgram_1_it + 1;
	TValue x, y = 0;			// hashwert des q-Grams
	x = hash(shape, qgram_1_it);
	++pos[x];
	TPosition i = 1;
	TPosition num_qgrams = length(text) - shapeSpan(shape) + 1;
	while( i < num_qgrams)
	{
		y = hashNext(shape, qgram_1_it, qgram_2_it, x);
		++qgram_1_it;
		++qgram_2_it;
		++pos[y];
		x = y;
		++i;
	}

	//cumulative sum
	TPosition len_pos = length(pos);
	i = 1;
	while(i<len_pos)
	{
		pos[i] += pos[i-1];
		++i;
	}
	
	//build index
	qgram_1_it = begin(text,Standard());
	qgram_2_it = qgram_1_it + 1;
	x = hash(shape, qgram_1_it);
	index[--pos[x]] = 0;
	i = 1;
	while(i < num_qgrams)
	{
		y = hashNext(shape, qgram_1_it, qgram_2_it, x);
		++qgram_1_it;
		++qgram_2_it;
		index[--pos[y]] = i;

		// für die nächste Runde
		x = y;
		++i;
	}


}







}




#endif //#ifndef SEQAN_HEADER_...
