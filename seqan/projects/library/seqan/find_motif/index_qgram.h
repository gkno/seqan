/*
 *  Created by Anne-Katrin Emde on 3.1.2007.
 */

#ifndef SEQAN_HEADER_INDEX_QGRAM_H
#define SEQAN_HEADER_INDEX_QGRAM_H

namespace SEQAN_NAMESPACE_MAIN
{

	struct QGram {};



///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// intPow, rechnet a hoch b aus (int)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//hmmmmmmmmmmmmmmmmmm?
inline int intPow(int a, int b)
{
	int c = 1;	
	for (int i = 0; i < b; ++i)
	{
		c *= a;
	}
	return c;
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hashfunktion, weist einem q-Gram eine Position zu (alphabetisch)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename TString, typename TSequenceIterator>
int 
hash(Shape<TString, SimpleShape> & shape, TSequenceIterator qgram_it)	
{

//	typename Size<Shape<TString,SimpleShape> >::Type const alp_size = ValueSize<Shape<TString, SimpleShape> >::VALUE;

	int pos = 0;			// resultierende Position des q-Grams
	int i = 0;
	int span = seedSize(shape);
	while(i < span)
	{
		pos = pos * ValueSize<Shape<TString, SimpleShape> >::VALUE + (int)*qgram_it;
		++qgram_it;
		++i;
	}


	return pos;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hash_next, weist einem q-Gram eine Position zu (alphabetisch) (bei exakten Seeds abgeleitete Berechnung,
// bei exakte gegappte Seeds Aufruf von hash)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename TSequenceIterator, typename TString>
int
hash_next(Shape<TString, SimpleShape> & shape, TSequenceIterator qgram_1_it, TSequenceIterator qgram_2_it, int & x)
{

//	typename Size<Shape<TString,SimpleShape> >::Type const alp_size = ValueSize<Shape<TString, SimpleShape> >::VALUE;
	int term = powTerm(shape);
	int span = seedSize(shape);
	int y;
    y = (x - (int)*qgram_1_it * term) * ValueSize<Shape<TString, SimpleShape> >::VALUE + (int)*(qgram_2_it+span);
	
	return y;
}



template <typename TString, typename TSequenceIterator>
int 
hash(Shape<TString, GappedShape1> & shape, TSequenceIterator qgram_it)	
{

//	typename Size<Shape<TString,GappedShape1> >::Type const alp_size = ValueSize<Shape<TString, GappedShape1> >::VALUE;

	int pos = 0;			// resultierende Position des q-Grams
	int * pshape = & shape[0];
	int i = 0;
	int span = seedSize(shape);
	while(i < span)
	{
		if(i == *pshape)
		{
			++pshape;
		}
		else
		{
			pos = pos * ValueSize<Shape<TString, SimpleShape> >::VALUE + (int)*qgram_it;
		}
		++i;
		++qgram_it;
		
	}

	return pos;
}



template <typename TSequenceIterator, typename TString>
int
hash_next(Shape<TString, GappedShape1> & shape, TSequenceIterator qgram_1_it, TSequenceIterator qgram_2_it, int & x)
{
	return hash(shape, qgram_2);
}



template <typename TString, typename TSequenceIterator>
int 
hash(Shape<TString, GappedShape2> & shape, TSequenceIterator qgram_it)	
{

//	typename Size<Shape<TString,GappedShape2> >::Type const alp_size = ValueSize<Shape<TString, GappedShape2> >::VALUE;

	int pos = 0;			// resultierende Position des q-Grams

	int len_shape = shapeLen(shape);
	int j = 0;
	while (j < shape[0])
	{
		pos = pos * ValueSize<Shape<TString, SimpleShape> >::VALUE + (int)*qgram_it;
		++qgram_it;
		++j;
	}		
	int i = 1;
	while (i < len_shape)
	{
		qgram_it += shape[i];
		//cout << shape[i]<<"\n";
		++i;
		j = 0;
		while (j < shape[i])
		{
			pos = pos * ValueSize<Shape<TString, SimpleShape> >::VALUE + (int)*qgram_it;
		//	cout << (int)*qgram_it <<" ";
			++qgram_it;
			++j;
		}		
		++i;
	}

	return pos;
}


template <typename TSequenceIterator, typename TString>
int
hash_next(Shape<TString, GappedShape2> & shape, TSequenceIterator qgram_1_it, TSequenceIterator qgram_2_it, int & x)
{
	return hash(shape, qgram_2);
}

/*template <typename TString, typename TSequenceIterator>
int 
hash(Shape<TString, GappedShape3> & shape, TSequenceIterator qgram_it)	
{

//	typename Size<Shape<TString,GappedShape3> >::Type const alp_size = ValueSize<Shape<TString, GappedShape3> >::VALUE;

	int pos = 0;			// resultierende Position des q-Grams
	int * pshape = & shape[0];
	int * pshape_end = & shape[shapeLen(shape)-1]+1;
	
	while(pshape < pshape_end)
	{
		if(*pshape>0){
			qgram_it += *pshape;
			++pshape;
			continue;
		}
		pos = pos * ValueSize<Shape<TString, SimpleShape> >::VALUE + (int)*qgram_it;
		std::cout << "*qgram_it=" << *qgram_it << "\n";
		std::cout << "pos=" << pos << "\n";
		++pshape;
	}


	return pos;
}


template <typename TSequenceIterator, typename TString>
int
hash_next(Shape<TString, GappedShape3> & shape, TSequenceIterator & qgram_1, TSequenceIterator & qgram_2, int & x)
{
	return hash(shape, qgram_2);
}
*/




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// baut QGramIndex (array von zeiger auf listenelemente) auf dem Text auf 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
template <typename TSequence, typename TArray, typename TList, typename TString, typename TShapeType>
void
createQGramIndex(TSequence & text, 
			 Shape<TString, TShapeType> & shape,
			 TArray & pos, 
			 TList & index)
{
	typedef typename Position<TSequence >::Type TPosition;
	typedef typename Iterator<TSequence >::Type TIterator;
	typedef typename Size<TSequence>::Type TSize;
	TSize q = seedSize(shape);
	
	// Array der Groesse von pos. Jedes Element enth?t die Anzahl des Vorkommens des entsprechenden q-Garms 
	// in der Query. Dazu durchlaufen wir einmal die Query

	String<int> count_ar;
	resize(count_ar, length(pos));
	typename Iterator<String<int> >::Type count_ar_it = begin(count_ar);
	typename Iterator<String<int> >::Type count_ar_end = end(count_ar);
	while(count_ar_it != count_ar_end)
	{
		*count_ar_it = 0;
		++count_ar_it;
	}
	

	typename TIterator qgram_1_it = begin(text);
	typename TIterator qgram_2_it = qgram_1_it + 1;

	// erstes q-Gram
	int x, y;			// hashwert des q-Grams
	x = hash(shape, qgram_1_it);
	++count_ar[x];
	TPosition i = 1;
	TPosition num_qgrams = length(text) - q + 1;
	while( i < num_qgrams)
	{
		y = hash_next(shape, qgram_1_it++, qgram_2_it++, x);
		++count_ar[y];
		x = y;
		++i;
	}

	// Wir durchlaufen das Array pos und setzen die Zeiger auf die richtige Position in der Liste,
	// je 1 nach dem letzten q-gram des blocks (abhaengig von # des entspr. q-Grams (count_ar[i])	
	int s = 0;
	int len_pos = length(pos);
	i = 0;
	while(i<len_pos)
	{
		s += count_ar[i];
		pos[i] = & index[s];
		++i;
	}
	
	// Wir durchlaufen die Query zum zweiten Mal und tragen in der Liste die Positionen der q-Grams in der
	// richtigen Stelle ein. Dazu schieben wir den Zeiger je 1 nach vorne und fuehren dort die Position ein
	// (damit es beim mehrfachen Vorkommen eines q-Grams auch klappt)
	// Am Ende zeigen die Zeiger auf die Position des ersten entsprechenden q-Grams in der Query.

	// erstes q-Gram
	qgram_1_it = begin(text);
	qgram_2_it = qgram_1_it + 1;

	x = hash(shape, qgram_1_it);
	*(--pos[x]) = 0;
	
	// alle anderen, leiten sich voneinander ab
	i = 1;
	while(i < num_qgrams)
	{
		y = hash_next(shape, qgram_1_it++, qgram_2_it++, x);
		*(--pos[y]) = i;

		// f? die n?hste Runde
		x = y;
		++i;
	}


}

//template <typename TSequence, typename TArray, typename TList, typename TString, typename TShapeType>
//void
//createQGramIndex2(TSequence & text, 
//			 Shape<TString, TShapeType> & shape,
//			 TArray & pos, 
//			 TList & index)
//{
//	typedef Position<TSequence >::Type TPosition;
//	typedef typename Size<TSequence>::Type TSize;
//	TSequence qgram_1, qgram_2;		
//	TSize q = seedSize(shape);
//	
//	// Array der Groesse von pos. Jedes Element enthaelt die Anzahl des Vorkommens des entsprechenden q-Garms 
//	// in der Query. Dazu durchlaufen wir einmal die Query
//
//	String<int> count_ar;
//	resize(count_ar, length(pos));
//	Iterator<String<int> >::Type count_ar_it = begin(count_ar);
//	Iterator<String<int> >::Type count_ar_end = end(count_ar);
//	while(count_ar_it != count_ar_end)
//	{
//		*count_ar_it = 0;
//		++count_ar_it;
//	}
//	
//
//	// erstes q-Gram
//	int x, y;			// hashwert des q-Grams
//	qgram_1 = infix(text, 0, q);
//	x = hash(shape, qgram_1);
//	++count_ar[x];
//	TPosition i = 1;
//	TPosition num_qgrams = length(text) - q + 1;
//	while( i < num_qgrams)
//	{
//		qgram_2 = infix(text, i, i+q);
//		y = hash_next(shape, qgram_1, qgram_2, x);
//		++count_ar[y];
//		qgram_1 = qgram_2;
//		x = y;
//		++i;
//	}
//
//	TPosition text_len = length(text);
//	TPostition B_end = text_len/2; 	
//	TPostition A_end = B_end/2; 	
//	TPostition C_end = B_end + A_end; 	
//	TPostition D_end = text_len; 	
//	
//	String<Pair<int,int> > bin_ar;
//	int bin_ar_size = 2000;
//	resize(bin_ar, bin_ar_size);//hier muss ne gescheite zahl hin
//	typedef typename Iterator<String<Pair<int,int> > >::Type BinIterator;
//
//	BinIterator bin_a_begin = begin(bin_ar);
//	BinIterator bin_b_begin = &bin_ar[bin_ar_size/4];
//	BinIterator bin_c_begin = &bin_ar[bin_ar_size/2];
//	BinIterator bin_d_begin = &bin_ar[bin_ar_size/2 + bin_ar_size/4];
//	
//	
//	
//	// Wir durchlaufen das Array pos und setzen die Zeiger auf die richtige Position in der Liste,
//	// je 1 nach dem letzten q-gram des blocks (abh?gig von # des entspr. q-Grams (count_ar[i])	
//	int s = 0;
//	int len_pos = length(pos);
//	i = 1;
////	pos[0] = &index[0];
//	while(i<len_pos)
//	{
//		count_ar[i] += count_ar[i-1];    
////		pos[i] = & index[count_ar[i-1]];
//		++i;
//	}
//	
//	// erstes q-Gram
//	qgram_1 = infix(text, 0, q);
//	x = hash(shape, qgram_1);
//	//*(--pos[x]) = 0;
//	
//	// alle anderen, leiten sich voneinander ab
//	i = 1;
//	while(i < num_qgrams)
//	{
//		qgram_2 = infix(text, i, i + q);
//		y = hash_next(shape, qgram_1, qgram_2, x);
//	//	*(--pos[y]) = i;
//
//		// f? die n?hste Runde
//		qgram_1 = qgram_2;
//		x = y;
//		++i;
//	}
//
//
//}


}

#endif //#ifndef SEQAN_HEADER_...
