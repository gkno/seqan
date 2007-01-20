/*
 *  Created by Anne-Katrin Emde on 3.1.2007.
 */

#ifndef SEQAN_HEADER_INDEX_QGRAM_H
#define SEQAN_HEADER_INDEX_QGRAM_H

namespace SEQAN_NAMESPACE_MAIN
{

	struct QGram1 {};
	struct QGram2 {};
	struct QGram3 {};



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



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// baut QGramIndex (array von zeiger auf listenelemente) auf dem Text auf 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
template <typename TSequence, typename TArray, typename TList, typename TString, typename TShapeType>
void
createQGramIndex(TSequence & text, 
			 Shape<TString, TShapeType> & shape,
			 QGram1 const &,
			 TArray & pos, 
			 TList & index)
{
	typedef typename Position<TSequence >::Type TPosition;
	typedef typename Iterator<TSequence >::Type TIterator;
	typedef typename Size<TSequence>::Type TSize;
	TSize q = seedSize(shape);
	
	// Array der Größe von pos. Jedes Element enthält die Anzahl des Vorkommens des entsprechenden q-Garms 
	// in der Query. Dazu durchlaufen wir einmal die Query

	Iterator<TArray>::Type pos_it = begin(pos);
	Iterator<TArray>::Type pos_end = end(pos);
	while(pos_it != pos_end)
	{
		*pos_it = 0;
		++pos_it;
	}

	typename TIterator qgram_1_it = begin(text);
	typename TIterator qgram_2_it = qgram_1_it + 1;

	// erstes q-Gram
	int x, y = 0;			// hashwert des q-Grams
	x = hash(shape, qgram_1_it);
	++pos[x];
	TPosition i = 1;
	TPosition num_qgrams = length(text) - q + 1;
	while( i < num_qgrams)
	{
		y = hash_next(shape, qgram_1_it, qgram_2_it, x);
		++qgram_1_it;
		++qgram_2_it;
		++pos[y];
		x = y;
		++i;
	}

	// Wir durchlaufen das Array pos und setzen die Zeiger auf die richtige Position in der Liste,
	// je 1 nach dem letzten q-gram des blocks (abhängig von # des entspr. q-Grams (pos[i])	
	int len_pos = length(pos);
	i = 1;
	while(i<len_pos)
	{
		pos[i] += pos[i-1];
		++i;
	}
	
	// Wir durchlaufen die Query zum zweiten Mal und tragen in der Liste die Positionen der q-Grams in der
	// richtigen Stelle ein. Dazu schieben wir den Zeiger je 1 nach vorne und fügen dort die Position ein
	// (damit es beim mehrfachen Vorkommen eines q-Grams auch klappt)
	// Am Ende zeigen die Zeiger auf die Position des ersten entsprechenden q-Grams in der Query.

	// erstes q-Gram
	qgram_1_it = begin(text);
	qgram_2_it = qgram_1_it + 1;

	x = hash(shape, qgram_1_it);
	index[--pos[x]] = 0;
	
	// alle anderen, leiten sich voneinander ab
	i = 1;
	while(i < num_qgrams)
	{
		y = hash_next(shape, qgram_1_it, qgram_2_it, x);
		++qgram_1_it;
		++qgram_2_it;
		index[--pos[y]] = i;

		// für die nächste Runde
		x = y;
		++i;
	}


}


template <typename TList, typename TIterator>
void _copy (TList & index,
		   TIterator & bin_x_begin, TIterator & bin_x_it)
{
//	TBin = String<Pair<int,int> >
	while(bin_x_it>bin_x_begin)
	{
		--bin_x_it;
		index[(*bin_x_it).i2] = (*bin_x_it).i1;
	}

}



template <typename TSequence, typename TArray, typename TList, typename TString, typename TShapeType>
void
createQGramIndex(TSequence & text, 
			 Shape<TString, TShapeType> & shape,
			 QGram2 const &,
			 TArray & pos, 
			 TList & index)
{
	typedef typename Position<TSequence >::Type TPosition;
	typedef typename Position<TSequence >::Type TPosPosition;
	typedef typename Iterator<TSequence >::Type TIterator;
	typedef typename Size<TSequence>::Type TSize;
	TSize q = seedSize(shape);
	
	// Array der Größe von pos. Jedes Element enthält die Anzahl des Vorkommens des entsprechenden q-Garms 
	// in der Query. Dazu durchlaufen wir einmal die Query

	Iterator<TArray>::Type pos_it = begin(pos);
	Iterator<TArray>::Type pos_end = end(pos);
	while(pos_it != pos_end)
	{
		*pos_it = 0;
		++pos_it;
	}
	
	// erstes q-Gram
	int x, y;			// hashwert des q-Grams

	TIterator qgram_1_it = begin(text);
	TIterator qgram_2_it = qgram_1_it + 1;
	x = hash(shape, qgram_1_it);
	++pos[x];
	TPosition i = 1;
	TPosition num_qgrams = length(text) - q + 1;
	while( i < num_qgrams)
	{
		y = hash_next(shape, qgram_1_it, qgram_2_it, x);
		++qgram_1_it;
		++qgram_2_it;
		++pos[y];
		x = y;
		++i;
	}

	int len_pos = length(pos);
	i = 1;
	while(i<len_pos)
	{
		pos[i] += pos[i-1];    
		++i;
	}

	TPosition text_len = length(text);
	TPosition const block_size = 4096; 	
	TPosition const num_blocks = (text_len/block_size) + 1; 	
	
	String<String<Pair<TPosition,TPosPosition> > > bin_ar;
	TPosition const bin_ar_size = 4096*2;
	TPosition const bin_ar_bin_size = bin_ar_size/num_blocks;
	resize(bin_ar, num_blocks);
	for(TPosition k = 0; k < num_blocks; ++k)
		resize(bin_ar[k], bin_ar_bin_size);
	typedef typename Iterator<String<Pair<TPosition,TPosPosition> > >::Type BinIterator;
	//typedef typename Iterator<String<Pair<TPosition,TPosPosition> > >::Type BinIterator;

	String<BinIterator> bin_its;
	resize(bin_its,num_blocks);
	for(TPosition k = 0; k < num_blocks; ++k)
		bin_its[k] = begin(bin_ar[k]);
	
	String<BinIterator> bin_ends;
	resize(bin_ends,num_blocks);
	for(TPosition k = 0; k < num_blocks; ++k)
		bin_ends[k] = end(bin_ar[k]);

	// Wir durchlaufen das Array pos und setzen die Zeiger auf die richtige Position in der Liste,
	// je 1 nach dem letzten q-gram des blocks (abhängig von # des entspr. q-Grams (pos[i])	
	qgram_1_it = begin(text);
	qgram_2_it = qgram_1_it + 1;
	x = hash(shape, qgram_1_it);
	--pos[x];
	
	TPosition bin_x = pos[x]/block_size;

	*bin_its[bin_x] = Pair<TPosition,TPosPosition>(0,pos[x]);
	++bin_its[bin_x];

	// alle anderen, leiten sich voneinander ab
	i = 1;
	while(i < num_qgrams)
	{
		y = hash_next(shape, qgram_1_it, qgram_2_it, x);
		++qgram_1_it;
		++qgram_2_it;
		//find the right temporary bin
		--pos[y];
		bin_x = pos[y]/block_size;
		*bin_its[bin_x] = Pair<TPosition,TPosPosition>(i,pos[y]);
		++bin_its[bin_x];
		if(bin_its[bin_x] == bin_ends[bin_x])
		{
			_copy(index,begin(bin_ar[bin_x]),bin_its[bin_x]);
		}
		// für die nächste Runde
		x = y;
		++i;

	}

	for(TPosition k = 0; k < num_blocks; ++k)
		_copy(index,begin(bin_ar[k]),bin_its[k]);


	}






template <typename TSequence, typename TArray, typename TList, typename TString, typename TShapeType>
void
createQGramIndex(TSequence & text, 
			 Shape<TString, TShapeType> & shape,
			 QGram3 const &,
			 TArray & pos, 
			 TList & index)
{
	typedef typename Position<TSequence >::Type TPosition;
	typedef typename Iterator<TSequence >::Type TIterator;
	typedef typename Size<TSequence>::Type TSize;
	TSize q = seedSize(shape);
	
	// Array der Größe von pos. Jedes Element enthält die Anzahl des Vorkommens des entsprechenden q-Garms 
	// in der Query. Dazu durchlaufen wir einmal die Query

	Iterator<TArray>::Type pos_it = begin(pos);
	Iterator<TArray>::Type pos_end = end(pos);
	while(pos_it != pos_end)
	{
		*pos_it = 0;
		++pos_it;
	}
	
	// erstes q-Gram
	int x, y;			// hashwert des q-Grams

	TIterator qgram_1_it = begin(text);
	TIterator qgram_2_it = qgram_1_it + 1;
	x = hash(shape, qgram_1_it);
	++pos[x];
	TPosition i = 1;
	TPosition num_qgrams = length(text) - q + 1;
	while( i < num_qgrams)
	{
		y = hash_next(shape, qgram_1_it, qgram_2_it, x);
		++qgram_1_it;
		++qgram_2_it;
		++pos[y];
		x = y;
		++i;
	}

	TPosition text_len = length(text);
	TPosition B_end = text_len/2; 	
	TPosition A_end = B_end/2; 	
	TPosition C_end = B_end + A_end; 	
	TPosition D_end = text_len; 	
	
	String<Pair<int,int> > bin_ar;
	int const bin_ar_size = 4096;
	int const bin_ar_bin_size = bin_ar_size/4;
	resize(bin_ar, bin_ar_size);//hier muss ne gescheite zahl hin
	typedef typename Iterator<String<Pair<int,int> > >::Type BinIterator;

	BinIterator bin_a_begin = begin(bin_ar);
	BinIterator bin_b_begin = bin_a_begin + bin_ar_bin_size;
	BinIterator bin_c_begin = bin_b_begin + bin_ar_bin_size;
	BinIterator bin_d_begin = bin_c_begin + bin_ar_bin_size;
	BinIterator bin_ar_end = end(bin_ar);
	

	BinIterator bin_a_it = bin_a_begin;
	BinIterator bin_b_it = bin_b_begin;
	BinIterator bin_c_it = bin_c_begin;
	BinIterator bin_d_it = bin_d_begin;
	
	

	// Wir durchlaufen das Array pos und setzen die Zeiger auf die richtige Position in der Liste,
	// je 1 nach dem letzten q-gram des blocks (abhängig von # des entspr. q-Grams (pos[i])	
	int s = 0;
	int len_pos = length(pos);
	i = 1;
	while(i<len_pos)
	{
		pos[i] += pos[i-1];    
		++i;
	}

	qgram_1_it = begin(text);
	qgram_2_it = qgram_1_it + 1;
	x = hash(shape, qgram_1_it);
	--pos[x];
	
	int bin_x = pos[x]/bin_ar_bin_size;

	if(bin_x == 0)
	{
		*bin_a_it = Pair<int,int>(0,pos[x]);
		++bin_a_it;
	}
	else{
		if(bin_x == 1)
		{
			*bin_b_it = Pair<int,int>(0,pos[x]);
			++bin_b_it;	
		}
		else{
			if(bin_x == 2)
			{
				*bin_c_it = Pair<int,int>(0,pos[x]);
				++bin_c_it;	
			}
			else
			{
				*bin_d_it = Pair<int,int>(0,pos[x]);
				++bin_d_it;
			}
		}
	}
	
	int a,b,c,d;
	a=b=c=d=0;
	// alle anderen, leiten sich voneinander ab
	i = 1;
	while(i < num_qgrams)
	{
		y = hash_next(shape, qgram_1_it, qgram_2_it, x);
		++qgram_1_it;
		++qgram_2_it;
		//find the right temporary bin
		--pos[y];
		bin_x = pos[y]/bin_ar_bin_size;
		if(bin_x == 0)
		{
			*bin_a_it = Pair<int,int>(i,pos[y]);
			++bin_a_it;
			if(bin_a_it==bin_b_begin)
			{
				_copy(index,bin_a_begin,bin_a_it);
				++a;
			}
		}
		else
		{
			if(bin_x == 1)
			{
				*bin_b_it = Pair<int,int>(i,pos[y]);
				++bin_b_it;
				if(bin_b_it==bin_c_begin)
				{
					_copy(index,bin_b_begin,bin_b_it);
					b++;
				}
			}
			else
			{
				if(bin_x == 2)
				{
					*bin_c_it = Pair<int,int>(i,pos[y]);
					++bin_c_it;
					if(bin_c_it==bin_d_begin)
					{
						_copy(index,bin_c_begin,bin_c_it);
						++c;
					}			
				}
				else
				{
					*bin_d_it = Pair<int,int>(i,pos[y]);
					++bin_d_it;
					if(bin_d_it==bin_ar_end)
					{
						_copy(index,bin_d_begin,bin_d_it);
						d++;
					}
				}
			}
		}
		// für die nächste Runde
		x = y;
		++i;
	}


	_copy(index,bin_a_begin,bin_a_it);
	_copy(index,bin_b_begin,bin_b_it);
	_copy(index,bin_c_begin,bin_c_it);
	_copy(index,bin_d_begin,bin_d_it);



	}



}




#endif //#ifndef SEQAN_HEADER_...
