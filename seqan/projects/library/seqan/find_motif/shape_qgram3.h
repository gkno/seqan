#ifndef SEQAN_HEADER_SHAPE_QGRAM3_H
#define SEQAN_HEADER_SHAPE_QGRAM3_H

namespace SEQAN_NAMESPACE_MAIN
{

struct GappedShape3;


// Objekt Shape, klassiefiziert die Seedart. Speichert die Anzahl und die genauen Stellen innerhalb eines q-Grams,
// die bei der Seed-Suche irrelevant sind.
// Weiterhin speichert es die Alpabhet-groesse, die Groesse des q-Grams q und berechnet damit einen fuer die 
// Berechnung des hash-Wertes des entsprechenden q-Grams benoetigten term.
// Ist die Anzhal irrelevanter Stellen = 0, handelt es sich um exakte Seeds, sonst um exakte gegappte seeds.
template<typename TString>
class Shape<TString,GappedShape3>
{

public:

	int q;
	int num_gaps;
	int /*const*/ term;
	//shape_len ist die l?ge des shape vektors
	int shape_len;
	//[num_springen num_springen...]
	String<int> shape;
	
	Shape()
	{
	}
	// Konstruktor
	Shape(int num_gaps, int q, int alp_size, int shape_len):
	q(q),
	shape_len(shape_len),
	num_gaps(num_gaps),
	term(intPow(alp_size, q-1-num_gaps))
	{
		resize(shape, shape_len);
		for(int i = 0; i < shape_len; ++i)
			shape[i] = 1;
	}

	// Kopierkonstruktor
	Shape(Shape const & rhs):
	q(rhs.q),
	num_gaps(rhs.num_gaps),
	shape_len(rhs.shape_len),
	term(rhs.term),
	shape(rhs.shape)
	{
	}	


	int & operator[] (int offset){ return shape[offset];}



	~Shape(){}
	
};

// anne katrin
template <typename TString>
void
setShape(Shape<TString,GappedShape3> & me,int num_gaps, int q, int alp_size, int shape_len)
{
	me.num_gaps = num_gaps;
	me.q = q;
	me.shape_len = shape_len;
	me.term = intPow(alp_size, q-1-num_gaps);
	resize(me.shape, shape_len);
	for (int i = 0; i < shape_len; ++i)
		me.shape[i] = 1;
	me.shape[0] = 0;

}

// my version
/*template <typename TString>
void
setShape(Shape<TString,GappedShape3> & me,int num_gaps, int q, int alp_size, int shape_len)
{
	me.num_gaps = num_gaps;
	me.q = q;
	me.shape_len = shape_len;
	me.term = intPow(alp_size, q-1-num_gaps);
	resize(me.shape, shape_len);
	for (int i = 0; i < shape_len; ++i)
		me.shape[i] = 0;
}*/



template <typename TString>
inline int &
numGaps(Shape<TString, GappedShape3>  & me)
{
	return me.num_gaps;
}


template <typename TString>
inline int &
shapeLen(Shape<TString, GappedShape3>  & me)
{
	return me.shape_len;
}

// my version
/*template <typename TString, typename TSequenceIterator>
int 
hash(Shape<TString, GappedShape3> & shape, TSequenceIterator qgram_it)	
{

	typedef Shape<TString, GappedShape3> TShape;
	typedef typename Size<TShape>::Type TValue;

	typedef typename Position<TString>::Type TPosition;
	TPosition pos = 0;			// resultierende Position des q-Grams

	typename Iterator< String<TValue> >::Type iter, iter_end;
	iter = begin(shape.shape);
	iter_end = end(shape.shape);
	while(iter!=iter_end)
	{
		std::cout << "value in shape: " << *iter << "\n";
		int int_letter = (int)*(qgram_it+(*iter));
		std::cout << "int letter: " << int_letter << "\n";
		pos = 
			pos * static_cast<int>(ValueSize< TShape >::VALUE) + int_letter;
		++qgram_it;
		++iter;
	}*/



	/*TValue * pshape = begin(shape.shape);
	TValue * pshape_end = end(shape.shape);
	
	std::cout << "shape.q=" << shape.q << "\n";
	TSequenceIterator qgram_it_end = qgram_it+shape.q+1;
	while(qgram_it!=qgram_it_end)
	{
		std::cout << "*qgram_it: " << *qgram_it << "\n";
		qgram_it += *pshape;
		pos = pos * static_cast<int>(ValueSize< TShape >::VALUE) + (int)*qgram_it;
		++qgram_it;
		++pshape;
	}*/

	/*return pos;

}*/

// anne
template <typename TString, typename TSequenceIterator>
int 
hash(Shape<TString, GappedShape3> & shape, TSequenceIterator qgram_it)	
{

	typedef typename Size<Shape<TString,GappedShape3> >::Type TValue;
	typedef typename Position<TString>::Type TPosition;
	TPosition pos = 0;			// resultierende Position des q-Grams
	TValue * pshape = begin(shape.shape);
	TValue * pshape_end = end(shape.shape);//
			
	while(pshape < pshape_end)
	{
		qgram_it += *pshape;
		pos = pos * ValueSize<Shape<TString, SimpleShape> >::VALUE + (int)*qgram_it;
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

// my version
/*template <typename TString>
void 
stringToShape(String<char> & shape_string,Shape<TString,GappedShape3> & shape)
{
	int count_gaps = 0;
	int shape_len = 0;
	typename Iterator<String<char> >::Type string_it, string_end;
	string_it = begin(shape_string);
	string_end = end(shape_string);

	while(string_it!=string_end)
	{
		if(*string_it!='_')
		{
			++shape_len;
			++string_it;
		}
		else
		{
			while(*string_it=='_')
			{
				++count_gaps;
				++string_it;
			}
			++shape_len;
		}
	}

    setShape(shape,count_gaps, length(shape_string),
		     ValueSize<Shape<TString,GappedShape3> >::VALUE,
			 shape_len);

	string_it = begin(shape_string);
	int i = 0;
	while(string_it!=string_end)
	{
		if(*string_it!='_')
		{
			shape[i] = 0;
			++string_it;
		}
		else
		{
			while(*string_it=='_')
			{
				++shape[i];
				++string_it;
			}
		}
		++i;
	}
}*/

/*template <typename TString>
void 
stringToShape(String<char> & shape_string,Shape<TString,GappedShape3> & shape)
{
	int count_gaps = 0;
	int shape_len = 0;
	typename Iterator<String<char> >::Type string_it, string_it_last, string_end;
	string_it = begin(shape_string);
	string_it_last = string_it;
	++string_it;
	string_end = end(shape_string);
	while(string_it != string_end) 
	{
		if(*string_it!='_')
		{
			if(*string_it_last=='_')
			{
				++count_gaps;;
			}
			++shape_len;
			string_it_last = string_it;
			++string_it;
		}
		else
		{
			if(*string_it_last!='_')
			{
				++shape_len;
				++string_it_last;
				++string_it;
				while(*string_it == '_')
				{
					++count_gaps;
					++string_it;
				}
			}
			else
			{
				++count_gaps;
				++string_it_last;
				++string_it;
				while(*string_it == '_')
				{
					++count_gaps;
					++string_it;
				}
			}
			++count_gaps;
			++shape_len;
			string_it_last = string_it;
			++string_it;
		}
	}

	if(*string_it_last!='_')
	{
		++shape_len;
	}
	else
	{
		--string_it_last;
		if(*string_it_last!='_')
		{
			++shape_len;
			++count_gaps;
		}
		else
		{
			++count_gaps;
		}
	}

	std::cout << "#GAPS=" << count_gaps << "\n";
	std::cout << "shape_len=" << shape_len << "\n";

	setShape(shape, count_gaps, length(shape_string),
		     ValueSize<Shape<TString,GappedShape3> >::VALUE,
			 shape_len);

	int i = 0;
	string_it = begin(shape_string);
	while(string_it!=end(shape_string))
	{
		if(*string_it!='_')
		{
			shape[i] = 0;
			++string_it;
		}
		else
		{
			while(*string_it=='_')
			{
				++shape[i];
				++string_it;
			}
		}
		++i;
	}
}*/

// korrigiert
template <typename TString>
void 
stringToShape(String<char> & shape_string,Shape<TString,GappedShape3> & shape)
{
	int count_gaps = 0;
	int shape_len = 0;
	typename Iterator<String<char> >::Type string_it, string_end;
	string_it = begin(shape_string);
	string_end = end(shape_string);

	while(string_it!=string_end)
	{
		if(*string_it!='_')
		{
			++string_it;
		}
		else
		{
			while(*string_it=='_')
			{
				++count_gaps;
				++string_it;
			}
		}
	}

	shape_len = length(shape_string) - count_gaps;
	setShape(shape,count_gaps,length(shape_string),ValueSize<Shape<TString,GappedShape3> >::VALUE,shape_len);
	string_it = begin(shape_string);
	int j = 0;
	while(string_it < string_end && j < shape_len) 
	{
		if(*string_it != '_')
		{
			++j;
		}
		else
		{
			++shape[j];
		}
		++string_it;
	}
}

// anne-katrin
/*template <typename TString>
void 
stringToShape(String<char> & shape_string,Shape<TString,GappedShape3> & shape)
{
	int count_gaps = 0;
	int shape_len = 0;
	typename Iterator<String<char> >::Type string_it, string_it_last, string_end;
	string_it = begin(shape_string);

	string_it_last = string_it;
	++string_it;
	string_end = end(shape_string);
	while(string_it < string_end) 
	{
		if(*string_it == '_')
				++count_gaps;
		++string_it_last;
		++string_it;
	}

	shape_len = length(shape_string) - count_gaps;
	setShape(shape,count_gaps,length(shape_string),ValueSize<Shape<TString,GappedShape3> >::VALUE,shape_len);
	string_it = begin(shape_string);
	int j = 0;
	while(string_it < string_end && j < shape_len) 
	{
		if(*string_it != '_')
		{
			++j;
		}
		else
		{
			++shape[j];
		}
		++string_it;
	}
	//shape[shape_len-1] = 0;

}*/


}	// namespace seqan

#endif