/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id: markov_model.h 2009-03-10 utro@math.unipa.it srombo@deis.unical.it $
 ==========================================================================*/

//////////////////////////////////////////////////////////////////////////////

#ifndef SEQAN_HEADER_STATISTICAL_INDEX_MARKOV_MODEL_H
#define SEQAN_HEADER_STATISTICAL_INDEX_MARKOV_MODEL_H

#include<seqan/align.h>
#include <seqan/index.h>

namespace seqan//SEQAN_NAMESPACE_MAIN
{
/*
with matrix calss matrices
*/
/**
.Class.MarkovModel:
..summary:Gives a suitable representation of a Marcov Chain.
..cat:MarkovModel
..signature:MarkovModel<TAlphabet,[TFloat, TSpec]>
..param.TAlphabet:The alphabet type
..param.TFloat:The type of the exploited arrays
..param.TSpec:The MarkovModel type
.Memvar.MarkovModel#order:
..class:Class.MarkovModel
..summary:The MarkovModel order
..type:nolink:int
.Memvar.MarkovModel#transition:
..class:Class.MarkovModel
..summary:The transition matrix
...type:Class.Matrix
.Memvar.MarkovModel#stationaryDistribution:
..class:Class.MarkovModel
..summary:The vector of character distribution
..type:Class.String
.Memfunc.MarkovModel#MarkovModel:
..class:Class.MarkovModel
..summary:Constructor
..signature:MarkovModel(order_)
..param.order_:The order of the MarkovModel.
.Memfunc.MarkovModel#build:
..class:Class.MarkovModel
..summary:Given a training set, computes the transition matrix, the character stationary distributions and the auxiliary information that give raise to an instance of MarkovModel
..signature:build(strings)
..param.strings:The training set.
...type:Class.StringSet

.Memfunc.MarkovModel#set:
..class:Class.MarkovModel
..summary: Given e transition matrix, sets it as transition matrix of the MarkovModel and computes (if it is not available) the  vector of character distributions and the auxiliary information
..signature:set(iTransition)
..signature:set(iTransition, iStationaryDistribution)
..param.iTransition:The transition matrix.
...type:Class.Matrix
..param.iStationaryDistribution:The vector of character distributions.
...type:Class.String

.Memfunc.MarkovModel#emittedProbability:
..class:Class.MarkovModel
..summary:Computes the probability that a string (or a set of strings) is emitted by the MarkovModel.
..signature:emittedProbability(string)
..signature:emittedProbability(stringSet)
..param.string:The string whose emission probability has to be computed.
...type:Class.String
..param.stringSet:The set of strings whose emission probability has to be computed.
...type:Class.StringSet
..returns:A TFloat representing the emission probability.

.Memfunc.MarkovModel#write:
..class:Class.MarkovModel
..summary: Stores an instance of MarkovModel on a file
..signature:write(file)
..param.file:The file on which storing the MarkovModel.

.Memfunc.MarkovModel#read:
..class:Class.MarkovModel
..summary: Loads an instance of MarkovModel from a file
..signature:read(file)
..param.file:The file from which loading the MarkovModel.
*/
template <typename TAlphabet, typename TFloat = double, typename TSpec = Default>
class MarkovModel
{

public:

	//Definition of matrix types, could be set to two dimensional as soon as implemented
	typedef Matrix<TFloat, 2> TMatrix;
	//typedef String<TFloat> TMatrix;
	typedef String<TFloat> TVector;

	unsigned int order;
	TMatrix transition;
	TVector stationaryDistribution;
	//The following matrices are only for internal use of the class
	TMatrix _q;			//QPP^(morder-1)
	TMatrix _qppp;		//QPP^(morder-1)
	TMatrix _qppqpp;	//QPPQP^(morder-1)


	MarkovModel(unsigned int order_):
	order(order_)
	{

		unsigned int const alphabet_size = ValueSize<TAlphabet>::VALUE;
		//unsigned int const column_size = (unsigned int) std::pow((double) alphabet_size, (int) order);
		unsigned int column_size = (unsigned int) std::pow((double) alphabet_size, (int) order);

		//Special case of order 0 marko model:
		//A Bernoulli scheme is a special case of a Markov chain where the transition probability matrix has identical rows, which means that the next state is even independent of the current state (in addition to being independent of the past states).
		if(order==0)
		{
			column_size=alphabet_size;

		}
		//resize the matrix
		setLength(transition, 0, column_size);
		setLength(transition, 1, column_size);

		fill(transition,(TFloat) 0);

		clear(stationaryDistribution);
		resize(stationaryDistribution, column_size);
	}


	///////////////////////////////////////////////////////////////
	///// BUILD THE MODEL
	///////////////////////////////////////////////////////////////
	typedef String<Dna> TText;
	typedef Size<TText>::Type TSize;
	typedef StringSet<TText > TStringSet;

	//template <typename TAlphabet>
	void build(StringSet<String<TAlphabet> > const &strings)
	{

		typedef String<TAlphabet> TText;
		//typedef Size<TText>::Type TSize;
		typedef StringSet<TText > TStringSet;

		typedef Index<TStringSet, Index_QGram<SimpleShape> > TIndex;
		typedef typename Fibre<TIndex, QGram_Dir>::Type TDir;
		typedef typename Iterator<TDir, Standard>::Type TIter;
		unsigned int const alphabet_size = ValueSize<TAlphabet>::VALUE;
		//unsigned int const column_size = (unsigned int) std::pow((double) alphabet_size, (int) order);
		unsigned int column_size = (unsigned int) std::pow((double) alphabet_size, (int) order);

		//Special case of order 0 marko model:
		//A Bernoulli scheme is a special case of a Markov chain where the transition probability matrix has identical rows, which means that the next state is even independent of the current state (in addition to being independent of the past states).
		if(order==0)
		{
			column_size=alphabet_size;

		}

		TIndex ind(strings);
		resize(indexShape(ind), order + 1);
		indexRequire(ind, QGram_SADir());

		TIter itBegin = begin(indexDir(ind), Standard());
		TIter itEnd = end(indexDir(ind), Standard()) - 1;

		//Frequency of all q-grams for a markov model of order q-1 to calculate transition probabilitiesof (q-1) gram to next (q-1)gram
		String<TAlphabet> qgram;
		Shape<TAlphabet, SimpleShape> orderShape;
		resize(orderShape, order);
		if(order==0)
		{
			resize(orderShape, order+1);
		}

		//count transition of sequences according to model order
		for(TIter itA = itBegin; itA != itEnd; ++itA)
		{
			unhash(qgram, itA - itBegin, weight(indexShape(ind)));
			//std::cout<<"\n"<<qgram<<"\n"<<(*(itA+1) - *itA)<<"\n"<<hash(orderShape, begin(qgram));
			//old for the array: value(transition, hash(orderShape, begin(qgram)) * column_size + hash(orderShape, begin(qgram) + 1)) = *(itA+1) - *itA;

			if(order==0)
			{
				for(unsigned int row=0;row<column_size;++row)
				{
					value(transition, row,hash(orderShape, begin(qgram))) = *(itA+1) - *itA;

				}
			}
			else
			{
				//new for the matrix
				value(transition, hash(orderShape, begin(qgram)),hash(orderShape, begin(qgram) + 1)) = *(itA+1) - *itA;

			}
		}
		//std::cout<<"\n"<<transition<<"\n\n";

		//normalization, rows have to sum up to 1
		for(unsigned int row = 0; row < column_size;++row)
		{
			TFloat sum = 0;
			for(unsigned int col = 0; col < column_size;++col)
			{
				sum += value(transition, row,col);	//sum of the rows
			}
			if (sum != 0)
			{
				for(unsigned int col = 0; col < column_size;++col)
				{
					value(transition, row, col) /= sum;	//normalize by dividing by sum of rows
				}
			}
		}
//std::cout<<transition<<"\n\n";
		//Special case of order 0 marko model:
		//A Bernoulli scheme is a special case of a Markov chain where the transition probability matrix has identical rows, which means that the next state is even independent of the current state (in addition to being independent of the past states).
		if(order==0)
		{
			order=1;
		}
		//----Calculation of stationary Distribution-----
		TMatrix temp = transition;
		//initialise a variable t representing a good threshold to estimate the vector
		//after multiplying t times the transition matrix with itself
		unsigned int t=6;
		for (unsigned int i=0; i<t; i++)
		{
			temp=temp*temp;
		}

		for (unsigned int i=0; i<column_size; i++){
			value(stationaryDistribution,i)=value(temp, 0,i);

		}
		//std::cout<<temp<<"\n\n";

		//for(unsigned int i=0;i<length(stationaryDistribution);++i){std::cout<<stationaryDistribution[i]<<"\t";}
		//!is commented since I dont use it and it makes everything very slow for k>3 or 4 
		//!call ensureAuxMatrices(markovModel);
		//_computeAuxiliaryMatrices();


	}


	///////////////////////////////////////////////////////////////
	///// EMITTEDPROBABILITY
	///////////////////////////////////////////////////////////////


	template <typename TString, typename TSetSpec>
	TFloat emittedProbability(StringSet<TString, TSetSpec > const &stringSet)
	{
		TFloat p = 0;

		for(unsigned int i=0; i<length(stringSet); i++)
		{
			p+= emittedProbability(stringSet[i]);
		}

		return p;
	}

	template <typename TString>
    TFloat emittedProbability(TString const &string)
	{
		Shape<TAlphabet, SimpleShape> orderShape;
		resize(orderShape, order);

		int row = hash(orderShape,begin(string));
		TFloat p = value(stationaryDistribution,row);

		for(unsigned int i=1; i<(length(string)-order+1); i++)
		{
			int column=hash(orderShape,begin(string)+i);
			p*=value(transition,row,column);

			row = column;
		}

		return p;
	}


	///////////////////////////////////////////////////////////////
	///// SET THE MODEL
	///////////////////////////////////////////////////////////////


	void set(Matrix<TFloat,2> &iTransition)
	{
		unsigned int const alphabet_size = ValueSize<TAlphabet>::VALUE;
		unsigned int const column_size = (unsigned int) std::pow((double) alphabet_size, (int) order);

		transition = iTransition;

		Matrix<TFloat,2> temp = transition;
		//initialise a variable t representing a good threshold to estimate the vector
		//after multiplying e times the transition matrix with itself
		unsigned int t=6;
		for (unsigned int i=0; i<t; i++){
			temp=temp*temp;
		}

		for (unsigned int i=0; i<column_size; i++){
			value(stationaryDistribution,i)=value(temp, 0,i);	//Set stationary distribution to row of transition^6
		}

		//_computeAuxiliaryMatrices();
	}



	void set(Matrix<TFloat,2> &iTransition, String<TFloat> &iStationaryDistribution)
	{
		transition = iTransition;

		stationaryDistribution = iStationaryDistribution;

		//_computeAuxiliaryMatrices();
	}


	///////////////////////////////////////////////////////////////
	///// WRITE
	///////////////////////////////////////////////////////////////


	void write(FILE *file)
	{
		ensureAuxMatrices(*this);
		unsigned int const alphabet_size = ValueSize<TAlphabet>::VALUE;
		unsigned int const column_size = (unsigned int) std::pow((double) alphabet_size, (int) order);

		//write the transition matrix
		for(unsigned int row=0; row<column_size; row++){
			for(unsigned int col=0; col<column_size; col++){
			  fprintf(file,"%f ",value(transition, row,col));
			}
			fprintf(file,"\n");
		}
		//write the stationary distribution vector
		for(unsigned int row=0; row<column_size; row++){
			  fprintf(file,"%f ",value(stationaryDistribution, row));
			}
		fprintf(file,"\n");

		if(length(_q)){
			//write the auxiliary matrix
			for(unsigned int row=0; row<column_size; row++){
				for(unsigned int col=0; col<column_size; col++){
					fprintf(file,"%f ",value(_q, row,col));
				}
				fprintf(file,"\n");
			}

			for(unsigned int row=0; row<column_size; row++){
				for(unsigned int col=0; col<column_size; col++){
				  fprintf(file,"%f ",value(_qppp,row,col));
				}
				fprintf(file,"\n");
			}

			for(unsigned int row=0; row<column_size; row++){
				for(unsigned int col=0; col<column_size; col++){
				fprintf(file,"%f ",value(_qppqpp, row,col));
				}
				fprintf(file,"\n");
			}
		}
	}


	///////////////////////////////////////////////////////////////
	///// READ
	///////////////////////////////////////////////////////////////

	void read(FILE *file)
	{
		unsigned int const alphabet_size = ValueSize<TAlphabet>::VALUE;
		unsigned int const column_size = (unsigned int) std::pow((double) alphabet_size, (int) order);

		//read the transition matrix
		for(unsigned int row=0; row<column_size; row++)
		{
			for(unsigned int col=0; col<column_size; col++)
			{
			  fscanf(file,"%lf ", & value(transition, row,col));
			}
			fscanf(file,"\n");
		}
		//read the stationary distribution vector
		for(unsigned int row=0; row<column_size; row++)
		{
			  fscanf(file,"%lf ",&value(stationaryDistribution, row));
		}
		fscanf(file,"\n");

		if (!feof(file))
		{
			setLength(_q, 0, column_size);
			setLength(_q, 1, column_size);
			resize(_q);

			//read the auxiliary matrix
			for(unsigned int row=0; row<column_size; row++)
			{
				for(unsigned int col=0; col<column_size; col++)
				{
					fscanf(file,"%lf ",&value(_q, row,col));
				}
				fscanf(file,"\n");
			}
			setLength(_qppp, 0, column_size);
			setLength(_qppp, 1, column_size);
			resize(_qppp);
			for(unsigned int row=0; row<column_size; row++){
				for(unsigned int col=0; col<column_size; col++){
				  fscanf(file,"%lf ",&value(_qppp, row,col));
				}
				fscanf(file,"\n");
			}
			setLength(_qppqpp, 0, column_size);
			setLength(_qppqpp, 1, column_size);
			resize(_qppqpp);
			for(unsigned int row=0; row<column_size; row++){
				for(unsigned int col=0; col<column_size; col++){
					fscanf(file,"%lf ",&value(_qppqpp, row,col));
				}
				fscanf(file,"\n");
			}
		}
	}


	/////////////////////////////////////////////////////////////////////////////
	///// COMPUTE THE AUXILIARY MATRICES FOR THE VARIANCE AND Z-SCORE COMPUTATION
	/////////////////////////////////////////////////////////////////////////////

	/*
		.Memfunc.MarkovModel#_computeAuxiliaryMatrices:
		..class:Class.MarkovModel
		..summary:Computes the auxiliary information for statistical indices computation
		..signature:_computeAuxiliaryMatrices()
	*/

	void _computeAuxiliaryMatrices()
	{
	//std::cout<<"auxMat\n";
		//clear(_q);
		//clear(_qppp);
		//clear(_qppqpp);

		unsigned int const alphabet_size = ValueSize<TAlphabet>::VALUE;
		unsigned int const column_size = (unsigned int) std::pow((double) alphabet_size, (int) order);
		TMatrix I;
		TMatrix Ip;

		//resize the matrices
		setLength(I, 0, column_size);
		setLength(I, 1, column_size);
		fill(I, 0.0);

		setLength(Ip, 0, column_size);
		setLength(Ip, 1, column_size);
		fill(Ip, 0.0);

		for(unsigned int i=0; i<column_size; i++){
			value(I,i,i)=1.0;
			 for (unsigned int j=0; j<column_size; j++)
			    value(Ip,i,j)=value(stationaryDistribution,j);
		}


		_q=transition-I;
		_q=_q+Ip;
		_q=_computeInverseMatrix(_q);//works for simple non singular matrices, others not checked
		//original code: _q=_computeInverseMatrix(_matricialSum(_matricialDifference(transition, I), Ip));

		_qppp=_q*transition;
		_qppqpp = _qppp*transition;
		_qppqpp = _qppqpp*_q;
		for(unsigned int i=1; i<order; i++)
		{
			_qppp=_qppp*transition;
			_qppqpp=_qppqpp*transition;
		}
	}


};

//////////////////////////////////////////////////////////////////////////////
// _computeInverseMatrix
//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeInverseMatrix:
..summary:Computes the inverse matrix of a given matrix
..signature:_computeInverseMatrix(matrix,n)
..param.matrix:The matrix in input.
...type:String<TAlphabet, TSpec1>&
..param.n:The number of columns of the matrix.
...type:nolink:unsigned int
..returns:The inverse matrix of the matrix.
*/

template <typename TValue>
Matrix<TValue,2> _computeInverseMatrix(Matrix<TValue,2> &matrix)
{
	typedef Matrix<TValue,2> TMatrix;
	unsigned int n = length(matrix,0);
	TMatrix result;
	//resize the matrix
	setLength(result, 0, n);
	setLength(result, 1, n);
	fill(result, 0.0);

	//copy the matrix in result, since the procedure is in-place
	TMatrix tmp = matrix;

	//lu decomposition of a in-place
	String<TValue> indx=_ludcmp(tmp);

	String<TValue> col;

	unsigned int i;

	// inverse by columns
	for (unsigned int j=0; j<n; j++)
	{
		fill(col,n,0);
		if(j>0)
		{
			for(i=0; i<n; i++)
			{
				value(col,i)=0;
			}
		}
		value(col,j) = 1;

		_lubksb(tmp,indx,col);

		for (i=0; i<n; i++)
		{
			value(result, i,j)= value(col,i);
		}

	}

	return result;
}

/*
 *
 * LU decomposition
 *
 */

#define TINY 1.0e-20
template <typename TValue>
String<TValue> _ludcmp(Matrix<TValue,2> &result)
{
	int n = length(result,0);
	int i, imax, j, k,d;
	double big,dum,sum,temp;
	String<TValue> vv;
	fill(vv, n, 1.0);


	d = 1;
	for (i=1; i<=n; i++)
	{
		big = 0.0;
		for (j=1; j<=n; j++)
		{
			if ((temp=fabs(value(result, i-1,j-1)))>big)
			{
				big = temp;
			}
		}
		if (big==0.0)
		{
			std::cout<<"Singular matrix in routine ludcmp" << std::endl;
			exit(1);
		}

		value(vv,i-1) = 1.0/big;
	}
	String<TValue> indx;
	resize(indx,n);

	for (j=1; j<=n; j++)
	{
		for (i=1; i<j; i++)
		{
			sum = value(result,i-1,j-1);
			for (k=1; k<i; k++)
			{
				sum -= value(result, i-1,k-1)*value(result, k-1,j-1);
			}
			value(result, i-1,j-1) = sum;
		}
		big = 0.0;
		for (i=j; i<=n; i++)
		{
			sum = value(result,i-1,j-1);
			for (k=1; k<j; k++)
			{
				sum -= value(result, i-1,k-1)*value(result, k-1,j-1);
			}
			value(result, i-1,j-1) = sum;
			if ((dum = value(vv, i-1)*fabs(sum))>=big)
			{
				big = dum;
				imax = i;
			}
		}
		if (j != imax)
		{
			for (k=1; k<=n; k++)
			{
				dum = value(result,imax-1,k-1);
				value(result, imax-1,k-1) = value(result, j-1,k-1);
				value(result,j-1,k-1) = dum;
			}
			d = -(d);
			value(vv, imax-1)=value(vv,j-1);
		}

		value(indx, j-1) = imax;

		if (value(result, j-1,j-1) == 0.0)
		{
			value(result, j-1,j-1) = TINY;
		}
		if (j!=n)
		{
			dum = 1.0/(value(result,j-1,j-1));
			for (i=j+1; i<=n; i++)
			{
				value(result, i-1,j-1) *= dum;
			}
		}
  }

 return indx;

}


template <typename TValue>
void _lubksb(Matrix<TValue,2> &a, String<TValue> &indx, String<TValue>  &b)
{
	int n =length(a,0);	//Number of columns in matrix a
	int i, ii=0,ip,j;
	double sum;

	for (i=1; i<=n; i++)
	{
		ip = value(indx,i-1);
		sum = value(b,ip-1);
		value(b,ip-1) = value(b,i-1);
		if (ii)
		{
			for (j=ii;j<=i-1;j++)
			{
				sum -=value(a,i-1,j-1)*value(b,j-1);
			}
		}
		else
		if (sum)
		{
			ii=i;
		}
		value(b,i-1) = sum;
	}
	for (i=n; i>=1; i--)
	{
		sum = value(b,i-1);
		for (j=i+1; j<=n; j++)
		{
			sum -= value(a,i-1,j-1)*value(b,j-1);
		}
		value(b,i-1) = sum/value(a,i-1,i-1);
	}
}


//////////////////////////////////////////////////////////////////////////////
// Interface
//////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TFloat, typename TSpec>
void buildMarkovModel(MarkovModel<TAlphabet, TFloat, TSpec> &mm,
		StringSet<String<TAlphabet > > &stringSet)
{

	mm.build(stringSet);

}

template <typename TAlphabet, typename TFloat, typename TSpec>
void buildMarkovModel(MarkovModel<TAlphabet, TFloat, TSpec> &mm,
		StringSet<String<TAlphabet > > const &stringSet)
{

	mm.build(stringSet);

}


///////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TFloat, typename TSpec>
void setMarkovModel(MarkovModel<TAlphabet, TFloat, TSpec> & mm,
					Matrix<TFloat,2> &iTransition)
{
	mm.set(iTransition);
}



template <typename TAlphabet, typename TFloat, typename TSpec>
void setMarkovModel(MarkovModel<TAlphabet, TFloat, TSpec> & mm,
					Matrix<TFloat,2> &iTransition,
					String<TFloat> &iStationaryDistribution)
{
	mm.set(iTransition, iStationaryDistribution);
}

///////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TFloat, typename TSpec>
TFloat emittedProbability(MarkovModel<TAlphabet, TFloat, TSpec> & mm,
		StringSet<String<TAlphabet > > & stringSet)
{
	return mm.emittedProbability(stringSet);
}

template <typename TAlphabet, typename TFloat, typename TSpec>
TFloat emittedProbability(MarkovModel<TAlphabet, TFloat, TSpec> & mm,
						  String<TAlphabet> &string)
{
	return mm.emittedProbability(string);
}

//const
template <typename TAlphabet, typename TFloat, typename TSpec, typename TString, typename TSetSpec>
TFloat emittedProbability(MarkovModel<TAlphabet, TFloat, TSpec> & mm,
		StringSet<TString, TSetSpec> const &stringSet)
{
	return mm.emittedProbability(stringSet);
}

template <typename TAlphabet, typename TFloat, typename TSpec, typename TString>
TFloat emittedProbability(MarkovModel<TAlphabet, TFloat, TSpec> & mm,
						  TString const &string)
{
	return mm.emittedProbability(string);
}

///////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TFloat, typename TSpec>
void write(FILE *file,
		   MarkovModel<TAlphabet, TFloat, TSpec> & mm )
{
	mm.write(file);
}

///////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TFloat, typename TSpec>
void read(FILE *file,
		  MarkovModel<TAlphabet, TFloat, TSpec> & mm )
{
	mm.read(file);
}
//////////////////////////////////////////////////////////////////////////////
template <typename TAlphabet, typename TFloat, typename TSpec>
void ensureAuxMatrices(MarkovModel<TAlphabet, TFloat, TSpec> & mm )
{

	if(empty(mm._q)){
		mm._computeAuxiliaryMatrices();
	}
}


}

#endif //#ifndef SEQAN_HEADER_...
