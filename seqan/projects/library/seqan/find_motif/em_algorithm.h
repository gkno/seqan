#ifndef SEQAN_HEADER_EM_ALGORITHM_H
#define SEQAN_HEADER_EM_ALGORITHM_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

/**
.Function.em:
..summary:Represents the EM algorithm as used by MEME.
..cat:Motif Finding
..signature:em(profile,dataset_start,t,l,oops_model)
..signature:em(profile,dataset_start,t,l,gamma,zoops_model)
..signature:em(profile,dataset_start,t,l,lambda,tcm_model)
..param.profile:The  $profile$ object which is a set of @Class.FrequencyDistribution|frequency distributions@.
...remarks:$profile$ is of type $String<$ @Class.FrequencyDistribution@ $>$
..param.dataset_start:An iterator pointing to the first input sequence of a given dataset.
...type:Concept.Iterator
..param.t:The number of input sequences.
..param.l:The size of the motif.
..param.oops_model:The oops_model object.
...type:Tag.OOPS
..param.zoops_model:The zoops_model object.
...type:Tag.ZOOPS
..param.tcm_model:The tcm_model object.
...type:Tag.TCM
..param.gamma:The probability of sequence having a motif occurence.
..param.lambda:The probability of starting a motif occurence 
...remarks:$lambda$ is calculated by dividing $gamma$ by the length of the corresponding sequence.
..remarks:This version of EM is used in the MEME program of Bailey and Elkan. It is a Bayesian
          variant of the basic EM which allows multiple occurrences of a motif in any sequence and can 
		  therefore be performed on sequences of one of the model types @Tag.OOPS@, @Tag.ZOOPS@ and 
		  @Tag.TCM@. We use the EM algorithm of MEME for the refinement step of PROJECTION.
*/

//////////////////////////////////////////////////////////////////////////////
//	OOPS model
//////////////////////////////////////////////////////////////////////////////

template<typename TProfile, typename TIter, typename TType>
double
em(TProfile & profile, 
   TIter dataset_start, 
   TType const & t,
   TType const & l,
   OOPS const & oops)
{
	// matrix w - allocate space (memory)
	TType row_size = t;
	double ** matrix_w = new double*[row_size];
	for(TType pos=0; pos<row_size; ++pos)
	{
		TType col_size = length(*(dataset_start+pos))-l+1;
		matrix_w[pos] = new double[col_size];
	}

	// E-step: compute matrix z and joint log likelihood
	double log_likelihood = static_cast<double>(0);
	_computeEStep(matrix_w,log_likelihood,profile,dataset_start,t,l,oops);

	// M-step: refine profile
	_computeMStep(profile,dataset_start,matrix_w,t,l,oops);

	return log_likelihood;
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeEStep:
..summary:Represents the E-step of the EM algorithm for OOPS models.
..cat:Motif Finding
..signature:_computeEStep(matrix_w,joint_log_likelihood,profile,dataset_start,t,l,oops_model)
..param.matrix_w:The matrix_w object.
...remarks:w(i,j) is the probability that the motif pattern 
            starts at position j in the i-th input sequence.
..param.joint_log_likelihood:The joint log-likelihood.
..param.profile:The profile object which is a set of frequency distributions.
...type:Class.String
....signature:String<TFrequencyDistribution>
..param.dataset_start:n iterator pointing to the first input sequence of a given dataset.
..param.t:The number of input sequences.
..param.l:The size of the motif.
..param.oops_model:The oops_model object.
...type:Tag.OOPS
..remarks:The joint log likelihood is not computed in the M-step as MEME does, but rather in the E-step
          of the algorithm.
*/

template<typename TMatrix, typename TProfile, typename TIter, typename TType>
void
_computeEStep(TMatrix & matrix_w, 
	double & joint_log_likelihood,
	TProfile & profile, 
	TIter dataset_start, //start iterator of dataset 
	TType const & t,
	TType const & l,
	OOPS const & oops)
{
	typedef typename Value<TProfile>::Type TFrequencyDist;
	typedef typename Position<TProfile>::Type TPos;
	TProfile log_profile = profile;
	for(TPos i=0; i<length(log_profile); ++i)
	{
		logarithmize(log_profile[i]); 
	}

	//compute matrix w
	for(TType i=0; i<t; ++i)
	{
		TType m = length(*(dataset_start+i))-l+1;
		double motif_prob = static_cast<double>(0);
		TFrequencyDist letter_freq;
		for(TType j=0; j<m; ++j)
		{
			double sum_of_log_probs = static_cast<double>(0);
			if(j==0)
			{
				//motif region
				for(TPos h=1; h<=l; ++h)
				{
					motif_prob +=
						log_profile[h][(int)(*(dataset_start+i))[j+h-1]];
				}
				sum_of_log_probs += motif_prob;

				//background region
				for(TType h=0; h<m-1; ++h)
				{
					sum_of_log_probs +=
						log_profile[0][(int)(*(dataset_start+i))[j+l+h]];
				}
			}
			else
			{
				sum_of_log_probs = matrix_w[i][j-1]-motif_prob;
				motif_prob = 0;
				//motif region
				for(TPos h=1; h<=l; ++h)
				{
					motif_prob +=
						log_profile[h][(int)(*(dataset_start+i))[j+h-1]];
				}
				sum_of_log_probs += motif_prob;
				sum_of_log_probs -= log_profile[0][(int)(*(dataset_start+i))[j+l-1]];
				sum_of_log_probs += log_profile[0][(int)(*(dataset_start+i))[j-1]];
			}
			matrix_w[i][j] = sum_of_log_probs;
		}
	
		double log_of_sums = matrix_w[i][0];
		for(TType j=1; j<m; ++j)
		{
			if( (matrix_w[i][j]-log_of_sums)>DBL_MIN_EXP )
			{
				//log_of_sums = 
				//	log_of_sums+log(static_cast<double>(1)+exp(matrix_w[i][j]-log_of_sums));
				log_of_sums += 
					log(static_cast<double>(1)+exp(matrix_w[i][j]-log_of_sums));
			}
		}

		for(TType j=0; j<m; ++j)
		{
			double prob = matrix_w[i][j];
			matrix_w[i][j] = matrix_w[i][j]-log_of_sums;
			joint_log_likelihood += matrix_w[i][j]*prob;
			//exp(w_ij)=e^(w_ij) -> e is the base of the natural logarithm
			matrix_w[i][j] = exp(matrix_w[i][j]);
		}
		//sum_{i=1}^{t}log(1/m_i)
		joint_log_likelihood += 
			log(static_cast<double>(1))-log(static_cast<double>(m));
	}
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeStep_M:
..summary:Represents the M-step of the EM algorithm for OOPS models.
          Refines the background and motif component.
..cat:Motif Finding
..signature:_computeStep_M(profile,dataset_start,matrix_w,t,l,oops_model)
..param.profile:The profile object which is a set of frequency distributions.
...type:Class.String
....signature:String<TFrequencyDistribution>
..param.dataset_start:n iterator pointing to the first input sequence of a given dataset.
..param.matrix_w:The matrix_w object.
...remarks:w(i,j) is the probability that the motif pattern 
            starts at position j in the i-th input sequence.
..param.t:The number of input sequences.
..param.l:The size of the motif.
..param.oops_model:The oops_model object.
...type:Tag.OOPS
..remarks:The joint log likelihood is computed in the E-step of the algorithm.
*/

template<typename TProfile, typename TIter, typename TMatrix, typename TType>
void
_computeMStep(TProfile & profile, 
	TIter dataset_start,
	TMatrix const & matrix_w,
	TType const & t,
	TType const & l,
	OOPS const & oops)
{
	typedef typename Value<TProfile>::Type TFrequencyDist;
	typedef typename Value<TFrequencyDist>::Type TValue;
	typedef typename Position<TProfile>::Type TPos;

	// total_counts_of_letters is used for the computing of the background frequency
	TFrequencyDist total_counts_of_letters; //<-c
	absFreqOfLettersInSetOfSeqs(total_counts_of_letters,dataset_start+0, dataset_start+t);

	// refine motif component
	for(TPos h=1; h<=l; ++h)
	{
		TFrequencyDist letter_freq;//<-p'_h
		for(TType i=0; i<t; ++i)
		{
			TType m = length(*(dataset_start+i))-l+1;
			for(TType j=0; j<m; ++j)
			{
				TFrequencyDist frequency;
	            frequency[(int) (*(dataset_start+i))[j+h-1]] = matrix_w[i][j];
				letter_freq += frequency;
			}
		}
		profile[h] = letter_freq;
		total_counts_of_letters-=letter_freq;//c->p'_0
	}

	// refine background component
	profile[0] = total_counts_of_letters;
	
	// addPseudocount (if necessary) & normalize
	double epsilon = 0.1;
	normalize(profile, Pseudocount<TValue, CMode>(epsilon));

	// matrix w - deallocate space (memory)
	for(TType pos=0; pos<t; ++pos)
	{
		delete[] matrix_w[pos];
	}
	delete[] matrix_w;
}

//////////////////////////////////////////////////////////////////////////////
//	ZOOPS model
//////////////////////////////////////////////////////////////////////////////

// gamma=(1/t)*sum(i,1,t,Qi), Qi=sum(j,1,m,zij) (i=1,...,t)

template<typename TProfile, typename TIter, typename TType>
double
em(TProfile & profile, 
   TIter dataset_start, 
   TType const & t,
   TType const & l,
   double & gamma,
   ZOOPS const & zoops)
{
	// matrix w - allocate space (memory)
	TType row_size = t;
	double ** matrix_w = new double*[row_size];
	for(TType pos=0; pos<row_size; ++pos)
	{
		TType col_size = length(*(dataset_start+pos))-l+1;
		matrix_w[pos] = new double[col_size];
	}

	// compute matrix w and joint log likelihood, refine gamma
	double log_likelihood = static_cast<double>(0);
	_computeEStep(matrix_w,log_likelihood,profile,dataset_start,gamma,t,l,zoops);

	// refine profile
	_computeMStep(profile,dataset_start,matrix_w,t,l,zoops);

	return log_likelihood;
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeEStep:
..summary:Represents the E-step of the EM algorithm for ZOOPS models.
..cat:Motif Finding
..signature:_computeEStep(matrix_w,joint_log_likelihood,profile,dataset_start,gamma,t,l,zoops_model)
..param.matrix_w:The matrix_w object.
...remarks:w(i,j) is the probability that the motif pattern 
            starts at position j in the i-th input sequence.
..param.joint_log_likelihood:The joint log-likelihood.
..param.profile:The profile object which is a set of frequency distributions.
...type:Class.String
....signature:String<TFrequencyDistribution>
..param.dataset_start:n iterator pointing to the first input sequence of a given dataset.
..param.gamma:The probability of sequence having a motif occurence.
...type:$double$
..param.t:The number of input sequences.
..param.l:The size of the motif.
..param.zoops_model:The oops_model object.
...type:Tag.ZOOPS
..remarks:The joint log likelihood is not computed in the M-step as MEME does, but rather in the E-step
          of the algorithm.
..remarks:The parameter gamma is reestimated during the E-step to save calculation time.
*/

template<typename TMatrix, typename TProfile, typename TIter, typename TType>
void
_computeEStep(TMatrix & matrix_w,
			  double & joint_log_likelihood,
			  TProfile & profile,
			  TIter dataset_start,
			  double & gamma,
			  TType const & t,
			  TType const & l,
			  ZOOPS const & zoops)
{
	typedef typename Value<TProfile>::Type TFrequencyDist;
	typedef typename Position<TProfile>::Type TPos;
	double lambda_i = static_cast<double>(0); 
	double sum_of_Q_i = static_cast<double>(0);

	TProfile log_profile = profile;
	for(TPos i=0; i<length(log_profile); ++i)
	{
		logarithmize(log_profile[i]); 
	}

	//compute matrix w
	for(TType i=0; i<t; ++i)
	{
		TType m = length(*(dataset_start+i))-l+1;
		lambda_i = gamma/static_cast<double>(m);
		TFrequencyDist t_i;
		absFreqOfLettersInSeq(t_i, begin(*(dataset_start+i)), end(*(dataset_start+i)));
		double log_of_seq_prob_having_no_motif = 
			std::inner_product(begin(t_i), end(t_i), 
			                   begin(log_profile[0]), static_cast<double>(0));

		double motif_prob = static_cast<double>(0);
		//double sum_of_Q_i = static_cast<double>(0);
		TFrequencyDist letter_freq;
		for(TType j=0; j<m; ++j)
		{
			double sum_of_log_probs = static_cast<double>(0);
			if(j==0)
			{
				//motif region
				for(TPos h=1; h<=l; ++h)
				{
					motif_prob +=
						log_profile[h][(int)(*(dataset_start+i))[j+h-1]];
				}
				sum_of_log_probs += motif_prob;

				//background region
				for(TType h=0; h<m-1; ++h)
				{
					sum_of_log_probs +=
						log_profile[0][(int)(*(dataset_start+i))[j+l+h]];
				}
			}
			else
			{
				sum_of_log_probs = matrix_w[i][j-1]-motif_prob;
				motif_prob = 0;
				//motif region
				for(TPos h=1; h<=l; ++h)
				{
					motif_prob +=
						log_profile[h][(int)(*(dataset_start+i))[j+h-1]];
				}
				sum_of_log_probs += motif_prob;
				sum_of_log_probs -= log_profile[0][(int)(*(dataset_start+i))[j+l-1]];
				sum_of_log_probs += log_profile[0][(int)(*(dataset_start+i))[j-1]];
			}
			matrix_w[i][j] = sum_of_log_probs;
		}

		//compute subtrahend - log(x+y) = log(x)+log(1+exp(log(y)-log(x)))
		double log_of_sums = matrix_w[i][0];
		for(TType j=1; j<m; ++j)
		{
			if( (matrix_w[i][j]-log_of_sums)>DBL_MIN_EXP )
			{
				log_of_sums +=
					log(static_cast<double>(1)+exp(matrix_w[i][j]-log_of_sums));
			}
		}

		for(TType j=0; j<m; ++j)
		{
			//double prob = exp(matrix_w[i][j]);
			double prob = matrix_w[i][j];
			matrix_w[i][j] += log(lambda_i);
			double exponent =
				log_of_sums+log(lambda_i)-log_of_seq_prob_having_no_motif
			   -log(static_cast<double>(1)-gamma);
			if( exponent>DBL_MIN_EXP )
			{
				matrix_w[i][j] -=
					(log_of_seq_prob_having_no_motif+log(static_cast<double>(1)-gamma)
					+log(1+exp(exponent)));
			}
			else
			{
				matrix_w[i][j] -= 
					(log_of_seq_prob_having_no_motif+log(static_cast<double>(1)-gamma));
			}
			matrix_w[i][j] = exp(matrix_w[i][j]);
			joint_log_likelihood += matrix_w[i][j]*prob;
		}

		double Q_i = 
			std::accumulate(matrix_w[i]+0,matrix_w[i]+m,static_cast<double>(0));
		sum_of_Q_i+=Q_i;
		joint_log_likelihood +=
			((1-Q_i)*log_of_seq_prob_having_no_motif
		   +(Q_i*log(gamma/static_cast<double>(m))
		   +(1-Q_i)*log(static_cast<double>(1)-gamma)));
		/*joint_log_likelihood +=
			((1-Q_i)*exp(log_of_seq_prob_having_no_motif)
		   +(Q_i*log(gamma/static_cast<double>(m))
		   +(1-Q_i)*log(static_cast<double>(1)-gamma)));*/
	}
	// refine value of gamma
	gamma =  sum_of_Q_i*(static_cast<double>(1)/static_cast<double>(t));
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeStep_M:
..summary:Represents the M-step of the EM algorithm for ZOOPS models.
          Refines the background and motif component.
..cat:Motif Finding
..signature:_computeStep_M(profile,dataset_start,matrix_w,t,l,zoops_model)
..param.profile:The profile object which is a set of frequency distributions.
...type:Class.String
....signature:String<TFrequencyDistribution>
..param.dataset_start:n iterator pointing to the first input sequence of a given dataset.
..param.matrix_w:The matrix_w object.
...remarks:w(i,j) is the probability that the motif pattern 
            starts at position j in the i-th input sequence.
..param.t:The number of input sequences.
..param.l:The size of the motif.
..param.zoops_model:The zoops_model object.
...type:Tag.ZOOPS
..remarks:The joint log likelihood is computed in the E-step of the algorithm.
*/

template<typename TProfile, typename TIter, typename TMatrix, typename TType>
void
_computeMStep(TProfile & profile,
			  TIter dataset_start,
			  TMatrix const & matrix_w,
			  TType const & t,
	          TType const & l,
			  ZOOPS const & zoops)
{
	typedef typename Value<TProfile>::Type TFrequencyDist;
	typedef typename Value<TFrequencyDist>::Type TValue;
	typedef typename Position<TProfile>::Type TPos;

    // total_counts_of_letters is used for the refinement of the background frequency
	TFrequencyDist total_counts_of_letters;
	absFreqOfLettersInSetOfSeqs(total_counts_of_letters, dataset_start+0, dataset_start+t); 

	// refine motif component 
	for(TPos h=1; h<=l; ++h)
	{
		TFrequencyDist letter_freq;
		for(TType i=0; i<t; ++i)
		{
			TType m = length(*(dataset_start+i))-l+1;
			for(TType j=0; j<m; ++j)
			{
				TFrequencyDist frequency;
	            frequency[(int)(*(dataset_start+i))[j+h-1]] = matrix_w[i][j];
				letter_freq += frequency;
			}
		}
		profile[h] = letter_freq;
		total_counts_of_letters-=letter_freq; // refine background component
	}
	profile[0] = total_counts_of_letters;

	// addPseudocount (if necessary) & normalize
	double epsilon = 0.1;
	normalize(profile, Pseudocount<TValue, CMode>(epsilon));

	// matrix w - deallocate space (memory)
	for(TType pos=0; pos<t; ++pos)
	{
		delete[] matrix_w[pos];
	}
	delete[] matrix_w;
}

//////////////////////////////////////////////////////////////////////////////
//	TCM model
//////////////////////////////////////////////////////////////////////////////

/*
..remarks:Dataset X is converted into a new pseudo-dataset consisting of all the width-l
          overlapping subsequences by running a window of width-l along each sequence Xi 
          and writing down the string contained in the window. Then the new dataset is 
          modeled as though it were generated by a two-component mixture model where each 
          component generates a string of width-l.
..remarks:Xij:=a width-l pseudo-sequence (=[Xij,Xij+1,...,Xij+l-1])
..remarks:lambda=gamma/m
*/

template<typename TProfile, typename TIter, typename TType>
double
em(TProfile & profile,
   TIter dataset_start,
   TType const & t,
   TType const & l,
   double & lambda,
   TCM const & tcm)
{
	// matrix z - allocate space (memory)
	TType row_size = t;
	double ** matrix_w = new double*[row_size];
	for(TType pos=0; pos<row_size; ++pos)
	{
		TType col_size = length(*(dataset_start+pos))-l+1;
		matrix_w[pos] = new double[col_size];
	}

	// E-step: compute matrix w and the joint log likelihood
	double log_likelihood = static_cast<double>(0);
	_computeEStep(matrix_w,log_likelihood,profile,dataset_start,lambda,t,l,tcm);

	// M-step: refine profile
	_computeMStep(profile,dataset_start,matrix_w,t,l,tcm);

	return log_likelihood;
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeEStep:
..summary:Represents the E-step of the EM algorithm for TCM models.
..cat:Motif Finding
..signature:_computeEStep(matrix_w,joint_log_likelihood,profile,dataset_start,lambda,t,l,tcm_model)
..param.matrix_w:The matrix_w object.
...remarks:w(i,j) is the probability that the motif pattern 
            starts at position j in the i-th input sequence.
..param.joint_log_likelihood:The joint log-likelihood.
..param.profile:The profile object which is a set of frequency distributions.
...type:Class.String
....signature:String<TFrequencyDistribution>
..param.dataset_start:n iterator pointing to the first input sequence of a given dataset.
..param.lambda:The probability of starting a motif occurence 
...type:$double$
..param.t:The number of input sequences.
..param.l:The size of the motif.
..param.tcm_model:The tcm_model object.
...type:Tag.TCM
..remarks:The joint log likelihood is not computed in the M-step as MEME does, but rather in the E-step
          of the algorithm.
..remarks:The parameter lambda is reestimated during the E-step to save calculation time.
*/

template<typename TMatrix, typename TProfile, typename TIter, typename TType>
void
_computeEStep(TMatrix & matrix_w,
			  double & joint_log_likelihood,
			  TProfile & profile, 
			  TIter dataset_start,
			  double & lambda,
			  TType const & t,
			  TType const & l,
			  TCM const & tcm)
{
	typedef typename Value<TProfile>::Type TFrequencyDist;
	typedef typename Position<TProfile>::Type TPos;
	TProfile log_profile = profile;
	for(TPos i=0; i<length(log_profile); ++i)
	{
		logarithmize(log_profile[i]); 
	}
	double prev_lambda = lambda;
	lambda = static_cast<double>(0);

	//compute matrix z
	for(TType i=0; i<t; ++i)
	{
		TType m = length(*(dataset_start+i))-l+1;
		double prev_log_prob_given_theta0 = static_cast<double>(0);
		double Q_i = static_cast<double>(0);
		for(TType j=0; j<m; ++j)
		{
			double log_prob_given_theta1 = static_cast<double>(0);
			double log_prob_given_theta0 = static_cast<double>(0);
			if(j==0)
			{
				for(TPos h=1; h<=l; ++h)
				{
					int letter = (int) (*(dataset_start+i))[j+h-1];
					//assume that subsequence was generated by the motif component
					log_prob_given_theta1 +=
						log_profile[h][letter];

					//assume that subsequence was generated by the background component
					log_prob_given_theta0 +=
						log_profile[0][letter];
				}
				prev_log_prob_given_theta0 = log_prob_given_theta0;
			}
			else
			{
				for(TPos h=1; h<=l; ++h)
				{
					int letter = (int) (*(dataset_start+i))[j+h-1];
					//assume that subsequence was generated by the motif component
					log_prob_given_theta1 +=
						log_profile[h][letter];
				}
				log_prob_given_theta0 = 
					prev_log_prob_given_theta0-log_profile[0][(int) (*(dataset_start+i))[j-1]]
				  +log_profile[0][(int) (*(dataset_start+i))[j+l-1]];
				prev_log_prob_given_theta0 = log_prob_given_theta0;
			}
			double minuend = log_prob_given_theta1+log(prev_lambda);
			double subtrahend = static_cast<double>(0);
			double exponent = minuend-log_prob_given_theta0-log(static_cast<double>(1)-prev_lambda);
			if(exponent>DBL_MIN_EXP)
			{
				subtrahend =
					log_prob_given_theta0+log(static_cast<double>(1)-prev_lambda)
				   +log(static_cast<double>(1)+exp(exponent));
			}
			else
			{
				subtrahend =
					log_prob_given_theta0+log(static_cast<double>(1)-prev_lambda);
			}

			matrix_w[i][j] = exp(minuend-subtrahend);
			Q_i += matrix_w[i][j];
			joint_log_likelihood +=
				((static_cast<double>(1)-matrix_w[i][j])*log_prob_given_theta0
			    +matrix_w[i][j]*log_prob_given_theta1
			    +(static_cast<double>(1)-matrix_w[i][j])*log(static_cast<double>(1)-prev_lambda)
			    +matrix_w[i][j]*log(prev_lambda));
		}
		lambda += Q_i/static_cast<double>(m);

	}
	lambda = lambda/static_cast<double>(t);

	// apply a smoothing step to reduce the degree to which any two overlapping 
	// subsequences can both be assigned to the motif component
	_smoothingStep(matrix_w, dataset_start, t, l);
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeStep_M:
..summary:Represents the M-step of the EM algorithm for TCM models.
          Refines the background and motif component.
..cat:Motif Finding
..signature:_computeStep_M(profile,dataset_start,matrix_w,t,l,tcm_model)
..param.profile:The profile object which is a set of frequency distributions.
...type:Class.String
....signature:String<TFrequencyDistribution>
..param.dataset_start:n iterator pointing to the first input sequence of a given dataset.
..param.matrix_w:The matrix_w object.
...remarks:w(i,j) is the probability that the motif pattern 
            starts at position j in the i-th input sequence.
..param.t:The number of input sequences.
..param.l:The size of the motif.
..param.tcm_model:The tcm_model object.
...type:Tag.TCM
..remarks:The joint log likelihood is computed in the E-step of the algorithm.
*/

template<typename TProfile, typename TIter, typename TMatrix, typename TType>
void
_computeMStep(TProfile & profile,
			  TIter dataset_start,
			  TMatrix const & matrix_w,
			  TType const & t,
	          TType const & l,
			  TCM const & tcm)
{
	typedef typename Value<TProfile>::Type TFrequencyDist;
	typedef typename Value<TFrequencyDist>::Type TValue;
	typedef typename Position<TProfile>::Type TPos;

    // total_counts_of_letters is used for the refinement of the background frequency
	TFrequencyDist total_counts_of_letters;
	absFreqOfLettersInSetOfSeqs(total_counts_of_letters, dataset_start+0, dataset_start+t); 

	// refine motif component 
	for(TPos h=1; h<=l; ++h)
	{
		TFrequencyDist letter_freq;
		for(TType i=0; i<t; ++i)
		{
			TType m = length(*(dataset_start+i))-l+1;
			for(TType j=0; j<m; ++j)
			{
				TFrequencyDist frequency;
	            frequency[(int)(*(dataset_start+i))[j+h-1]] = matrix_w[i][j];
				letter_freq += frequency;
			}
		}
		profile[h] = letter_freq;
		total_counts_of_letters-=letter_freq; // refine background component
	}

	// refine background component
	TFrequencyDist freq;
	for(TType i=0; i<t; ++i)
	{
		TType m = length(*(dataset_start+i))-l+1;
		for(TType j=0; j<m; ++j)
		{
			for(TType h=1; h<=l; ++h)
			{
	            freq[(int)(*(dataset_start+i))[j+h-1]] += (static_cast<double>(1)-matrix_w[i][j]);
			}
		}
	}
	profile[0] = freq;

	// addPseudocount (if necessary) & normalize
	double epsilon = 0.1;
	normalize(profile, Pseudocount<TValue, CMode>(epsilon));

	// matrix w - deallocate space (memory)
	for(TType pos=0; pos<t; ++pos)
	{
		delete[] matrix_w[pos];
	}
	delete[] matrix_w;
}

//////////////////////////////////////////////////////////////////////////////
//Subfunctions
//////////////////////////////////////////////////////////////////////////////

/*
.Function._smoothingStep:
..summary:Applies a smoothing step to reduce the degree to which any two overlapping 
          subsequences can both be assigned to the motif component.
         (We do not want the model to predict that two overlapping substrings are both 
          motif occurences.)
..cat:Motif Finding
..signature:_smoothingStep(matrix_w,dataset_start,t,l)
..param.matrix_w:The matrix_w object.
...remarks:w(i,j) is the probability that the motif pattern 
            starts at position j in the i-th input sequence.
..param.dataset_start:An iterator pointing to the first input sequence of a given dataset.
..param.t:The number of input sequences.
..param.l:The size of the motif.
..remarks:Function is used by the EM algorithm for TCM models.
*/

template<typename TMatrix, typename TIter, typename TType>
void
_smoothingStep(TMatrix & matrix_w,
			   TIter dataset_start,
			   TType const & t,
			   TType const & l)
{
	for(TType i=0; i<t; ++i)
	{
		TType m = length(*(dataset_start+i))-l+1;
		for(TType k=0; k<l; ++k)
		{
			for(TType j=k; j<m; ++j)
			{
				// w_sum:=sum of Zij in the current window of size l (l:=l)
				double w_sum = 
					static_cast<double>(std::accumulate(matrix_w[i]+j, 
										matrix_w[i]+j+l, 
										static_cast<double>(0)));

				// w_max:=the largest Zij in the current window
				double w_max =
					static_cast<double>(*std::max_element(matrix_w[i]+j, 
														  matrix_w[i]+j+l));

				if(w_sum>1.0)
				{
					double factor = (static_cast<double>(1)-w_max)/(w_sum-w_max);
					for(TType h=j; h<j+l; ++h)
					{
						if(matrix_w[i][h]!=w_max)
						{
							matrix_w[i][h] *= factor;
						}
					}
				}
			}
		}
	}
}

} // SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...