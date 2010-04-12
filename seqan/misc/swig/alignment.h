#include <cstring>
#include <seqan.h>
#include <seqan/align.h>


namespace seqan
{

	#ifdef SWIG

	%typemap(in,numinputs=0)int& scoreMatrix(int temp)
	{
		$1=&temp;
	}
	%typemap(argout)int& scoreMatrix
	{    
		PyObject *o, *o2, *o3;
    o = PyInt_FromLong(*$1);
    if ((!$result) || ($result == Py_None))
		{
			$result = o;
    }
		else
		{
      if(!PyTuple_Check($result))
			 {
            PyObject *o2 = $result;
            $result = PyTuple_New(1);
            PyTuple_SetItem($result,0,o2);
        }
				o3 = PyTuple_New(1);
        PyTuple_SetItem(o3,0,o);
        o2 = $result;
        $result = PySequence_Concat(o2,o3);
        Py_DECREF(o2);
        Py_DECREF(o3);
    }
	}

	%typemap(in,numinputs=0) char** str1(char* temp) 
	{
		$1 = (char **) malloc((strlen(temp)+1)*sizeof(char *));
	}
	%typemap(in,numinputs=0) char** str2(char* temp) 
	{
		  $1 = (char **) malloc((strlen(temp)+1)*sizeof(char *));
	}

	%typemap(argout)char** str1
	{
		$result=SWIG_Python_AppendOutput($result, PyString_FromStringAndSize(*$1,strlen(*$1)));
	}
	%typemap(argout)char** str2
	{
		$result=SWIG_Python_AppendOutput($result,PyString_FromString(*$1));
	}

	%typemap(freearg) char ** 
	{
		free((char *) $1);
	}

	%typemap(in,numinputs=0) int** scoreMatrix(int* temp)
	{
		  $1 = (int **) malloc((sizeof(temp)+1)*sizeof(int *));
	}
	%typemap(in,numinputs=0) int& aSize(int temp)
	{
		  $1 = &temp;
	}

	%typemap(argout) (int** scoreMatrix, int& aSize)
	{
  	int i;
  	$result = PyList_New(*$2);
  	for (i = 0; i < *$2; i++)
	 {
    	PyObject *o =PyInt_FromLong(*(*$1+i));
    	PyList_SetItem($result,i,o);
  	}
	}


	%typemap(freearg) int ** 
	{
		free((int *) $1);
	}


	%typemap(in) (int len, int *value)
	{
		int i;
		if (!PyList_Check($input))
		{
		  PyErr_SetString(PyExc_ValueError, "Expecting a list");
		  return NULL;
		}
		$1 = PyList_Size($input);
		$2 = (int*) malloc(($1+1)*sizeof(int));
		for (i = 0; i < $1; i++)
		{
		  PyObject *s = PyList_GetItem($input,i);
		  if (!PyInt_Check(s))
			{
		      free($2);
		      PyErr_SetString(PyExc_ValueError, "List items must be integers");
		      return NULL;
		  }
		  $2[i] = (int)PyInt_AS_LONG(s);
		}
		$2[i] = 0;
	}

	%typemap(freearg) (int len, int *value) 
	{
		 if ($2) free($2);
	}

	%typemap(default) unsigned int matrix
	{
		 $1 = 0;
	}
	%typemap(default) unsigned int dia1
	{
		 $1 = 0;
	}
	%typemap(default) unsigned int dia2 
	{
		 $1 = 0;
	}


	%define DOCSTRING_
	"""getAminoAcidScoreMatrix(\"X\")->[long integer]; \n Returns a list containing long integer, which represent the specific scoring matrices X.\n A selectable value for X could be:\n Blosum30 \n Blosum45 \n Blosum62 \n Blosum80 \n Pam40 \n Pam120 \n Pam200 \n Pam250 \n Vtml200 \n Example: scoreMatrix=getAminoAcidScoreMatrix(\"Blosum30\");"""
	%enddef

	%define DOCSTRING2
	"printAlignment(AlignmentObject)->[Strings];\n Returns a list containing three strings.\n The first string illustrate the alignment for both sequence in one representation.\n List entry two and three contains the alignment sequence separately.\n"
	%enddef

	%feature("autodoc", DOCSTRING_) getAminoAcidScoreMatrix;

	%feature("autodoc", DOCSTRING2) printAlignment;

	#endif

	void getAminoAcidScoreMatrix(char* str,int** scoreMatrix, int& aSize)
	{
		if(std::strcmp(str,"Blosum30")==0)
		{
			typedef Score<int, ScoreMatrix<AminoAcid, _Blosum30> > TScore;
			int const* matrix=_ScoringMatrixData<int ,AminoAcid, _Blosum30 >::getData();
			*scoreMatrix = new int[TScore::TAB_SIZE];
			arrayCopy(matrix, matrix + TScore::TAB_SIZE, *scoreMatrix);
			aSize=TScore::TAB_SIZE;
		}

		else if (std::strcmp(str,"Blosum45")==0)
		{
			typedef Score<int, ScoreMatrix<AminoAcid, _Blosum45> > TScore;
			int const* matrix=_ScoringMatrixData<int ,AminoAcid, _Blosum45 >::getData();
			*scoreMatrix = new int[TScore::TAB_SIZE];
			arrayCopy(matrix, matrix + TScore::TAB_SIZE, *scoreMatrix);
			aSize=TScore::TAB_SIZE;
		}
		else if(std::strcmp(str,"Blosum62")==0)
		{
			typedef Score<int, ScoreMatrix<AminoAcid, _Blosum62> > TScore;
			int const* matrix=_ScoringMatrixData<int ,AminoAcid, _Blosum62 >::getData();
			*scoreMatrix = new int[TScore::TAB_SIZE];
			arrayCopy(matrix, matrix + TScore::TAB_SIZE, *scoreMatrix);
			aSize=TScore::TAB_SIZE;
		}
		else if(std::strcmp(str,"Blosum80")==0)
		{

			typedef Score<int, ScoreMatrix<AminoAcid, _Blosum80> > TScore;
			int const* matrix=_ScoringMatrixData<int ,AminoAcid, _Blosum80 >::getData();
			*scoreMatrix = new int[TScore::TAB_SIZE];
			arrayCopy(matrix, matrix + TScore::TAB_SIZE, *scoreMatrix);
			aSize=TScore::TAB_SIZE;

		}
		else if(std::strcmp(str,"Pam40")==0)
		{
	
			typedef Score<int, ScoreMatrix<AminoAcid, _Pam40> > TScore;
			int const* matrix=_ScoringMatrixData<int ,AminoAcid, _Pam40 >::getData();
			*scoreMatrix = new int[TScore::TAB_SIZE];
			arrayCopy(matrix, matrix + TScore::TAB_SIZE, *scoreMatrix);
			aSize=TScore::TAB_SIZE;
		}
		else if(std::strcmp(str,"Pam120")==0)
		{
			typedef Score<int, ScoreMatrix<AminoAcid, _Pam120> > TScore;
			int const* matrix=_ScoringMatrixData<int ,AminoAcid, _Pam120 >::getData();
			*scoreMatrix = new int[TScore::TAB_SIZE];
			arrayCopy(matrix, matrix + TScore::TAB_SIZE, *scoreMatrix);
			aSize=TScore::TAB_SIZE;
		}
		else if(std::strcmp(str,"Pam200")==0)
		{
			typedef Score<int, ScoreMatrix<AminoAcid, _Pam200> > TScore;
			int const* matrix=_ScoringMatrixData<int ,AminoAcid, _Pam200 >::getData();
			*scoreMatrix = new int[TScore::TAB_SIZE];
			arrayCopy(matrix, matrix + TScore::TAB_SIZE, *scoreMatrix);
			aSize=TScore::TAB_SIZE;
		}
		else if(std::strcmp(str,"Pam250")==0)
		{
			typedef Score<int, ScoreMatrix<AminoAcid, _Pam250> > TScore;
			int const* matrix=_ScoringMatrixData<int ,AminoAcid, _Pam250 >::getData();
			*scoreMatrix = new int[TScore::TAB_SIZE];
			arrayCopy(matrix, matrix + TScore::TAB_SIZE, *scoreMatrix);
			aSize=TScore::TAB_SIZE;
		}
		else if(std::strcmp(str,"Vtml200")==0)
		{
			typedef Score<int, ScoreMatrix<AminoAcid, _Vtml200> > TScore;
			int const* matrix=_ScoringMatrixData<int ,AminoAcid, _Vtml200 >::getData();
			*scoreMatrix = new int[TScore::TAB_SIZE];
			arrayCopy(matrix, matrix + TScore::TAB_SIZE, *scoreMatrix);
			aSize=TScore::TAB_SIZE;
		}
	}



	template<typename Type,typename TSpec, typename TAlignConfig>
	int _alignment(seqan::Align<seqan::String<Type> >*align,seqan::Score<int,TSpec>& score,TAlignConfig& ac,unsigned int dia1, unsigned int dia2,char *algo)
	{
		if(std::strcmp(algo,"NeedlemanWunsch")==0) return globalAlignment(*align,score,ac,NeedlemanWunsch());
		else if(std::strcmp(algo,"Gotoh")==0) return globalAlignment(*align,score,ac,Gotoh());
		else if(std::strcmp(algo,"BandedNeedlemanWunsch")==0) return  globalAlignment(stringSet(*align),score,ac,dia1,dia2,BandedNeedlemanWunsch());
		else if(std::strcmp(algo,"BandedGotoh")==0) return globalAlignment(stringSet(*align),score,ac,dia1,dia2,BandedGotoh());
	}
	 

	template<typename Type, typename TSpec>
	int _alignment(seqan::Align<seqan::String<Type> >*align,seqan::Score<int,TSpec>& score,unsigned int ac,char * algo,unsigned int dia1, unsigned int dia2)
	{
	int alignmentScore=0;
			switch(ac)
			{
				case 0: {
									AlignConfig<false,false,false,false>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);
									break;
								}
				case 1: {
									AlignConfig<false,false,false,true>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);
									break;
								}
				case 2: {
									AlignConfig<false,false,true,false>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);
									break;
								}
				case 3: {
									AlignConfig<false,false,true,true>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 4: {
									AlignConfig<false,true,false,false>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 5: {
									AlignConfig<false,true,false,true>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 6: {
									AlignConfig<false,true,true,false>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 7: {
									AlignConfig<false,true,true,true>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 8: {
									AlignConfig<true,false,false,false>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 9: {
									AlignConfig<true,false,false,true>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 10: {
									AlignConfig<true,false,true,false>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 11: {
									AlignConfig<true,false,true,true>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 12: {
									AlignConfig<true,true,false,false>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 13: {
									AlignConfig<true,true,false,true>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 14: {
									AlignConfig<true,true,true,false>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 15: {
									AlignConfig<true,true,true,true>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
			}

	return alignmentScore;
	}

//allgemeine Wrapper funktion fürs Alignment
	template< typename Type >
	seqan::Align<seqan::String<Type> >* alignSequence(char* seq1, char*seq2, int len, int *value, char* algo, unsigned int matrix, unsigned int dia1, unsigned int dia2, int& scoreMatrix)
	{
		Align<String<Type> >* align = new Align<String<Type> >();
		
		String<Type> _seq1=seq1; //ist nötig!! <-python strings are immutable!! swig does not change this behaviour
		String<Type> _seq2=seq2;
		
		resize(rows(*align), 2);
		assignSource(row(*align, 0), _seq1);
		assignSource(row(*align, 1), _seq2);

		if(len >=4 && ValueSize<Type>::VALUE * ValueSize<Type>::VALUE == len && std::strcmp(algo,"MyersHirschberg")!=0 && std::strcmp(algo,"Hirschberg")!=0 && std::strcmp(algo,"SmithWaterman")!=0)
		{
			Score<int, ScoreMatrix<Type> > score;
 			arrayCopy(value, value + len, score.data_tab);
			scoreMatrix=_alignment(align,score,matrix,algo,dia1,dia2);
		}
		else if(len == 4)
		{   
			Score<int, Simple > score(*(value),*(value+1),*(value+2),*(value+3));
			if(std::strcmp(algo,"MyersHirschberg")!=0 && std::strcmp(algo,"Hirschberg")!=0&&std::strcmp(algo,"SmithWaterman")!=0 )
				{ 
					scoreMatrix=_alignment(align,score,matrix,algo,dia1,dia2);
				}
				else if(std::strcmp(algo,"MyersHirschberg")==0)
				{
					scoreMatrix=globalAlignment(*align,score,MyersHirschberg());
				}
				else if(std::strcmp(algo,"Hirschberg")==0)
				{
					scoreMatrix=globalAlignment(*align,score,Hirschberg());		
				}
				else if(std::strcmp(algo,"SmithWaterman")==0)
				{
					scoreMatrix=localAlignment(*align,score, SmithWaterman());
				}
		}	
//std::cout << *align << std::endl;
		return align;
	} 

	template<typename Type>
	char* printAlignment(seqan::Align<seqan::String<Type> >&align,char **str1,char **str2)
	{
		typedef Align<seqan::String<Type> > const TAlign;
		typedef typename Row<TAlign>::Type TRow;
		typedef typename Position<typename Rows<TAlign>::Type>::Type TRowsPosition;
		typedef typename Position<TAlign>::Type TPosition;
		TRowsPosition row_count = length(rows(align));
		TPosition begin_ = beginPosition(cols(align));
		TPosition end_ = endPosition(cols(align));
		unsigned int baseCount=0;
		unsigned int leftSpace=6;
		std::string pAlignment;
		std::string seq1Alignment,seq2Alignment;
		while(begin_ < end_) 
		{
			unsigned int windowSize_ = 50;
			if ((begin_ + windowSize_)>end_) windowSize_=end_ - begin_;
			for(TRowsPosition i=0;i<2*row_count-1;++i)
			{
				for(unsigned int j = 0;j<leftSpace+2;++j) {
	//std::cout << " " << strlen(pAlignment) << std::endl; //strcat(pAlignment,"\n");
	pAlignment+=" "; }
				if ((i % 2)==0) 
				{
					TRow& row_ = row(align, i/2);
					typedef typename Iterator<typename Row<TAlign>::Type const, Standard>::Type TIter;
					TIter begin1_ = iter(row_, begin_);
					TIter end1_ = iter(row_, begin_ + windowSize_);
					for (; begin1_ != end1_; ++begin1_) 
					{
						if (isGap(begin1_))
						{
							 //strcat(pAlignment, "-");
							pAlignment+="-";
							if(i==0) {
	//						strcat(seq1Alignment, "-");
								seq1Alignment+="-";
							}
							else if(i==1)seq2Alignment+="-";   // strcat(seq2Alignment, "-");
						}
						else
						{
							char b=convert<char>(*begin1_);
							char* t=(char*)b;
							//strcat(pAlignment, &b );
							pAlignment+=b;
							if(i==0) seq1Alignment+=b;//strcat(seq1Alignment,&b);					
							else if(i==2)seq2Alignment +=b;//strcat(seq2Alignment,&b);					
						}
					}
				}
				else
				{
					for(unsigned int j = 0;j<windowSize_;++j)
					{
						if ((!isGap(row(align, (i-1)/2), begin_+j)) &&(!isGap(row(align, (i+1)/2), begin_+j)) && (row(align, (i-1)/2)[begin_+j]==row(align, (i+1)/2)[begin_+j]))
						{
							pAlignment+="|";//strcat(pAlignment,"|");
						}
						else
					 {
						pAlignment+=" ";//strcat(pAlignment," ");
						}
					} 
				}
				pAlignment+="\n";//strcat(pAlignment,"\n");
			}
			pAlignment+="\n";//strcat(pAlignment,"\n");
			begin_+=50;
		}
	pAlignment+="\n";//	strcpy(pAlignment,"\n");


	*str1 = new char[seq1Alignment.length()+1];
	*str2 = new char[seq2Alignment.length()+1];
	char* pp=new char[pAlignment.length()+1];
	strcpy(*str1,seq1Alignment.c_str());
	strcpy(*str2,seq2Alignment.c_str());
	strcpy(pp,pAlignment.c_str());
	return pp;
	}
/*
	char* printAlignment(seqan::Align<seqan::String<seqan::Dna5> >&align,char **str1,char **str2)
	{
		return _printAlignment(align,str1,str2);
	}

	char* printAlignment(seqan::Align<seqan::String<seqan::AminoAcid> >*align,char **str1,char **str2)
	{
		
		return _printAlignment(*align,str1,str2);
	}

	char* printAlignment(seqan::Align<seqan::String<char> >&align,char **str1,char **str2)
	{
		return _printAlignment(align,str1,str2);
	}
*/
}
