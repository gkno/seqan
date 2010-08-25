Fiona: v.expected.value.pvalue.several_error_corrections.iterations.different_reads_length
date :23.07.2010
TODO: Third method - with power
----------------
Necessaty files:
----------------
	Makefile 
	Fiona.cpp
----------------
Using libraries:
----------------
	SeqAn downloanding and instruction for instalation from http://www.seqan.de/

---------------------
Compile by executing:
--------------------- 
	make suffix

---------
Clean by:
---------
	make clean

------	
Using:
------ 
	./Fiona [options] [input_reads_file] [corrected_reads_file] 

--------------------
Optional parameters: 
--------------------
      	-f $1: choose the correction method. We propose two methods - by default the method is with pvalue estimation at the base of the Poisson distribution. To switch to another method it must give the strictness value ($1), which is float value necessaire to estimate the confidential intervall. 
	-p $1: float value fixing the p-value. By default 0.0001. If -f is choice, don't use -p. 
	-i $1: positive integer giving the number of times to run the correction method, by default $1 is equal to 3.
	-l $1 $2: positives integers for specify the levels between, which the search will be made. $1 is the top level after which the searching will be done. $2 is the value for the down level until the search will continue. By defalt $1 is estimate as 2 levels after the log4 of the total number of input reads and $2 is 10 levels after the top level. When there are errors in the two ends of the read, small than the log ration, it is better to beggin several levels before. More the level is high, more the time of execution is shorter. 
	-g $1: giving a float value for the length of the genome. By default $1 is estimate at the first level of the tree, before begging the correction.
	-m $1: positive integer for the number of mismatch that are accepted in the length of each read. By default $1 is set to 1. After that at each iteration this value is decreasing. 

--------------------	
Required Parameters:
--------------------
[input_read_file] - Path and name of the input file. 
	The input file with reads can be in FASTA format or just one reads by line separate or without a line. 
    The reads in this version are of equal length and contain only the letters {A,C,G,T}. 
	
[corrected_reads_file] - Path and name of output file containing all reads from [input_read_file]
	 with a part of them detect by FIONA as erroneous and corrected in the better manner.

	 

