/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de 
  ===========================================================================
  Copyright (C) 2010
  
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.
  
  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
  
  ===========================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ===========================================================================
  Tests for the SeqAn model store, I/O functionality.
  ===========================================================================*/

#include <seqan/basic.h>  // For test functionality.
#include <seqan/store.h>  // Header under test.

using namespace seqan;

/*
SEQAN_DEFINE_TEST(test_store_io_writeme) 
{
	FragmentStore<> store;
    char buffer[1023];
    strcpy(buffer, SEQAN_PATH_TO_PROJECTS());
    strcat(buffer, "/projects/tests/store/ex1.sam");
    
	std::ifstream samFile(buffer);
	SEQAN_ASSERT_TRUE(samFile);
	read(samFile, store, SAM());
	
    strcpy(buffer, SEQAN_PATH_TO_PROJECTS());
    strcat(buffer, "/projects/tests/store/ex1.fa");
    
	loadContigs(store, buffer);


    strcpy(buffer, SEQAN_PATH_TO_PROJECTS());
    strcat(buffer, "/projects/tests/store/ex1.sam.copy");
    
	std::ofstream samFileOut(buffer);
	SEQAN_ASSERT_TRUE(samFileOut);
	write(samFileOut, store, SAM());
}
*/
SEQAN_DEFINE_TEST(test_store_io_sam) 
{
	FragmentStore<> store;
    char buffer[1023];
    strcpy(buffer, SEQAN_PATH_TO_PROJECTS());
    strcat(buffer, "/projects/tests/store/ex1.sam.copy");
	MultiSeqFile sam1;
	open(sam1.concat, buffer);
	split(sam1, Raw());
    
	{
		// read reference SAM from file
		std::ifstream samFile(buffer);
		SEQAN_ASSERT_TRUE(samFile);
		read(samFile, store, SAM());
	}
	
    strcpy(buffer, SEQAN_PATH_TO_PROJECTS());
    strcat(buffer, "/projects/tests/store/ex1.fa");
    
	loadContigs(store, buffer);

    strcpy(buffer, SEQAN_TEMP_FILENAME());
	{
		// write SAM to temp file
		std::ofstream samFileOut(buffer);
		SEQAN_ASSERT_TRUE(samFileOut);
		write(samFileOut, store, SAM());
	}
	
	MultiSeqFile sam2;
	open(sam2.concat, buffer);
	split(sam2, Raw());
	
	SEQAN_ASSERT_TRUE(!empty(sam1));
	SEQAN_ASSERT_TRUE(!empty(sam2));
	for (unsigned i = 0; i < length(sam1); ++i)
	{
		if (sam1[i] != sam2[i])
		{
			std::cout << "    \t" << sam1[i] << std::endl;
			std::cout << " != \t" << sam2[i] << std::endl;
			SEQAN_ASSERT_FAIL("Files differ in line %d.", i);
		}
	}
}
