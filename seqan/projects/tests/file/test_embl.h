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
 Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
 Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
 ==========================================================================*/

#ifndef TESTS_FILE_TEST_FILE_EMBL_H_
#define TESTS_FILE_TEST_FILE_EMBL_H_

using namespace std;
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////
SEQAN_DEFINE_TEST(test_file_embl_file)
{
    char buffer[1023];
    strcpy(buffer, SEQAN_PATH_TO_PROJECTS());
    strcat(buffer, "/projects/tests/file/takifugu_scl_embl.txt");

	std::fstream strm; 
	strm.open(buffer, ios_base::in | ios_base::binary);

	String<char> line;
	String<char> feature_line;

	readLineType(strm, feature_line, "FT", Embl());
	//cout << feature_line << "\n";

	int count = 0;
	int next_pos = readFeature(feature_line, 0, line, "exon", Embl());
	while(next_pos != 0)
	{
		++count;
	//  cout << line << "\n";
		next_pos = readFeature(feature_line, next_pos, line, "exon", Embl());
	}

	SEQAN_ASSERT_EQ(count, 3);
}


SEQAN_DEFINE_TEST(test_file_embl_meta)
{
    char buffer[1023];
    strcpy(buffer, SEQAN_PATH_TO_PROJECTS());
    strcat(buffer, "/projects/tests/file/takifugu_scl_embl.txt");
    
	std::fstream strm; 
	strm.open(buffer, ios_base::in | ios_base::binary);

	String<char> line;

	String<char> feature_line;
	String<char> meta;
	readMeta(strm,meta,Embl());
	
	readLineType(meta, line, "KW", Embl());
	SEQAN_ASSERT_EQ(line, "SCL gene.");

	readLineType(meta, line, "RX", Embl());
	SEQAN_ASSERT_EQ(infix(line,0,28), "DOI; 10.1073/pnas.101532998.");
	SEQAN_ASSERT_TRUE(length(line) == 46u || length(line) == 47u);

	clear(line);
	readLineType(meta, feature_line, "FT", Embl());
	//cout << feature_line << "\n";

	int count = 0;
	int next_pos = readFeature(feature_line, 0, line, "CDS", Embl());
	while(next_pos != 0)
	{
		++count;
	//  cout << line << "\n";
		next_pos = readFeature(feature_line, next_pos, line, "CDS", Embl());
	}

	SEQAN_ASSERT_EQ(count, 1);
}

#endif  // TESTS_FILE_TEST_FILE_EMBL_H_
