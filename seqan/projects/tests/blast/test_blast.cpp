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
  Author: @@Your Name@@ <@@Your Email@@>
  ===========================================================================
  @@Description of what is tested here@@
  ===========================================================================
*/

// @@Lines and sections with two at-chars are comments that should either
// be removed or replaced by tests defined by you.@@


#include <seqan/file.h>
#include <seqan/blast.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <time.h>
#include <cstdio>
#include <vector>
#include <time.h>
#include <string>
#include <ctime>


// @@Create one header with tests for each of your headers under test.@@
#include "test_blast_parsing.h"

SEQAN_BEGIN_TESTSUITE(test_blast) {
    // Call tests.
	SEQAN_CALL_TEST(test_blast_store_report_int);
	SEQAN_CALL_TEST(test_blast_store_report_basic_int);
	SEQAN_CALL_TEST(test_blast_parsing_int_blastn);
	SEQAN_CALL_TEST(test_blast_parsing_int_blastp);	
	SEQAN_CALL_TEST(test_blast_parsing_basic_int_blastn);
    SEQAN_CALL_TEST(test_blast_parsing_basic_int_blastp);


    // Verify checkpoints.
    // @@
    // Remove this line and accordingly add a call for each of your headers
    // under test.  When the SeqAn testing mode is enabled, each call to
    // SEQAN_CHECKPOINTS will be registered with the testing system.
    // Verification of checkpoints means that we will check for each registered
    // checkpoint to be hit.
    // @@
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/blast/blast_base.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/blast/blast_parsing.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/blast/blast_base.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/blast/blast_report.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/blast/blast_hit.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/blast/blast_hsp.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/blast/blast_stream_report.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/blast/blast_stream_hit.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/blast/blast_iterator.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/blast/blast_hit_iterator.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/blast/blast_hsp_iterator.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/blast/blast_stream_hit_iterator.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/blast/blast_stream_hsp_iterator.h");
	
}
SEQAN_END_TESTSUITE
