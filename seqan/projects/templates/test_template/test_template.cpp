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

#include <seqan/basic.h>  // @@Includes testing infrastructure.@@
#include <seqan/file.h>   // @@Required to print strings in tests.@@

// @@Create one header with tests for each of your headers under test.@@
#include "test_template_strings.h"
#include "test_template_others.h"


SEQAN_BEGIN_TESTSUITE(test_template) {
    // Call tests.
    SEQAN_CALL_TEST(test_template_strings_example1);
    SEQAN_CALL_TEST(test_template_others_example1);

    // Verify checkpoints.
    // @@
    // Remove this line and accordingly add a call for each of your headers
    // under test.  When the SeqAn testing mode is enabled, each call to
    // SEQAN_CHECKPOINTS will be registered with the testing system.
    // Verification of checkpoints means that we will check for each registered
    // checkpoint to be hit.
    // @@
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence/lexical.h");
}
SEQAN_END_TESTSUITE
