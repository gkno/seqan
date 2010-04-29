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
  @@By its name, this would test the header template/others.h.@@
  ===========================================================================
*/

// @@Lines and sections with two at-chars are comments that should either
// be removed or replaced by tests defined by you.@@

// @@Replace symbol accordingly to your file name!@@
#ifndef TEST_TEMPLATE_TEST_TEMPLATE_OTHERS_H_
#define TEST_TEMPLATE_TEST_TEMPLATE_OTHERS_H_

#include <seqan/basic.h>     // @@For the testing infrastructure.@@
#include <seqan/sequence.h>  // @@Replace with your header under test.@@


// @@This is an empty example test.  Replace it with your own test.@@
SEQAN_DEFINE_TEST(test_template_others_example1) {
    using namespace seqan;

    // Just an empty test, show assertion with message.
    SEQAN_ASSERT_FAIL("Failure message!");
    // You will not see the following message.
    SEQAN_ASSERT_FAIL("Only one failure is reported per test!");
}

// @@Replace symbol accordingly to your file name!@@
#endif  // TEST_TEMPLATE_TEST_TEMPLATE_OTHERS_H_
