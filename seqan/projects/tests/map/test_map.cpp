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

#include <seqan/basic.h>

#include "test_map_map.h"
#include "test_map_sumlist.h"


SEQAN_BEGIN_TESTSUITE(test_map) {
    // Call Tests.
    SEQAN_CALL_TEST(test_map_map);
    SEQAN_CALL_TEST(test_map_sumlist);
}
SEQAN_END_TESTSUITE