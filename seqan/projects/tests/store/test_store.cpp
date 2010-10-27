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
  Tests for the SeqAn model store.
  ===========================================================================*/

#define SEQAN_ENABLE_CHECKPOINTS 0

#include <seqan/basic.h>
#include "test_store_io.h"

SEQAN_BEGIN_TESTSUITE(test_store) {
    SEQAN_CALL_TEST(test_store_io_sam);
}
SEQAN_END_TESTSUITE

