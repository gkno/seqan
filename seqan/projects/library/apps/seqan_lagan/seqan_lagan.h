/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de 
 ===========================================================================
  Copyright (C) 2010
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ===========================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
 ===========================================================================
  Global stuff for SeqAn::LAGAN.
 ===========================================================================*/

#ifndef SEQAN_LAGAN_H_
#define SEQAN_LAGAN_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

// ===========================================================================
// Tags, Enums, Classes, Specializations
// ===========================================================================

struct _Global;
typedef Tag<_Global> Global;

template <typename TSpec>
struct Options;

enum GlobalCommand
{
    COMMAND_LAGAN,
    COMMAND_CLASSIC_DP,
    COMMAND_EVALUATE
};

template <>
struct Options<Global>
{
    GlobalCommand command;
};

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

#endif  // SEQAN_LAGAN_H_
