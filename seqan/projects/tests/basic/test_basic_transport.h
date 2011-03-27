// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================

#ifndef TESTS_BASIC_TEST_BASIC_TRANSPORT_H_
#define TESTS_BASIC_TEST_BASIC_TRANSPORT_H_

struct MoveObj
{
	mutable int data_dat;

	MoveObj(int dat = 0): data_dat(dat) {}
	MoveObj(MoveObj const & other_): data_dat(other_.data_dat) 
	{ 
		other_.data_dat = 0; 
	}
	MoveObj const & operator = (MoveObj const & other_) const
	{ 
		data_dat = other_.data_dat; 
		other_.data_dat = 0; 
		return *this;
	}
	~MoveObj() {}
};

SEQAN_DEFINE_TEST(test_basic_transport)
{
	MoveObj o1(10);
	MoveObj o2;
	move(o2, o1);
	SEQAN_ASSERT_EQ(o1.data_dat, 0);
	SEQAN_ASSERT_EQ(o2.data_dat, 10);

	MoveObj const o3;
	move(o3, o2);
	SEQAN_ASSERT_EQ(o2.data_dat, 0);
	SEQAN_ASSERT_EQ(o3.data_dat, 10);

	MoveObj const o4;
	move(o4, o3);
	SEQAN_ASSERT_EQ(o3.data_dat, 0);
	SEQAN_ASSERT_EQ(o4.data_dat, 10);

	move(o1, o4);
	SEQAN_ASSERT_EQ(o4.data_dat, 0);
	SEQAN_ASSERT_EQ(o1.data_dat, 10);
}

#endif  // #ifndef TESTS_BASIC_TEST_BASIC_TRANSPORT_H_
