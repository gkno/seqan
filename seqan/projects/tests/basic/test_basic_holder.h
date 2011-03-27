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
// Author: Manuel Holtgrewe <manuel.holtgrew@fu-berlin.de>
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Tests for the Holder classes.
// ==========================================================================

#ifndef TESTS_BASIC_TEST_BASIC_HOLDER_H_
#define TESTS_BASIC_TEST_BASIC_HOLDER_H_

#include <seqan/basic.h>

using namespace seqan;

//////////////////////////////////////////////////////////////////////////////
// Helper class for constructor/destructor call counting.
// This is needed for test of Holder.

struct CtorDtorCounter
{
	static int static_ctors;
	static int static_dtors;

	int data_value;

	CtorDtorCounter():
		data_value(0)
	{
		++static_ctors;
	}
	CtorDtorCounter(CtorDtorCounter const & other_):
		data_value(other_.data_value)
	{
		++static_ctors;
	}
	~CtorDtorCounter()
	{
		++static_dtors;
	}
	CtorDtorCounter & operator = (CtorDtorCounter const & other_)
	{
		data_value = other_.data_value;
		return *this;
	}
};

int CtorDtorCounter::static_ctors = 0;
int CtorDtorCounter::static_dtors = 0;

//____________________________________________________________________________
// test for holder class

SEQAN_DEFINE_TEST(test_basic_holder)
{
	{
//ctors
		Holder<CtorDtorCounter> ho1;
		SEQAN_ASSERT(empty(ho1));
		SEQAN_ASSERT(!dependent(ho1));

		create(ho1);
		SEQAN_ASSERT(!empty(ho1));
		SEQAN_ASSERT(!dependent(ho1));

		Holder<CtorDtorCounter> ho2(ho1);

		Holder<CtorDtorCounter> ho3(value(ho1));
		SEQAN_ASSERT(!dependent(ho1));


//create
		CtorDtorCounter rco1;
		create(ho3, rco1);

//setValue
		setValue(ho3, rco1);
		SEQAN_ASSERT(dependent(ho3));

		rco1.data_value = 10;
		create(ho3);
		SEQAN_ASSERT_EQ(value(ho3).data_value, 10);

		CtorDtorCounter rco2;
		rco2.data_value = 20;

//operator = (value) => assignValue
		ho2 = rco2;
		SEQAN_ASSERT_EQ(value(ho2).data_value, 20);
		SEQAN_ASSERT(!dependent(ho2));

		rco2.data_value = 30;
		SEQAN_ASSERT_EQ(value(ho2).data_value, 20);

//operator = (holder) => assign
		setValue(ho1, rco1);
		ho1 = ho2;
		SEQAN_ASSERT_EQ(value(ho1).data_value, 20);

//clear
		clear(ho3);
		SEQAN_ASSERT(empty(ho3));

		assign(ho2, ho3);
		SEQAN_ASSERT(empty(ho2));

//conversion operator
		rco1 = ho1;
		SEQAN_ASSERT_EQ(rco1.data_value, 20);

//moveValue
		moveValue(ho1, rco2);
		SEQAN_ASSERT_EQ(rco1.data_value, 30);

	}

	SEQAN_ASSERT_EQ(CtorDtorCounter::static_ctors, CtorDtorCounter::static_dtors);

//test default implementations of addRef and releaseRef

//test const object holders
/*
	typedef char Bla[100];
	Holder<Bla const> cho1 = "test";*/
}

#endif  // #ifndef TESTS_BASIC_TEST_BASIC_HOLDER_H_
