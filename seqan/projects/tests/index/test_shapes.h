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
 Author: David Weese <david.weese@fu-berlin.de>
 ==========================================================================*/

#ifndef TESTS_INDEX_TEST_SHAPES_H
#define TESTS_INDEX_TEST_SHAPES_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

template <typename TShape1, typename TShape2>
void testShape(TShape1 shape1, TShape2 shape2, bool dump)
{
	if (dump) std::cout << std::endl;
	              // 012345678901234 len=15
	DnaString dna = "CGGTACGTAAGTTAG";
	DnaString dna1, dna2;

	TShape1 shape1b(shape1);
	TShape2 shape2b(shape2);

	SEQAN_ASSERT_EQ(length(shape1), length(shape2));
	SEQAN_ASSERT_EQ(weight(shape1), weight(shape2));
	
	Iterator<DnaString>::Type it = begin(dna);	
	unsigned H1b, H2b;
	for (int i = length(dna); i >= 0; --i)
	{
		{
			unsigned H1 = hash(shape1, it, i);
			unsigned H2 = hash(shape2, it, i);
			if (i >= (int)length(shape1) && i != (int)length(dna))
			{
				H1b = hashNext(shape1b, it);
				H2b = hashNext(shape2b, it);
			} else {
				H1b = hash(shape1b, it, i);
				H2b = hash(shape2b, it, i);
			}

			unhash(dna1, H1, weight(shape1));
			unhash(dna2, H2, weight(shape2));
			if (dump) std::cout << std::dec << i << "\t" << std::hex << H1 << " " << H2 << ' ' << H1b << ' ' << H2b <<"\t" << dna1 << " " << dna2 << "\t" << (H1 == H2 && H1 == H1b && H1b == H2b);

			SEQAN_ASSERT_EQ(H1, H2);
			SEQAN_ASSERT_EQ(H1, H1b);
			SEQAN_ASSERT_EQ(H1b, H2b);

			if (i >= (int)length(shape1)) 
			{
				unsigned H1c = hash(shape1, it);
				unsigned H2c = hash(shape2, it);
				if (dump) std::cout << " " << (H1c == H2c && H1 == H1c);
			}
		}
		if (dump) std::cout << "\t";
		{
			unsigned H1 = hashUpper(shape1, it, i);
			unsigned H2 = hashUpper(shape2, it, i);
			
			unhash(dna1, H1, weight(shape1));
			unhash(dna2, H2, weight(shape2));
			if (dump) std::cout << std::dec << i << "\t" << std::hex << H1 << " " << H2 << "\t" << dna1 << " " << dna2 << "\t" << (H1 == H2);

			SEQAN_ASSERT_EQ(H1, H2);
		}
		if (dump) std::cout << std::endl;
		++it;
	}
}

SEQAN_DEFINE_TEST(testShapes)
{
	Shape<Dna, SimpleShape > shapeA(6);
	testShape(shapeA, Shape<Dna, UngappedShape<6> >(), false);
	
	                   // 012345678  len=9
	CharString pattern = "11100110100";
	Shape<Dna, GenericShape> shapeB(pattern);
	testShape(shapeB, Shape<Dna, GappedShape<HardwiredShape<1,1,3,1,2> > >(), false);

	pattern = "11110011";
	Shape<Dna, OneGappedShape> shapeC(pattern);
	testShape(shapeC, Shape<Dna, GappedShape<HardwiredShape<1,1,1,3,1> > >(), false);
}

//////////////////////////////////////////////////////////////////////////////


} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
