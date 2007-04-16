#include <iostream>
#include <fstream>
#include <typeinfo>
#include <time.h>
#include <cstdio>
#include <vector>
#include <time.h>

#define SEQAN_DEBUG
#define SEQAN_TEST
#define TEST_PATH "projects/tests/index/"

#include <seqan/index.h>
#include <seqan/file.h>

using namespace std;
using namespace seqan;

/*
void testGappedShapes()
{
	String<char> shape_string = "__x_xxx_xxx_x";
	Shape<Dna,GappedShape> shape1;
	stringToShape(shape1, shape_string);
	Shape<Dna,GappedShape> shape2 = Shape<Dna,GappedShape>(shape1);

	for(unsigned i = 0; i < shape1.shape_len; ++i)
        SEQAN_ASSERT(shape1[i] == shape2[i]);
	SEQAN_ASSERT(shape1.num_gaps == shape2.num_gaps);
	SEQAN_ASSERT(shape1.span == shape2.span);
	SEQAN_ASSERT(shape1.shape_len == shape2.shape_len);
	SEQAN_ASSERT(length(shape1) == length(shape2));
	SEQAN_ASSERT(shapeCountBlanks(shape1) == shapeCountBlanks(shape2));

	Shape<Dna,GappedShape> shape3 = Shape<Dna,GappedShape>(5, 13);
	shape3[0]=2;
	shape3[1]=2;
	shape3[2]=1;
	shape3[3]=1;
	shape3[4]=2;
	shape3[5]=1;
	shape3[6]=1;
	shape3[7]=2;
	for(int i = 0; i < 8; ++i)
        SEQAN_ASSERT(shape1[i] == shape3[i]);


}
*/

void testUngappedShapes()
{
	String<char> shape_string = "xxxx";
	Shape<Dna,SimpleShape> shape1;
	stringToShape(shape1, shape_string);
	Shape<Dna,SimpleShape> shape2 = Shape<Dna,SimpleShape>(shape1);

	SEQAN_ASSERT(shape1.span == shape2.span);
	SEQAN_ASSERT(shape1.leftFactor == shape2.leftFactor);
	SEQAN_ASSERT(length(shape1) == length(shape2));
	SEQAN_ASSERT(shapeCountBlanks(shape1) == shapeCountBlanks(shape2));
	

	Shape<Dna,SimpleShape> shape3 = Shape<Dna,SimpleShape>(4);
	SEQAN_ASSERT(shape3.leftFactor == 64);
	SEQAN_ASSERT(shape1.leftFactor == shape3.leftFactor);


}

/*
void testQGramIndexSchnell()
{
	clock_t start, finish;
	double duration;

	String<Dna> text;
	fstream strm_t;
	strm_t.open(TEST_PATH "fasta.txt", ios_base::in);
	read(strm_t, text, Fasta());
	String<Dna> next;
	resize(next,length(text));
	for(int i = 0; i < 1; ++i)							// datei zu ende?
	{
		arrayCopyForward(begin(text),end(text),begin(next));
		append(text, next);
	}
	strm_t.close();

	String<char> shape_string = "x___x__x____xxx__xx";
	int q = length(shape_string);
	Shape<Dna,GappedShape> shape;
	stringToShape(shape, shape_string);

	typedef Position<String<Dna> >::Type TPosition;
	String<TPosition> index;
	resize(index, length(text) - q + 2);	
	
	String<TPosition> pos;	
	int pos_size = _intPow(ValueSize<Dna>::VALUE, q - shapeCountBlanks(shape));
	pos_size += 1;	
	resize(pos, pos_size);

	start = clock();
	createQGramIndex(index, pos, text, shape);
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	//cout << "\nQGramIndex bauen dauert: " << duration << " Sekunden.\n\n";
	
	
}
*/
/*
void testGappedQGramIndex()
{
	String<Dna> text = "CTGAACCCTAAACCCT";
	String<char> shape_string = "x_x";
	int q = length(shape_string);
	Shape<Dna,GappedShape> shape;
	stringToShape(shape, shape_string);

	typedef Position<String<Dna> >::Type TPosition;
	String<TPosition> index;
	resize(index, length(text) - q + 2);
	
	String<TPosition> pos;
    int pos_size = _intPow(ValueSize<Dna>::VALUE, q - shapeCountBlanks(shape));
	pos_size += 1;	
	resize(pos, pos_size);

	createQGramIndex(index, pos, text, shape);
	
	SEQAN_ASSERT(pos[0] == 0);
	SEQAN_ASSERT(pos[1] == 1);
	SEQAN_ASSERT(pos[2] == 5);
	SEQAN_ASSERT(pos[3] == 5);
	SEQAN_ASSERT(pos[4] == 5);
	SEQAN_ASSERT(pos[5] == 6);
	SEQAN_ASSERT(pos[6] == 8);
	SEQAN_ASSERT(pos[7] == 9);
	SEQAN_ASSERT(pos[8] == 11);
	SEQAN_ASSERT(pos[9] == 12);
	SEQAN_ASSERT(pos[10] == 12);
	SEQAN_ASSERT(pos[11] == 12);
	SEQAN_ASSERT(pos[12] == 12);
	SEQAN_ASSERT(pos[13] == 14);
	SEQAN_ASSERT(pos[14] == 14);
	SEQAN_ASSERT(pos[15] == 14);
	SEQAN_ASSERT(pos[16] == 14);

	SEQAN_ASSERT(index[0] == 9);
	SEQAN_ASSERT(index[1] == 11);
	SEQAN_ASSERT(index[2] == 10);
	SEQAN_ASSERT(index[3] == 4);
	SEQAN_ASSERT(index[4] == 3);
	SEQAN_ASSERT(index[5] == 7);
	SEQAN_ASSERT(index[6] == 12);
	SEQAN_ASSERT(index[7] == 5);
	SEQAN_ASSERT(index[8] == 0);
	SEQAN_ASSERT(index[9] == 13);
	SEQAN_ASSERT(index[10] == 6);
	SEQAN_ASSERT(index[11] == 2);
	SEQAN_ASSERT(index[12] == 8);
	SEQAN_ASSERT(index[13] == 1);
	
}
*/
void testUngappedQGramIndex()
{
	String<Dna> text = "CTGAACCCTAAACCCT";
	String<char> shape_string = "xx";
	int q = length(shape_string);
	Shape<Dna,SimpleShape> shape;
	stringToShape(shape, shape_string);

	typedef Position<String<Dna> >::Type TPosition;
	String<TPosition> index;
	resize(index, length(text) - q + 1);
	
	String<TPosition> pos;
    int pos_size = _intPow(ValueSize<Dna>::VALUE, q - shapeCountBlanks(shape));
	pos_size += 1;	
	resize(pos, pos_size);

	createQGramIndex(index, pos, text, shape);
	
	
	SEQAN_ASSERT(pos[0] == 0);
	SEQAN_ASSERT(pos[1] == 3);
	SEQAN_ASSERT(pos[2] == 5);
	SEQAN_ASSERT(pos[3] == 5);
	SEQAN_ASSERT(pos[4] == 5);
	SEQAN_ASSERT(pos[5] == 5);
	SEQAN_ASSERT(pos[6] == 9);
	SEQAN_ASSERT(pos[7] == 9);
	SEQAN_ASSERT(pos[8] == 12);
	SEQAN_ASSERT(pos[9] == 13);
	SEQAN_ASSERT(pos[10] == 13);
	SEQAN_ASSERT(pos[11] == 13);
	SEQAN_ASSERT(pos[12] == 13);
	SEQAN_ASSERT(pos[13] == 14);
	SEQAN_ASSERT(pos[14] == 14);
	SEQAN_ASSERT(pos[15] == 15);

	SEQAN_ASSERT(index[0] == 3);
	SEQAN_ASSERT(index[1] == 9);
	SEQAN_ASSERT(index[2] == 10);
	SEQAN_ASSERT(index[3] == 4);
	SEQAN_ASSERT(index[4] == 11);
	SEQAN_ASSERT(index[5] == 5);
	SEQAN_ASSERT(index[6] == 6);
	SEQAN_ASSERT(index[7] == 12);
	SEQAN_ASSERT(index[8] == 13);
	SEQAN_ASSERT(index[9] == 0);
	SEQAN_ASSERT(index[10] == 7);
	SEQAN_ASSERT(index[11] == 14);
	SEQAN_ASSERT(index[12] == 2);
	SEQAN_ASSERT(index[13] == 8);
	SEQAN_ASSERT(index[14] == 1);
}


void Main_TestQGram()
{
	SEQAN_TREPORT("TEST QGRAM BEGIN")

//	testGappedQGramIndex();
	testUngappedQGramIndex();
//	testQGramIndexSchnell();
//	testGappedShapes();
	testUngappedShapes();

	debug::verifyCheckpoints("projects/library/seqan/index/index_qgram.h");
	debug::verifyCheckpoints("projects/library/seqan/index/shape_base.h");
	debug::verifyCheckpoints("projects/library/seqan/index/shape_gapped.h");


	SEQAN_TREPORT("TEST QGRAM END")
}
