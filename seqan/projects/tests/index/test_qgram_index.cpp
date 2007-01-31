/*
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



template <typename TShapeSpec, typename TConstruct>
void testQGramIndexSchnell3()
{
	clock_t start, finish;
	double duration;

	String<Dna> text;
	fstream strm_t;
	//strm_t.open("z:/testdb.fasta", ios_base::in);
	strm_t.open(TEST_PATH "fasta.txt", ios_base::in);
	read(strm_t, text, Fasta());
	strm_t.close();
	cout << length(text)<<"<-textlen\n";

	String<char> shape_aussehen = "xxx_x__xxx_x__xxx_x";
	int q = length(shape_aussehen);
	Shape<Dna,TShapeSpec> shape;
	stringToShape(shape_aussehen, shape);

	typedef Position<String<Dna> >::Type TPosition;
	String<TPosition> index;
	resize(index, length(text) - q + 2);	// # q-Grams in der Query +1, damit der letzte Zeiger ein Nachfolgezeiger hat um zu stoppen
	
	String<TPosition> pos;
    //int pos_size = intPow((int)ValueSize<Dna>::VALUE, q - numGaps(shape));
	int pos_size = (int) pow((float)ValueSize<Dna>::VALUE, q - numGaps(shape));
	pos_size += 1;		// damit der letzte Zeiger ein Nachfolgerzeiger hat, um zu wissen wann er stoppen muss
	resize(pos, pos_size);
	for(int i = 0; i < pos_size; ++i)
		pos[i] = 0;

	start = clock();
	TConstruct const construct;
	createQGramIndex(text, shape, construct, pos, index);
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "\nQGramIndex bauen dauert: " << duration << " Sekunden.\n\n";
	
	
}

template <typename TShapeSpec, typename TConstruct>
void testQGramIndexSchnell2()
{
	clock_t start, finish;
	double duration;

	String<Dna> text;
	fstream strm_t;
	//strm_t.open("z:/testdb1.fasta", ios_base::in);
	strm_t.open(TEST_PATH "fasta.txt", ios_base::in);
	read(strm_t, text, Fasta());
	strm_t.close();
	cout << length(text)<<"<-textlen\n";

	String<char> shape_aussehen = "x___x_xxx";
//	String<char> shape_aussehen = "x_x_x_x_x_x_x_x";
	int q = length(shape_aussehen);
	Shape<Dna,TShapeSpec> shape;
	stringToShape(shape_aussehen, shape);

	typedef Position<String<Dna> >::Type TPosition;
	String<TPosition> index;
	resize(index, length(text) - q + 2);	// # q-Grams in der Query +1, damit der letzte Zeiger ein Nachfolgezeiger hat um zu stoppen
	
	String<TPosition> pos;
    int pos_size = (int)pow((float)ValueSize<Dna>::VALUE, q - numGaps(shape));
	pos_size += 1;		// damit der letzte Zeiger ein Nachfolgerzeiger hat, um zu wissen wann er stoppen muss
	resize(pos, pos_size);
	for(int i = 0; i < pos_size; ++i)
		pos[i] = 0;

	start = clock();
	TConstruct const construct;
	createQGramIndex(text, shape, construct, pos, index);
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "\nQGramIndex bauen dauert: " << duration << " Sekunden.\n\n";
	
	
}


template <typename TShapeSpec, typename TConstruct>
void testQGramIndex()
{
	String<Dna> text = "CTGAACCCTAAACCCT";
	String<char> shape_aussehen = "x_x";
	int q = length(shape_aussehen);
	Shape<Dna,TShapeSpec> shape;
	stringToShape(shape_aussehen, shape);

	typedef Position<String<Dna> >::Type TPosition;
	String<TPosition> index;
	resize(index, length(text) - q + 2);	// # q-Grams in der Query +1, damit der letzte Zeiger ein Nachfolgezeiger hat um zu stoppen
	
	String<TPosition> pos;
    int pos_size = (int)pow((float)ValueSize<Dna>::VALUE, q - numGaps(shape));
	pos_size += 1;		// damit der letzte Zeiger ein Nachfolgerzeiger hat, um zu wissen wann er stoppen muss
	resize(pos, pos_size);
	for(int i = 0; i < pos_size; ++i)
		pos[i] = 0;

	TConstruct const construct;
	createQGramIndex(text, shape, construct, pos, index);
	
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


void Main_TestQGram()
{
	SEQAN_TREPORT("TEST QGRAM BEGIN")

		testQGramIndex<GappedShape3,QGram2>();
//		testQGramIndexSchnell3<GappedShape3,QGram2>();
	//	testQGramIndexSchnell2<GappedShape1,QGram1>();
	SEQAN_TREPORT("TEST QGRAM END")
}
*/
