// FRAGMENT(solution)
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
    typedef String<char> TSequence;
    typedef Align<TSequence,ArrayGaps> TAlign;
    typedef Row<TAlign>::Type TRow;
    typedef Iterator<TRow>::Type TRowIterator;

    TSequence seq1 = "ACGTCACCTC";
    TSequence seq2 = "ACGGGCCTATC";

    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align,0),seq1);
    assignSource(row(align,1),seq2);

    TRow & row1 = row(align,0);
    TRow & row2 = row(align,1);

    insertGaps(row1,2,2);
    insertGap(row1,7);
    insertGaps(row2,9,2);


    TRowIterator itRow1 = begin(row1);
    TRowIterator itEndRow1 = end(row1);
    TRowIterator itRow2 = begin(row2);
    int gapCount = 0;
    for(;itRow1 != itEndRow1; ++itRow1, ++itRow2)
    {
        if(isGap(itRow1))
        {
            gapCount += countGaps(itRow1);
            itRow1 += countGaps(itRow1);
        }
        if(isGap(itRow2))
        {
            gapCount += countGaps(itRow2);
            itRow1 += countGaps(itRow2);
        }
    }
    ::std::cout << "Number of gaps: " << gapCount << ::std::endl;
}
