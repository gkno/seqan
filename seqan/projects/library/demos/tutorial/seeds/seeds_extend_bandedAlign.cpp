//FRAGMENT(includes)
#include <seqan/seeds.h>

using namespace seqan;

// FRAGMENT(writeSeed)
template<typename TSeed, typename TSeq>
void
writeSeed(TSeed & seed, TSeq const & seq0, TSeq const & seq1) {
    std::cout << "Seed from position " << leftPosition(seed, 0);
    std::cout << " to " << rightPosition(seed, 0) << ": ";
    std::cout << infix(seq0, leftPosition(seed, 0), rightPosition(seed, 0)+1) << std::endl;
    std::cout << "Seed from position " << leftPosition(seed, 1);
    std::cout << " to " << rightPosition(seed, 1) << ": ";
    std::cout << infix(seq1, leftPosition(seed, 1), rightPosition(seed, 1)+1) << std::endl;
}

//FRAGMENT(main)
int main() {
    typedef Seed<> TSeed;
    
    DnaString seq0 = "ATCATCAGTTATACTTTACCCAGGC";
    DnaString seq1 = "ATTCAGCATACTTTCCATGAAGC";

    TSeed seed(10, 7, 7);
    writeSeed(seed, seq0, seq1);    

// FRAGMENT(extension)
    typedef int TScore;
    Score<TScore, Simple> scoreMatrix(1, -1, -1);
    TScore scoreDropOff = 1;
    extendSeed(seed, scoreDropOff, scoreMatrix, seq0, seq1, 2, GappedXDrop());

    std::cout << std::endl << "After extension:" << std::endl;
    writeSeed(seed, seq0, seq1);  

// FRAGMENT(banded-alignment)
    Align<DnaString> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seq0);
    assignSource(row(align, 1), seq1);
    
    std::cout << std::endl << "Banded Alignment:" << std::endl;
    std::cout << "Score: " << bandedAlignment(align, seed, 2, scoreMatrix) << std::endl;
    std::cout << align;

// FRAGMENT(global-alignment)
    Align<DnaString> globalAlign;
    resize(rows(globalAlign), 2);
    assignSource(row(globalAlign, 0), seq0);
    assignSource(row(globalAlign, 1), seq1);

    std::cout << std::endl << "Global Alignment:" << std::endl;
    std::cout << "Score: " << globalAlignment(globalAlign, scoreMatrix) << std::endl;
    std::cout << globalAlign;

    return 0;
}
