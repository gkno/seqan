// FRAGMENT(includes)
#include <iostream>
#include <seqan/sequence.h>

using namespace seqan;

// FRAGMENT(printFunction1)
template <typename TText>
void printAlign(TText const & genomeFragment, TText const & read)
{
        std::cout <<  "Alignment " << std::endl;
        std::cout << "  genome : ";
        for (unsigned j = 0; j < length(genomeFragment); ++j){
            std::cout << genomeFragment[j];
        }  
        std::cout << std::endl;
        std::cout << "  read   : ";
        for (unsigned j = 0; j < length(read); ++j){
            std::cout << read[j];
        }     
        std::cout << std::endl;
}
// FRAGMENT(printFunction2)
template <typename TInfix, typename TText>
void printAlign(TInfix const & genomeFragment, TText const & read)
{
        std::cout <<  "Alignment " << std::endl;
        std::cout << "  genome : ";
        for (unsigned j = 0; j < length(genomeFragment); ++j){
            std::cout << genomeFragment[j];
        }  
        std::cout << std::endl;
        std::cout << "  read   : ";
        for (unsigned j = 0; j < length(read); ++j){
            std::cout << read[j];
        }     
        std::cout << std::endl;
}


int main()
{
    Dna5String genome = "ATGGTTTCAACGTAATGCTGAACATGTCGCGTGTACTGACTATGCATGCATGACTG";

    String<Dna5String> readList;  
    resize(readList, 3);
    readList[0] = "GGTTTCGACG";
    readList[1] = "AAGATGTCGC";
    readList[2] = "TATGCATGAT";

    String<unsigned> beginPositions;
    resize(beginPositions, 3);
    beginPositions[0] = 2;
    beginPositions[1] = 20;
    beginPositions[2] = 42;

// FRAGMENT(initInfix)
    Infix<Dna5String>::Type inf;    

    for (unsigned i = 0; i < length(readList); ++i){
        
        unsigned begin = beginPositions[i];
        Dna5String read = readList[i];
// FRAGMENT(useInfix)
        unsigned end = begin + length(read);
        inf = infix(genome, begin, end);
        printAlign(inf, read);
    }
    return 0;
}
