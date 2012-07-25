// FRAGMENT(includes)
#include <iostream>
#include <seqan/sequence.h>

using namespace seqan;

// FRAGMENT(printFunction1)
// Function to print simple alignment between two sequences with the same length
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
    // We have given a genome sequence
    Dna5String genome = "ATGGTTTCAACGTAATGCTGAACATGTCGCGT";
    // A read sequence
    Dna5String read = "TGGTNTCA";
    // And the begin position of a given alignment between the read and the genome
    unsigned beginPosition = 1;
// FRAGMENT(initInfix)
    Infix<Dna5String>::Type inf;       
// FRAGMENT(useInfix)
    unsigned end = begin + length(read);
    inf = infix(genome, begin, end);
    printAlign(inf, read);
    return 0;
}
