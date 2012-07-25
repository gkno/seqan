#include <iostream>
#include <seqan/sequence.h>

using namespace seqan;

// Function to print simple alignment between two sequences with the same length
// .. for two sequences of the same type
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
// .. for one infix and one other sequence type
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

    // Create Infix of type Dna5String
    Infix<Dna5String>::Type inf;       
    unsigned endPosition = begin + length(read);
    // Get the corresponding infix sequence of genome
    inf = infix(genome, begin, endPosition);
    printAlign(inf, read);
    return 0;
}
