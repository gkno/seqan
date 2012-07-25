#include <iostream>
#include <seqan/sequence.h>

using namespace seqan;

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


int main()
{
    // We have given a genome sequence
    Dna5String genome = "ATGGTTTCAACGTAATGCTGAACATGTCGCGT";
    // A read sequence
    Dna5String read = "TGGTNTCA";
    // And the begin position of a given alignment between the read and the genome
    unsigned beginPosition = 1;

    Dna5String genomeFragment;       
    // We have to create a copy of the corresponding fragment of the genome, where the read aligns to
    for (unsigned i = 0; i < length(read); ++i){
        appendValue(genomeFragment, genome[beginPosition+i]);
    }
    // Call of our function to print the simple alignment
    printAlign(genomeFragment, read);
  
    return 0;
}
