#include <iostream>
#include <seqan/sequence.h>

using namespace seqan;


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
    Dna5String genome = "ATGGTTTCAACGTAATGCTGAACATGTCGCGTGTACTGACTATGCATGCATGACTG";

    String<Dna5String> readList;  
    resize(readList, 3);
    readList[0] = "GGTTTCGACG";
    readList[1] = "AAGATGTCGC";
    readList[2] = "TATGCATGAT";

    // Below we create a list, containing the begin positions
    // for each read of its given alignment to the genome
    String<unsigned> beginPositions;
    resize(beginPositions, 3);
    beginPositions[0] = 2;
    beginPositions[1] = 20;
    beginPositions[2] = 42;

    Dna5String genomeFragment;
    // For each read out of the list we want to call the printAlign function 
    // to print the alignment between the read and the genome
    for (unsigned i = 0; i < length(readList); ++i){
        Dna5String read = readList[i];
        unsigned begin = beginPositions[i];
        clear(genomeFragment);
        // We have to create a copy of the corresponding fragment of the genome, where the read align to
        for (unsigned j = 0; j < length(read); ++j){
            appendValue(genomeFragment, genome[begin+j]);
        }
        printAlign(genomeFragment, read);
    }

    return 0;
}
