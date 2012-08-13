//Adjust the code to handle Strings of type Dna, Dna5 and AminoAcid. Print the text on screen and observe how it changes depending on the type of the string.

#include <iostream>
#include <seqan/file.h>
#include <seqan/sequence.h>

int computeScore(seqan::String<seqan::Dna> subText, seqan::String<seqan::Dna> pattern)
{
    int localScore = 0;
    for (unsigned i = 0; i < seqan::length(pattern); ++i)
        if (subText[i] == pattern[i])
            ++localScore;
    
    return localScore;
}

int main()
{
    seqan::String<seqan::Dna> text = "This is an awesome tutorial to get to now the basic principles of SeqAn!";
    seqan::String<seqan::Dna> pattern = "tutorial";

    std::cout << text << std::endl;

    seqan::String<int> score;
    seqan::resize(score, seqan::length(text), 0);

    for (unsigned i = 0; i < seqan::length(text) - seqan::length(pattern) + 1; ++i)
        score[i] = computeScore(seqan::infix(text, i, i + seqan::length(pattern)), pattern);

    for (unsigned i = 0; i < seqan::length(score); ++i)
        std::cout << score[i] << " ";
    std::cout << std::endl;

    return 0;
}
