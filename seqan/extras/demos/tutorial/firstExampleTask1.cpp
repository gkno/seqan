//Adjust the code to handle Strings of type Dna, Dna5 and AminoAcid. Print the text on screen and observe how it changes depending on the type of the string.

#include <iostream>
#include <seqan/sequence.h>

int computeLocalScore(seqan::String<char> subText, seqan::String<char> pattern)
{
    int localScore = 0;
    for (unsigned i = 0; i < seqan::length(pattern); ++i)
        if (subText[i] == pattern[i])
            ++localScore;
    
    return localScore;
}

seqan::String<int> computeScore(seqan::String<char> text, seqan::String<char> pattern)
{
    seqan::String<int> score;
    seqan::resize(score, seqan::length(text), 0);

    for (unsigned i = 0; i < seqan::length(text) - seqan::length(pattern) + 1; ++i)
        score[i] = computeLocalScore(infix(text, i, i + seqan::length(pattern)), pattern);
    
    return score;
}

void print(seqan::String<int> score)
{
    for (unsigned i = 0; i < seqan::length(score); ++i)
        std::cout << score[i] << " ";
    std::cout << std::endl;
}

int main()
{
    seqan::String<char> text = "This is an awesome tutorial to get to now the basic principles of SeqAn!";
    seqan::String<char> pattern = "tutorial";
    seqan::String<int> score = computeScore(text, pattern);

    print(score);

    return 0;
}
