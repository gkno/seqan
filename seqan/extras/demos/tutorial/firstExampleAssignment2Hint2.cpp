//Adjust the code to handle Strings of type Dna, Dna5 and AminoAcid. Print the text on screen and observe how it changes depending on the type of the string.

#include <iostream>
#include <seqan/sequence.h>

int computeScore(seqan::String<char> subText, seqan::String<char> pattern)
{
    int localScore = 0;
    for (unsigned i = 0; i < seqan::length(pattern); ++i)
        if (subText[i] == pattern[i])
            ++localScore;
    
    return localScore;
}

template <typename TText>
void print(TText const & text)
{
    std::cout << text << std::endl;
}

int main()
{
    seqan::String<char> text = "This is an awesome tutorial to get to now the basic principles of SeqAn!";
    seqan::String<char> pattern = "tutorial";

    print(text);

    seqan::String<int> score;
    seqan::resize(score, seqan::length(text), 0);

    for (unsigned i = 0; i < seqan::length(text) - seqan::length(pattern) + 1; ++i)
        score[i] = computeScore(seqan::infix(text, i, i + seqan::length(pattern)), pattern);

    print(score);

    return 0;
}
