// Write a print function that print a text on the screen. If the type of the string is int include a space between the values. Do NOT change existing code.
// Hint: Start by implementing the generic version.

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

void print(seqan::String<int> const & text)
{
    for (unsigned i = 0; i < seqan::length(text); ++i)
        std::cout << text[i] << " ";
    std::cout << std::endl;
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

