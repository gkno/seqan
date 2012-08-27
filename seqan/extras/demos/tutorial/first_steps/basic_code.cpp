// FRAGMENT(all)
#include <iostream>
#include <seqan/file.h>
#include <seqan/sequence.h>

int main()
{
    // Initialization
    seqan::String<char> text = "This is an awesome tutorial to get to know SeqAn!";
    seqan::String<char> pattern = "tutorial";

    seqan::String<int> score;
    resize(score, seqan::length(text));

    // Computation of the similarities
    // Iteration over the text (outer loop)
    for (unsigned i = 0; i < seqan::length(text) - seqan::length(pattern) + 1; ++i)
    {
        int localScore = 0;
        // Iteration over the pattern for character comparison
        for (unsigned j = 0; j < seqan::length(pattern); ++j)
        {
            if (text[i + j] == pattern[j])
                ++localScore;
        }
        score[i] = localScore;
    }

    // Printing the result
    for (unsigned i = 0; i < seqan::length(score); ++i)
        std::cout << score[i] << " ";
    std::cout << std::endl;

    return 0;
}
