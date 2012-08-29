// ==========================================================================
//                                mini_bowtie
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#include <iostream>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/index.h>
#include <seqan/index_fm.h>
#include <seqan/store.h>

using namespace seqan;

template <typename TIter, typename TStringSet>
void search(TIter & it, TStringSet const & pattern)
{
    typedef typename Iterator<TStringSet const, Standard>::Type TPatternIter;

    for (TPatternIter patternIt = begin(pattern, Standard()); patternIt != end(pattern, Standard()); ++patternIt)
    {
        // exact search on pattern half
        std::cerr << "Pattern: " << value(patternIt) << std::endl;
        unsigned startApproxSearch = length(value(patternIt)) / 2;
        if (goDown(it, infix(value(patternIt), startApproxSearch + 1, length(value(patternIt)))))
        {
            for (unsigned i = startApproxSearch; ; --i)
            {
                Dna character = getValue(patternIt)[i];
                for (Dna5 c = MinValue<Dna>::VALUE; c < +ValueSize<Dna>::VALUE; ++c)
                {
                    if (c != character)
                    {
                        TIter localIt = it;
                        if (goDown(localIt, c))
                        {
                            if (goDown(localIt, infix(value(patternIt), 0, i)))
                            {
                               //HIT 
                            }
                        }
                    }
                }
                if (!goDown(it, character))
                    break;
                else if (i == 0)
                {
                    // HIT 
                    break;
                }
            }
        }
        goRoot(it);
    }
}

int main(int argc, char *argv[]) 
{
    typedef String<Dna5> TString;
    typedef StringSet<String<Dna5> > TStringSet;
    typedef Index<StringSet<TString>, FMIndex<> > TIndex;
    typedef Iterator<TIndex, TopDown<> >::Type TIter;

    // 0) Handle command line arguments.
    if (argc < 3) {
        std::cerr << "Invalid number of arguments." << std::endl
                  << "USAGE: minimapper GENOME.fasta READS.fasta OUT.sam" << std::endl;
        return 1;
    }
    // 1) Load contigs and reads.
    FragmentStore<> fragStore;
    if (!loadContigs(fragStore, argv[1])) return 1;
    if (!loadReads(fragStore, argv[2])) return 1;

    StringSet<TString> text;
    for (unsigned i = 0; i < length(fragStore.contigStore); ++i)
        appendValue(text, fragStore.contigStore[i].seq);
        
    TIndex fmIndex(text);
    TIter it(fmIndex);
    search(it, fragStore.readSeqStore);
    clear(fmIndex);
    clear(it);

    reverse(text);
    reverse(fragStore.readSeqStore);

    fmIndex = TIndex(text);
    it = TIter(fmIndex);
    search(it, fragStore.readSeqStore);
    clear(fmIndex);
    clear(it);

    reverse(text);
    reverse(fragStore.readSeqStore);
    std::ofstream samFile(argv[3], std::ios_base::out);
    write(samFile, fragStore, Sam());

    return 0;
}
