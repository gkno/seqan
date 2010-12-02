//FRAGMENT(header)
/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de 
  ===========================================================================
  Copyright (C) 2010
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
  ===========================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ===========================================================================
  Minimapper -- Minimal read mapping in SeqAn.
  ===========================================================================
  This file contains code for a minimal read mapper with heavy restrictions.
  The restrictions are explained in the tutorial chapter, together with
  suggestions on how to extend this code.
  ===========================================================================*/
//FRAGMENT(includes)
#include <cstdio>
#include <iostream>

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <seqan/misc/misc_cmdparser.h>

using namespace seqan;

//FRAGMENT(typedefs)
// Some typedefs.
typedef FragmentStore<>::TReadSeqStore              TReadSeqStore;
typedef Value<TReadSeqStore>::Type                  TReadSeq;
typedef FragmentStore<>::TContigStore               TContigStore;
typedef Value<TContigStore>::Type                   TContigStoreElement;
typedef TContigStoreElement::TContigSeq             TContigSeq;
typedef Index<TReadSeqStore, Index_QGram<Shape<Dna, UngappedShape<11> >, OpenAddressing> > TIndex;
typedef Pattern<TIndex, Swift<SwiftSemiGlobal> >    TPattern;
typedef Finder<TContigSeq, Swift<SwiftSemiGlobal> > TFinder;
typedef FragmentStore<>::TAlignedReadStore          TAlignedReadStore;
typedef Value<TAlignedReadStore>::Type              TAlignedRead;

// const double EPSILON = 0.08;



struct Options
{
    double  epsilon;
    bool    rev_comp;

    Options()
    {
        epsilon = 0.08;
        rev_comp = false;
    }
    
};





//FRAGMENT(main-input)
int main(int argc, const char *argv[]) {
    // 0) Handle command line arguments.
    if (argc < 3) {
        std::cerr << "Invalid number of arguments." << std::endl
                  << "USAGE: minimapper GENOME.fasta READS.fasta OUT.sam" << std::endl;
        return 1;
    }

    CommandLineOption eps = CommandLineOption("e",
                                              "epsilon",
                                              1,
                                              "maximum error rate",
                                              OptionType::Double,
                                              0.08);

    CommandLineOption revcomp = CommandLineOption("rc",
                                              "reverse-complement",
                                              1,
                                              "align against both strands",
                                              OptionType::Bool,
                                              false);

    CommandLineParser parser;
    addOption(parser, eps);
    addOption(parser, revcomp);
    parse(parser, argc, argv);

    Options MyOptions;
    getOptionValueShort(parser, "e", MyOptions.epsilon);

    // 1) Load contigs and reads.
    FragmentStore<> fragStore;
    if (!loadContigs(fragStore, argv[1])) return 1;
    if (!loadReads(fragStore, argv[2])) return 1;

//     if (MyOptions.rev_comp)
//     {
//         int orglength = length(fragStore.contigStore);
//         resize(fragStore.contigStore, orglength *2);
//         for (int i=0; i < orglength; ++i)
//             assign(fragStore.contigStore[i+orglength],
//                       DnaStringReverseComplement(fragStore.contigStore[i]));
//     }
            

    
//FRAGMENT(pattern-finder)
    // 2) Build an index over all reads and a SWIFT pattern over this index.
    TIndex index(fragStore.readSeqStore);
    TPattern pattern(index);

    // String<TVerifyPattern> patterns
    // resize()
    // for ... setHost()
    
//FRAGMENT(swift)
    // 3) Enumerate all epsilon matches.
    for (unsigned i = 0; i < length(fragStore.contigStore); ++i) {
        TFinder finder(fragStore.contigStore[i].seq);
        // Verify Finder defined here
        while (find(finder, pattern, MyOptions.epsilon)) {
//FRAGMENT(verification)
            // Verify match.
            // clear(verifyFinder)
            // pattern = patterns[i]
            Finder<TContigSeq> verifyFinder(fragStore.contigStore[i].seq);
            setPosition(verifyFinder, beginPosition(finder));
//             Pattern<TReadSeq, HammingSimple> verifyPattern(fragStore.readSeqStore[position(pattern).i1]);
            Pattern<TReadSeq, Myers<> > verifyPattern(fragStore.readSeqStore[position(pattern).i1]);
            setScoreLimit(verifyPattern,
                          -floor(
                              MyOptions.epsilon
                              * length(fragStore.readSeqStore[position(pattern).i1])));
                              
            while (find(verifyFinder, verifyPattern)
                    && position(verifyFinder) < endPosition(infix(finder))
                  )
            {
                bool b = findBegin(verifyFinder, verifyPattern, getScore(verifyPattern));
                SEQAN_ASSERT_TRUE(b);
                TAlignedRead match(length(fragStore.alignedReadStore),
                                   position(pattern).i1,
                                   i,
                                   beginPosition(verifyFinder),
                                   endPosition(verifyFinder));
                std::cout << "BEGIN: " << beginPosition(verifyFinder)
                          << " END: " << endPosition(verifyFinder) << std::endl; 
                
                appendValue(fragStore.alignedReadStore, match);
//                 std::cout << fragStore.readSeqStore[position(pattern).i1] << std::endl;
//                 std::cout << infix(verifyFinder) << std::endl;
            }
        }
    }

    convertMatchesToGlobalAlignment(fragStore, Score<int>(0, -999, -1001, -1000), True());
//FRAGMENT(main-output)
    // 4) Write out SAM file.
    std::ofstream samFile(argv[3], std::ios_base::out);
    write(samFile, fragStore, SAM());

    return 0;
}
