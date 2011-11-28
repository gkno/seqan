// ==========================================================================
//                  ALF - Alignment free sequence comparison
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
//
// ==========================================================================
// Author: Jonathan Goeke <goeke@molgen.mpg.de>
// ==========================================================================

#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/misc/edit_environment.h>
#include <seqan/misc/misc_cmdparser.h>
#include <seqan/alignment_free.h>

#include "alf.h"

using namespace seqan;
using namespace std;

void printHelp()
{
    cerr << endl << "Usage: alf -i <InputMultiFasta.fa> -o <outputFile.txt>" << endl;
    cerr << endl << "Advanced: alf --method [N2/D2/D2z/D2Star] --kmerSize [integer] --reverseComplement [min/max/mean/bothStrands] --mismatches [0/1] --mismatchWeight [double] --inputFile <InputMultiFasta.fa> --outputFile <outputFile.txt>" << endl;
    cerr << endl << "Options:" << endl;
    cerr << "  -m,   --method        \t\t" << "Method to use, either N2, D2,D2Star (background model estimation on both sequences), D2z" << endl;
    cerr << "  -i,   --inputFile     \t\t" << "Name of the multifasta inputfile" << endl;
    cerr << "  -o,   --outputFile    \t\t" << "Name of file to which the tab delimited matrix with pairwise scores will be written" << endl;
    cerr << "  -k,   --kmerSize      \t\t" << "Size of the kmers" << endl;
    cerr << "  -mo,  --bgModelOrder  \t\t" << "Order of background markov model" << endl;
    cerr << "  -rc,  --reverseComplement\t" << "N2 only, default: only input strand (alternative: mean,min,max or 'bothStrands' to score both strands simultaneously)" << endl;
    cerr << "  -n,  --norm\t" << "normalize scores, only N2 [default TRUE] and D2Star [default FALSE]" << endl;
    cerr << "  -mm,  --mismatches\t" << "0 or 1:calculate N2 using the kmer-neighbourhood with one mismatch" << endl;
    cerr << "  -mmw,  --mismatchWeight\t" << "weight of counts for words with mismatches" << endl;
    cerr << "  -kwf, --kmerWeightsFile\t\t" << "Print kmerWeights to this file, N2 and D2star only" << endl;
    cerr << "  -v, --verbose\t\t" << " Print details on progress to the screen" << endl;
    cerr << "  -rt, --runningTime\t\t" << "Print running time into this file" << endl;

    cerr << "  -h,   --help          \t\t" << "This screen" << endl;

}

int main(int argc, char * argv[])
{
    std::cout << "\n******************************************" << endl;
    std::cout << "***                 ALF                ***" << endl;
    std::cout << "*** Alignment free sequence comparison ***" << endl;
    std::cout << "******************************************" << endl << endl;
    // std::cout.precision(25);
    // Declare all parameters
    string method = "N2";
    int kmerSize = 4;
    int bgModelOrder = 1;
    string outFile = "";
    string inFile = "";
    string revCom = "";
    bool norm = true;

    bool verbose = false;
    string kmerWeightsFile = "";
    unsigned mismatches = 0;
    double mismatchWeight = 0.1;

    double startTime= sysTime();
    string runTimeFile = "";
    
    // Read in comman line arguments
    for (int arg = 1; arg < argc; arg++)
    {
        if (argv[arg][0] == '-')
        {
            // parse options
            if (strcmp(argv[arg], "-m") == 0 || strcmp(argv[arg], "--method") == 0)
            {
                ++arg;
                method = argv[arg];
                std::cout << "\n--method: "; // Display option
                for (unsigned int i = 0; i < length(method); ++i)
                {
                    std::cout << method[i];
                }
                continue;
            }
            if (strcmp(argv[arg], "-i") == 0 || strcmp(argv[arg], "--inputFile") == 0)
            {
                ++arg;
                inFile = argv[arg];
                std::cout << "\n--inputFile: "; // Display option
                for (unsigned int i = 0; i < length(inFile); ++i)
                {
                    std::cout << inFile[i];
                }
                continue;
            }
            if (strcmp(argv[arg], "-o") == 0 || strcmp(argv[arg], "--outputFile") == 0)
            {
                ++arg;
                outFile = argv[arg];
                std::cout << "\n--outputFile: ";
                for (unsigned int i = 0; i < length(outFile); ++i)
                {
                    std::cout << outFile[i];
                }
                continue;
            }
            if (strcmp(argv[arg], "-k") == 0 || strcmp(argv[arg], "--kmerSize") == 0)
            {
                ++arg;
                kmerSize = atoi(argv[arg]);
                std::cout << "\n--kmerSize: " << (int) (kmerSize - 0);
                continue;
            }
            if (strcmp(argv[arg], "-mo") == 0 || strcmp(argv[arg], "--bgModelOrder") == 0)
            {
                ++arg;
                bgModelOrder = atoi(argv[arg]);
                std::cout << "\n--bgModelOrder: " << bgModelOrder;
                continue;
            }
            if (strcmp(argv[arg], "-rc") == 0 || strcmp(argv[arg], "--reverseComplement") == 0)
            {
                ++arg;
                revCom = argv[arg];
                std::cout << "\n--reverseComplement: ";
                for (unsigned int i = 0; i < length(revCom); ++i)
                {
                    std::cout << revCom[i];
                }
                continue;
            }
            if (strcmp(argv[arg], "-n") == 0 || strcmp(argv[arg], "--norm") == 0)
            {
                ++arg;
                if (strcmp(argv[arg], "true") == 0)
                {
                    norm = true;
                }
                else
                {
                    norm = false;
                }
                std::cout << "\n--norm: " << norm;
                continue;
            }
            if (strcmp(argv[arg], "-mm") == 0 || strcmp(argv[arg], "--mismatches") == 0)
            {
                ++arg;
                mismatches = atof(argv[arg]);
                std::cout << "\n--mismatches: " << mismatches;
                continue;
            }
            if (strcmp(argv[arg], "-mmw") == 0 || strcmp(argv[arg], "--mismatchWeight") == 0)
            {
                ++arg;
                mismatchWeight = atof(argv[arg]);
                std::cout << "\n--mismatchWeight: " << mismatchWeight;
                continue;
            }
            if (strcmp(argv[arg], "-kwf") == 0 || strcmp(argv[arg], "--kmerWeightsFile") == 0)
            {
                ++arg;
                kmerWeightsFile = argv[arg];
                std::cout << "\n--kmerWeightsFile: ";
                for (unsigned int i = 0; i < length(kmerWeightsFile); ++i)
                {
                    std::cout << kmerWeightsFile[i];
                }
                continue;
            }
            if (strcmp(argv[arg], "-v") == 0 || strcmp(argv[arg], "--verbose") == 0)
            {
                ++arg;
                if (strcmp(argv[arg], "true") == 0)
                {
                    verbose = true;
                }
                else
                {
                    verbose = false;
                }
                std::cout << "\n--verbose: " << verbose;
                continue;
            }
            if (strcmp(argv[arg], "-rt") == 0 || strcmp(argv[arg], "--runningTime") == 0)
            {
		++arg;
		runTimeFile = argv[arg];
		std::cout << "\n--runningTime: ";
                for (unsigned int i = 0; i < length(runTimeFile); ++i)
                {
                    std::cout << runTimeFile[i];
                }
		continue;
            }
            if (strcmp(argv[arg], "-h") == 0 || strcmp(argv[arg], "--help") == 0)
            {
                printHelp();
                return 0;
            }

        }

    }
    std::cout << "\n";
    // Definition of type DNA string sets
    typedef String<Dna5> TText;
    typedef Size<TText>::Type TSize;
    typedef StringSet<TText> TStringSet;

    // Definition of mxn two-dimensional matrix
    typedef Matrix<double, 2> TMatrix;

    TText seqIID1 =
    	"TAGGTTTTCCGAAAAGGTAGCAACTTTACGTGATCAAACCTCTGACGGGGTTTTCCCCGTCGAAATTGGGTG"
        "TTTCTTGTCTTGTTCTCACTTGGGGCATCTCCGTCAAGCCAAGAAAGTGCTCCCTGGATTCTGTTGCTAACG"
        "AGTCTCCTCTGCATTCCTGCTTGACTGATTGGGCGGACGGGGTGTCCACCTGACGCTGAGTATCGCCGTCAC"
        "GGTGCCACATGTCTTATCTATTCAGGGATCAGAATTTATTCAGGAAATCAGGAGATGCTACACTTGGGTTAT"
        "CGAAGCTCCTTCCAAGGCGTAGCAAGGGCGACTGAGCGCGTAAGCTCTAGATCTCCTCGTGTTGCAACTACA"
        "CGCGCGGGTCACTCGAAACACATAGTATGAACTTAACGACTGCTCGTACTGAACAATGCTGAGGCAGAAGAT"
        "CGCAGACCAGGCATCCCACTGCTTGAAAAAACTATNNNNCTACCCGCCTTTTTATTATCTCATCAGATCAAG";
    TText seqIID2 =
        "ACCGACGATTAGCTTTGTCCGAGTTACAACGGTTCAATAATACAAAGGATGGCATAAACCCATTTGTGTGAA"
        "AGTGCCCATCACATTATGATTCTGTCTACTATGGTTAATTCCCAATATACTCTCGAAAAGAGGGTATGCTCC"
        "CACGGCCATTTACGTCACTAAAAGATAAGATTGCTCAAANNNNNNNNNACTGCCAACTTGCTGGTAGCTTCA"
        "GGGGTTGTCCACAGCGGGGGGTCGTATGCCTTTGTGGTATACCTTACTAGCCGCGCCATGGTGCCTAAGAAT"
        "GAAGTAAAACAATTGATGTGAGACTCGACAGCCAGGCTTCGCGCTAAGGACGCAAAGAAATTCCCTACATCA"
        "GACGGCCGCGNNNAACGATGCTATCGGTTAGGACATTGTGCCCTAGTATGTACATGCCTAATACAATTGGAT"
        "CAAACGTTATTCCCACACACGGGTAGAAGAACNNNNATTACCCGTAGGCACTCCCCGATTCAAGTAGCCGCG";

    TText seqM1_1 =
        "TGAATGCTATGAGCTTTTTTATGAGTGGACTAATGGGGGTAAAAATTTTAAGGTCAGTTTTTCAAAAAGGTA"
        "TCCCATCAAGATGTAAGCTGAGTATTGACGCTGAAAATCAAGGGCAGAGTGATGGTGTTTTGAACCCTTGTG"
        "ATCATGGGCCATAAAATAGTTTCTTTCTTTCTAGGTAAATTAACCTGCTACACCCTTACCCAACTTTGATCA"
        "CTTGACCTTTGCTGGATTATCCATTGATTTGTTGGCCCCTGACCCTTCTGGGTATTATCTTGATGAAAGAAG"
        "TTTAACGGTGTCAGCTACTTTCAAATATTTAACAGTAACACTGAATTGCAAGGTTTTTTGTAAATGAGAAAG"
        "GTGGAAGACTGTGATAGCAGACAGACTTCTACACTCGTTAGGCCAAGGTGTGACTTCTAGAAAAGTGTTAAG"
        "GTCAGATAGAAATGCCAACAGAAAGGACATGATTTTAGTTGGTGATTGCCTATACCAATTTTTCCTTC";
    TText seqM1_2 =
        "TCCCATGTCTTCCTGAGTTGTTGAGTGCTTGTTGCTCATACCTATTGAAAGATTGGGACTGTTCCCCCCCTT"
        "TTTTTTAAGGGATTAAACGAGATCTCGTTTACTTACGGCATTAGTCACTGATCTAGCCTCAATTGGTGCTTT"
        "CAAGCATTAGGGTTGATGAGTGTTTCATGCTAACTACTACAACCTTTCTTTTGTAGCGCCATTTTGACATGA"
        "CCACAATTTCCACCCTGAAAAGGCCAATGTACCCAACTGTGAATCCTTCCATGCAGGCCATCTATCTTACCT"
        "ATATTGGGGAGTGTAGTGAAATGGTTCTTTCCCTGCTGAGTAGGCCAAGCAATTTACAACAACCCAGTGCCA"
        "ATAGGTAAGGAGGTGAGAAATAGAGGGGGGGAACACAGTTGGAATTTCATATCATGTCTGAATGCTGTTTTA"
        "TGCTTCAGCAGGAATGTGGGTTATCTCCAGTACCCGGCCTGATGATTCTAATTCTTTAATTTTTTTGG";

    TMatrix myMatrix;       // myMatrix stores pairwise kmerScores

    TStringSet mySequenceSet;       // mySequenceSet stores all sequences from the multi-fasta file

    if (inFile.compare(""))     // read in file
    {
        MultiSeqFile multiSeqFile;
        open(multiSeqFile.concat, inFile.c_str(), OPEN_RDONLY);
        AutoSeqFormat format;
        guessFormat(multiSeqFile.concat, format);
        split(multiSeqFile, format);
        unsigned seqCount = length(multiSeqFile);
        StringSet<String<Dna5Q> > seqs;
        StringSet<CharString> seqIDs;
        reserve(mySequenceSet, seqCount, Exact());
        reserve(seqIDs, seqCount, Exact());
        String<Dna5Q> seq;
        TText sequenceTMP;
        CharString qual;
        CharString id;
        for (unsigned i = 0; i < seqCount; ++i)
        {
            assignSeq(sequenceTMP, multiSeqFile[i], format);        // read sequence
            assignSeqId(id, multiSeqFile[i], format);       // read sequence id


            appendValue(mySequenceSet, sequenceTMP, Generous());
            appendValue(seqIDs, id, Generous());
        }
    }
    else        // if no file is given some standard strings are used
    {
        printHelp();
        appendValue(mySequenceSet,seqIID1);
        appendValue(mySequenceSet,seqIID2);
        // appendValue(mySequenceSet,seqIID3);
        cout << "\n Example Run:\n";
        //appendValue(mySequenceSet, seqM1_1);
        //appendValue(mySequenceSet, seqM1_2);

    }
    //std::cout.precision(3);
    if (!method.compare("D2"))
    {
        //-----D2-----
        AFScore<D2> myScoreD2(kmerSize, verbose);
	//myScoreD2.verbose=verbose;
        alignmentFreeComparison(myMatrix, mySequenceSet, myScoreD2);
    }
    else if (!method.compare("D2z"))
    {
        //-----D2z-----
        AFScore<D2z> myScoreD2z(kmerSize, bgModelOrder, verbose);
        alignmentFreeComparison(myMatrix, mySequenceSet, myScoreD2z);
    }
    else if (!method.compare("D2Star"))
    {
        //-----D2star-----
       // String<char> file;
        //file = kmerWeightsFile;
        AFScore<D2Star> myScoreD2Star(kmerSize, bgModelOrder, verbose);
        alignmentFreeComparison(myMatrix, mySequenceSet, myScoreD2Star);
    }
    else if (!method.compare("N2"))
    {
        //-----N2-----
        String<char> file;
        file = kmerWeightsFile;
        AFScore<N2> myScoreN2(kmerSize, bgModelOrder, revCom, mismatches, mismatchWeight, file, verbose);
        alignmentFreeComparison(myMatrix, mySequenceSet, myScoreN2);
    }

    if (outFile.compare(""))
    {
        ofstream myfile;
        myfile.open(outFile.c_str());
        myfile << myMatrix;
        myfile.close();
    }
    else
    {
        cout << "\n" << myMatrix;
    }
    if (runTimeFile.compare(""))
    {
        ofstream myfile;
        myfile.open(runTimeFile.c_str());
        double endTime= sysTime();
        myfile <<(endTime-startTime);
        myfile.close();
    }

    return 0;
}