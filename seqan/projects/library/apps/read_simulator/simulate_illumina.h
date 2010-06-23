#ifndef SIMULATE_ILLUMINA_H_
#define SIMULATE_ILLUMINA_H_

#include <numeric>

#include <seqan/store.h>
#include <seqan/sequence_journal.h>

#include "read_simulator.h"

using namespace seqan;


// Error probabilities for Illumina reads of length 36, from
// Anne-Katrin's interpol36subsBiases.dat.
const static double ERROR_PROBABILITIES_36[36] = {
    0.005,
    0.003671429,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0037,
    0.004142857,
    0.004585714,
    0.005028571,
    0.005471429,
    0.005914286,
    0.006,
    0.0066,
    0.007242857,
    0.008057143,
    0.009257143,
    0.01071429,
    0.01492857,
    0.01271429,
    0.0141,
    0.01827143,
    0.02057143,
    0.02305714,
    0.0267,
    0.03037143,
    0.03534286,
    0.038
};


// Error probabilities for Illumina reads of length 50, from
// Anne-Katrin's interpol50subsBiases.dat.
const static double ERROR_PROBABILITIES_50[50] = {
    0.005,
    0.004051020,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.003510204,
    0.003826531,
    0.004142857,
    0.004459184,
    0.00477551,
    0.005091837,
    0.005408163,
    0.00572449,
    0.006,
    0.006,
    0.006346939,
    0.006979592,
    0.007306122,
    0.007867347,
    0.008816327,
    0.009510204,
    0.01071429,
    0.01387755,
    0.01397959,
    0.01258163,
    0.01384694,
    0.01618367,
    0.01966327,
    0.02057143,
    0.02191837,
    0.02476531,
    0.02714286,
    0.02961224,
    0.03340816,
    0.03610204,
    0.038
};


// Error probabilities for Illumina reads of length 100, from
// Anne-Katrin's interpol100subsBiases.dat.
const static double ERROR_PROBABILITIES_100[100] = {
    0.005,
    0.004530303,
    0.004060606,
    0.003590909,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.0035,
    0.003606061,
    0.003762626,
    0.003919192,
    0.004075758,
    0.004232323,
    0.004388889,
    0.004545455,
    0.00470202,
    0.004858586,
    0.005015152,
    0.005171717,
    0.005328283,
    0.005484848,
    0.005641414,
    0.00579798,
    0.005954545,
    0.006,
    0.006,
    0.006,
    0.006161616,
    0.006474747,
    0.006787879,
    0.007050505,
    0.00720707,
    0.007363636,
    0.007560606,
    0.008030303,
    0.0085,
    0.008969697,
    0.00929293,
    0.00960606,
    0.009919192,
    0.01116162,
    0.01272727,
    0.01429293,
    0.01457071,
    0.01378788,
    0.01300505,
    0.01272222,
    0.01334848,
    0.01397475,
    0.01477778,
    0.0165,
    0.01822222,
    0.01994444,
    0.02030303,
    0.02061616,
    0.02092929,
    0.02209091,
    0.0235,
    0.02490909,
    0.02613636,
    0.02723232,
    0.02832828,
    0.02972727,
    0.03160606,
    0.03348485,
    0.03518182,
    0.03612121,
    0.03706061,
    0.038
};


struct _IlluminaReads;
typedef Tag<_IlluminaReads> IlluminaReads;


struct IlluminaOptions {
    // true iff help is to be shown.
    bool showHelp;
    // true iff verbosity is enabled.
    bool verbose;
    // true if very verbose is enabled.
    bool veryVerbose;

    // Basic Read Simulation Parameters.

    // Seed to use for the random number generator.
    unsigned seed;
    // Length of the reads to simulate.
    unsigned readLength;
    // Number of reads (pairs) to simulate.
    unsigned numReads;
    // true iff a random reference sequence is to be used.
    bool useRandomSequence;
    // Length of random sequence to be simulated.
    unsigned randomSourceLength;
    // Maximal number of errors per read.
    unsigned maxErrorsPerRead;
    // true iff only the forward strand is to be simulated.
    bool onlyForward;
    // true iff only the reverse strand is to be simulated.
    bool onlyReverse;
    // The output file.  Defaults to REFERENCE-FILE.reads.fastq,
    // possibly with a ".1" or ".2" before the ".fastq" if mate pairs
    // are simulated.
    CharString outputFile;
    // Path to the SAM file to generate.  Defaults to fastq file name
    // with suffix ".sam"
    CharString samFile;
    // true iff qualities are to be simulated.
    bool simulateQualities;

    // Mate-Pair Related Options.

    // true iff generating mate pairs is enabled.
    bool generateMatePairs;
    // Mate-pair library length.
    unsigned libraryLength;
    // Mate-pair library length tolerance.
    unsigned libraryLengthError;

    // Sequencing Error Probability Options.

    // Path to the mismatch error distribution file.  If empty,
    // built-ins (available for n=36, 50, 100) will be used.
    CharString errorDistributionFile;
    // Probability of an insertion.
    double probabilityInsert;
    // Probability of a deletion.
    double probabilityDelete;
    // Standard error in %*100 around error probability for quality simulation.
    double qualityErrorFactor;

    // Source Sequence Repeat Parameters.

    // TODO(holtgrew): Mostly interesting for randomly generated source sequences.

    // Haplotype parameters.

    // Number of haplotypes to generated.  All are generated from the input genome.
    unsigned numHaplotypes;
    // SNP rate.
    double haplotypeSnpRate;
    // Indel rate.
    double haplotypeIndelRate;
    // Smallest number of indels.
    unsigned haplotypeIndelRangeMin;
    // Largest number of indels.
    unsigned haplotypeIndelRangeMax;

    IlluminaOptions()
            : showHelp(false),
              verbose(true),
              veryVerbose(true),
              seed(0),
              readLength(36),
              numReads(1000),
              useRandomSequence(false),
              randomSourceLength(1000*1000),
              maxErrorsPerRead(4),
              onlyForward(false),
              onlyReverse(false),
              outputFile(""),
              samFile(""),
              simulateQualities(false),
              generateMatePairs(false),
              libraryLength(1000),
              libraryLengthError(100),
              errorDistributionFile(""),
              probabilityInsert(0.01),
              probabilityDelete(0.01),
              qualityErrorFactor(0.1),
              numHaplotypes(1),
              haplotypeSnpRate(0.01),
              haplotypeIndelRate(0.01),
              haplotypeIndelRangeMin(4),
              haplotypeIndelRangeMax(6)
    {}
};


template <typename TStream>
TStream & operator<<(TStream & stream, IlluminaOptions const & options) {
    stream << "illumina-options {" << std::endl
           << "  seed:                   " << options.seed << std::endl
           << "  readLength:             " << options.readLength << std::endl
           << "  numReads:               " << options.numReads << std::endl
           << "  useRandomSequence:      " << (options.useRandomSequence ? "true" : "false") << std::endl
           << "  randomSourceLength:     " << options.randomSourceLength << std::endl
           << "  maxErrorsPerRead:       " << options.maxErrorsPerRead << std::endl
           << "  onlyForward:            " << (options.onlyForward ? "true" : "false") << std::endl
           << "  onlyReverse:            " << (options.onlyReverse ? "true" : "false") << std::endl
           << "  outputFile:             \"" << options.outputFile << "\"" <<std::endl
           << "  samFile:                \"" << options.samFile << "\"" <<std::endl
           << "  simulateQualities:      " << (options.simulateQualities ? "true" : "false") << std::endl
           << "  generateMatePairs:      " << (options.generateMatePairs ? "true" : "false") << std::endl
           << "  libraryLength:          " << options.libraryLength << std::endl
           << "  libraryLengthError:     " << options.libraryLengthError << std::endl
           << "  errorDistributionFile:  \"" << options.errorDistributionFile << "\"" <<std::endl
           << "  probabilityInsert:      " << options.probabilityInsert << std::endl
           << "  probabilityDelete:      " << options.probabilityDelete << std::endl
           << "  qualityErrorFactor:     " << options.qualityErrorFactor << std::endl
           << "  numHaplotypes:          " << options.numHaplotypes << std::endl
           << "  haplotypeSnpRate:       " << options.haplotypeSnpRate << std::endl
           << "  haplotypeIndelRate:     " << options.haplotypeIndelRate << std::endl
           << "  haplotypeIndelRangeMin: " << options.haplotypeIndelRangeMin << std::endl
           << "  haplotypeIndelRangeMax: " << options.haplotypeIndelRangeMax << std::endl
           << "}" << std::endl;
    return stream;
}


// Instructions for simulating a read in the Illumina model.
// TODO(holtgrew): For 454 reads, the edit string must specify what is inserted.
struct ReadSimulationInstruction {
    unsigned haplotype;
    unsigned contigId;
    bool isForward;
    size_t beginPos;
    size_t endPos;
    // Number of characters added/removed to the string by indels.
    unsigned delCount;
    unsigned insCount;
    String<ErrorType> editString;
    String<int> qualities;

    ReadSimulationInstruction() : delCount(0), insCount(0) {}
};


template <typename TStream>
TStream & operator<<(TStream & stream, ReadSimulationInstruction const & inst) {
    stream << "(haplotype=" << inst.haplotype << ", contigId=" << inst.contigId << ", isForward=" << inst.isForward << ", beginPos=" << inst.beginPos << ", endPos=" << inst.endPos << ", insCount=" << inst.insCount << ", delCount=" << inst.delCount << ", editString=";
    for (unsigned i = 0; i < length(inst.editString); ++i) {
        stream << "MEID"[inst.editString[i]];
    }
    stream << ", qualities=[";
    for (unsigned i = 0; i < length(inst.qualities); ++i) {
        if (i != 0)
            stream << ", ";
        stream << inst.qualities[i];
    }
    stream << "])";
    return stream;
}


void buildEditString(ReadSimulationInstruction & inst, unsigned readLength, String<double> const & errorProbabilities) {
    SEQAN_ASSERT_EQ(readLength * 4, length(errorProbabilities));
    clear(inst.editString);
    reserve(inst.editString, static_cast<size_t>(1.2 * readLength), Generous());
    inst.delCount = 0;
    inst.insCount = 0;

    for (unsigned i = 0; i < readLength; /*NOP*/) {
        double x = mtRandDouble();
        double pMatch    = errorProbabilities[i * 4 + ERROR_TYPE_MATCH];
        double pMismatch = errorProbabilities[i * 4 + ERROR_TYPE_MISMATCH];
        double pInsert   = errorProbabilities[i * 4 + ERROR_TYPE_INSERT];
        if (x < pMatch) {
//             std::cout << "match" << std::endl;
            // match
            i += 1;
            appendValue(inst.editString, ERROR_TYPE_MATCH);
        } else if (x < pMatch + pMismatch) {
//             std::cout << "mismatch" << std::endl;
            // mismatch
            i += 1;
            appendValue(inst.editString, ERROR_TYPE_MISMATCH);
        } else if (x < pMatch + pMismatch + pInsert) {
//             std::cout << "insert" << std::endl;
            // insert
            if (length(inst.editString) > 0 && back(inst.editString == ERROR_TYPE_DELETE)) {
                inst.delCount -= 1;
                eraseBack(inst.editString);
            } else {
                i += 1;
                inst.insCount += 1;
                appendValue(inst.editString, ERROR_TYPE_INSERT);
            }
        } else {
//             std::cout << "delete" << std::endl;
            // Decrement string size, do not add a delete if string is
            // too short, possibly remove insert from edit string.
            if (length(inst.editString) > 0) {
                if (back(inst.editString == ERROR_TYPE_INSERT)) {
                    i -= 1;
                    inst.insCount -= 1;
                    eraseBack(inst.editString);
                } else {
                    inst.delCount += 1;
                    appendValue(inst.editString, ERROR_TYPE_DELETE);
                }
            }
        }
    }
//     for (unsigned i = 0; i < length(inst.editString); ++i)
//         std::cout << " " << inst.editString[i];
//     std::cout << std::endl;
//     std::cout << "del " << inst.delCount << " ins " << inst.insCount << std::endl;
    SEQAN_ASSERT_EQ(readLength, length(inst.editString) - inst.delCount);
}


void buildQualityValues(ReadSimulationInstruction & inst, String<double> const & errorProbabilities, IlluminaOptions const & options) {
    SEQAN_ASSERT_GT(length(inst.editString), 0u);
    clear(inst.qualities);
    fill(inst.qualities, length(inst.editString), 0, Exact());
    String<double> tmp;
    fill(tmp, inst.endPos - inst.beginPos + inst.delCount, 0, Exact());

    // TODO(holtgrew): Quality computation is HIGHLY bogus.
    //
    // We derive the read qualities from the error probabilities and
    // the actual match/mismatch/indel event at a point.  For a match,
    // we compute the quality as the phred score of x where x is
    // picked around options.qualityErrorFactor around the error
    // probability for this location.  For mismatches, we do the same
    // and half the resulting score.  Matches and mismatches are
    // handled in the first round.
    //
    // In a second round, we compute the scores for indels as half the
    // arithmetic mean of the quality values of the neighbouring
    // bases.
    //
    // First round.
    for (unsigned i = 0, j = 0; i < length(inst.editString); i++) {
        SEQAN_ASSERT_LEQ(j, inst.endPos - inst.beginPos + inst.delCount);
        if (inst.editString[i] == ERROR_TYPE_MATCH || inst.editString[i] == ERROR_TYPE_MISMATCH) {
            double p = 1 - errorProbabilities[j * 4 + ERROR_TYPE_MATCH];
            double delta = mtRandDouble() * 2 * options.qualityErrorFactor * p;
            double x = p - options.qualityErrorFactor * p + delta;
            int score = -10 * std::log10(x);
            if (inst.editString[i] == ERROR_TYPE_MISMATCH)
                score = std::ceil(score / 2.0);
            // Store scores.
            inst.qualities[i] = score;
            tmp[j] = score;
            j += 1;
        }
    }

    // Second round.
    for (unsigned i = 0, j = 0; i < length(inst.editString); i++) {
        if (inst.editString[i] == ERROR_TYPE_INSERT || inst.editString[i] == ERROR_TYPE_DELETE) {
            if (j == 0) {
                // Indel at beginning.
                inst.qualities[i] = tmp[0];
            } else if (j == inst.endPos - inst.beginPos - 1) {
                // Indel at end.
                inst.qualities[i] = back(tmp);
            } else {
                // Indel in center.
                inst.qualities[i] = 0.25 * (tmp[j] + tmp[j + 1]);
            }
        } else {
            j += 1;
        }
    }
}


/* Indels change the size of a sequence.  We use Strings of DeltaEntry
   structs to efficiently convert from the original position
 */
struct DeltaEntry {
    size_t pos;
    long int delta;
};


bool operator<(DeltaEntry const & a, DeltaEntry const & b) {
    return a.pos < b.pos;
}


void buildHaplotype(StringSet<SequenceJournal<String<Dna5>, SortedArray> > & haplotype,
                    FragmentStore<MyFragmentStoreConfig> & fragmentStore,
                    IlluminaOptions const & options) {
    resize(haplotype, length(fragmentStore.contigStore), Exact());
    String<Dna5> buffer;
    reserve(buffer, options.haplotypeIndelRangeMax);

//     std::cout << back(haplotype)._journalEntries << std::endl;

    for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i) {
        clear(haplotype[i]);
        setHost(haplotype[i], fragmentStore.contigStore[i].seq);
        String<Dna5> const & contig = fragmentStore.contigStore[i].seq;
        SequenceJournal<String<Dna5>, SortedArray> & haplotypeContig = haplotype[i];

        // j is position in original sequence, k is position in haplotype
        for (size_t j = 0, k = 0; j < length(contig);) {
            if (k >= length(haplotypeContig))
                std::cout << haplotypeContig._journalEntries << std::endl;
            SEQAN_ASSERT_LT_MSG(k, length(haplotypeContig), "k = %lu, length(contig) == %lu", k, length(contig));

            double x = mtRandDouble();
            if (x < options.haplotypeSnpRate) {
                // SNP
                Dna5 c = Dna5(4 * mtRandDouble());
                if (c == contig[j])
                    c = Dna5(ordValue(c) + 1);
                assignValue(haplotypeContig, k, c);
                j += 1;
                k += 1;
            } else if (x < options.haplotypeSnpRate + options.haplotypeIndelRate) {
                // Indel of random length.
                unsigned rangeLen = options.haplotypeIndelRangeMax - options.haplotypeIndelRangeMin;
                unsigned indelLen = options.haplotypeIndelRangeMin + static_cast<unsigned>(mtRandDouble() * rangeLen);
                if (mtRandDouble() < 0.5) {
                    // Insertion.
                    clear(buffer);
                    for (unsigned ii = 0; ii < indelLen; ++ii)
                        appendValue(buffer, Dna5(5 * mtRandDouble()));
                    insert(haplotypeContig, k, buffer);
                    k += indelLen;
                } else {
                    // Deletion.
                    erase(haplotypeContig, k, k + indelLen);
                    j += indelLen;
                }
            } else {
                // Match.
                j += 1;
                k += 1;
            }
        }
        
//         std::cout << "contig = " << contig << std::endl;
//         std::cout << "back(haplotype) = " << back(haplotype) << std::endl;
//         std::cout << "haplotype contig = " << haplotypeContig << std::endl;
    }
//     std::cout << back(haplotype)._journalEntries << std::endl;
}


template <typename TString>
void applySimulationInstructions(TString & read, ReadSimulationInstruction const & inst, IlluminaOptions const & options, IlluminaReads const &)
{
    typedef typename Value<TString>::Type TAlphabet;

    SEQAN_ASSERT_EQ(length(inst.qualities), length(inst.editString));
    
    TString tmp;
    reserve(tmp, length(read) + inst.insCount - inst.delCount);
    unsigned j = 0;
//     std::cout << "length(read) == " << length(read) << std::endl;
//     std::cout << "beg = " << inst.beginPos << std::endl;
//     std::cout << "end = " << inst.endPos << std::endl;
//     std::cout << "#   = " << inst.endPos - inst.beginPos << std::endl;
//     std::cout << "insCount = " << inst.insCount << std::endl;
//     std::cout << "delCount = " << inst.delCount << std::endl;
//     for (unsigned i = 0; i < length(inst.editString); ++i)
//         std::cout << "i == " << i << "  edit string[i] = " << inst.editString[i] << "  len = " << length(inst.editString) << std::endl;
    for (unsigned i = 0; i < length(inst.editString); ++i) {
        SEQAN_ASSERT_LEQ(j, i);

        TAlphabet c;
        switch (inst.editString[i]) {
            case ERROR_TYPE_MATCH:
                SEQAN_ASSERT_LT_MSG(j, length(read), "i = %u", i);
                appendValue(tmp, read[j]);
                assignQualityValue(back(tmp), inst.qualities[i]);
                j += 1;
                break;
            case ERROR_TYPE_MISMATCH:
                c = TAlphabet(mtRandDouble() * (ValueSize<TAlphabet>::VALUE - 1));
                SEQAN_ASSERT_LT_MSG(j, length(read), "i = %u", i);
                if (c == read[j])
                    c = TAlphabet(ordValue(c) + 1);
                appendValue(tmp, c);
                assignQualityValue(back(tmp), inst.qualities[i]);
                j += 1;
                break;
            case ERROR_TYPE_INSERT:
                appendValue(tmp, TAlphabet(mtRandDouble() * ValueSize<TAlphabet>::VALUE));
                assignQualityValue(back(tmp), inst.qualities[i]);
                break;
            case ERROR_TYPE_DELETE:
                j += 1;
                break;
            default:
                SEQAN_ASSERT_FAIL("Invalid error type.");
        }
    }
    SEQAN_ASSERT_EQ(j, length(read));
    SEQAN_ASSERT_GEQ(length(tmp), options.readLength);
//     std::cerr << "read == " << read << std::endl << "tmp  == " << tmp << std::endl;
    resize(tmp, options.readLength, Exact());
//     std::cerr << "read == " << read << std::endl << "tmp  == " << tmp << " (truncated)" << std::endl;
    move(read, tmp);
}

// Build a read simulation instructions for a haplotype.
//
// pick a contig, probability is proportional to the length
// pick a start position, end position = start position + read length
// pick whether to match on the forward or reverse strand
// simulate edit string
// build quality values
//   Juliane Dohm, MPI Solexa quality...
// possibly adjust mate if left read has insert at the beginning or right read has insert at the right
int buildReadSimulationInstruction(
        String<ReadSimulationInstruction> & instructions,
        unsigned const & haplotypeId,
        StringSet<SequenceJournal<String<Dna5>, SortedArray> > const & haplotype,
        String<double> const & relativeContigLengths,
        String<double> const & errorProbabilities,
        IlluminaOptions const & options)
{
    ReadSimulationInstruction inst;
    inst.haplotype = haplotypeId;

    // We have to retry simulation if the mate pair did not fit in.
    bool invalid = false;
    do {
		clear(instructions);
        invalid = false;  // By default, we do not want to repeat.
        // Pick contig id, probability is proportional to the length.
        double x = mtRandDouble();
        for (unsigned i = 0; i < length(relativeContigLengths); ++i) {
            if (x < relativeContigLengths[i]) {
                inst.contigId = i -1;
                break;
            }
        }
        // Pick whether on forward or reverse strand.
        if (options.onlyForward)
            inst.isForward = true;
        else if (options.onlyReverse)
            inst.isForward = false;
        else
            inst.isForward = mtRandDouble() < 0.5 ? true : false;
        // Pick a start position.
        if (length(haplotype[inst.contigId]) < options.readLength) {
            std::cerr << "ERROR: haplotype (== " << length(haplotype[inst.contigId]) << ") < read length!" << std::endl;
            return 1;
        }
        inst.beginPos = mtRand() % (length(haplotype[inst.contigId]) - options.readLength);
        inst.endPos = inst.beginPos + options.readLength;
        // Build edit string.
        buildEditString(inst, options.readLength, errorProbabilities);
        // Build quality values.
        buildQualityValues(inst, errorProbabilities, options);
        // If the number of reads does not equal the number of inserts
        // then we have to adjust the read positions.
        if (inst.delCount != inst.insCount) {
            int delta = static_cast<int>(inst.delCount) - static_cast<int>(inst.insCount);
            //                 std::cout << __LINE__ << " delta == " << delta << std::endl;
            inst.endPos += delta;
            if (inst.endPos > length(haplotype[inst.contigId])) {
                delta = inst.endPos - length(haplotype[inst.contigId]);
                //                     std::cout << __LINE__ << " delta == " << delta << std::endl;
                inst.endPos -= delta;
                inst.beginPos -= delta;
            }
            SEQAN_ASSERT_EQ(inst.endPos - inst.beginPos + inst.insCount - inst.delCount,
                            options.readLength);
        }
        // Append read to result list.
        appendValue(instructions, inst);

        // Maybe create a mate for this read.
        if (options.generateMatePairs) {
            // Pick a library length.
            // TODO(holtgrew): Library length should really be a random variable.
            size_t libraryLength = options.libraryLength;  // TODO(holtgrew): std distributed around mean with error from options
            // Compute start and end position.
            inst.endPos = inst.beginPos + options.readLength + libraryLength;
            inst.beginPos = inst.endPos - options.readLength;
            // Verify that the mate fits right of the originally simulated read.
            size_t contigLength = length(haplotype[inst.contigId]);
            if (inst.beginPos > contigLength || inst.endPos > contigLength) {
                // Mate did not fit!  Remove previously added read and set
                // invalid to true so we repeat this simulation.
                SEQAN_ASSERT_GT(length(instructions), 0u);
                eraseBack(instructions);
                invalid = true;
                if (options.verbose)
                    std::cerr << "INFO: Mate did not fit! Repeating..." << std::endl;
            }
            // Build edit string.
            buildEditString(inst, options.readLength, errorProbabilities);
            // If there were != inserts and deletes then we have to adjust
            // the mate reads' boundaries but keep the library length the
            // same.
            if (inst.delCount != inst.insCount) {
                int delta = static_cast<int>(inst.delCount) - static_cast<int>(inst.insCount);
                //                     std::cout << __LINE__ << " delta == " << delta << std::endl;
                inst.endPos += delta;
                if (inst.endPos > length(haplotype[inst.contigId])) {
                    delta = inst.endPos - length(haplotype[inst.contigId]);
                    //                         std::cout << __LINE__ << " delta == " << delta << std::endl;
                    inst.endPos -= delta;
                    inst.beginPos -= delta;
                }
            }
            // Make sure to keep the library length.
            if (inst.endPos - back(instructions).beginPos != options.libraryLength) {
                int delta = static_cast<__int64>(inst.endPos) - static_cast<__int64>(back(instructions).beginPos) - static_cast<__int64>(options.libraryLength);
                //                     std::cout << __LINE__ << " delta == " << delta << std::endl;
                inst.beginPos -= delta;
                inst.endPos -= delta;
            }
            //                 if (inst.endPos - inst.beginPos + inst.insCount - inst.delCount != options.readLength) {
            //                     for (unsigned i = 0; i < length(inst.editString); ++i)
            //                         std::cout << " " << inst.editString[i];
            //                     std::cout << std::endl;
            //                 }
            SEQAN_ASSERT_EQ(inst.endPos - inst.beginPos + inst.insCount - inst.delCount,
                            options.readLength);
            // Build quality values.
            buildQualityValues(inst, errorProbabilities, options);
            // Append read to result list.
            appendValue(instructions, inst);
        }
    } while (invalid);
	
	if (options.generateMatePairs)
		SEQAN_ASSERT_EQ(length(instructions), 2u);
	else
		SEQAN_ASSERT_EQ(length(instructions), 1u);
	
    return 0;
}

/**
..param.fragmentStore:FragmentStore with the contigs and where to write reads to.
..param.options:Options for the simulation.
..param.errorDistribution:Error distribution, indexed by pos * 4 + ERROR_TYPE_{MATCH,MISMATCH,INSERT,DELETE}.
..param.tag:Tag for specifying reads to simulate.
*/
int simulateReadsMain(FragmentStore<MyFragmentStoreConfig> & fragmentStore,
                      IlluminaOptions const & options,
                      String<double> const & errorProbabilities,
                      IlluminaReads const &) {
    typedef FragmentStore<MyFragmentStoreConfig> TFragmentStore;
    typedef Value<TFragmentStore::TMatePairStore>::Type TMatePairStoreElement;
    
    if (options.verbose)
        std::cerr << "Simulating reads..." << std::endl;

    typedef Position<CharString>::Type TPos;

    // First, we randomly pick the haplotype for each read to be simulated.
    String<unsigned> haplotypeIds;
    reserve(haplotypeIds, options.numReads);
    for (size_t i = 0; i < options.numReads; ++i)
        appendValue(haplotypeIds, mtRand() % options.numHaplotypes);

    // We do not build all haplotypes at once since this could cost a
    // lot of memory.
    //
    // for each haplotype id
    //   simulate haplotype
    //   for each simulation instruction for this haplotype:
    //     build simulated read
    reserve(fragmentStore.readSeqStore, options.numReads, Exact());
    reserve(fragmentStore.readNameStore, options.numReads, Exact());
    char readName[1024];
    char outFileName[151];
    snprintf(outFileName, 150, "%s", toCString(options.outputFile));
    String<bool> flipped;
    for (unsigned haplotypeId = 0; haplotypeId < options.numHaplotypes; ++haplotypeId) {
        std::cerr << "Simulating for haplotype #" << haplotypeId << "..." << std::endl;
        StringSet<SequenceJournal<String<Dna5>, SortedArray> > haplotypeContigs;
        buildHaplotype(haplotypeContigs, fragmentStore, options);

        // Build partial sums over relative contig lengths so we can pick the contigs later on.
        size_t totalLength = 0;
        for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i)
            totalLength += length(fragmentStore.contigStore[i].seq);
        String<double> relativeContigLengths;
        resize(relativeContigLengths, length(fragmentStore.contigStore) + 1, Exact());
        front(relativeContigLengths) = 0.0;
        for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i) {
            double l = static_cast<double>(length(fragmentStore.contigStore[i].seq));
            relativeContigLengths[i + 1] = l / totalLength;
        }
        std::partial_sum(begin(relativeContigLengths), end(relativeContigLengths), begin(relativeContigLengths));
        back(relativeContigLengths) = 1.0;

        // Simulate the reads...
        if (options.verbose)
            std::cerr << "Simulating reads for haplotype #" << haplotypeId << "..." << std::endl;

//         std::cerr << "Journal: " << haplotypeContigs[0]._journalEntries << std::endl;

        for (unsigned j = 0; j < length(haplotypeIds); ++j) {
            if (haplotypeIds[j] != haplotypeId)
                continue;  // Guard against instructions on wrong haplotype.

            // Build simulation instructions.
            String<ReadSimulationInstruction> instructions;
            int res = buildReadSimulationInstruction(instructions, haplotypeId, haplotypeContigs, relativeContigLengths, errorProbabilities, options);
            if (res != 0)
                return res;

            int previousMateNum = 0;
            for (unsigned k = 0; k < length(instructions); ++k) {
                ReadSimulationInstruction const & inst = instructions[k];
                // Apply simulation instructions.
                SEQAN_ASSERT_EQ(length(fragmentStore.readSeqStore), length(fragmentStore.readNameStore));
                // Cut out segment from haplotype.
                String<Dna5Q> read = infix(haplotypeContigs[inst.contigId], inst.beginPos, inst.endPos);
                applySimulationInstructions(read, inst, options, IlluminaReads());
                // Append read sequence to read seq store and mate pair to read
                // name store.  This also yields the read id.  We will generate
                // and append the read name below, depending on the read id.
                unsigned readId;
                if (options.generateMatePairs)
                    readId = appendRead(fragmentStore, read, length(fragmentStore.matePairStore));
                else
                    readId = appendRead(fragmentStore, read);

                // Get expected begin/end position in the original sequence.  If we decide to flip this read later on, we will modify the WIT store.
                TPos origBeginPos = virtualToHostPosition(haplotypeContigs[inst.contigId], inst.beginPos);
                TPos origEndPos = virtualToHostPosition(haplotypeContigs[inst.contigId], inst.endPos);

                // Print info about read and haplotype.
                if (options.veryVerbose) {
                    std::cout << ",-- Read #" << readId << std::endl
                              << "| original infix:  " << infix(fragmentStore.contigStore[inst.contigId].seq, origBeginPos, origEndPos) << std::endl
                              << "| haplotype infix: " << infix(haplotypeContigs[inst.contigId], inst.beginPos, inst.endPos) << std::endl
                              << "| read:            " << read << std::endl
                              << "`-- " << std::endl;
                }

                // Generate read name.
                if (options.generateMatePairs) {
                    // Generate the mate num \in {1, 2}, randomly but consistent so two entries belonging together have different nums.  This also decides about the flipping.
                    int mateNum = 3 - previousMateNum;
                    if (readId % 2 == 0) {
                        mateNum = mtRand() % 2 + 1;
						SEQAN_ASSERT_GEQ(mateNum, 1);
						SEQAN_ASSERT_LEQ(mateNum, 2);
                        previousMateNum = mateNum;
                        appendValue(flipped, mateNum == 2);
                    } else {
						SEQAN_ASSERT_EQ(flipped[readId - 1], mateNum == 1);
                        appendValue(flipped, mateNum == 1);
                    }
                    sprintf(readName, "%s.%09u/%d contig=%s haploid=%u length=%lu orig_begin=%lu orig_end=%lu edit_string=", outFileName, readId / 2, mateNum, toCString(fragmentStore.contigNameStore[inst.contigId]), haplotypeId, length(read), origBeginPos, origEndPos);
                } else {
                    sprintf(readName, "%s.%09u contig=%s haploid=%u length=%lu orig_begin=%lu orig_end=%lu edit_string=", outFileName, readId, toCString(fragmentStore.contigNameStore[inst.contigId]), haplotypeId, length(read), origBeginPos, origEndPos);
                }
                for (unsigned i = 0; i < length(inst.editString); ++i) {
                    char buffer[2] = "*";
                    buffer[0] = "MEID"[static_cast<int>(inst.editString[i])];
                    strcat(readName, buffer);
                }
                appendValue(fragmentStore.readNameStore, readName);

                // Tentatively add matches to aligned read store.  We will
                // maybe flip begin and end position below in the "flipping and
                // reordering" step and convert the matches to a global
                // alignment in the "convertMatchesToGlobalAlignment" call.
                if (options.generateMatePairs)
                    appendAlignedRead(fragmentStore, readId, inst.contigId, origBeginPos, origEndPos, length(fragmentStore.matePairStore));
                else
                    appendAlignedRead(fragmentStore, readId, inst.contigId, origBeginPos, origEndPos);

                // Perform flipping and reordering.
                if (options.generateMatePairs) {
                    if (readId % 2 == 1) {  // Only flip and append mate pair info after simulating second mate.
                        // Append mate pair element to fragment store's mate pair store.
                        TMatePairStoreElement matePair;
                        matePair.readId[0] = readId - 1 + flipped[readId];
                        matePair.readId[1] = readId - flipped[readId];
                        appendValue(fragmentStore.matePairStore, matePair);

                        // The first mate always comes from the forward strand.
                        append(fragmentStore.readNameStore[readId - 1], " strand=forward");
                        // The second read always comes from the reverse strand.
                        reverseComplementInPlace(fragmentStore.readSeqStore[readId]);
                        // Note: readId is also last index of aligned read store because we only have one alignment per read!
                        std::swap(fragmentStore.alignedReadStore[readId].beginPos, fragmentStore.alignedReadStore[readId].endPos);
                        append(fragmentStore.readNameStore[readId], " strand=reverse");

                        // Maybe, we write out the mates in different order, i.e. flip the entries in the stores.
                        if (flipped[readId]) {
                            SEQAN_ASSERT_TRUE(flipped[readId - 1]);
                            std::swap(fragmentStore.readSeqStore[readId - 1], fragmentStore.readSeqStore[readId]);
                            std::swap(fragmentStore.readNameStore[readId - 1], fragmentStore.readNameStore[readId]);
                            std::swap(fragmentStore.alignedReadStore[readId - 1], fragmentStore.alignedReadStore[readId]);
                            std::cout << "flipped" << std::endl;
                            std::cout << fragmentStore.readNameStore[readId - 1] << std::endl << fragmentStore.readNameStore[readId] << std::endl;
                            std::cout << fragmentStore.alignedReadStore[readId - 1] << std::endl << fragmentStore.alignedReadStore[readId] << std::endl;
                        } else {
                            SEQAN_ASSERT_NOT(flipped[readId - 1]);
                            std::cout << "not flipped" << std::endl;
                            std::cout << fragmentStore.readNameStore[readId - 1] << std::endl << fragmentStore.readNameStore[readId] << std::endl;
                            std::cout << fragmentStore.alignedReadStore[readId - 1] << std::endl << fragmentStore.alignedReadStore[readId] << std::endl;
                        }
                    }
                } else {
                    if (mtRandDouble() < 0.5) {
                        reverseComplementInPlace(back(fragmentStore.readSeqStore));
                        append(back(fragmentStore.readNameStore), " strand=reverse");
                        // Note: readId is also last index of aligned read store because we only have one alignment per read!
                        std::swap(fragmentStore.alignedReadStore[readId].beginPos, fragmentStore.alignedReadStore[readId].endPos);
                    } else {
                        append(back(fragmentStore.readNameStore), " strand=forward");
                    }
                }
            }
			if (options.generateMatePairs) {
                // When generating mate pairs, an even number of reads is generated in each step.
 				SEQAN_ASSERT_EQ(length(fragmentStore.alignedReadStore) % 2, 0u);
				SEQAN_ASSERT_EQ(length(fragmentStore.readNameStore) % 2, 0u);
				SEQAN_ASSERT_EQ(length(fragmentStore.readSeqStore) % 2, 0u);
			}
		}
    }

    // Last but not least, convert the matches collected before to a global alignment.
    convertMatchesToGlobalAlignment(fragmentStore, Score<int, EditDistance>());
    
    if (options.verbose)
        std::cerr << "Simulated " << length(fragmentStore.readSeqStore) << " reads" << std::endl;

    return 0;
}


int simulateReads(IlluminaOptions options, CharString referenceFilename, IlluminaReads const &) {
    // Print options.
    std::cerr << options;
    std::cerr << "reference file: " << referenceFilename << std::endl;

    // Seed RNG.
    std::srand(options.seed);
    mtRandInit(false);

    // Load or randomly generate the reference sequence.
    FragmentStore<MyFragmentStoreConfig> fragmentStore;
    if (options.useRandomSequence) {
        // TODO(holtgrew): Make setting random output file name possible.
        referenceFilename = "random.fasta";
        std::cerr << "Generating random sequence of length " << options.randomSourceLength
                  << " to file \"" << referenceFilename << "\"." << std::endl;
        int ret = writeRandomSequence(options.randomSourceLength, referenceFilename);
        if (ret != 0)
            return ret;
    }
    // Generate output file name if necessary.
    if (options.outputFile == "") {
        options.outputFile = referenceFilename;
        append(options.outputFile, ".fastq");
    }
    if (options.samFile == "") {
        options.samFile = options.outputFile;
        append(options.samFile, ".sam");
    }
    std::cerr << "Loading reference sequence from \"" << referenceFilename << "\"" << std::endl;
    if (!loadContigs(fragmentStore, referenceFilename)) {
        std::cerr << "Could not load reference sequence from " << referenceFilename << std::endl;
        return 1;
    }

    // Load error probabilities from file or set from built-in data.
    String<double> mismatchProbabilities;
    resize(mismatchProbabilities, options.readLength, Exact());
    if (options.errorDistributionFile == "") {
        if (options.readLength == 36) {
            for (unsigned i = 0; i < 36; ++i)
                mismatchProbabilities[i] = ERROR_PROBABILITIES_36[i];
        } else if (options.readLength == 50) {
            for (unsigned i = 0; i < 50; ++i)
                mismatchProbabilities[i] = ERROR_PROBABILITIES_50[i];
        } else if (options.readLength == 100) {
            for (unsigned i = 0; i < 100; ++i)
                mismatchProbabilities[i] = ERROR_PROBABILITIES_100[i];
        } else {
            std::cerr << "WARNING: No error distribution file given and n != 36, 50, 100.  Using uniform distribution." << std::endl;
            double sum = 0;
            if (options.readLength > 1) {
                for (unsigned i = 0; i < options.readLength - 1; ++i) {
                    sum += 1.0 / options.readLength;
                    mismatchProbabilities[i] = 1.0 / options.readLength;
                }
            }
            mismatchProbabilities[options.readLength - 1] = 1 - sum;
        }
    } else {
        std::cerr << "Cannot load error distribution from file yet!" << std::endl;
        return 1;
    }
	// Prepare log error distribution;
    String<double> errorDistribution;
	resize(errorDistribution, 4 * options.readLength, Exact());
	// Log probs for seeing 1s at positions 0...optionMaxN-1.
	double remainingProb = 1.0 - options.probabilityInsert - options.probabilityDelete;
    for (unsigned j = 0; j < options.readLength; ++j) {
		errorDistribution[j * 4 + ERROR_TYPE_MISMATCH] = mismatchProbabilities[j];
		errorDistribution[j * 4 + ERROR_TYPE_INSERT]   = options.probabilityInsert;
		errorDistribution[j * 4 + ERROR_TYPE_DELETE]   = options.probabilityDelete;
		errorDistribution[j * 4 + ERROR_TYPE_MATCH]    = remainingProb - mismatchProbabilities[j];
	}

    // Kick off the read generation.
    int ret = simulateReadsMain(fragmentStore, options, errorDistribution, IlluminaReads());
    if (ret != 0)
        return ret;

    // Write out results.
    if (options.generateMatePairs) {
        // Build filename with '.1.' infix.
        CharString mateFilename = options.outputFile;
        size_t dotPos = 0;
        for (size_t i = 0; i < length(mateFilename); ++i)
            if (mateFilename[i] == '.')
                dotPos = i;
        infix(mateFilename, dotPos, dotPos + 1) = ".1.";
        // Write out first mates.
        std::cerr << "Writing resulting reads to \"" << mateFilename << "\" mates/1" << std::endl;
        StringSet<String<Dna5Q>, Dependent<> > reads;
        StringSet<CharString, Dependent<> > readNames;
        for (size_t i = 0; i < length(fragmentStore.readNameStore); i += 2) {
            appendValue(readNames, fragmentStore.readNameStore[i]);
            appendValue(reads, fragmentStore.readSeqStore[i]);
        }
        {
            std::fstream fstrm(toCString(mateFilename), std::ios_base::out);
            if (not fstrm.is_open()) {
                std::cerr << "Could not open out file \"" << mateFilename << "\"" << std::endl;
                return 1;
            }
            write(fstrm, readNames, reads, Fastq());
        }
        // Build filename with '.2.' infix.
        infix(mateFilename, dotPos + 1, dotPos + 2) = "2";
        // Write out second mates.
        std::cerr << "Writing resulting reads to \"" << mateFilename << "\" mates/2" << std::endl;
        clear(reads);
        clear(readNames);
        for (size_t i = 1; i < length(fragmentStore.readNameStore); i += 2) {
            appendValue(readNames, fragmentStore.readNameStore[i]);
            appendValue(reads, fragmentStore.readSeqStore[i]);
        }
        {
            std::fstream fstrm(toCString(mateFilename), std::ios_base::out);
            if (not fstrm.is_open()) {
                std::cerr << "Could not open out file \"" << mateFilename << "\"" << std::endl;
                return 1;
            }
            write(fstrm, readNames, reads, Fastq());
        }
    } else {
        std::cerr << "Writing resulting reads to \"" << options.outputFile << "\"" << std::endl;
        std::fstream fstrm(toCString(options.outputFile), std::ios_base::out);
        if (not fstrm.is_open()) {
            std::cerr << "Could not open out file \"" << options.outputFile << "\"" << std::endl;
            return 1;
        }
        write(fstrm, fragmentStore.readNameStore, fragmentStore.readSeqStore, Fastq());
    }
    std::cerr << "Writing SAM file \"" << options.samFile << "\"" << std::endl;
    {
        std::fstream fstrm(toCString(options.samFile), std::ios_base::out);
        if (not fstrm.is_open()) {
            std::cerr << "Could not open SAM file \"" << options.samFile << "\"" << std::endl;
            return 1;
        }
        write(fstrm, fragmentStore, SAM());
    }
    return 0;
}

#endif  // SIMULATE_ILLUMINA_H_
