#ifndef SIMULATE_ILLUMINA_H_
#define SIMULATE_ILLUMINA_H_

#include <numeric>

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

    // Basic Read Simulation Parameters.

    // Seed to use for the random number generator.
    unsigned seed;
    // Length of the reads to simulate.
    unsigned readLength;
    // Number of reads to simulate.
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
              seed(0),
              readLength(36),
              numReads(1000),
              useRandomSequence(false),
              randomSourceLength(1000*1000),
              maxErrorsPerRead(4),
              onlyForward(false),
              onlyReverse(false),
              outputFile(""),
              simulateQualities(false),
              generateMatePairs(true),
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
           << "  simulateQualities:      " << (options.simulateQualities ? "true" : "false") << std::endl
           << "  generateMatePairs:      " << (options.generateMatePairs ? "true" : "false") << std::endl
           << "  libraryLength:          " << options.libraryLength << std::endl
           << "  libraryLengthError:     " << options.libraryLengthError << std::endl
           << "  errorDistributionFile:  \"" << options.errorDistributionFile << "\"" <<std::endl
           << "  probabilityInsert:      " << options.probabilityInsert << std::endl
           << "  probabilityDelete:      " << options.probabilityDelete << std::endl
           << "  qualityErrorFactor: " << options.qualityErrorFactor << std::endl
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


void buildEditString(ReadSimulationInstruction & inst, String<double> const & errorProbabilities) {
    size_t sourceLength = inst.endPos - inst.beginPos;  // Length of original data in genome.
    SEQAN_ASSERT_EQ_MSG(length(errorProbabilities), 4 * sourceLength, "Error probability vector len must be == 4 * read len.");
    clear(inst.editString);
    reserve(inst.editString, static_cast<size_t>(1.2 * sourceLength), Generous());
    inst.delCount = 0;
    inst.insCount = 0;

    for (unsigned i = 0; i < sourceLength; /*NOP*/) {
        double x = mtRandDouble();
        double pMatch    = errorProbabilities[i * 4 + ERROR_TYPE_MATCH];
        double pMismatch = errorProbabilities[i * 4 + ERROR_TYPE_MISMATCH];
        double pInsert   = errorProbabilities[i * 4 + ERROR_TYPE_INSERT];
        if (x < pMatch) {
            // match
            appendValue(inst.editString, ERROR_TYPE_MATCH);
            i += 1;
        } else if (x < pMatch + pMismatch) {
            // mismatch
            appendValue(inst.editString, ERROR_TYPE_MISMATCH);
            i += 1;
        } else if (x < pMatch + pMismatch + pInsert) {
            // insert
            inst.insCount += 1;
            if (length(inst.editString) > 0 && back(inst.editString == ERROR_TYPE_DELETE)) {
                resize(inst.editString, length(inst.editString) - 1, Exact());
            } else {
                appendValue(inst.editString, ERROR_TYPE_INSERT);
            }
        } else {
            // Decrement string size, do not add a delete if string is too short.
            if (i == 0)
                continue;
            i -= 1;
            // delete
            if (length(inst.editString) > 0) {
                inst.delCount += 1;
                if (back(inst.editString == ERROR_TYPE_INSERT)) {
                    resize(inst.editString, length(inst.editString) - 1, Exact());
                } else {
                    appendValue(inst.editString, ERROR_TYPE_DELETE);
                }
            }
        }
    }
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


bool operator<(DeltaEntry const & a, DeltaEntry const & b) const {
    return a.pos < b.pos;
}


void buildHaplotype(StringSet<String<Dna5Q> > & haplotype,
                    String<String<std::pair<_t> >
                    FragmentStore<> & fragmentStore,
                    IlluminaOptions const & options) {
    clear(haplotype);
    for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i) {
        appendValue(haplotype, fragmentStore.contigStore[i].seq);
    }
}


/**
..param.fragmentStore:FragmentStore with the contigs and where to write reads to.
..param.options:Options for the simulation.
..param.errorDistribution:Error distribution, indexed by pos * 4 + ERROR_TYPE_{MATCH,MISMATCH,INSERT,DELETE}.
..param.tag:Tag for specifying reads to simulate.
*/
int simulateReadsMain(FragmentStore<> & fragmentStore,
                      IlluminaOptions const & options,
                      String<double> const & errorProbabilities,
                      IlluminaReads const &) {
    if (options.verbose)
        std::cerr << "Simulating reads..." << std::endl;

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

    // Build simulation instructions.
    //   for i = 0..n:
    //     pick a haplotype, uniform distribution
    //     pick a contig, probability is proportional to the length
    //     pick a start position, end position = start position + read length
    //     pick whether to match on the forward or reverse strand
    //     simulate edit string
    //     build quality values
    //       Juliane Dohm, MPI Solexa quality...
    //     possibly adjust mate if left read has insert at the beginning or right read has insert at the right
    String<ReadSimulationInstruction> simulationInstructions;
    reserve(simulationInstructions, options.numReads);
    for (size_t i = 0; i < options.numReads; ++i) {
        // We have to retry simulation if the mate pair did not fit in.
        bool invalid = false;
        do {
            invalid = false;  // By default, we do not want to repeat.
            ReadSimulationInstruction inst;
            // We only consider two haplotypes at the moment.
            inst.haplotype = mtRandDouble() < 0.5 ? 0 : 1;
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
            if (length(fragmentStore.contigStore[inst.contigId].seq) < options.readLength) {
                std::cerr << "ERROR: contig length (== " << length(fragmentStore.contigStore[inst.contigId].seq) << ") < read length!" << std::endl;
                return 1;
            }
            inst.beginPos = mtRand() % (length(fragmentStore.contigStore[inst.contigId].seq) - options.readLength);
            inst.endPos = inst.beginPos + options.readLength;
            // Build edit string.
            buildEditString(inst, errorProbabilities);
            // Build quality values.
            buildQualityValues(inst, errorProbabilities, options);
            // Append read to result list.
            appendValue(simulationInstructions, inst);

            // Maybe create a mate for this read.
            if (options.generateMatePairs) {
                // Pick a library length.
                // TODO(holtgrew): Library length should really be a random variable.
                size_t libraryLength = options.libraryLength;  // std distributed around mean with error from options
                // Compute start and end position.
                inst.beginPos += libraryLength;
                inst.endPos += libraryLength;
                // Verify that the mate fits right of the originally simulated read.
                size_t contigLength = length(fragmentStore.contigStore[inst.contigId].seq);
                if (inst.beginPos > contigLength || inst.endPos > contigLength) {
                    // Mate did not fit!  Remove previously added read and set
                    // invalid to true so we repeat this simulation.
                    SEQAN_ASSERT_GT(length(simulationInstructions), 0u);
                    resize(simulationInstructions, length(simulationInstructions) - 1, Exact());
                    invalid = true;
                    if (options.verbose)
                        std::cerr << "INFO: Mate did not fit! Repeating..." << std::endl;
                }
                // Build edit string.
                buildEditString(inst, errorProbabilities);
                // Build quality values.
                buildQualityValues(inst, errorProbabilities, options);
	            // Append read to result list.
    	        appendValue(simulationInstructions, inst);
            }
        } while (invalid);
    }

    // TODO(holtgrew): If mate pairs are to be simulated, change end of left one and beginning of right one to adjust for additional characters because of many deletions.

    // We do not build all haplotypes at once since this could cost a
    // lot of memory.
    //
    // for each haplotype id
    //   simulate haplotype
    //   for each simulation instruction for this haplotype:
    //     build simulated read
    StringSet<String<Dna5q> > readSeqs;
    resize(readSeqs, length(simulationInstructions), Exact());
    StringSet<CharString> readIds;
    resize(readIds, length(simulationInstructions), Exact());
    for (unsigned i = 0; i < options.numHaplotypes; ++i) {
        std::cerr << "Simulating for haplotype #" << i << "..." << std::endl;
        StringSet<String<Dna5Q> > haplotypeContigs;
        buildHaplotype(haplotypeContigs, fragmentStore, options);

        std::cerr << "Simulating reads for haplotype #" << i << "..." << std::endl;
        for (unsigned j = 0; j < length(simulationInstructions); ++j) {
            if (simulationInstructions[j].haplotype != i)
                continue;  // Guard against instructions on wrong haplotype.
            
        }
    }

    return 0;
}


int simulateReads(IlluminaOptions options, CharString referenceFilename, IlluminaReads const &) {
    // Print options.
    std::cerr << options;
    std::cerr << "reference file: " << referenceFilename << std::endl;

    // Generate output file name if necessary.
    if (options.outputFile == "") {
        options.outputFile = referenceFilename;
        append(options.outputFile, ".fastq");
    }

    // Seed RNG.
    std::srand(options.seed);
    mtRandInit(false);

    // Load or randomly generate the reference sequence.
    FragmentStore<> fragmentStore;
    if (options.useRandomSequence) {
        // TODO(holtgrew): Make setting random output file name possible.
        referenceFilename = "random.fasta";
        std::cerr << "Generating random sequence of length " << options.randomSourceLength
                  << " to file \"" << referenceFilename << "\"." << std::endl;
        int ret = writeRandomSequence(options.randomSourceLength, referenceFilename);
        if (ret != 0)
            return ret;
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
    // TODO(holtgrew): Care about writing optimal alignments to a WIT file.
    // TODO(holtgrew): Care about mate-paired reads.
    std::cerr << "Writing resulting reads to \"" << options.outputFile << "\"" << std::endl;
    std::fstream fstrm(toCString(options.outputFile), std::ios_base::out);
    if (not fstrm.is_open()) {
        std::cerr << "Could not open out file \"" << options.outputFile << "\"" << std::endl;
        return 1;
    }
    write(fstrm, fragmentStore.readNameStore, fragmentStore.readSeqStore, Fastq());

    return 0;
}


#endif  // SIMULATE_ILLUMINA_H_

