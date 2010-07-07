#ifndef SIMULATE_ILLUMINA_H_
#define SIMULATE_ILLUMINA_H_

#include <cmath>

#include <seqan/store.h>

#include "read_simulator.h"
#include "simulate_illumina_data.h"

using namespace seqan;

// ============================================================================
// Enums, Tags, Classes.
// ============================================================================

struct _IlluminaReads;
typedef Tag<_IlluminaReads> IlluminaReads;

template<>
struct Options<IlluminaReads> : public Options<Global>
{
    // Length of the reads to simulate.
    unsigned readLength;

    // Base Calling Error Model Parameters.

    // Path to the mismatch error distribution file.  If empty,
    // built-ins (available for n=36, 50, 100) will be used, or uniform
    // distribution is assumed.
    CharString errorDistributionFile;
    // Probability of an insertion.
    double probabilityInsert;
    // Probability of a deletion.
    double probabilityDelete;
    // Standard error in %*100 around error probability for quality simulation.
    double qualityErrorFactor;

    Options()
            : readLength(36),
              errorDistributionFile(""),
              probabilityInsert(0.01),
              probabilityDelete(0.01),
              qualityErrorFactor(0.1)
    {}
};

template<>
struct ModelParameters<IlluminaReads>
{
    String<double> mismatchProbabilities;
    String<double> errorDistribution;
};

// ============================================================================
// Metafunctions.
// ============================================================================

// ============================================================================
// Functions.
// ============================================================================

template <typename TStream>
TStream & operator<<(TStream & stream, Options<IlluminaReads> const & options) {
    stream << static_cast<Options<Global> >(options);
    stream << "illumina-options {" << std::endl
           << "  readLength:             " << options.readLength << std::endl
           << "  errorDistributionFile:  \"" << options.errorDistributionFile << "\"" <<std::endl
           << "  probabilityInsert:      " << options.probabilityInsert << std::endl
           << "  probabilityDelete:      " << options.probabilityDelete << std::endl
           << "  qualityErrorFactor:     " << options.qualityErrorFactor << std::endl
           << "}" << std::endl;
    return stream;
}

template <>
struct ReadSimulationInstruction<IlluminaReads> : public ReadSimulationInstruction<Global> {};

void setUpCommandLineParser(CommandLineParser & parser,
                            IlluminaReads const &)
{
    setUpCommandLineParser(parser);

    addSection(parser, "Illumina Read Lengths");

    addOption(parser, CommandLineOption("n",  "read-length", "The length of the reads to simulate.  Default: 36.", OptionType::Integer | OptionType::Label));
    addHelpLine(parser, "All resulting reads will have the same length.");

    addSection(parser, "Illumina Error Model");

    addOption(parser, CommandLineOption("d",  "error-distribution", "File containing mismatch qualities.  If left blank, defaults are used.  Defaults are available for n = 36, 50, 100.  Default: \"\".", OptionType::String));
    addOption(parser, CommandLineOption("pi", "prob-insert", "Probability of an insertion.  Default: 0.01.", OptionType::Double));
    addOption(parser, CommandLineOption("pd", "prob-delete", "Probability of a deletion.  Default: 0.01.", OptionType::Double));
}

int parseCommandLineAndCheckModelSpecific(Options<IlluminaReads> & options,
                                          CommandLineParser & parser)
{
    if (isSetLong(parser, "read-length"))
        getOptionValueLong(parser, "read-length", options.readLength);

    if (isSetLong(parser, "error-distribution"))
        getOptionValueLong(parser, "error-distribution", options.errorDistributionFile);
    if (isSetLong(parser, "prob-insert"))
        getOptionValueLong(parser, "prob-insert", options.probabilityInsert);
    if (isSetLong(parser, "prob-delete"))
        getOptionValueLong(parser, "prob-delete", options.probabilityDelete);

    return 0;
}

// Called in simulateReads() to setup model specific data between loading
// or generating the reference sequence and actually simulating the reads.
int simulateReadsSetupModelSpecificData(ModelParameters<IlluminaReads> & parameters,
                                         Options<IlluminaReads> const & options)
{
    // Load error probabilities from file or set from built-in data.
    String<double> & mismatchProbabilities = parameters.mismatchProbabilities;
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
    String<double> & errorDistribution = parameters.errorDistribution;
	resize(errorDistribution, 4 * options.readLength, Exact());
	// Log probs for seeing 1s at positions 0...optionMaxN-1.
	double remainingProb = 1.0 - options.probabilityInsert - options.probabilityDelete;
    for (unsigned j = 0; j < options.readLength; ++j) {
		errorDistribution[j * 4 + ERROR_TYPE_MISMATCH] = mismatchProbabilities[j];
		errorDistribution[j * 4 + ERROR_TYPE_INSERT]   = options.probabilityInsert;
		errorDistribution[j * 4 + ERROR_TYPE_DELETE]   = options.probabilityDelete;
		errorDistribution[j * 4 + ERROR_TYPE_MATCH]    = remainingProb - mismatchProbabilities[j];
	}

    return 0;
}

template <typename TRNG>
unsigned pickReadLength(TRNG const &, Options<IlluminaReads> const & options)
{
    return options.readLength;
}

template <typename TRNG, typename TContig>
void buildSimulationInstructions(ReadSimulationInstruction<IlluminaReads> & inst, TRNG & rng, unsigned readLength, TContig const & contig, ModelParameters<IlluminaReads> const & parameters, Options<IlluminaReads> const & options) {
    String<double> errorProbabilities = parameters.errorDistribution;
    
    SEQAN_ASSERT_EQ(readLength * 4, length(errorProbabilities));
    clear(inst.editString);
    reserve(inst.editString, static_cast<size_t>(1.2 * readLength), Generous());
    inst.delCount = 0;
    inst.insCount = 0;

    //
    // Build Edit String.
    //
    for (unsigned i = 0; i < readLength; /*NOP*/) {
        double x = pickRandomNumber(rng, PDF<Uniform<double> >(0, 1));
        double pMatch    = errorProbabilities[i * 4 + ERROR_TYPE_MATCH];
        double pMismatch = errorProbabilities[i * 4 + ERROR_TYPE_MISMATCH];
        double pInsert   = errorProbabilities[i * 4 + ERROR_TYPE_INSERT];
        if (x < pMatch) {
            // match
            i += 1;
            appendValue(inst.editString, ERROR_TYPE_MATCH);
        } else if (x < pMatch + pMismatch) {
            // mismatch
            i += 1;
            appendValue(inst.editString, ERROR_TYPE_MISMATCH);
        } else if (x < pMatch + pMismatch + pInsert) {
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
    SEQAN_ASSERT_EQ(readLength, length(inst.editString) - inst.delCount);


    //
    // Adjust Positions.
    //

    // If the number of reads does not equal the number of inserts
    // then we have to adjust the read positions.
    if (inst.delCount != inst.insCount) {
        int delta = static_cast<int>(inst.delCount) - static_cast<int>(inst.insCount);
        inst.endPos += delta;
        if (inst.endPos > length(contig)) {
            delta = inst.endPos - length(contig);
            inst.endPos -= delta;
            inst.beginPos -= delta;
        }
        SEQAN_ASSERT_EQ(inst.endPos - inst.beginPos + inst.insCount - inst.delCount,
                        readLength);
    }

    //
    // Quality Simulation.
    //

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
            double delta = pickRandomNumber(rng, PDF<Uniform<double> >(0, 1)) * 2 * options.qualityErrorFactor * p;
            double x = p - options.qualityErrorFactor * p + delta;
            int score = static_cast<int>(round(-10 * std::log10(x)));
            if (inst.editString[i] == ERROR_TYPE_MISMATCH)
                score = static_cast<int>(std::ceil(score / 2.0));
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
                inst.qualities[i] = static_cast<int>(round(tmp[0]));
            } else if (j == inst.endPos - inst.beginPos - 1) {
                // Indel at end.
                inst.qualities[i] = static_cast<int>(round(back(tmp)));
            } else {
                // Indel in center.
                inst.qualities[i] = static_cast<int>(round(0.25 * (tmp[j] + tmp[j + 1])));
            }
        } else {
            j += 1;
        }
    }
}


template <typename TRNG, typename TString>
void applySimulationInstructions(TString & read, TRNG & rng, ReadSimulationInstruction<IlluminaReads> const & inst, Options<IlluminaReads> const & options)
{
    typedef typename Value<TString>::Type TAlphabet;

    SEQAN_ASSERT_EQ(length(inst.qualities), length(inst.editString));
    
    TString tmp;
    reserve(tmp, length(read) + inst.insCount - inst.delCount);
    unsigned j = 0;
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
                c = TAlphabet(pickRandomNumber(rng, PDF<Uniform<double> >(0, 1)) * (ValueSize<TAlphabet>::VALUE - 1));
                SEQAN_ASSERT_LT_MSG(j, length(read), "i = %u", i);
                if (c == read[j])
                    c = TAlphabet(ordValue(c) + 1);
                appendValue(tmp, c);
                assignQualityValue(back(tmp), inst.qualities[i]);
                j += 1;
                break;
            case ERROR_TYPE_INSERT:
                appendValue(tmp, TAlphabet(pickRandomNumber(rng, PDF<Uniform<double> >(0, 1)) * ValueSize<TAlphabet>::VALUE));
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

    resize(tmp, options.readLength, Exact());
    move(read, tmp);
}

#endif  // SIMULATE_ILLUMINA_H_
