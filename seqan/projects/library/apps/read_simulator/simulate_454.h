#ifndef SIMULATE_454_H_
#define SIMULATE_454_H_

#include "simulate_454_base_calling.h"

// Maximal homopolymer length we will observe.
const unsigned MAX_HOMOPOLYMER_LEN = 40;

// ============================================================================
// Enums, Tags, Classes.
// ============================================================================

struct _LS454Reads;
typedef Tag<_LS454Reads> LS454Reads;

template <>
struct Options<LS454Reads> : public Options<Global>
{
    // Read Length Parameters.

    // Iff true, read lengths follow a uniform distribution, otherwise a
    // standard distribution will be used.
    bool readLengthUniform;

    // Average read length.
    double readLengthMean;

    // For standard distributed read lengths, this is the standard deviation,
    // for uniform read length the interval around the average to use for
    // picking the read lengths.
    double readLengthError;

    // Base Calling Error Model Parameters.

    // If set, $\sigma = k * \sqrt(r)$, otherwise $\sigma = k * r$ is used.
    bool sqrtInStdDev;

    // Proportionality factor for calculating standard deviation proportional
    // to sqrt(homopolymer length).
    double k;

    Options()
            : readLengthUniform(false),
              readLengthMean(400),
              readLengthError(40),
              sqrtInStdDev(false),
              k(0.15)
    {}
};

template <>
struct ReadSimulationInstruction<LS454Reads> : public ReadSimulationInstruction<Global> {
    // For each insertion in the edit string, this string provides the
    // nucleotide types to be inserted at this point.
    String<Dna5> insertionNucleotides;
};

template<>
struct ModelParameters<LS454Reads>
{
    ThresholdMatrix thresholdMatrix;
};

// ============================================================================
// Metafunctions.
// ============================================================================

// ============================================================================
// Functions.
// ============================================================================

template <typename TStream>
TStream & operator<<(TStream & stream, Options<LS454Reads> const & options) {
    stream << static_cast<Options<Global> >(options);
    stream << "454-options {" << std::endl
           << "  readLengthUniform:  " << options.readLengthUniform << std::endl
           << "  readLengthMean:     " << options.readLengthMean << std::endl
           << "  readLengthError:    " << options.readLengthError << std::endl
           << "  sqrtInStdDev:       " << options.sqrtInStdDev << std::endl
           << "  k:                  " << options.k << std::endl
           << "}" << std::endl;
    return stream;
}

void setUpCommandLineParser(CommandLineParser & parser,
                            LS454Reads const &)
{
    setUpCommandLineParser(parser);

    addSection(parser, "454 Read Length Parameters");

    addOption(parser, CommandLineOption("nu",  "read-length-uniform", "If set, the read lengths are simulated with a uniform distribution, with standard distribution otherwise.  Default: false.", OptionType::Bool));
    addOption(parser, CommandLineOption("nm",  "read-length-mean", "The mean of the read lengths.  Default: 400.", OptionType::Double));
    addOption(parser, CommandLineOption("ne",  "read-length-error", "The standard deviation (for standard distribution) and interval length (for uniform distribution) for the read length.  Default: 40.", OptionType::Double));

    addSection(parser, "454 Error Model Parameters");

    addOption(parser, CommandLineOption("sq",  "sqrt-in-std-dev", "If set, no square root is used in error calculation.  Default: Don't use sqrt.", OptionType::Bool));
    addOption(parser, CommandLineOption("k", "proportionality-factor", "Proportionality factor for calculating standard deviation proportional to sqrt(homopolymer length).  Default: 0.15", OptionType::Double));
}

int parseCommandLineAndCheckModelSpecific(Options<LS454Reads> & options,
                                          CommandLineParser & parser)
{
    if (isSetLong(parser, "read-length-uniform"))
        options.readLengthUniform = true;
    if (isSetLong(parser, "read-length-mean"))
        getOptionValueLong(parser, "read-length-mean", options.readLengthMean);
    if (isSetLong(parser, "read-length-error"))
        getOptionValueLong(parser, "read-length-error", options.readLengthError);

    if (isSetLong(parser, "sqrt-in-std-dev"))
        options.sqrtInStdDev = true;
    if (isSetLong(parser, "proportionality-factor"))
        getOptionValueLong(parser, "proportionality-factor", options.k);

    return 0;
}

// For 454 reads, we do not need model specific data (yet?).
int simulateReadsSetupModelSpecificData(ModelParameters<LS454Reads> & parameters,
                                        Options<LS454Reads> const & options)
{
    setK(parameters.thresholdMatrix, options.k);
    setUseSqrt(parameters.thresholdMatrix, options.sqrtInStdDev);
    return 0;
}

unsigned pickReadLength(Options<LS454Reads> const & options)
{
    if (options.readLengthUniform) {
        // Pick uniformly.
        double len = options.readLengthMean - options.readLengthError;
        len += mtRandDouble() * 2 * options.readLengthError;
        return static_cast<unsigned>(round(len));
    } else {
        // Pick from normal distribution.
        double len = normRand(options.readLengthMean, options.readLengthError);
        return static_cast<unsigned>(round(len));
    }
}

template <typename TContig>
void buildSimulationInstructions(ReadSimulationInstruction<LS454Reads> & inst, unsigned readLength, TContig const & contig, ModelParameters<LS454Reads> const & parameters, Options<LS454Reads> const & options)
{
//     std::cout << __FILE__ << ":" << __LINE__ << " -- inst == " << inst << std::endl;

    typedef Iterator<String<Dna5>, Standard>::Type TIterator;
    
    if (inst.endPos == inst.beginPos)
        return;

    //
    // Perform Flowcell Simulation.
    //
    reserve(inst.editString, readLength, Generous());
    clear(inst.editString);
    reserve(inst.qualities, readLength, Generous());
    clear(inst.qualities);

    // Get a copy of the haplotype region we are considering.
    String<Dna5> haplotypeInfix = infix(contig, inst.beginPos, inst.endPos);
    SEQAN_ASSERT_EQ(readLength, inst.endPos - inst.beginPos);

    // In the flow cell simulation, we will simulate light intensities which
    // will be stored in observedIntensities.
    String<double> observedIntensities;
    reserve(observedIntensities, 4 * readLength);
    String<Dna5> observedBases;
    // We also store the real homopolymer length.
    String<unsigned> realBaseCount;

    // Initialize information about the current homopolymer length.
    unsigned homopolymerLength = 0;
    Dna homopolymerType = haplotypeInfix[0];
    while (haplotypeInfix[homopolymerLength] == homopolymerType)
        ++homopolymerLength;

    // Simulate flowcell.
    for (unsigned i = 0, j = 0; i < readLength; ++j, j = j % 4) {  // i indicates first pos of current homopolymer, j indicates flow phase
        if (ordValue(homopolymerType) == j) {
            // Simulate positive flow observation.
            double l = homopolymerLength;
            double sigma = options.k * (options.sqrtInStdDev ? sqrt(l) : l);
            double intensity = _max(0.0, normRand(homopolymerLength, sigma));
            appendValue(observedIntensities, intensity);
            appendValue(realBaseCount, homopolymerLength);
            // Get begin pos and length of next homopolymer.
            i += homopolymerLength;
            homopolymerType = haplotypeInfix[i];
            homopolymerLength = 0;
            while (haplotypeInfix[i + homopolymerLength] == homopolymerType && ((i + homopolymerLength) < (inst.endPos - inst.beginPos)))
                ++homopolymerLength;
        } else {
            // Simulate negative flow observation.
            //
            // Constants taken from MetaSim paper which have it from the
            // original 454 publication.
//             static const double negativeFlowMean = 0.23;
//             static const double negativeFlowStdDev = 0.15;
//             double intensity = lognormRand(negativeFlowMean, negativeFlowStdDev);
            // TODO(holtgrew): Use something that makes more sense for now.
            double intensity = _max(0.0, normRand(0, 0.01));
            appendValue(observedIntensities, intensity);
            appendValue(realBaseCount, 0);
        }
//         std::cout << "observed == " << back(observedIntensities) << ", real == " << back(realBaseCount) << std::endl;
    }

//     std::cout << "infix == " << infix(contig, inst.beginPos, inst.endPos) << std::endl;
//     std::cout << "real base count == ";
//     for (unsigned y = 0; y < length(realBaseCount); ++y)
//         std::cout << realBaseCount[y] << " ";
//     std::cout << std::endl;

    inst.insCount = 0;
    inst.delCount = 0;

    // Call bases, from this build the edit string and maybe qualities.  We
    // only support the "inter" base calling method which was published by
    // the MetaSim authors in the PLOS paper.
    typedef Iterator<String<double>, Standard>::Type IntensitiesIterator;
    int i = 0;  // Flow round, Dna(i) gives base.
    for (IntensitiesIterator it = begin(observedIntensities); it != end(observedIntensities); ++it, ++i) {
        double threshold = getThreshold(parameters.thresholdMatrix, floor(*it), ceil(*it));
        unsigned calledBaseCount = static_cast<unsigned>(*it < threshold ? floor(*it) : ceil(*it));
        // Add any matches.
        unsigned j = 0;
        for (; j < _min(calledBaseCount, realBaseCount[i]); ++j)
            appendValue(inst.editString, ERROR_TYPE_MATCH);
        // Add insertions, if any.
        for (; j < calledBaseCount; ++j) {
            appendValue(inst.insertionNucleotides, Dna(i % 4));
            appendValue(inst.editString, ERROR_TYPE_INSERT);
            inst.insCount += 1;
        }
        // Add deletions, if any.
        for (; j < realBaseCount[i]; ++j) {
            appendValue(inst.editString, ERROR_TYPE_DELETE);
            inst.delCount += 1;
        }
        // Compute likelihood for calling the bases, given this intensity and the Phred score from this.
        double densitySum = 0;
        for (unsigned j = 0; j <= _max(4u, 2 * MAX_HOMOPOLYMER_LEN); ++j)  // Anecdotally through plot in maple: Enough to sum up to 4.
            densitySum += dispatchDensityFunction(parameters.thresholdMatrix, j, *it);
        double x = 0;  // Probability of seeing < (j+1) bases.
//         std::cout << "called base count " << calledBaseCount << std::endl;
        for (j = 0; j < calledBaseCount; ++j) {
            x += dispatchDensityFunction(parameters.thresholdMatrix, j == 0 ? 0.0001 : j, *it);
            unsigned phredScore = -static_cast<unsigned>(10 * ::std::log10(x / densitySum));
            appendValue(inst.qualities, phredScore);
//             std::cout << "j=" << j << " x=" << x << " densitySum=" << densitySum << " phredScore=" << phredScore << " log=" << ::std::log10(x / densitySum) << " p=" << x / densitySum << std::endl;
        }
    }

//     std::cout << __FILE__ << ":" << __LINE__ << " -- inst == " << inst << std::endl;
}

template <typename TString>
void applySimulationInstructions(TString & read, ReadSimulationInstruction<LS454Reads> const & inst, Options<LS454Reads> const & /*options*/)
{
//     std::cout << __FILE__ << ":" << __LINE__ << " -- length(read) == " << length(read) << std::endl;
    typedef typename Value<TString>::Type TAlphabet;

    SEQAN_ASSERT_EQ(length(inst.insertionNucleotides), inst.insCount);
    SEQAN_ASSERT_EQ(length(inst.qualities) + inst.delCount, length(inst.editString));
    SEQAN_ASSERT_EQ(length(read), length(inst.editString) - inst.insCount);
    
    TString tmp;
    reserve(tmp, length(read) + inst.insCount - inst.delCount);
    unsigned j = 0;  // Index in read
    unsigned k = 0;  // Index in inst.insertionNucleotides
    unsigned l = 0;  // Index in inst.qualities
    for (unsigned i = 0; i < length(inst.editString); ++i) {
        SEQAN_ASSERT_LEQ(j, i);

        TAlphabet c;
        switch (inst.editString[i]) {
            case ERROR_TYPE_MATCH:
                SEQAN_ASSERT_LT_MSG(j, length(read), "i = %u", i);
                appendValue(tmp, read[j]);
                assignQualityValue(back(tmp), inst.qualities[l++]);
                j += 1;
                break;
            case ERROR_TYPE_MISMATCH:
                SEQAN_ASSERT_FAIL("No mismatches should occur for 454 reads!");
                break;
            case ERROR_TYPE_INSERT:
                appendValue(tmp, inst.insertionNucleotides[k++]);
                assignQualityValue(back(tmp), inst.qualities[l++]);
                break;
            case ERROR_TYPE_DELETE:
                j += 1;
                break;
            default:
                SEQAN_ASSERT_FAIL("Invalid error type.");
        }
    }
    SEQAN_ASSERT_EQ(j, length(read));

    move(read, tmp);
}

#endif  // SIMULATE_454_H_
