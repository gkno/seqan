#ifndef READ_ANALYZER_READ_ANALYZER_H_
#define READ_ANALYZER_READ_ANALYZER_H_

using namespace seqan;

/*
.Class.ReadEvaluationResult
..summary:Stores base and quality counts for reads.
*/
struct ReadEvaluationResult
{
    String<size_t> baseCountOverall;  // arr[base]
    String<String<size_t> > baseCountPerPosition;  // arr[base][pos]
    String<String<String<size_t> > > qualityCountsPerPositionAndBase;  // arr[base][quality][pos]

    ReadEvaluationResult() {}
};

/**
.Function.setReadLength
..summary:Set the maximal read lenth for an EvaluationResult object.  This adjusts the internal buffer sizes and has to be called before using the object.
*/
void setReadLength(ReadEvaluationResult & result, unsigned readLength)
{
    clear(result.baseCountOverall);
    fill(result.baseCountOverall, 5, 0);

    clear(result.baseCountPerPosition);
    resize(result.baseCountPerPosition, 5);
    for (unsigned i = 0; i < 5; ++i)
        fill(result.baseCountPerPosition[i], readLength, 0);

    clear(result.qualityCountsPerPositionAndBase);
    resize(result.qualityCountsPerPositionAndBase, 5);
    for (unsigned i = 0; i < 5; ++i) {
        resize(result.qualityCountsPerPositionAndBase[i], 63);
        for (unsigned j = 0; j < 63; ++j) {
            fill(result.qualityCountsPerPositionAndBase[i][j], readLength, 0);
        }
    }
}

/**
.Function.countBaseWithQualityAtPosition
..summary:Tally the given (quality, base, position) combination in the EvaluationResult object.
*/
inline void countBaseWithQualityAtPosition(ReadEvaluationResult & result, Dna5Q base, size_t position)
{
    // Count overall bases.
    result.baseCountOverall[ordValue(base)] += 1;

    // Count overall bases per position.
    result.baseCountPerPosition[ordValue(base)][position] += 1;

    // Count overall quality per base per position.
    int b = ordValue(base);
    int quality = getQualityValue(base);
    result.qualityCountsPerPositionAndBase[b][quality][position] += 1;
}

/**
.Function.performReadEvaluation
..summary:Perform an evaluation of the read base counts and qualities.
*/
template <typename TFragmentStore>
void performReadEvaluation(ReadEvaluationResult & result, TFragmentStore & fragmentStore)
{
    typedef typename TFragmentStore::TReadSeqStore TReadSeqStore;
    typedef typename Iterator<TReadSeqStore, Standard>::Type TReadSeqStoreIterator;

    for (TReadSeqStoreIterator it = begin(fragmentStore.readSeqStore); it != end(fragmentStore.readSeqStore); ++it) {
        for (unsigned i = 0; i < length(*it); ++i)
            countBaseWithQualityAtPosition(result, (*it)[i], i);
    }
}

/**
.Function.printReadEvaluationResults
..summary:Print statistical metrics about the read evaluation results.
*/
void printReadEvaluationResults(ReadEvaluationResult const & result)
{
    std::cout << "Overall Base Frequencies" << std::endl;
    size_t sum = 0;
    for (unsigned i = 0; i < 5; ++i)
        sum += result.baseCountOverall[i];
    std::cout << "base   ratio      count" << std::endl;
    for (unsigned i = 0; i < 5; ++i) {
        std::cout << "   " << Dna5(i) << " ";
        printf("%6.2f%%  %9lu\n", 100.0 * result.baseCountOverall[i] / sum, result.baseCountOverall[i]);
    }

    std::cout << std::endl << std::endl << "Base Frequencies Per Position" << std::endl;
    std::cout << "position      A      C      G      T      N" << std::endl;
    for (unsigned i = 0; i < length(result.baseCountPerPosition[0]); ++i) {  // position
        size_t sum = 0;
        for (unsigned j = 0; j < 5; ++j)  // base
            sum += result.baseCountPerPosition[j][i];
        printf("    %4u %6.2f %6.2f %6.2f %6.2f %6.2f\n",
               i,
               100.0 * result.baseCountPerPosition[0][i] / sum,
               100.0 * result.baseCountPerPosition[1][i] / sum,
               100.0 * result.baseCountPerPosition[2][i] / sum,
               100.0 * result.baseCountPerPosition[3][i] / sum,
               100.0 * result.baseCountPerPosition[4][i] / sum);
    }

    std::cout << std::endl << std::endl << "Mean Quality/Std Dev Per Base Per Position" << std::endl;
    std::cout << "position      A   sd A      C   sd C      G   sd G      T   sd T      N   sd N      *   sd *" << std::endl;
    for (unsigned i = 0; i < length(result.baseCountPerPosition[0]); ++i) {  // position
        printf("    %4u", i);
        for (unsigned j = 0; j < 5; ++j) {  // base
            // Compute mean.
            size_t sum = 0;
            size_t count = 0;
            for (unsigned k = 0; k < 63; ++k) {  // qualities
                count += result.qualityCountsPerPositionAndBase[j][k][i];
                sum += k * result.qualityCountsPerPositionAndBase[j][k][i];
            }
            double mean = 1.0 * sum / count;
            if (count == 0) {
                printf(" %6s %6s", "-", "-");
            } else {
                // Compute standard deviation.
                double devSum = 0;
                for (unsigned k = 0; k < 63; ++k) {  // qualities
                    double x = k - mean;
                    devSum += result.qualityCountsPerPositionAndBase[j][k][i] * x * x;
                }
                double stdDev = sqrt(devSum / count);
                printf(" %6.2f %6.2f", mean, stdDev);
            }
        }
        {
            // Compute mean.
            size_t sum = 0;
            size_t count = 0;
            for (unsigned j = 0; j < 5; ++j) {  // base
                for (unsigned k = 0; k < 63; ++k) {  // qualities
                    count += result.qualityCountsPerPositionAndBase[j][k][i];
                    sum += k * result.qualityCountsPerPositionAndBase[j][k][i];
                }
            }
            double mean = 1.0 * sum / count;
            if (count == 0) {
                printf(" %6s %6s", "-", "-");
            } else {
                // Compute standard deviation.
                double devSum = 0;
                for (unsigned j = 0; j < 5; ++j) {  // base
                    for (unsigned k = 0; k < 63; ++k) {  // qualities
                        double x = k - mean;
                        devSum += result.qualityCountsPerPositionAndBase[j][k][i] * x * x;
                    }
                }
                double stdDev = sqrt(devSum / count);
                printf(" %6.2f %6.2f", mean, stdDev);
            }
        }
        printf("\n");
    }
}

/*
Read Alignment Statistics
  base position (0..read length-1)
  error kind {insert, delete, mismatch}
    mismatch type
 */

#endif  // READ_ANALYZER_READ_ANALYZER_H_
