/* Stress test for the sequence journal that randomly modifies a
   string and compares against modifications of the string.
 */
#undef SEQAN_ENABLE_DEBUG
#define SEQAN_ENABLE_DEBUG 1

#include <cstdlib>
#include <sstream>
#include <string>

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/misc/misc_random.h>
#include <seqan/sequence_journal.h>

using namespace seqan;

const unsigned SEED = 42;
const unsigned INITIAL_LENGTH = 10;//1000;//*1000;
const unsigned NUM_CHANGES = 1000*1000;
const unsigned MAX_INSERT = 10000;

#define RAND_CHAR() ('A' + mtRand() % ('Z' - 'A'))

int main(int, char **)
{
    // Initialize random number generators with a fixed seed.
    std::srand(SEED);
    mtRandInit(false);

    // Build random reference and host string.
    String<char> string;
    reserve(string, INITIAL_LENGTH);
    String<char> host;
    reserve(host, INITIAL_LENGTH);
    for (unsigned i = 0; i < INITIAL_LENGTH; ++i) {
        char c = RAND_CHAR();
        appendValue(string, c);
        appendValue(host, c);
    }

    // Output of initial sequences.
//     std::cout << "reference = " << string << std::endl;
//     std::cout << "host = " << host << std::endl;

    typedef SequenceJournal<String<char>, Unbalanced> TSequenceJournal;

    // Construct sequence journal on host.
    TSequenceJournal sequenceJournal(host);

//     unsigned nextId = 0;
//     std::cerr << "digraph {" << std::endl;

    // We will use a string stream to test the string result of tmp.
    {
        std::stringstream tmp;
        tmp << sequenceJournal;
//         SEQAN_ASSERT_EQ(string, tmp.str());
//         std::cout << "string = " << string << std::endl;
//         std::cout << "jrnld  = " << tmp.str() << std::endl;
//         std::cout << "  tree = " << sequenceJournal._journalTree << std::endl;
//         std::cout << "  orig = " << value(sequenceJournal._host) << std::endl;
//         std::cout << "  buff = " << sequenceJournal._insertionBuffer << std::endl;
//         journalTreeToDot(std::cerr, nextId, sequenceJournal._journalTree);
    }

    size_t expectedLength = length(string);

    for (unsigned i = 0; i < NUM_CHANGES; ++i) {
        std::cout << "i == " << i << std::endl;
        unsigned changeType = mtRand() % 3;
        if (changeType == 0) {  // edit
            if (length(string) == 0)
                continue;
            unsigned begin = 0;
            unsigned end = 0;
            while (begin == end) {
                begin = mtRand() % length(string);
                end = mtRand() % (length(string) + 1);
            }
            if (begin > end)
                std::swap(begin, end);            
            unsigned len = end - begin;
            String<char> buffer;
            reserve(buffer, len);
            for (unsigned i = 0; i < len; ++i)
                appendValue(buffer, RAND_CHAR());
            // Perform insert.
//             std::cout << "assignInfix(sequenceJournal, " << begin << ", " << end << ", \"" << buffer << "\")" << std::endl;
            std::cout << "assignInfix(sequenceJournal, " << begin << ", " << end << ", buffer, len(buffer) == " << length(buffer) << ")" << std::endl;
            infix(string, begin, end) = buffer;
//             std::cout << "pre assign infix " << length(sequenceJournal) << std::endl;
            assignInfix(sequenceJournal, begin, end, buffer);
//             std::cout << "post assign infix " << length(sequenceJournal) << std::endl;
        } else if (changeType == 1) {  // insert
            unsigned begin = 0;
            unsigned len = 0;
            while (len == 0) {
                if (length(string) == 0)
                    begin = 0;
                else
                    begin = mtRand() % length(string);
                len = mtRand() % MAX_INSERT + 1;
            }
            expectedLength += len;
            String<char> buffer;
            reserve(buffer, len);
            for (unsigned i = 0; i < len; ++i)
                appendValue(buffer, RAND_CHAR());
            // Perform insert.
//             std::cout << "insert(sequenceJournal, " << begin << ", \"" << buffer << "\")" << std::endl;
            std::cout << "insert(sequenceJournal, " << begin << ", buffer)" << std::endl;
            infix(string, begin, begin) = buffer;
            insert(sequenceJournal, begin, buffer);
        } else if (changeType == 2) {  // delete
            if (length(string) == 0)
                continue;
            unsigned begin = 0;
            unsigned end = 0;
            while (begin == end) {
                begin = mtRand() % length(string);
                end = mtRand() % (length(string) + 1);
            }
            if (begin > end)
                std::swap(begin, end);
            expectedLength -= (end - begin);
            // Perform erase.
//             std::stringstream tmp;
//             tmp << sequenceJournal;
//             std::cout << ",---" << std::endl;
//             std::cout << "| string = " << string << std::endl;
//             std::cout << "| jrnld  = " << tmp.str() << std::endl;
//             std::cout << "|   tree = " << sequenceJournal._journalTree << std::endl;
//             std::cout << "|   orig = " << value(sequenceJournal._host) << std::endl;
//             std::cout << "|   buff = " << sequenceJournal._insertionBuffer << std::endl;
            std::cout << "erase(sequenceJournal, " << begin << ", " << end << ")" << std::endl;
//             std::cout << "`---" << std::endl;
//             std::cout << sequenceJournal._journalTree << std::endl;
            erase(string, begin, end);
            erase(sequenceJournal, begin, end);
//             std::cout << sequenceJournal._journalTree << std::endl;
        } else {
            SEQAN_ASSERT_FAIL("Invalid change type.");
        }

        {
            // Check via stream operator<< into stringstream.
            std::stringstream tmp;
            tmp << sequenceJournal;
            SEQAN_ASSERT_EQ(expectedLength, length(tmp.str()));
            SEQAN_ASSERT_EQ(expectedLength, length(string));
            SEQAN_ASSERT_EQ(string, tmp.str());
//             std::cout << "string = " << string << std::endl << "tmp.str() = " << tmp.str() << std::endl;
            // Check via iterator on the journal string.
            std::string buffer;
            reserve(buffer, length(sequenceJournal));
            typedef Iterator<TSequenceJournal, Standard>::Type TIterator;
//             std::cout << sequenceJournal._journalTree << std::endl;
            for (TIterator it = begin(sequenceJournal), itend = end(sequenceJournal, Standard()); it != itend; ++it) {
                appendValue(buffer, *it);
            }
//             std::cout << "buffer = " << buffer << std::endl;
            SEQAN_ASSERT_EQ(expectedLength, length(buffer));
            SEQAN_ASSERT_EQ(buffer, tmp.str());
        }

        {
            // Check operator+ and operator+= on sequence journal iterators;
            typedef Iterator<TSequenceJournal, Standard>::Type TSequenceJournalIterator;
            typedef Iterator<String<char>, Standard>::Type TCharStringIterator;

            TSequenceJournalIterator sjIt = begin(sequenceJournal);
            TCharStringIterator csIt = begin(string);
            size_t remaining = length(sequenceJournal);

            while (remaining > 1) {
                SEQAN_ASSERT_TRUE(csIt != end(string) - 1);
                size_t len = mtRand() % (remaining + 1);
                remaining -= len;
                if (remaining == 0)
                    break;
                SEQAN_ASSERT_EQ(*(sjIt + len), *(csIt + len));
                sjIt += len;
                csIt += len;
                SEQAN_ASSERT_EQ(*sjIt, *csIt);
            }
        }
    }

//     std::cerr << "}" << std::endl;

    return 0;
}
