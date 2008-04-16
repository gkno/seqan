#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#define SEQAN_DEBUG_PIZZACHILI
#define SEQAN_DEBUG

#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/find.h>

using namespace std;
using namespace seqan;

//#define main pizzachili_test
//#include "../library/demos/index_pizzachili.cpp"
//#undef main

template <typename TIndexTrait, typename TChar>
struct TestHelper {
    typedef TIndexTrait index_trait;
    typedef TChar char_t;
    typedef String<char_t, PizzaChili<index_trait> > string_t;
    typedef Index<String<char_t>, PizzaChili<index_trait> > index_t;

    template <typename TStr>
    static index_t test_index_create(TStr const& str) {
        index_t idx(str);
        cout << "Index text (set by c'tor):" << endl;
        cout << indexText(idx) << endl;
        
        //indexText(idx) = str;
        setIndexText(idx, str);
        cout << "Index text (set explicitly):" << endl;
        cout << indexText(idx) << endl;

        return idx;
    }

    template <typename TStr>
    static void test_index_find(index_t& idx, TStr const& needle) {
        typedef Finder<index_t> finder_t;

        finder_t finder(idx);

        typedef vector<typename Position<finder_t>::Type> hits_t;
        typedef typename hits_t::const_iterator const_iter_t;
        hits_t hits;
        while (find(finder, needle))
            hits.push_back(position(finder));

        if (hits.empty())
            cout << "No matches found." << endl;
        else {
            cout << "Hits: " << hits.size() << endl;
            cout << "Position, Suffix" << endl;

            string_t found_text = indexText(idx);
            typename Size<string_t>::Type len = length(found_text);

            for (const_iter_t i = hits.begin(); i != hits.end(); ++i)
                cout << setw(8) << *i << ", " << infix(found_text, *i, len) << endl;
        }
    }

    static void test_index_save(index_t& idx, char const* filename) {
        if (save(idx, filename))
            cout << "Index successfully saved." << endl;
        else {
            cout << "ERROR: save(idx, \"" << filename << "\") failed." << endl;
            struct {} ex;
            throw ex;
        }
    }

    static void test_index_load(index_t& idx, char const* filename) {
        if (open(idx, filename))
            cout << "Index successfully loaded." << endl;
        else {
            cout << "ERROR: open(idx, \"" << filename << "\") failed." << endl;
            struct {} ex;
            throw ex;
        }
    }

    static void test_all() {
        try {
            {
                //string text = "This is the best test with a bast jest";
                DnaString text = "GATTACATAG";
                index_t idx = test_index_create(text);

                DnaString needle = "AT";
                test_index_find(idx, needle);
                test_index_save(idx, "indexdata");
            }

            {
                index_t idx;
                test_index_load(idx, "indexdata");

                //string_t needle = "est";
                DnaString needle = "AT";
                test_index_find(idx, needle);
            }
        }
        catch (...) {
            cout << "Aborting test run due to failure." << endl;
        }
    }
};

int main() {
    cout << "Test AF" << endl;
    TestHelper<PizzaChili_AF, Dna>::test_all();

    cout << endl << "Test CCSA" << endl;
    TestHelper<PizzaChili_CCSA, Dna>::test_all();

    cout << endl << "Test FM" << endl;
    TestHelper<PizzaChili_FM, Dna>::test_all();

    cout << endl << "Test RSA" << endl;
    TestHelper<PizzaChili_RSA, Dna>::test_all();

    cout << endl << "Test SA" << endl;
    TestHelper<PizzaChili_SA, Dna>::test_all();

    cout << endl <<"Text SADA" << endl;
    TestHelper<PizzaChili_SADA, Dna>::test_all();

    return 0;
    //test_all();
    //pizzachili_test();

    typedef Index<CharString, PizzaChili<PizzaChili_SA> > index_t;
    index_t idx;
    indexText(idx) = "Dies ist ein Test mit Bast";
    // Position:      0    5   10   15   20   25
    //cout << indexText(idx) << endl;

    Finder<index_t> finder(idx);
    String<char> needle = "st";
    while (find(finder, needle)) {
        //Position<Finder<index_t> >::Type pos = position(finder);
        //cout << pos << ":\t" << infix(indexText(idx), pos, length(indexText(idx))) << endl;
    }

    String<char, PizzaChili<PizzaChili_SA> > index_test = "Hey, this is a test.";
    cout << "index_test[13]: " << index_test[13] << endl;
    index_test[13] = 'b';
    //cout << index_test << endl;

    /*
    cout << infix(indexText(idx), begin(indexText(idx)), end(indexText(idx))) << endl;
    cout << "With positions:" << endl;
    cout << prefix(indexText(idx), 10) << endl;
    cout << suffix(indexText(idx), 10) << endl;
    cout << "With iterators:" << endl;
    cout << prefix(indexText(idx), begin(indexText(idx)) + 10) << endl;
    cout << suffix(indexText(idx), begin(indexText(idx)) + 10) << endl;
    */
}
