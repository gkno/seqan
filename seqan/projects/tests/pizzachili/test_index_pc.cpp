#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>

#define SEQAN_DEBUG_PIZZACHILI
#define SEQAN_DEBUG
#define SEQAN_TEST

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
        cout << "---------- Index text (set by c'tor):" << endl;
        cout << prefix(indexText(idx), 20) << endl;
        
        //indexText(idx) = str;
        setIndexText(idx, str);
        cout << "---------- Index text (set explicitly):" << endl;
        cout << prefix(indexText(idx), 20) << endl;

        return idx;
    }

    template <typename TStr>
    static size_t test_index_find(index_t& idx, TStr const& needle) {
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
            typename Size<TStr>::Type len = length(needle) + 10;

            for (const_iter_t i = hits.begin(); i != hits.end(); ++i)
                cout << setw(8) << *i << ", " << infix(found_text, *i, *i + len) << endl;
        }

        return hits.size();
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

    static void test_const_specs(index_t const& idx) {
        cout << "---------- Test with const fibre: " << endl;
        cout << prefix(indexText(idx), 20) << "..." << endl;
    }

    template <typename TStr>
    static void test_infix(TStr& str) {
        cout << "---------- Infix iterator test:" << endl;
        typename Iterator<TStr>::Type start = begin(str);
        cout << infix(str, start + 10, start + 20) << endl;
        cout << "Prefix iterator test:" << endl;
        cout << prefix(str, start + 10) << endl;
        cout << "Suffix iterator test:" << endl;
        cout << suffix(str, start + 10) << endl;

        cout << "---------- Infix position test:" << endl;
        cout << infix(str, 10, 20) << endl;
        cout << "Prefix position test:" << endl;
        cout << prefix(str, 10) << endl;
        cout << "Suffix iterator test:" << endl;
        cout << suffix(str, 10) << endl;
    }

    static void test_infix() {
        string_t str = "---------- This is the best test with a bast jest.";
        cout << "----- Non-const" << endl;
        test_infix(str);
        cout << "----- Const" << endl;
        string_t const& crstr = str;
        test_infix(crstr);
    }

    static void test_all() {
        try {
            {
                ifstream ifs("../example.txt");
                SEQAN_TASSERT(ifs);

                string text((istreambuf_iterator<char>(ifs)), istreambuf_iterator<char>());
                index_t idx = test_index_create(text);

                cout << "---------- Searching with arbitrary needle ..." << endl;
                string_t needle = "est";
                SEQAN_TASSERT(test_index_find(idx, needle) == 14);

                cout << "---------- Searching with char* needle ..." << endl;
                char const* pc_needle = "Ishmael";
                SEQAN_TASSERT(test_index_find(idx, pc_needle) == 2);

                cout << "---------- Saving:" << endl;
                test_index_save(idx, "indexdata");
                test_const_specs(idx);
            }

            {
                index_t idx;
                cout << "---------- Loading:" << endl;
                test_index_load(idx, "indexdata");
                cout << "---------- After loading: " << endl;
                cout << prefix(indexText(idx), 20) << endl;
            }
        }
        catch (...) {
            cout << "Aborting test run due to failure." << endl;
        }
    }
};

int main() {
    SEQAN_TREPORT("TEST BEGIN");

    cout << "Test AF" << endl;
    TestHelper<PizzaChili_AF, char>::test_all();

    cout << endl << "Test CCSA" << endl;
    TestHelper<PizzaChili_CCSA, char>::test_all();

    cout << endl << "Test FM" << endl;
    TestHelper<PizzaChili_FM, char>::test_all();

    cout << endl << "Test RSA" << endl;
    TestHelper<PizzaChili_RSA, char>::test_all();

    cout << endl << "Test SA" << endl;
    TestHelper<PizzaChili_SA, char>::test_all();

    cout << endl <<"Text SADA" << endl;
    TestHelper<PizzaChili_SADA, char>::test_all();

    cout << endl << "Other tests" << endl;
    TestHelper<PizzaChili_AF, char>::test_infix();

    debug::verifyCheckpoints(const_cast<char*>("projects/library/seqan/index/index_pizzachili.h"));
    debug::verifyCheckpoints(const_cast<char*>("projects/library/seqan/index/index_pizzachili_find.h"));
    debug::verifyCheckpoints(const_cast<char*>("projects/library/seqan/index/index_pizzachili_string.h"));

    SEQAN_TREPORT("TEST END");
}
