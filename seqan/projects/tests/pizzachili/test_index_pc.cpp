#include <iostream>

#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/find.h>
#include <seqan/index/index_pizzachili.h>
#include <seqan/index/index_pizzachili_find.h>

using namespace std;
using namespace seqan;

int main() {
    typedef Index<CharString, PizzaChili<> > index_t;
    index_t idx;
    indexText(idx) = "Dies ist ein Test mit Bast";
    // Position:      0    5   10   15   20   25
    cout << indexText(idx) << endl;

    Finder<index_t> finder(idx);
    String<char> needle = "st";
    while (find(finder, needle)) {
        Position<Finder<index_t> >::Type pos = position(finder);
        cout << pos << ":\t" << infix(indexText(idx), pos, length(indexText(idx))) << endl;
    }

    String<char, PizzaChili<> > index_test = "Hey, this is a test.";
    cout << "index_test[13]: " << index_test[13] << endl;
    index_test[13] = 'b';
    cout << index_test << endl;

    cout << infix(indexText(idx), begin(indexText(idx)), end(indexText(idx))) << endl;
    cout << "With positions:" << endl;
    cout << prefix(indexText(idx), 10) << endl;
    cout << suffix(indexText(idx), 10) << endl;
    cout << "With iterators:" << endl;
    cout << prefix(indexText(idx), begin(indexText(idx)) + 10) << endl;
    cout << suffix(indexText(idx), begin(indexText(idx)) + 10) << endl;
}