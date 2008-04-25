#include <iostream>
#include <seqan/index.h>

using namespace std;
using namespace seqan;

template <typename TSpec>
void testPizzaChili() {
///The following code creates a Pizza & Chili index and assigns it a text.
    typedef Index<String<char>, PizzaChili<TSpec> > index_t;
    index_t index_pc;
    indexText(index_pc) = "This is the best test with a bast jest.";

///Of course, we can access the text as usual:
    cout << indexText(index_pc) << endl;

///Now we may create a default finder and search for a substring. The index
///is only now created because its evaluation is lazy. Immediately after
///the index has been created, the $indexText$ is discarded to save memory.
///Notice that the results returned by the finder might not be in the order
///of their occurrence in the text.
    Finder<index_t> finder(index_pc);
    while (find(finder, "est"))
        cout << "Hit at position " << position(finder) << endl;

///We may query the text of the index. Notice that this returns a string
///without any real content. The string is lazily evaluated in order to
///save memory. Only the substring we are actually displaying will be
///loaded into memory.
/// $indexText(..)$ is a shortcut for $getFibre(.., PizzaChili_Text())$.
    typename Fibre<index_t, PizzaChili_Text>::Type text = indexText(index_pc);
    cout << "infix(text, 12, 21): " << infix(text, 12, 21) << endl;

///We can save the index structure on disc and load it again.
///Notice, however, that not all Pizza & Chili libraries support saving
///and loading at the moment. Please refer to the documentation of the
///different @Tag.Pizza & Chili Index Tags@ for details.
    save(index_pc, "pizzachili");
    index_t index2;
    open(index2, "pizzachili");
    cout << indexText(index2) << endl;
}

int main() {
    cout << "Test the alphabet-friendly FM index:" << endl;
    testPizzaChili<PizzaChili_AF>();
    cout << endl << "Test the compressed compact suffix array index:" << endl;
    testPizzaChili<PizzaChili_CCSA>();
    return 0;
}
