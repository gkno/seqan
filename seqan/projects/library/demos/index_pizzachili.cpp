#include <iostream>
#include <seqan/index.h>

using namespace std;
using namespace seqan;

int main ()
{
///The following code creates a Pizza & Chili index and assigns it the text
///$"tobeornottobe"$.
    typedef Index<String<char>, PizzaChili<> > index_t;
    index_t index_pc;
    indexText(index_pc) = "tobeornottobe";

///Now we may create a default finder and search for a substring. The index
///is only now created because its evaluation is lazy. Immediately after
///the index has been created, the $indexText$ is discarded to save memory.
    Finder<index_t> finder(index_pc);
    while (find(finder, "be"))
        cout << "Hit at " << position(finder) << endl;

///We may query the text of the index. Notice that this returns a string
///without any real content. The string is lazily evaluated in order to
///save memory. Only the substring we are actually displaying will be
///loaded into memory.
/// $indexText(..)$ is a shortcut for $getFibre(.., PizzaChili_Text())$.
    Fibre<index_t, PizzaChili_Text>::Type text = indexText(index_pc);
    cout << "infix(text, 2, 7): " << infix(text, 2, 7) << endl;

///We can save the index structure in disc and load it again.
    save(index_pc, "pizzachili");
    index_t index2;
    open(index2, "pizzachili");
    cout << indexText(index2) << endl;
}
