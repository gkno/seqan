// FRAGMENT(includes)
#include <iostream>
#include <seqan/store.h>
#include <seqan/misc/misc_svg.h>

using namespace seqan;

// FRAGMENT(import)
int main ()
{
    FragmentStore<> store;
    loadContigs(store, "ex1.fa");
    std::ifstream file("ex1.sam");
    read(file, store, Sam());

// FRAGMENT(ascii)
    AlignedReadLayout layout;
    layoutAlignment(layout, store);
    for (unsigned i = 0; i < length(store.contigStore); ++i)
        printAlignment(std::cout, Raw(), layout, store, i, 50, 150, 0, 30);

// FRAGMENT(svg)
    SVGFile svg("layout.svg");
    for (unsigned i = 0; i < length(store.contigStore); ++i)
        printAlignment(svg, Raw(), layout, store, i, 50, 150, 0, 30);

	return 0;
}
