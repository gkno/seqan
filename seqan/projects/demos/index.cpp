#include <seqan/index.h>
using namespace seqan;

int main ()
{
    Index< String<char> > index("tobeornottobe");
	Finder< Index< String<char> > > finder(index);

	std::cout << "hit at ";
	while (find(finder, "be"))
		std::cout << position(finder) << " ";
	std::cout << std::endl;

	return 0;
}
