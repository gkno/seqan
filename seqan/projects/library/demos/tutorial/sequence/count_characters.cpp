#include <iostream>
#include <seqan/sequence.h>

using namespace seqan;

void countOneMers(CharString const & str) {
	String<int> table;
	resize(table, 'z' - 'a' + 1);
	
	for (int i = 0; i < length(str); ++i) {
		table[str[i] - 'a'] += 1;
	}
	
	for (int i = 0; i < 'z' - 'a' + 1; ++i) {
		if (table[i] == 0)
			continue;
		std::cout << 'a' + i << " " << table[i] << std::endl;
	}
}

int main(int argc, char **argv) {
	std::cout << "String: helloworld" << std::endl;
	countOneMers("helloworld");
	
	std::cout << "String: mississippi" << std::endl;
	countOneMers("mississippi");
	
	std::cout << "String: banana" << std::endl;
	countOneMers("banana");
	
	return 0;
}
