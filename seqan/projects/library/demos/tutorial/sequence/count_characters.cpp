#include <iostream>
#include <seqan/sequence.h>

void countOneMers(const CharString &str) {
  String<int> table;
  resize(table, 'z' - 'a' + 1);

  for (int i = 0; i < length(str); ++i) {
    table[str[i] - 'a'] += 1;
  }

  for (int i = 0; i < 'z' - 'a' + 1; ++i) {
    std::cout << 'a' + i << " " << table[i] << std::endl;
  }
}

int main(int argc, char **argv) {
  std::cout << "String: hello world" << std::endl;
  countOneMers("hello world");
  std::cout << "String: mississippi" << std::endl;
  countOneMers("mississippi");
  std::cout << "String: banana" << std::endl;
  countOneMers("banana");
  return 0;
}
