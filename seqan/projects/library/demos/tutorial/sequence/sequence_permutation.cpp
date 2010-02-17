#include <iostream>

#include <seqan/file.h>
#include <seqan/sequence.h>

using namespace seqan;

// Helper function for printPermutations().
void printPermutationsRec(String<char> &current, int pos) {
  if (pos < length(current)) {
    for (char c = 'a'; c <= 'z'; ++c) {
      current[pos] = c;
      printPermutationsRec(current, pos + 1);
    }
  } else {
    std::cout << current << std::endl;;
  }
  
  current[pos] = 'a';
}


// Print all permutations of the alphabet {a, ..., z} of length len.
void printPermutations(int len) {
  String<char> current;
  // Build first value;
  for (int i = 0; i < len; ++i)
    current += "a";
  // Now, current == "a" * len.
  
  // To iterate: human -- to recurse: DIVINE.
  printPermutationsRec(current, 0);
}


int main(int argc, char **argv) {
  printPermutations(3);
  return 0;
}