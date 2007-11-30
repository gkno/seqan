#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>


int main()
{
	seqan::Peptide prot = "anypeptide";
	std::cout << length(prot) << std::endl;  //output: 10

	prot += "anewend";
	std::cout << prot << std::endl;          //ouput: "ANYPEPTIDEANEWEND"

	return 0;
}
