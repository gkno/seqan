#include <iostream>
#include <fstream>
#include <cstdio>

#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;
using namespace std;

int main()
{
///Create a fasta file:
/// Open a standard library stream for binary write. 
/// You can use both C++ iostreams or old style FILE pointer.
	FILE * fl = fopen("testfile.fa", "wb");
    write(fl, "aacagtattagaccactaggaccct", "a test file", Fasta());
	close (fl);
///Read a fasta file into a string:
/// Open a stram for binary read.
	fstream fstrm;
	fstrm.open("testfile.fa", ios_base::in | ios_base::binary);
	String<char> fasta_tag;
	String<Dna> fasta_seq;

	readMeta(fstrm, fasta_tag, Fasta());
	cout << fasta_tag << "\n";	//prints "a test file"

	read(fstrm, fasta_seq, Fasta());
	cout << fasta_seq << "\n";	//prints the sequence
	fstrm.close();
///Opens a file using a file reader string:
	String<Dna, FileReader<Fasta> > fr("testfile.fa");
	cout << fr << "\n";			//prints the sequence

}

