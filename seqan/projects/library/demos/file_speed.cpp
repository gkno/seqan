#define SEQAN_PROFILE

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/system.h>

using namespace seqan;
using namespace std;

const int blockSize = 1 << 12;
const int repeats = 1 << 17;

char block1[blockSize] = "This a test string";
char block2[blockSize];

template <typename TFile>
void testThroughput(const char *fileName)
{
	TFile myFile;
	typename aRequest<TFile>::Type req1, req2;

	if (!open(myFile, fileName, OPEN_WRONLY | OPEN_CREATE)) {
		cout << "Could not open for writing" << endl;
		return;
	}

	SEQAN_PROTIMESTART(iotime);

	for (unsigned i = 0; i < repeats; ++i) 
	{
		waitFor(req1);
		awriteAt(myFile, block1, blockSize,    2*i  * blockSize, req1);
		waitFor(req2);
		awriteAt(myFile, block2, blockSize, (2*i+1) * blockSize, req2);
	}
	waitFor(req1);
	waitFor(req2);

	cout << ((repeats*blockSize / (512.0 * 1024.0)) / SEQAN_PROTIMEDIFF(iotime));
	cout << " MB/s" << endl;

	close(myFile);
}

int main() 
{
	testThroughput< FILE* >("file_speedF.bin");
	testThroughput< File< Sync<> > >("file_speedS.bin");
	testThroughput< File< Async<> > >("file_speedA.bin");
	return 0;
}
