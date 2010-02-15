#include <iostream>
#include <seqan/basic.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
	typedef String<AminoAcid> TAminoAcidString;
	TAminoAcidString sourceSeq = "MQDRVKRPMNAFIVWSRDQRRKMALEN";
