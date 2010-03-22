//FRAGMENT(main)
#include <iostream>
#include <seqan/align.h>
#include <seqan/graph_msa.h>

using namespace seqan;

int main()
{
//FRAGMENT(init)
	typedef String<AminoAcid> TSequence;
	StringSet<TSequence> seq;
	appendValue(seq,"GKGDPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGE");
	appendValue(seq,"MQDRVKRPMNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRPRRKAKMLPK");
	appendValue(seq,"MKKLKKHPDFPKKPLTPYFRFFMEKRAKYAKLHPEMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDHPDLIQNAKK");
	appendValue(seq,"MHIKKPLNAFMLYMKEMRANVVAESTLKESAAINQILGRRWHALSREEQAKYYELARKERQLHMQLYPGWSARDNYGKKKKRKREK");

//FRAGMENT(alignment)
	Blosum62 sc(-1, -11);
	Graph<Alignment<StringSet<TSequence, Dependent<> > > > aliG(seq);
	globalMsaAlignment(aliG, sc);
	::std::cout << aliG << ::std::endl;

	return 0;
}
