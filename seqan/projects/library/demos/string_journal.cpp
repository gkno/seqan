#define SEQAN_TEST_SKEW3

#include <iostream>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/file.h>
#include <seqan/sequence_journal/string_journal_base.h>
#include "/home/dethecor/Resources/seqan/projects/tests/index/test_index_creation.h"

using namespace std;
using namespace seqan;

int main(){

   /*
   String< char > str_alloc = "tobeornottobe";

   Index< String< char > > testindex_alloc( str_alloc );
   std::cout << "\nrequire SA  . . ." << std::endl;
   indexCreate( testindex_alloc, ESA_SA(), Skew3() );

   std::cout << "\nrequire LCP . . ." << std::endl;
   indexCreate( testindex_alloc, ESA_LCP() );

   String< char, Journal< char, Alloc<>, Alloc<>, Sloppy > > str = "tobeornottobe";
   String< char > ins = "wakkawakka";



   Index< String< char, Journal< char, Alloc<>, Alloc<>, Sloppy > > > testindex( str );

   std::cout << "\nrequire SA  . . ." << std::endl;
   indexCreate( testindex, ESA_SA(), Skew7() );



   String< size_t, Journal< size_t, Alloc<>, Alloc<>, Sloppy > > fibre_sa = indexSA( testindex );
   std::cout << "SA (Skew7):" << std::endl;
   for( size_t i = 0; i < length( fibre_sa ); ++i ){
      std::cout << i << " : " << fibre_sa[i] << std::endl;
   }
   String< size_t > bar = indexSA( testindex_alloc );

   std::cout << "\nrequire LCP . . ." << std::endl;
   indexCreate( testindex, ESA_LCP() );

   StringSet< String< char, Journal< char, Alloc<>, Alloc<>, Sloppy > > > set;
   resize( set, 3 );

   String< char, Journal< char, Alloc<>, Alloc<>, Sloppy > > s1 = "tobeornottobefsagfcvjsdabcvikonhweroifhasjhcgiusagfdseahfsdzhfnoasdhfk";
   String< char, Journal< char, Alloc<>, Alloc<>, Sloppy > > s2 = "theplacetogojfsdhajhfsakudefdsjhihkkkwkwwppakhifkgsaoiuhhdkjdjfnklsuchiuhsuerfjsdfkjbnsdjfbskjfzazrjbdsmnfbslahbfls";
   String< char, Journal< char, Alloc<>, Alloc<>, Sloppy > > s3 = "yyyyyyaaaaaasdafhiiiejjjdfhduckjdjofiwmfgsdhipgsa";

   set[0] = s1;
   set[1] = s2;
   set[2] = s3;

   Iterator< StringSet< String< char, Journal< char, Alloc<>, Alloc<>, Sloppy > > > >::Type it_set = begin(set);
   while( it_set != end( set ) ){
    cout << *it_set << endl;
    ++it_set;
   }


   Index< StringSet< String< char, Journal< char, Alloc<>, Alloc<>, Sloppy > > > > myIndex( set );
   //Finder< Index< StringSet< String< char, Journal< char, Alloc<>, Alloc<>, Sloppy > > > > > myFinder( myIndex );

   indexRequire( myIndex, ESA_SA() );
   indexRequire( myIndex, ESA_LCP() );

   insert( 2, set[0], begin(ins), 6 );
   insert( 5, set[1], begin(ins), 4 );
   insert( 7, set[2], begin(ins), 2 );
   insert( 4, set[0], begin(ins), 1 );
   erase( set[0], 3, 5 );
   insert( 7, set[1], begin(ins), 3 );
   insert( 9, set[2], begin(ins), 5 );
   erase( set[1], 4, 8 );
   erase( set[0], 8 );
   erase( set[2], 0, 3 );

   it_set = begin(set);
   while( it_set != end( set ) ){
    cout << *it_set << endl;
    ++it_set;
    it_set += 0;
   }
   */

    StringSet< String< Dna, Journal< Dna, Alloc<>, Alloc<>, Sloppy > > > set;
    resize( set, 2 );

    fstream fstrm;
    fstrm.open("/home/dethecor/Resources/Sequences/dmel_pangolin_iso_E.fasta", ios_base::in | ios_base::binary);

    String< char > fasta_tag;

    readMeta(fstrm, fasta_tag, Fasta());
    cout << fasta_tag << "\n";	//prints "a test file"

    read(fstrm, set[0], Fasta());

    fstrm.close();

    fstrm.open("/home/dethecor/Resources/Sequences/dmel_pangolin_iso_G.fasta", ios_base::in | ios_base::binary);

    readMeta(fstrm, fasta_tag, Fasta());
    cout << fasta_tag << "\n";	//prints "a test file"

    read(fstrm, set[1], Fasta());

    fstrm.close();

    Index< StringSet< String< Dna, Journal< Dna, Alloc<>, Alloc<>, Sloppy > > > > myIndex( set );

    indexRequire( myIndex, ESA_SA() );
    indexRequire( myIndex, ESA_LCP() );
    timespec start,finish;

    String< Dna > ins = "gattaca";

    insert( 256, set[0], begin(ins), 2 );
    //insert( 512, set[1], begin(ins), 1 );

    cout << "Starting measurement!" << endl;

    clock_gettime(CLOCK_REALTIME, &start);

    synchronize_index( myIndex, set );

    clock_gettime(CLOCK_REALTIME, &finish);

    std::cout << "Execution time for synchronization" << std::endl << "\tFrom: " << start.tv_nsec << "ns\tTo: " << finish.tv_nsec << "ns"<< std::endl << "\tDifference: " << ( (double)( finish.tv_nsec - start.tv_nsec ) / 1000000.0 ) << "ms" << std::endl;

    String< SAValue< Index< StringSet< String< char, Journal< char, Alloc<>, Alloc<>, Sloppy > > > > >::Type, Alloc< > > ds_sa;

    String< char > foo = concat( set );
    resize( ds_sa, length( foo ) );

    clock_gettime(CLOCK_REALTIME, &start);

    createSuffixArray( ds_sa, foo, Shawarma< DeepShallow >() );

    clock_gettime(CLOCK_REALTIME, &finish);

    std::cout << "Execution time for rebuild with Shawarma< DeepShallow >" << std::endl << "\tFrom: " << start.tv_nsec << "ns\tTo: " << finish.tv_nsec << "ns"<< std::endl << "\tDifference: " << ( (double)( finish.tv_nsec - start.tv_nsec ) / 1000000.0 ) << "ms" << std::endl;

    return 0;
}

