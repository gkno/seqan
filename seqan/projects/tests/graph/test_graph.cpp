#include <iostream>
#include <fstream>
#include <typeinfo>
#include <time.h>
#include <cstdio>
#include <vector>
#include <time.h>
#include <string>

#define SEQAN_DEBUG
#define SEQAN_TEST

#define TEST_PATH "projects/tests/graph/"
#define LIB_PATH "projects/library/seqan/graph/"

#include <seqan/graph.h>

using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////

void Test_IdManager() {
//____________________________________________________________________________
// IdManager
	typedef Id<IdManager<> >::Type TIdType;
	IdManager<> idm;
	
	// Obtain Ids
	TIdType id=obtainId(idm);
	SEQAN_TASSERT(id == 0)
	SEQAN_TASSERT(idInUse(idm,0) == true)
	id=obtainId(idm);
	SEQAN_TASSERT(id == 1)
	id=obtainId(idm);
	SEQAN_TASSERT(id == 2)
	id=obtainId(idm);
	SEQAN_TASSERT(id == 3)
	id=obtainId(idm);
	SEQAN_TASSERT(id == 4)
	SEQAN_TASSERT(idInUse(idm,4) == true)
	SEQAN_TASSERT(getIdUpperBound(idm) == 5)
	SEQAN_TASSERT(idCount(idm) == 5)

	// Release Ids
	releaseId(idm,3);
	SEQAN_TASSERT(idInUse(idm,3) == false)
	releaseId(idm,1);
	SEQAN_TASSERT(idInUse(idm,1) == false)
	releaseId(idm,2);
	SEQAN_TASSERT(idInUse(idm,2) == false)
	SEQAN_TASSERT(idCount(idm) == 2)
	SEQAN_TASSERT(getIdUpperBound(idm) == 5)
	releaseId(idm,4);  // Now we can shrink id range
	SEQAN_TASSERT(idCount(idm) == 1)
	SEQAN_TASSERT(getIdUpperBound(idm) == 4)

	// Ids are reused
	id=obtainId(idm);
	id=obtainId(idm);
	id=obtainId(idm);
	id=obtainId(idm);
	SEQAN_TASSERT(getIdUpperBound(idm) == 5)
	releaseId(idm,3);
	SEQAN_TASSERT(getIdLowerBound(idm) == 0)
	releaseId(idm,0);
	SEQAN_TASSERT(idInUse(idm,0) == false)
	SEQAN_TASSERT(getIdLowerBound(idm) == 1)
	releaseId(idm,2);
	SEQAN_TASSERT(idInUse(idm,2) == false)
	id=obtainId(idm);
	SEQAN_TASSERT(id == 2)
	SEQAN_TASSERT(idInUse(idm,2) == true)
	
	// Check copy constructor and assignment operator
	IdManager<> idm2(idm);
	SEQAN_TASSERT(idCount(idm2) == 3)
	releaseAll(idm2);
	SEQAN_TASSERT(idCount(idm2) == 0)
	SEQAN_TASSERT(getIdUpperBound(idm2) == 0)
	SEQAN_TASSERT(getIdLowerBound(idm2) == 0)


//____________________________________________________________________________
// Dummy IdManager

	// Dummy IdManager
	typedef Id<IdManager<void> >::Type TIdType;
	IdManager<void> id_dummy;
	
	// Obtain Ids
	TIdType idd=obtainId(id_dummy);
	SEQAN_TASSERT(idd == 0)
	idd=obtainId(id_dummy); // Always zero
	SEQAN_TASSERT(idd == 0)
	SEQAN_TASSERT(idInUse(id_dummy, 1) == false) // Always false
	SEQAN_TASSERT(idInUse(id_dummy, 2) == false)
	idd=obtainId(id_dummy);
	SEQAN_TASSERT(idd == 0)
	idd=obtainId(id_dummy);
	SEQAN_TASSERT(idd == 0)
	idd=obtainId(id_dummy);
	SEQAN_TASSERT(idd == 0)
	SEQAN_TASSERT(getIdUpperBound(id_dummy) == 5) 
	SEQAN_TASSERT(getIdLowerBound(id_dummy) == 0)
	SEQAN_TASSERT(idCount(id_dummy) == 5) //But: Correct id count

	// Release Ids
	releaseId(id_dummy,3); 
	SEQAN_TASSERT(idCount(id_dummy) == 4)
	releaseId(id_dummy,1); 
	SEQAN_TASSERT(idCount(id_dummy) == 3)
	IdManager<void> id_dummy2(id_dummy);
	SEQAN_TASSERT(idCount(id_dummy2) == 3)
	idd=obtainId(id_dummy2);
	SEQAN_TASSERT(idd == 0)
	id_dummy = id_dummy2;
	SEQAN_TASSERT(idCount(id_dummy) == 4)
	releaseAll(id_dummy);
	SEQAN_TASSERT(idCount(id_dummy) == 0)
	SEQAN_TASSERT(getIdUpperBound(id_dummy) == 0)
	SEQAN_TASSERT(getIdLowerBound(id_dummy) == 0)
}

//////////////////////////////////////////////////////////////////////////////

void Test_EdgeStump() {
//____________________________________________________________________________
// Test all cargoless EdgeStumps in a list
	// No cargo, list, source, id
	EdgeStump<void, true, true, true> es1;
	_assignId(&es1, 4);
	SEQAN_TASSERT(_getId(&es1) == 4)
	assignTarget(&es1, 5);
	SEQAN_TASSERT(getTarget(&es1) == 5)
	target(&es1)=7;
	SEQAN_TASSERT(getTarget(&es1) == 7)
	assignSource(&es1, 20);
	SEQAN_TASSERT(getSource(&es1) == 20)
	source(&es1)=30;
	SEQAN_TASSERT(getSource(&es1) == 30)
	assignCargo(&es1, 15);  //Does nothing, no cargo -> pointer to 0
	SEQAN_TASSERT(getCargo(&es1) == (void*) 0)
	SEQAN_TASSERT(cargo(&es1) == (void*) 0)
	SEQAN_TASSERT(_getId(&es1) == 4)
	EdgeStump<void, true, true, true> const es1_const(es1);
	SEQAN_TASSERT(getCargo(&es1_const) == (void*) 0)
	SEQAN_TASSERT(cargo(&es1_const) == (void*) 0)
	SEQAN_TASSERT(getTarget(&es1_const) == 7)
	SEQAN_TASSERT(target(&es1_const) == 7)
	SEQAN_TASSERT(source(&es1_const) == 30)
	SEQAN_TASSERT(getSource(&es1_const) == 30)
	SEQAN_TASSERT(_getId(&es1_const) == 4)
	EdgeStump<void, true, true, true> es11(es1);
	EdgeStump<void, true, true, true> es12(es1);
	nextT(&es1) = &es11;
	SEQAN_TASSERT(getNextT(&es1) == &es11)
	assignNextT(&es1, &es12);
	SEQAN_TASSERT(getNextT(&es1) == &es12)
	nextS(&es1) = &es11;
	SEQAN_TASSERT(getNextS(&es1) == &es11)
	SEQAN_TASSERT(getNextT(&es1) == &es12)
	assignNextS(&es1, &es12);
	SEQAN_TASSERT(getNextS(&es1) == &es12)

	// No cargo, list, source, no id
	EdgeStump<void, true, true, false> es2;
	_assignId(&es2, 4);
	SEQAN_TASSERT(_getId(&es2) == 0) // No id, always 0
	assignTarget(&es2, 5);
	SEQAN_TASSERT(getTarget(&es2) == 5)
	target(&es2)=7;
	SEQAN_TASSERT(getTarget(&es2) == 7)
	assignSource(&es2, 20);
	SEQAN_TASSERT(getSource(&es2) == 20)
	source(&es2)=30;
	SEQAN_TASSERT(getSource(&es2) == 30)
	assignCargo(&es2, 15);  //Does nothing, no cargo -> pointer to 0
	SEQAN_TASSERT(getCargo(&es2) == (void*) 0)
	SEQAN_TASSERT(cargo(&es2) == (void*) 0)
	SEQAN_TASSERT(_getId(&es2) == 0)
	EdgeStump<void, true, true, false> const es2_const(es2);
	SEQAN_TASSERT(getCargo(&es2_const) == (void*) 0)
	SEQAN_TASSERT(cargo(&es2_const) == (void*) 0)
	SEQAN_TASSERT(getTarget(&es2_const) == 7)
	SEQAN_TASSERT(target(&es2_const) == 7)
	SEQAN_TASSERT(source(&es2_const) == 30)
	SEQAN_TASSERT(getSource(&es2_const) == 30)
	SEQAN_TASSERT(_getId(&es2_const) == 0)
	EdgeStump<void, true, true, false> es21(es2);
	EdgeStump<void, true, true, false> es22(es2);
	nextT(&es2) = &es21;
	SEQAN_TASSERT(getNextT(&es2) == &es21)
	assignNextT(&es2, &es22);
	SEQAN_TASSERT(getNextT(&es2) == &es22)
	nextS(&es2) = &es21;
	SEQAN_TASSERT(getNextS(&es2) == &es21)
	SEQAN_TASSERT(getNextT(&es2) == &es22)
	assignNextS(&es2, &es22);
	SEQAN_TASSERT(getNextS(&es2) == &es22)

	// No cargo, list, no source, id
	EdgeStump<void, true, false, true> es3;
	_assignId(&es3, 4);
	SEQAN_TASSERT(_getId(&es3) == 4)
	assignTarget(&es3, 5);
	SEQAN_TASSERT(getTarget(&es3) == 5)
	target(&es3)=7;
	SEQAN_TASSERT(getTarget(&es3) == 7)
	assignSource(&es3, 20);
	SEQAN_TASSERT(getSource(&es3) == 0)  // No source
	SEQAN_TASSERT(source(&es3) == 0) // No source
	assignCargo(&es3, 15);  //Does nothing, no cargo -> pointer to 0
	SEQAN_TASSERT(getCargo(&es3) == (void*) 0)
	SEQAN_TASSERT(cargo(&es3) == (void*) 0)
	SEQAN_TASSERT(_getId(&es3) == 4)
	EdgeStump<void, true, false, true> const es3_const(es3);
	SEQAN_TASSERT(getCargo(&es3_const) == (void*) 0)
	SEQAN_TASSERT(cargo(&es3_const) == (void*) 0)
	SEQAN_TASSERT(getTarget(&es3_const) == 7)
	SEQAN_TASSERT(target(&es3_const) == 7)
	SEQAN_TASSERT(getSource(&es3_const) == 0)
	SEQAN_TASSERT(source(&es3_const) == 0)
	SEQAN_TASSERT(_getId(&es3_const) == 4)
	EdgeStump<void, true, false, true> es31(es3);
	EdgeStump<void, true, false, true> es32(es3);
	nextT(&es3) = &es31;
	SEQAN_TASSERT(getNextT(&es3) == &es31)
	assignNextT(&es3, &es32);
	SEQAN_TASSERT(getNextT(&es3) == &es32)
	assignNextS(&es3, &es32); // No source
	SEQAN_TASSERT(nextS(&es3) == 0)
	SEQAN_TASSERT(getNextS(&es3) == 0)
	SEQAN_TASSERT(getNextT(&es3) == &es32)

	// No cargo, list, no source, no id
	EdgeStump<void, true, false, false> es4;
	_assignId(&es4, 4);
	SEQAN_TASSERT(_getId(&es4) == 0) // No id
	assignTarget(&es4, 5);
	SEQAN_TASSERT(getTarget(&es4) == 5)
	target(&es4)=7;
	SEQAN_TASSERT(getTarget(&es4) == 7)
	assignSource(&es4, 20);
	SEQAN_TASSERT(getSource(&es4) == 0) // No source
	SEQAN_TASSERT(source(&es4) == 0)
	assignCargo(&es4, 15);  //Does nothing, no cargo -> pointer to 0
	SEQAN_TASSERT(getCargo(&es4) == (void*) 0)
	SEQAN_TASSERT(cargo(&es4) == (void*) 0)
	SEQAN_TASSERT(_getId(&es4) == 0)
	EdgeStump<void, true, false, false> const es4_const(es4);
	SEQAN_TASSERT(getCargo(&es4_const) == (void*) 0)
	SEQAN_TASSERT(cargo(&es4_const) == (void*) 0)
	SEQAN_TASSERT(getTarget(&es4_const) == 7)
	SEQAN_TASSERT(target(&es4_const) == 7)
	SEQAN_TASSERT(getSource(&es4_const) == 0)
	SEQAN_TASSERT(source(&es4_const) == 0)
	SEQAN_TASSERT(_getId(&es4_const) == 0)
	EdgeStump<void, true, false, false> es41(es4);
	EdgeStump<void, true, false, false> es42(es4);
	nextT(&es4) = &es41;
	SEQAN_TASSERT(getNextT(&es4) == &es41)
	assignNextT(&es4, &es42);
	SEQAN_TASSERT(getNextT(&es4) == &es42)
	assignNextS(&es4, &es42);
	SEQAN_TASSERT(getNextS(&es4) == 0)
	SEQAN_TASSERT(nextS(&es4) == 0)

//____________________________________________________________________________
// Test all EdgeStumps with cargo in a list
	// Cargo, list, source, id
	EdgeStump<unsigned int, true, true, true> es5;
	_assignId(&es5, 4);
	SEQAN_TASSERT(_getId(&es5) == 4)
	assignTarget(&es5, 5);
	SEQAN_TASSERT(getTarget(&es5) == 5)
	target(&es5)=7;
	SEQAN_TASSERT(getTarget(&es5) == 7)
	assignSource(&es5, 20);
	SEQAN_TASSERT(getSource(&es5) == 20)
	source(&es5)=30;
	SEQAN_TASSERT(getSource(&es5) == 30)
	assignCargo(&es5, 15); 
	SEQAN_TASSERT(getCargo(&es5) == 15)
	SEQAN_TASSERT(cargo(&es5) == 15)
	SEQAN_TASSERT(_getId(&es5) == 4)
	EdgeStump<unsigned int, true, true, true> const es5_const(es5);
	SEQAN_TASSERT(getCargo(&es5_const) == 15)
	SEQAN_TASSERT(cargo(&es5_const) == 15)
	SEQAN_TASSERT(getTarget(&es5_const) == 7)
	SEQAN_TASSERT(target(&es5_const) == 7)
	EdgeStump<unsigned int, true, true, true> es51(es5);
	EdgeStump<unsigned int, true, true, true> es52(es5);
	nextT(&es5) = &es51;
	SEQAN_TASSERT(getNextT(&es5) == &es51)
	assignNextT(&es5, &es52);
	SEQAN_TASSERT(getNextT(&es5) == &es52)
	nextS(&es5) = &es51;
	SEQAN_TASSERT(getNextS(&es5) == &es51)
	SEQAN_TASSERT(getNextT(&es5) == &es52)
	assignNextS(&es5, &es52);
	SEQAN_TASSERT(getNextS(&es5) == &es52)

	// Cargo, list, source, no id
	EdgeStump<unsigned int, true, true, false> es6;
	_assignId(&es6, 4);
	SEQAN_TASSERT(_getId(&es6) == 0) // No id, always 0
	assignTarget(&es6, 5);
	SEQAN_TASSERT(getTarget(&es6) == 5)
	target(&es6)=7;
	SEQAN_TASSERT(getTarget(&es6) == 7)
	assignSource(&es6, 20);
	SEQAN_TASSERT(getSource(&es6) == 20)
	source(&es6)=30;
	SEQAN_TASSERT(getSource(&es6) == 30)
	assignCargo(&es6, 15); 
	SEQAN_TASSERT(getCargo(&es6) == 15)
	SEQAN_TASSERT(cargo(&es6) == 15)
	SEQAN_TASSERT(_getId(&es6) == 0)
	EdgeStump<unsigned int, true, true, false> const es6_const(es6);
	SEQAN_TASSERT(getCargo(&es6_const) == 15)
	SEQAN_TASSERT(cargo(&es6_const) == 15)
	SEQAN_TASSERT(getTarget(&es6_const) == 7)
	SEQAN_TASSERT(target(&es6_const) == 7)
	EdgeStump<unsigned int, true, true, false> es61(es6);
	EdgeStump<unsigned int, true, true, false> es62(es6);
	nextT(&es6) = &es61;
	SEQAN_TASSERT(getNextT(&es6) == &es61)
	assignNextT(&es6, &es62);
	SEQAN_TASSERT(getNextT(&es6) == &es62)
	nextS(&es6) = &es61;
	SEQAN_TASSERT(getNextS(&es6) == &es61)
	SEQAN_TASSERT(getNextT(&es6) == &es62)
	assignNextS(&es6, &es62);
	SEQAN_TASSERT(getNextS(&es6) == &es62)

	// Cargo, list, no source, id
	EdgeStump<unsigned int, true, false, true> es7;
	_assignId(&es7, 4);
	SEQAN_TASSERT(_getId(&es7) == 4)
	assignTarget(&es7, 5);
	SEQAN_TASSERT(getTarget(&es7) == 5)
	target(&es7)=7;
	SEQAN_TASSERT(getTarget(&es7) == 7)
	assignSource(&es7, 20);
	SEQAN_TASSERT(getSource(&es7) == 0)  // No source
	SEQAN_TASSERT(source(&es7) == 0) // No source
	assignCargo(&es7, 15); 
	SEQAN_TASSERT(getCargo(&es7) == 15)
	SEQAN_TASSERT(cargo(&es7) == 15)
	SEQAN_TASSERT(_getId(&es7) == 4)
	EdgeStump<unsigned int, true, false, true> const es7_const(es7);
	SEQAN_TASSERT(getCargo(&es7_const) == 15)
	SEQAN_TASSERT(cargo(&es7_const) == 15)
	SEQAN_TASSERT(getTarget(&es7_const) == 7)
	SEQAN_TASSERT(target(&es7_const) == 7)
	EdgeStump<unsigned int, true, false, true> es71(es7);
	EdgeStump<unsigned int, true, false, true> es72(es7);
	nextT(&es7) = &es71;
	SEQAN_TASSERT(getNextT(&es7) == &es71)
	assignNextT(&es7, &es72);
	SEQAN_TASSERT(getNextT(&es7) == &es72)
	assignNextS(&es7, &es72); // No source
	SEQAN_TASSERT(nextS(&es7) == 0)
	SEQAN_TASSERT(getNextS(&es7) == 0)
	SEQAN_TASSERT(getNextT(&es7) == &es72)

	// Cargo, list, no source, no id
	EdgeStump<unsigned int, true, false, false> es8;
	_assignId(&es8, 4);
	SEQAN_TASSERT(_getId(&es8) == 0) // No id
	assignTarget(&es8, 5);
	SEQAN_TASSERT(getTarget(&es8) == 5)
	target(&es8)=7;
	SEQAN_TASSERT(getTarget(&es8) == 7)
	assignSource(&es8, 20);
	SEQAN_TASSERT(getSource(&es8) == 0) // No source
	SEQAN_TASSERT(source(&es8) == 0)
	assignCargo(&es8, 15); 
	SEQAN_TASSERT(getCargo(&es8) == 15)
	SEQAN_TASSERT(cargo(&es8) == 15)
	SEQAN_TASSERT(_getId(&es8) == 0)
	EdgeStump<unsigned int, true, false, false> const es8_const(es8);
	SEQAN_TASSERT(getCargo(&es8_const) == 15)
	SEQAN_TASSERT(cargo(&es8_const) == 15)
	SEQAN_TASSERT(getTarget(&es8_const) == 7)
	SEQAN_TASSERT(target(&es8_const) == 7)
	EdgeStump<unsigned int, true, false, false> es81(es8);
	EdgeStump<unsigned int, true, false, false> es82(es8);
	nextT(&es8) = &es81;
	SEQAN_TASSERT(getNextT(&es8) == &es81)
	assignNextT(&es8, &es82);
	SEQAN_TASSERT(getNextT(&es8) == &es82)
	assignNextS(&es8, &es82);
	SEQAN_TASSERT(getNextS(&es8) == 0)
	SEQAN_TASSERT(nextS(&es8) == 0)

//____________________________________________________________________________
// Test all cargoless EdgeStumps in an array
	// No cargo, no list, source, id
	EdgeStump<void, false, true, true> es9;
	_assignId(&es9, 4);
	SEQAN_TASSERT(_getId(&es9) == 4)
	assignTarget(&es9, 5);
	SEQAN_TASSERT(getTarget(&es9) == 5)
	target(&es9)=7;
	SEQAN_TASSERT(getTarget(&es9) == 7)
	assignSource(&es9, 20);
	SEQAN_TASSERT(getSource(&es9) == 20)
	source(&es9)=30;
	SEQAN_TASSERT(getSource(&es9) == 30)
	assignCargo(&es9, 15);  //Does nothing, no cargo -> pointer to 0
	SEQAN_TASSERT(getCargo(&es9) == (void*) 0)
	SEQAN_TASSERT(cargo(&es9) == (void*) 0)
	SEQAN_TASSERT(_getId(&es9) == 4)
	EdgeStump<void, false, true, true> const es9_const(es9);
	SEQAN_TASSERT(getCargo(&es9_const) == (void*) 0)
	SEQAN_TASSERT(cargo(&es9_const) == (void*) 0)
	SEQAN_TASSERT(getTarget(&es9_const) == 7)
	SEQAN_TASSERT(target(&es9_const) == 7)

	// No cargo, no list, source, no id
	EdgeStump<void, false, true, false> es10;
	_assignId(&es10, 4);
	SEQAN_TASSERT(_getId(&es10) == 0) // No id
	assignTarget(&es10, 5);
	SEQAN_TASSERT(getTarget(&es10) == 5)
	target(&es10)=7;
	SEQAN_TASSERT(getTarget(&es10) == 7)
	assignSource(&es10, 20);
	SEQAN_TASSERT(getSource(&es10) == 20)
	source(&es10)=30;
	SEQAN_TASSERT(getSource(&es10) == 30)
	assignCargo(&es10, 15);  //Does nothing, no cargo -> pointer to 0
	SEQAN_TASSERT(getCargo(&es10) == (void*) 0)
	SEQAN_TASSERT(cargo(&es10) == (void*) 0)
	SEQAN_TASSERT(_getId(&es10) == 0)
	EdgeStump<void, false, true, false> const es10_const(es10);
	SEQAN_TASSERT(getCargo(&es10_const) == (void*) 0)
	SEQAN_TASSERT(cargo(&es10_const) == (void*) 0)
	SEQAN_TASSERT(getTarget(&es10_const) == 7)
	SEQAN_TASSERT(target(&es10_const) == 7)

	// No cargo, no list, no source, id
	EdgeStump<void, false, false, true> es_11;
	_assignId(&es_11, 4);
	SEQAN_TASSERT(_getId(&es_11) == 4)
	assignTarget(&es_11, 5);
	SEQAN_TASSERT(getTarget(&es_11) == 5)
	target(&es_11)=7;
	SEQAN_TASSERT(getTarget(&es_11) == 7)
	assignSource(&es_11, 20);
	SEQAN_TASSERT(getSource(&es_11) == 0) // No source
	SEQAN_TASSERT(source(&es_11) == 0)
	assignCargo(&es_11, 15);  //Does nothing, no cargo -> pointer to 0
	SEQAN_TASSERT(getCargo(&es_11) == (void*) 0)
	SEQAN_TASSERT(cargo(&es_11) == (void*) 0)
	SEQAN_TASSERT(_getId(&es_11) == 4)
	EdgeStump<void, false, false, true> const es_11_const(es_11);
	SEQAN_TASSERT(getCargo(&es_11_const) == (void*) 0)
	SEQAN_TASSERT(cargo(&es_11_const) == (void*) 0)
	SEQAN_TASSERT(getTarget(&es_11_const) == 7)
	SEQAN_TASSERT(target(&es_11_const) == 7)

	
	// No cargo, no list, no source, no id
	EdgeStump<void, false, false, false> es_12;
	_assignId(&es_12, 4);
	SEQAN_TASSERT(_getId(&es_12) == 0)
	assignTarget(&es_12, 5);
	SEQAN_TASSERT(getTarget(&es_12) == 5)
	target(&es_12)=7;
	SEQAN_TASSERT(getTarget(&es_12) == 7)
	assignSource(&es_12, 20);
	SEQAN_TASSERT(getSource(&es_12) == 0) // No source
	SEQAN_TASSERT(source(&es_12) == 0)
	assignCargo(&es_12, 15);  //Does nothing, no cargo -> pointer to 0
	SEQAN_TASSERT(getCargo(&es_12) == (void*) 0)
	SEQAN_TASSERT(cargo(&es_12) == (void*) 0)
	SEQAN_TASSERT(_getId(&es_12) == 0)
	EdgeStump<void, false, false, false> const es_12_const(es_12);
	SEQAN_TASSERT(getCargo(&es_12_const) == (void*) 0)
	SEQAN_TASSERT(cargo(&es_12_const) == (void*) 0)
	SEQAN_TASSERT(getTarget(&es_12_const) == 7)
	SEQAN_TASSERT(target(&es_12_const) == 7)

//____________________________________________________________________________
// Test all EdgeStumps with cargo in an array
	// cargo, no list, source, id
	EdgeStump<unsigned int, false, true, true> es_13;
	_assignId(&es_13, 4);
	SEQAN_TASSERT(_getId(&es_13) == 4)
	assignTarget(&es_13, 5);
	SEQAN_TASSERT(getTarget(&es_13) == 5)
	target(&es_13)=7;
	SEQAN_TASSERT(getTarget(&es_13) == 7)
	assignSource(&es_13, 20);
	SEQAN_TASSERT(getSource(&es_13) == 20)
	source(&es_13)=30;
	SEQAN_TASSERT(getSource(&es_13) == 30)
	assignCargo(&es_13, 15);  
	SEQAN_TASSERT(getCargo(&es_13) == 15)
	SEQAN_TASSERT(cargo(&es_13) == 15)
	SEQAN_TASSERT(_getId(&es_13) == 4)
	EdgeStump<unsigned int, false, true, true> const es_13_const(es_13);
	SEQAN_TASSERT(getCargo(&es_13_const) == 15)
	SEQAN_TASSERT(cargo(&es_13_const) == 15)
	SEQAN_TASSERT(getTarget(&es_13_const) == 7)
	SEQAN_TASSERT(target(&es_13_const) == 7)

	// cargo, no list, source, no id
	EdgeStump<unsigned int, false, true, false> es_14;
	_assignId(&es_14, 4);
	SEQAN_TASSERT(_getId(&es_14) == 0) // No id
	assignTarget(&es_14, 5);
	SEQAN_TASSERT(getTarget(&es_14) == 5)
	target(&es_14)=7;
	SEQAN_TASSERT(getTarget(&es_14) == 7)
	assignSource(&es_14, 20);
	SEQAN_TASSERT(getSource(&es_14) == 20)
	source(&es_14)=30;
	SEQAN_TASSERT(getSource(&es_14) == 30)
	assignCargo(&es_14, 15); 
	SEQAN_TASSERT(getCargo(&es_14) == 15)
	SEQAN_TASSERT(cargo(&es_14) == 15)
	SEQAN_TASSERT(_getId(&es_14) == 0)
	EdgeStump<unsigned int, false, true, false> const es_14_const(es_14);
	SEQAN_TASSERT(getCargo(&es_14_const) == 15)
	SEQAN_TASSERT(cargo(&es_14_const) == 15)
	SEQAN_TASSERT(getTarget(&es_14_const) == 7)
	SEQAN_TASSERT(target(&es_14_const) == 7)

	// cargo, no list, no source, id
	EdgeStump<unsigned int, false, false, true> es_15;
	_assignId(&es_15, 4);
	SEQAN_TASSERT(_getId(&es_15) == 4)
	assignTarget(&es_15, 5);
	SEQAN_TASSERT(getTarget(&es_15) == 5)
	target(&es_15)=7;
	SEQAN_TASSERT(getTarget(&es_15) == 7)
	assignSource(&es_15, 20);
	SEQAN_TASSERT(getSource(&es_15) == 0) // No source
	SEQAN_TASSERT(source(&es_15) == 0)
	assignCargo(&es_15, 15); 
	SEQAN_TASSERT(getCargo(&es_15) == 15)
	SEQAN_TASSERT(cargo(&es_15) == 15)
	SEQAN_TASSERT(_getId(&es_15) == 4)
	EdgeStump<unsigned int, false, false, true> const es_15_const(es_15);
	SEQAN_TASSERT(getCargo(&es_15_const) == 15)
	SEQAN_TASSERT(cargo(&es_15_const) == 15)
	SEQAN_TASSERT(getTarget(&es_15_const) == 7)
	SEQAN_TASSERT(target(&es_15_const) == 7)

	// cargo, no list, no source, no id
	EdgeStump<unsigned int, false, false, false> es_16;
	_assignId(&es_16, 4);
	SEQAN_TASSERT(_getId(&es_16) == 0)
	assignTarget(&es_16, 5);
	SEQAN_TASSERT(getTarget(&es_16) == 5)
	target(&es_16)=7;
	SEQAN_TASSERT(getTarget(&es_16) == 7)
	assignSource(&es_16, 20);
	SEQAN_TASSERT(getSource(&es_16) == 0) // No source
	SEQAN_TASSERT(source(&es_16) == 0)
	assignCargo(&es_16, 15);  
	SEQAN_TASSERT(getCargo(&es_16) == 15)
	SEQAN_TASSERT(cargo(&es_16) == 15)
	SEQAN_TASSERT(_getId(&es_16) == 0)
	EdgeStump<unsigned int, false, false, false> const es_16_const(es_16);
	SEQAN_TASSERT(getCargo(&es_16_const) == 15)
	SEQAN_TASSERT(cargo(&es_16_const) == 15)
	SEQAN_TASSERT(getTarget(&es_16_const) == 7)
	SEQAN_TASSERT(target(&es_16_const) == 7)
}


//////////////////////////////////////////////////////////////////////////////

void Test_Directed() {
//____________________________________________________________________________
// Graph without edge cargo but with edge ids

	typedef Graph<> StandardGraph;
	typedef VertexDescriptor<StandardGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<StandardGraph>::Type TEdgeDescriptor;
	
	StandardGraph g;
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	SEQAN_TASSERT(empty(g) == true)

	// Add vertex
	TVertexDescriptor v0 = addVertex(g);
	SEQAN_TASSERT(v0 == 0)
	SEQAN_TASSERT(outDegree(g, v0) == 0)	
	SEQAN_TASSERT(inDegree(g, 0) == 0)
	SEQAN_TASSERT(degree(g, 0) == 0)
	SEQAN_TASSERT(numVertices(g) == 1)
	SEQAN_TASSERT(empty(g) == false)
	
	// Add edge
	TEdgeDescriptor e1 =addEdge(g,v0,v0);
	SEQAN_TASSERT(findEdge(g, v0, v0) == e1)
	SEQAN_TASSERT(_getVertexString(g)[0] == e1)
	SEQAN_TASSERT(getIdUpperBound(_getVertexIdManager(g)) == 1)
	SEQAN_TASSERT(getIdUpperBound(_getEdgeIdManager(g)) == 1)
	SEQAN_TASSERT(targetVertex(g, e1) == 0)
	SEQAN_TASSERT(sourceVertex(g, e1) == 0)  //Expensive in standard graph!
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(outDegree(g, v0) == 1)	
	SEQAN_TASSERT(inDegree(g, v0) == 1)
	SEQAN_TASSERT(degree(g, v0) == 2)
	
	// Add further edges and vertices
	TVertexDescriptor v1 = addVertex(g);
	TEdgeDescriptor e2 =addEdge(g,0,1);
	SEQAN_TASSERT(v1 == 1)
	SEQAN_TASSERT(numVertices(g) == 2)
	SEQAN_TASSERT(targetVertex(g, e2) == 1)
	SEQAN_TASSERT(sourceVertex(g, e2) == 0)
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(outDegree(g, v0) == 2)	
	SEQAN_TASSERT(inDegree(g, 1) == 1)
	SEQAN_TASSERT(inDegree(g, 0) == 1)	
	SEQAN_TASSERT(degree(g, 0) == 3)
		
	// Add more vertices and edges
	addVertex(g);  //2
	TVertexDescriptor v3 = addVertex(g);  //3
	addVertex(g);  //4
	addEdge(g,3,4);
	TEdgeDescriptor my_edge = addEdge(g,3,1);
	addEdge(g,3,0);
	SEQAN_TASSERT(v3 == 3)
	SEQAN_TASSERT(numVertices(g) == 5)
	SEQAN_TASSERT(targetVertex(g, e2) == 1)
	SEQAN_TASSERT(sourceVertex(g, e2) == 0)
	SEQAN_TASSERT(targetVertex(g, my_edge) == 1)
	SEQAN_TASSERT(sourceVertex(g, my_edge) == 3)
	SEQAN_TASSERT(numEdges(g) == 5)
	SEQAN_TASSERT(outDegree(g, v3) == 3)	
	
	// Graph drawing
	removeEdge(g,0,0); // ToDo: Drawing of self edges
	addEdge(g,4,3);
	// Raw output
	// std::cout << g << ::std::endl;
	// File output
	fstream strm;
	strm.open(TEST_PATH "my_graph.dot", ios_base::out | ios_base::trunc);
	write(strm,g,DotDrawing());
	strm.close();
	// File read
	StandardGraph gTmp;
	strm.open(TEST_PATH "my_graph.dot", ios_base::in);
	read(strm,gTmp,DotDrawing());
	strm.close();

	removeEdge(g,4,3);
	addEdge(g,0,0);

	// Remove edges
	removeEdge(g,my_edge);
	removeEdge(g,0,1);
	SEQAN_TASSERT(numEdges(g) == 3)

	// Remove vertices 
	TEdgeDescriptor e3 = addEdge(g,3,3);
	addEdge(g,1,3);
	addEdge(g,0,3);
	addEdge(g,0,4);
	SEQAN_TASSERT(outDegree(g, 0) == 3)
	SEQAN_TASSERT(outDegree(g, 1) == 1)
	SEQAN_TASSERT(targetVertex(g, e3) == 3)
	SEQAN_TASSERT(sourceVertex(g, e3) == 3)
	removeVertex(g, v3);
	SEQAN_TASSERT(outDegree(g, 0) == 2)
	SEQAN_TASSERT(outDegree(g, 1) == 0)
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 2)

	// Clear graph
	clearEdges(g);
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 0)
	addEdge(g,2,0);
	addEdge(g,4,1);
	clearVertices(g);
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	addVertex(g);addVertex(g);addVertex(g);
	addVertex(g);addVertex(g);
	addEdge(g,2,0);
	addEdge(g,4,1);
	clear(g);
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	addVertex(g);addVertex(g);addVertex(g);
	addVertex(g);addVertex(g);
	addEdge(g,2,0);
	addEdge(g,4,1);
	addEdge(g,4,2);
	removeVertex(g,3);
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(outDegree(g, 4) == 2)
	SEQAN_TASSERT(inDegree(g, 4) == 0)

	// Transpose
	transpose(g); 
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(outDegree(g, 4) == 0)
	SEQAN_TASSERT(inDegree(g, 4) == 2)
	StandardGraph g_copy(g);
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 0)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)
	addVertex(g_copy);
	addEdge(g_copy,3,0);
	g_copy = g;
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 0)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)
	//Copies the graph and transposes just the copy
	transpose(g,g_copy);  // g does not change!
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 2)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 0)
	removeVertex(g,0);



	// Adjacency matrix
	Matrix<unsigned int> mat;
	getAdjacencyMatrix(g, mat);
	unsigned int len = length(mat, 0);
	SEQAN_TASSERT(getValue(mat, 1*len+4) == 1)
	SEQAN_TASSERT(getValue(mat, 2*len+4) == 1)
	SEQAN_TASSERT(getValue(mat, 2*len+2) == 0)

//____________________________________________________________________________
//Graph with edge cargo and edge ids
	typedef Pair<char, int> TPair;
	typedef Directed<TPair> TEdges;
	typedef VertexDescriptor<Graph<TEdges> >::Type TVertexDescriptor2;
	typedef EdgeDescriptor<Graph<TEdges> >::Type TEdgeDescriptor2;

	Graph<TEdges> g2;
	SEQAN_TASSERT(numVertices(g2) == 0)
	SEQAN_TASSERT(numEdges(g2) == 0)
	TVertexDescriptor2 ver0 = addVertex(g2);
	SEQAN_TASSERT(ver0 == 0)
	SEQAN_TASSERT(numVertices(g2) == 1)
	TVertexDescriptor2 ver1 = addVertex(g2);
	SEQAN_TASSERT(ver1 == 1)
	SEQAN_TASSERT(numVertices(g2) == 2)
	TEdgeDescriptor2 ed1 =addEdge(g2,ver0,ver0, TPair('a',3));
	TEdgeDescriptor2 ed2 =addEdge(g2,0,1);
	SEQAN_TASSERT((getCargo(ed1)).i1 == 'a')
	SEQAN_TASSERT((getCargo(ed1)).i2 == 3)
	SEQAN_TASSERT(targetVertex(g2, ed1) == v0)
	SEQAN_TASSERT(targetVertex(g2, ed1) == 0)
	SEQAN_TASSERT(sourceVertex(g2, ed1) == 0)
	SEQAN_TASSERT(targetVertex(g2, ed2) == 1)
	SEQAN_TASSERT(numEdges(g2) == 2)
	assignCargo(ed2, TPair('b',4));
	SEQAN_TASSERT((getCargo(ed1)).i1 == 'a')
	SEQAN_TASSERT((getCargo(ed1)).i2 == 3)
	SEQAN_TASSERT((getCargo(ed2)).i1 == 'b')
	SEQAN_TASSERT((getCargo(ed2)).i2 == 4)
	cargo(ed1) = TPair('c',1);
	SEQAN_TASSERT((getCargo(ed1)).i1 == 'c')
	SEQAN_TASSERT((getCargo(ed1)).i2 == 1)
	addVertex(g2);
	addVertex(g2);
	addVertex(g2);
	TEdgeDescriptor2 ed4 =addEdge(g2,1,4);
	cargo(ed4) = TPair('z',100);
	removeVertex(g2, 2);
	Graph<TEdges> g2_copy(g2);
	SEQAN_TASSERT(numVertices(g2_copy) == 4)
	SEQAN_TASSERT(numEdges(g2_copy) == 3)
	clearEdges(g2_copy);
	SEQAN_TASSERT(numVertices(g2_copy) == 4)
	SEQAN_TASSERT(numEdges(g2_copy) == 0)
	clearVertices(g2_copy);
	SEQAN_TASSERT(numVertices(g2_copy) == 0)
	addVertex(g2_copy);addVertex(g2_copy);
	addEdge(g2_copy,0,1);
	clear(g2_copy);
	SEQAN_TASSERT(numVertices(g2_copy) == 0)
	addVertex(g2_copy);addVertex(g2_copy);
	addEdge(g2_copy,0,1);
	SEQAN_TASSERT(numEdges(g2) == 3)
	SEQAN_TASSERT(outDegree(g2, 0) == 2)
	SEQAN_TASSERT(inDegree(g2, 0) == 1)
	transpose(g2, g2_copy);
	SEQAN_TASSERT(outDegree(g2_copy, 0) == 1)
	SEQAN_TASSERT(inDegree(g2_copy, 0) == 2)
	SEQAN_TASSERT(numEdges(g2_copy) == 3)
	TEdgeDescriptor2 edgCargo = addEdge(g2, 0, 0, TPair('m',3));
	SEQAN_TASSERT((getCargo(edgCargo)).i1 == 'm')
	SEQAN_TASSERT((getCargo(edgCargo)).i2 == 3)

//____________________________________________________________________________
//Graph without edge cargo and without edge ids
	typedef Directed<void, WithoutEdgeId> TEdges3;
	typedef VertexDescriptor<Graph<TEdges3> >::Type TVertexDescriptor3;
	typedef EdgeDescriptor<Graph<TEdges3> >::Type TEdgeDescriptor3;

	Graph<TEdges3> g3;
	addVertex(g3);addVertex(g3);addVertex(g3);
	addVertex(g3);addVertex(g3);
	addEdge(g3,1,4);
	SEQAN_TASSERT(numVertices(g3) == 5)
	SEQAN_TASSERT(numEdges(g3) == 1)
	TEdgeDescriptor3 edge3 = addEdge(g3,0,4);
	SEQAN_TASSERT(_getId(edge3) == 0);
	SEQAN_TASSERT(getCargo(edge3) == (void*) 0);
	addEdge(g3,0,2);
	addEdge(g3,0,0);
	removeEdge(g3,0,4);
	removeEdge(g3,0,2);
	SEQAN_TASSERT(numEdges(g3) == 2)

		
	clear(g);
	TVertexDescriptor edges[] = {0,1, 1,2};
	unsigned int numEdg = 2;
	std::string nameEd[] = {"ar", "ae"};
	addEdges(g, edges, numEdg);
	String<std::string> edMap;
	resizeEdgeMap(g, edMap, nameEd);
	SEQAN_TASSERT(getProperty(edMap, 0) == "ar")
	SEQAN_TASSERT(getProperty(edMap, 1) == "ae")
}

//////////////////////////////////////////////////////////////////////////////

void Test_Undirected() {
//____________________________________________________________________________
// Graph without edge cargo but with edge ids

	typedef Graph<Undirected<void> > StandardGraph;
	typedef VertexDescriptor<StandardGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<StandardGraph>::Type TEdgeDescriptor;
	
	StandardGraph g;
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	SEQAN_TASSERT(empty(g) == true)

	// Add vertex
	TVertexDescriptor v0 = addVertex(g);
	SEQAN_TASSERT(v0 == 0)
	SEQAN_TASSERT(outDegree(g, v0) == 0)	
	SEQAN_TASSERT(inDegree(g, 0) == 0)
	SEQAN_TASSERT(degree(g, 0) == 0)
	SEQAN_TASSERT(numVertices(g) == 1)
	SEQAN_TASSERT(empty(g) == false)
	
	// Add edge
	// TEdgeDescriptor e1 =addEdge(g,v0,v0);  // Self edges are not allowed in undirected graphs
	TVertexDescriptor v1 = addVertex(g);
	TEdgeDescriptor e =addEdge(g,0,1);
	SEQAN_TASSERT(findEdge(g, 0, 1) == e)
	SEQAN_TASSERT(_getVertexString(g)[0] == e)
	SEQAN_TASSERT(getIdUpperBound(_getVertexIdManager(g)) == 2)
	SEQAN_TASSERT(getIdUpperBound(_getEdgeIdManager(g)) == 1)
	SEQAN_TASSERT(v1 == 1)
	SEQAN_TASSERT(numVertices(g) == 2)
	SEQAN_TASSERT(targetVertex(g, e) == 1)
	SEQAN_TASSERT(sourceVertex(g, e) == 0)
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(outDegree(g, v0) == 1)	
	SEQAN_TASSERT(inDegree(g, 1) == 1)
	SEQAN_TASSERT(inDegree(g, 0) == 1)	
	SEQAN_TASSERT(degree(g, 0) == 1)

	// Add more vertices and edges
	addVertex(g);  //2
	TVertexDescriptor v3 = addVertex(g);  //3
	addVertex(g);  //4
	addEdge(g,3,4);
	TEdgeDescriptor my_edge = addEdge(g,3,1);
	addEdge(g,3,0);
	SEQAN_TASSERT(v3 == 3)
	SEQAN_TASSERT(numVertices(g) == 5)
	SEQAN_TASSERT(targetVertex(g, my_edge) == 3)
	SEQAN_TASSERT(sourceVertex(g, my_edge) == 1)
	SEQAN_TASSERT(numEdges(g) == 4)
	SEQAN_TASSERT(outDegree(g, v3) == 3)
	SEQAN_TASSERT(inDegree(g, v3) == 3)
	SEQAN_TASSERT(degree(g, v3) == 3)

	// Graph drawing
	// Raw output
	// std::cout << g << ::std::endl;
	// File output
	fstream strm;
	strm.open(TEST_PATH "my_undirected_graph.dot", ios_base::out | ios_base::trunc);
	write(strm,g,DotDrawing());
	strm.close();
	// File read
	StandardGraph gTmp;
	strm.open(TEST_PATH "my_undirected_graph.dot", ios_base::in);
	read(strm,gTmp,DotDrawing());
	strm.close();

	// Remove edges
	removeEdge(g,my_edge);
	removeEdge(g,0,1);
	SEQAN_TASSERT(numEdges(g) == 2)

	
	// Remove vertices 
	addVertex(g);  //5
	addEdge(g,5,2);
	addEdge(g,2,3);
	addEdge(g,1,3);
	addEdge(g,1,4);
	SEQAN_TASSERT(outDegree(g, 3) == 4)
	SEQAN_TASSERT(outDegree(g, 4) == 2)
	removeVertex(g, v3);
	SEQAN_TASSERT(outDegree(g, 4) == 1)
	SEQAN_TASSERT(outDegree(g, 0) == 0)
	SEQAN_TASSERT(numVertices(g) == 5)
	SEQAN_TASSERT(numEdges(g) == 2)

	// Clear graph
	clearEdges(g);
	SEQAN_TASSERT(numVertices(g) == 5)
	SEQAN_TASSERT(numEdges(g) == 0)
	addEdge(g,2,0);
	addEdge(g,4,1);
	clearVertices(g);
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	addVertex(g);addVertex(g);addVertex(g);
	addVertex(g);addVertex(g);
	addEdge(g,2,0);
	addEdge(g,4,1);
	clear(g);
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	addVertex(g);addVertex(g);addVertex(g);
	addVertex(g);addVertex(g);
	addEdge(g,2,0);
	addEdge(g,4,1);
	addEdge(g,4,2);
	removeVertex(g,3);
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(outDegree(g, 4) == 2)
	SEQAN_TASSERT(inDegree(g, 4) == 2)

	// Transpose
	transpose(g); 
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(outDegree(g, 4) == 2)
	SEQAN_TASSERT(inDegree(g, 4) == 2)
	StandardGraph g_copy(g);
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 2)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)
	addVertex(g_copy);
	addEdge(g_copy,3,0);
	g_copy = g;
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 2)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)
	//Copies the graph and transposes just the copy
	transpose(g,g_copy);  // g does not change!
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 2)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)

	// Adjacency matrix
	Matrix<unsigned int> mat;
	getAdjacencyMatrix(g, mat);
	unsigned int len = getIdUpperBound(g.data_id_managerV);
	SEQAN_TASSERT(getValue(mat,0*len+2) == 1)
	SEQAN_TASSERT(getValue(mat,3*len+2) == 0)
	SEQAN_TASSERT(getValue(mat,0*len+2) == getValue(mat,2*len+0))
	SEQAN_TASSERT(getValue(mat,1*len+4) == getValue(mat,4*len+1))
	SEQAN_TASSERT(getValue(mat,2*len+4) == getValue(mat,4*len+2))

//____________________________________________________________________________
//Graph with edge cargo and edge ids
	typedef Pair<char, int> TPair;
	typedef Undirected<TPair> TEdges;
	typedef VertexDescriptor<Graph<TEdges> >::Type TVertexDescriptor2;
	typedef EdgeDescriptor<Graph<TEdges> >::Type TEdgeDescriptor2;

	Graph<TEdges> g2;
	SEQAN_TASSERT(numVertices(g2) == 0)
	SEQAN_TASSERT(numEdges(g2) == 0)
	TVertexDescriptor2 ver0 = addVertex(g2);
	SEQAN_TASSERT(ver0 == 0)
	SEQAN_TASSERT(numVertices(g2) == 1)
	TVertexDescriptor2 ver1 = addVertex(g2);
	SEQAN_TASSERT(ver1 == 1)
	SEQAN_TASSERT(numVertices(g2) == 2)
	TEdgeDescriptor2 ed1 =addEdge(g2,0,1);
	SEQAN_TASSERT(targetVertex(g2, ed1) == 1)
	SEQAN_TASSERT(sourceVertex(g2, ed1) == 0)
	SEQAN_TASSERT(numEdges(g2) == 1)
	assignCargo(ed1, TPair('a',3));
	SEQAN_TASSERT((getCargo(ed1)).i1 == 'a')
	SEQAN_TASSERT((getCargo(ed1)).i2 == 3)
	cargo(ed1) = TPair('c',1);
	SEQAN_TASSERT((getCargo(ed1)).i1 == 'c')
	SEQAN_TASSERT((getCargo(ed1)).i2 == 1)
	addVertex(g2);
	addVertex(g2);
	addVertex(g2);
	TEdgeDescriptor2 ed4 =addEdge(g2,1,4);
	cargo(ed4) = TPair('z',100);
	removeVertex(g2, 2);
	Graph<TEdges> g2_copy(g2);
	SEQAN_TASSERT(numVertices(g2_copy) == 4)
	SEQAN_TASSERT(numEdges(g2_copy) == 2)
	clearEdges(g2_copy);
	SEQAN_TASSERT(numVertices(g2_copy) == 4)
	SEQAN_TASSERT(numEdges(g2_copy) == 0)
	clearVertices(g2_copy);
	SEQAN_TASSERT(numVertices(g2_copy) == 0)
	addVertex(g2_copy);addVertex(g2_copy);
	addEdge(g2_copy,0,1);
	clear(g2_copy);
	SEQAN_TASSERT(numVertices(g2_copy) == 0)
	addVertex(g2_copy);addVertex(g2_copy);
	addEdge(g2_copy,0,1);
	transpose(g2, g2_copy);
	SEQAN_TASSERT(outDegree(g2_copy, 0) == 1)
	SEQAN_TASSERT(inDegree(g2_copy, 0) == 1)
	SEQAN_TASSERT(numEdges(g2_copy) == 2)
	TEdgeDescriptor2 edgCargo = addEdge(g2, 0, 3, TPair('m',3));
	SEQAN_TASSERT((getCargo(edgCargo)).i1 == 'm')
	SEQAN_TASSERT((getCargo(edgCargo)).i2 == 3)

//____________________________________________________________________________
//Graph without edge cargo and without edge ids
	typedef Undirected<void, WithoutEdgeId> TEdges3;
	typedef VertexDescriptor<Graph<TEdges3> >::Type TVertexDescriptor3;
	typedef EdgeDescriptor<Graph<TEdges3> >::Type TEdgeDescriptor3;

	Graph<TEdges3> g3;
	addVertex(g3);addVertex(g3);addVertex(g3);
	addVertex(g3);addVertex(g3);
	addEdge(g3,1,4);
	SEQAN_TASSERT(numVertices(g3) == 5)
	SEQAN_TASSERT(numEdges(g3) == 1)
	TEdgeDescriptor3 edge3 = addEdge(g3,0,4);
	SEQAN_TASSERT(_getId(edge3) == 0);
	SEQAN_TASSERT(getCargo(edge3) == (void*) 0);
	addEdge(g3,0,2);
	addEdge(g3,0,1);
	removeEdge(g3,0,4);
	removeEdge(g3,0,2);
	SEQAN_TASSERT(numEdges(g3) == 2)
	removeInEdges(g3,1);
	SEQAN_TASSERT(numEdges(g3) == 0)
	

//____________________________________________________________________________
// Undirected graph iterators
	typedef Graph<Undirected<> > TGraphIter;
	typedef VertexDescriptor<TGraphIter>::Type TVertexDescriptorIter;
	typedef EdgeDescriptor<TGraphIter>::Type TEdgeDescriptorIter;
	
	TGraphIter gIter;
	addVertex(gIter);addVertex(gIter);addVertex(gIter);addVertex(gIter);
	addVertex(gIter);addVertex(gIter);addVertex(gIter);addVertex(gIter);
	removeVertex(gIter,0);
	removeVertex(gIter,5);
	addEdge(gIter,2,7);
	addEdge(gIter,2,3);
	addEdge(gIter,2,4);
	addEdge(gIter,4,3);
	addEdge(gIter,3,6);
	addEdge(gIter,4,6);

	typedef Iterator<TGraphIter, OutEdgeIterator>::Type TOutEdgeIterator;
	TOutEdgeIterator itOutEdge(gIter,3);
	// Both ways are fast for undirected graphs
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itOutEdge))==3)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itOutEdge))==6)
	SEQAN_TASSERT(sourceVertex(itOutEdge)==3)
	SEQAN_TASSERT(targetVertex(itOutEdge)==6)
	SEQAN_TASSERT(sourceVertex(gIter, value(itOutEdge))==3)
	SEQAN_TASSERT(targetVertex(gIter, *itOutEdge)==6)
	SEQAN_TASSERT(atEnd(itOutEdge)==false)
	SEQAN_TASSERT(atBegin(itOutEdge)==true)
	goNext(itOutEdge);
	SEQAN_TASSERT(atEnd(itOutEdge)==false)
	SEQAN_TASSERT(atBegin(itOutEdge)==false)
	SEQAN_TASSERT(sourceVertex(itOutEdge)==3)
	SEQAN_TASSERT(targetVertex(itOutEdge)==4)
	++itOutEdge;
	itOutEdge++;
	SEQAN_TASSERT(atEnd(itOutEdge)==true)
	SEQAN_TASSERT(atBegin(itOutEdge)==false)
	goPrevious(itOutEdge);
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itOutEdge))==2)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itOutEdge))==3)
	--itOutEdge;
	itOutEdge--; 
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itOutEdge))==3)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itOutEdge))==6)
	itOutEdge--;
	itOutEdge--;
	SEQAN_TASSERT(atBegin(itOutEdge)==true)
	TOutEdgeIterator itEdge2(itOutEdge);
	TOutEdgeIterator itEdge3;
	itEdge3 = itOutEdge;
	SEQAN_TASSERT(itOutEdge == itEdge2)
	SEQAN_TASSERT(itEdge2 == itEdge3)
	goEnd(itOutEdge);
	SEQAN_TASSERT(itEdge2 != itOutEdge)
	goEnd(itEdge2);
	SEQAN_TASSERT(itEdge2 == itOutEdge)
	goBegin(itEdge2);
	SEQAN_TASSERT(atBegin(itEdge2)==true)
	SEQAN_TASSERT(&gIter == &hostGraph(itOutEdge))

	
	typedef Iterator<TGraphIter, EdgeIterator>::Type TEdgeIterator;
	TEdgeIterator itEdge(gIter);
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==2)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==4)
	SEQAN_TASSERT(atBegin(itEdge)==true)
	SEQAN_TASSERT(atEnd(itEdge)==false)
	goNext(itEdge);
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==2)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==3)
	SEQAN_TASSERT(atBegin(itEdge)==false)
	SEQAN_TASSERT(atEnd(itEdge)==false)
	goNext(itEdge);
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==2)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==7)
	++itEdge;
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==3)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==6)
	itEdge++;
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==3)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==4)
	goNext(itEdge);
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==4)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==6)
	goNext(itEdge);
	SEQAN_TASSERT(atBegin(itEdge)==false)
	SEQAN_TASSERT(atEnd(itEdge)==true)
	goPrevious(itEdge);	
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==4)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==6)
	--itEdge;	
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==3)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==4)
	itEdge--;	
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==3)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==6)
	goPrevious(itEdge);	
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==2)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==7)
	goPrevious(itEdge);	
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==2)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==3)
	goPrevious(itEdge);	
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==2)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==4)
	SEQAN_TASSERT(atBegin(itEdge)==true)
	SEQAN_TASSERT(atEnd(itEdge)==false)
}


//////////////////////////////////////////////////////////////////////////////

void Test_Automaton() {
//____________________________________________________________________________
// Standard automaton: No edge cargo

	typedef Graph<Automaton<Dna> > StandardAutomaton;
	typedef VertexDescriptor<StandardAutomaton>::Type TVertexDescriptor;
	typedef EdgeDescriptor<StandardAutomaton>::Type TEdgeDescriptor;
	
	StandardAutomaton g;
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	SEQAN_TASSERT(empty(g) == true)

	// Add vertex
	createRoot(g);
	TVertexDescriptor v0 = getRoot(g);
	SEQAN_TASSERT(v0 == 0)
	SEQAN_TASSERT(outDegree(g, v0) == 0)	
	SEQAN_TASSERT(inDegree(g, 0) == 0)
	SEQAN_TASSERT(degree(g, 0) == 0)
	SEQAN_TASSERT(numVertices(g) == 1)
	SEQAN_TASSERT(empty(g) == false)

	// Add edge
	TEdgeDescriptor e1 =addEdge(g,v0,v0,'a');
	SEQAN_TASSERT(findEdge(g, 0, 'a') == e1)
	SEQAN_TASSERT(&_getVertexString(g)[0].data_edge[0] == e1)
	SEQAN_TASSERT(getIdUpperBound(_getVertexIdManager(g)) == 1)
	SEQAN_TASSERT(getIdUpperBound(_getEdgeIdManager(g)) == 1)
	SEQAN_TASSERT(_getId(e1) == 0)
	SEQAN_TASSERT(_getId(e1) == 0)
	SEQAN_TASSERT(targetVertex(g, e1) == 0)
	SEQAN_TASSERT(sourceVertex(g, e1) == 0) 
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(outDegree(g, v0) == 1)	
	SEQAN_TASSERT(inDegree(g, v0) == 1)
	SEQAN_TASSERT(degree(g, v0) == 2)
	
	// Add further edges and vertices
	TVertexDescriptor v1 = addVertex(g);
	TEdgeDescriptor e2 =addEdge(g,0,1,'g');
	SEQAN_TASSERT(_getId(e2) == 1)
	SEQAN_TASSERT(v1 == 1)
	SEQAN_TASSERT(numVertices(g) == 2)
	SEQAN_TASSERT(targetVertex(g, e2) == 1)
	SEQAN_TASSERT(sourceVertex(g, e2) == 0)
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(outDegree(g, v0) == 2)	
	SEQAN_TASSERT(inDegree(g, 1) == 1)
	SEQAN_TASSERT(inDegree(g, 0) == 1)	
	SEQAN_TASSERT(degree(g, 0) == 3)

	// Add more vertices and edges
	addVertex(g);  //2
	TVertexDescriptor v3 = addVertex(g);  //3
	addVertex(g);  //4
	addEdge(g,3,4,'g');
	TEdgeDescriptor my_edge = addEdge(g,3,1,'c');
	SEQAN_TASSERT(_getId(my_edge) == 3)
	addEdge(g,3,0,'t');
	SEQAN_TASSERT(v3 == 3)
	SEQAN_TASSERT(numVertices(g) == 5)
	SEQAN_TASSERT(targetVertex(g, e2) == 1)
	SEQAN_TASSERT(sourceVertex(g, e2) == 0)
	SEQAN_TASSERT(targetVertex(g, my_edge) == 1)
	SEQAN_TASSERT(sourceVertex(g, my_edge) == 3)
	SEQAN_TASSERT(numEdges(g) == 5)
	SEQAN_TASSERT(outDegree(g, v3) == 3)	

	// Output
	// Raw output
	// std::cout << g << ::std::endl;
	// File output
	fstream strm;
	strm.open(TEST_PATH "my_automaton.dot", ios_base::out | ios_base::trunc);
	write(strm,g,DotDrawing());
	strm.close();
	// File read
	StandardAutomaton gTmp;
	strm.open(TEST_PATH "my_automaton.dot", ios_base::in);
	read(strm,gTmp,DotDrawing());
	strm.close();

	// Remove edges
	removeEdge(g,3,1,'c');
	removeEdge(g,0,1,'g');
	SEQAN_TASSERT(numEdges(g) == 3)

	// Remove vertices 
	TEdgeDescriptor e3 = addEdge(g,3,3,'a');
	addEdge(g,1,3,'a');
	addEdge(g,0,3,'c');
	addEdge(g,0,4,'t');
	SEQAN_TASSERT(outDegree(g, 0) == 3)
	SEQAN_TASSERT(outDegree(g, 1) == 1)
	SEQAN_TASSERT(targetVertex(g, e3) == 3)
	SEQAN_TASSERT(sourceVertex(g, e3) == 3)
	SEQAN_TASSERT(numEdges(g) == 7)
	removeVertex(g, v3);
	SEQAN_TASSERT(outDegree(g, 0) == 2)
	SEQAN_TASSERT(outDegree(g, 1) == 0)
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 2)

	// Clear graph
	clearEdges(g);
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 0)
	addEdge(g,2,0,'a');
	addEdge(g,4,1,'c');
	clearVertices(g);
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	addVertex(g);addVertex(g);addVertex(g);
	addVertex(g);addVertex(g);
	addEdge(g,2,0,'t');
	addEdge(g,4,1,'g');
	clear(g);
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	addVertex(g);addVertex(g);addVertex(g);
	addVertex(g);addVertex(g);
	addEdge(g,2,0,'c');
	addEdge(g,4,1,'g');
	addEdge(g,4,2,'t');
	removeVertex(g,3);
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(outDegree(g, 4) == 2)
	SEQAN_TASSERT(inDegree(g, 4) == 0)

	//Transposes the graph in-place
	transpose(g); 
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(outDegree(g, 4) == 0)
	SEQAN_TASSERT(inDegree(g, 4) == 2)
	StandardAutomaton g_copy(g);
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 0)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)
	addVertex(g_copy);
	addEdge(g_copy,3,0,'a');
	g_copy = g;
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 0)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)
	//Copies the graph and transposes just the copy
	transpose(g,g_copy);  // g does not change!
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 2)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 0)
	removeVertex(g,0);

	// Adjacency matrix
	Matrix<unsigned int> mat;
	getAdjacencyMatrix(g, mat);
	unsigned int len = length(mat,0);
	SEQAN_TASSERT(getValue(mat,1*len+4) == 1)
	SEQAN_TASSERT(getValue(mat,2*len+4) == 1)
	SEQAN_TASSERT(getValue(mat,0*len+2) == 0)

	// Test iterators
	typedef Iterator<StandardAutomaton, VertexIterator>::Type TVertexIterator;
	TVertexIterator itVert(g);
	SEQAN_TASSERT(getValue(itVert) == 1)
	++itVert;
	SEQAN_TASSERT(getValue(itVert) == 2)
	itVert++;
	SEQAN_TASSERT(getValue(itVert) == 4)
	goNext(itVert);
	SEQAN_TASSERT(atEnd(itVert) == true)

	addEdge(g,1,2,'T');
	typedef Iterator<StandardAutomaton, OutEdgeIterator>::Type TOutEdgeIterator;
	TOutEdgeIterator itEdge(g,1);
	// Slow
	SEQAN_TASSERT(sourceVertex(g, getValue(itEdge))==1)
	SEQAN_TASSERT(targetVertex(g, getValue(itEdge))==4)
	// Fast
	SEQAN_TASSERT(sourceVertex(itEdge)==1)
	SEQAN_TASSERT(targetVertex(itEdge)==4)
	
	SEQAN_TASSERT(sourceVertex(g, value(itEdge))==1)
	SEQAN_TASSERT(targetVertex(g, *itEdge)==4)
	SEQAN_TASSERT(atEnd(itEdge)==false)
	SEQAN_TASSERT(atBegin(itEdge)==true)
	goNext(itEdge);
	SEQAN_TASSERT(sourceVertex(itEdge)==1)
	SEQAN_TASSERT(targetVertex(itEdge)==2)
	++itEdge;
	itEdge++;
	SEQAN_TASSERT(atEnd(itEdge)==true)
	SEQAN_TASSERT(atBegin(itEdge)==false)
	goPrevious(itEdge);
	SEQAN_TASSERT(sourceVertex(g, getValue(itEdge))==1)
	SEQAN_TASSERT(targetVertex(g, getValue(itEdge))==2)
	--itEdge;
	itEdge++; 
	SEQAN_TASSERT(sourceVertex(g, getValue(itEdge))==1)
	SEQAN_TASSERT(targetVertex(g, getValue(itEdge))==2)
	itEdge--;
	itEdge--;
	SEQAN_TASSERT(atBegin(itEdge)==true)
	TOutEdgeIterator itEdge2(itEdge);
	TOutEdgeIterator itEdge3;
	itEdge3 = itEdge;
	SEQAN_TASSERT(itEdge == itEdge2)
	SEQAN_TASSERT(itEdge2 == itEdge3)
	goEnd(itEdge);
	SEQAN_TASSERT(itEdge2 != itEdge)
	goEnd(itEdge2);
	SEQAN_TASSERT(itEdge2 == itEdge)
	goBegin(itEdge2);
	SEQAN_TASSERT(atBegin(itEdge2)==true)
	SEQAN_TASSERT(&g == &hostGraph(itEdge))

	typedef Iterator<StandardAutomaton, EdgeIterator>::Type TEdgeIterator;
	TEdgeIterator itEd(g);
	SEQAN_TASSERT(sourceVertex(g, getValue(itEd))==1)
	SEQAN_TASSERT(targetVertex(g, getValue(itEd))==4)
	SEQAN_TASSERT(sourceVertex(g, value(itEd))==1)
	SEQAN_TASSERT(targetVertex(g, *itEd)==4)
	SEQAN_TASSERT(atEnd(itEd)==false)
	SEQAN_TASSERT(atBegin(itEd)==true)
	goNext(itEd);
	SEQAN_TASSERT(sourceVertex(g, getValue(itEd))==1)
	SEQAN_TASSERT(targetVertex(g, getValue(itEd))==2)
	SEQAN_TASSERT(atEnd(itEd)==false)
	SEQAN_TASSERT(atBegin(itEd)==false)
	++itEd;
	SEQAN_TASSERT(atEnd(itEd)==false)
	SEQAN_TASSERT(atBegin(itEd)==false)
	// Slow
	SEQAN_TASSERT(sourceVertex(g, getValue(itEd))==2)
	SEQAN_TASSERT(targetVertex(g, getValue(itEd))==4)
	// Fast
	SEQAN_TASSERT(sourceVertex(itEd)==2)
	SEQAN_TASSERT(targetVertex(itEd)==4)
	itEd++;
	itEd++;
	SEQAN_TASSERT(atEnd(itEd)==true)
	SEQAN_TASSERT(atBegin(itEd)==false)
	goPrevious(itEd);
	SEQAN_TASSERT(sourceVertex(g, getValue(itEd))==2)
	SEQAN_TASSERT(targetVertex(g, getValue(itEd))==4)
	--itEd;
	SEQAN_TASSERT(sourceVertex(g, getValue(itEd))==1)
	SEQAN_TASSERT(targetVertex(g, getValue(itEd))==2)
	TEdgeIterator itEd2(g);
	TEdgeIterator itEd3;
	goBegin(itEd);
	itEd3 = itEd;
	SEQAN_TASSERT(itEd == itEd2)
	SEQAN_TASSERT(itEd2 == itEd3)
	goEnd(itEd);
	SEQAN_TASSERT(itEd2 != itEd)
	goEnd(itEd2);
	SEQAN_TASSERT(itEd2 == itEd)
	goBegin(itEd2);
	SEQAN_TASSERT(itEd2 != itEd)
	SEQAN_TASSERT(&hostGraph(itEd) == &g)

	typedef Iterator<StandardAutomaton, AdjacencyIterator>::Type TAdjacencyIterator;
	TAdjacencyIterator itAd(g,1);
	SEQAN_TASSERT(getValue(itAd) == 4)
	SEQAN_TASSERT(&hostGraph(itAd) == &g)
	SEQAN_TASSERT(value(itAd) == 4)
	SEQAN_TASSERT(*itAd == 4)
	SEQAN_TASSERT(atEnd(itAd)==false)
	SEQAN_TASSERT(atBegin(itAd)==true)
	goNext(itAd);
	SEQAN_TASSERT(getValue(itAd)==2)
	SEQAN_TASSERT(atEnd(itAd)==false)
	SEQAN_TASSERT(atBegin(itAd)==false)
	++itAd;
	SEQAN_TASSERT(atEnd(itAd)==true)
	SEQAN_TASSERT(atBegin(itAd)==false)
	goBegin(itAd);
	itAd++;
	itAd++;
	itAd++;
	SEQAN_TASSERT(atEnd(itAd)==true)
	SEQAN_TASSERT(atBegin(itAd)==false)
	goPrevious(itAd);
	SEQAN_TASSERT(getValue(itAd)==2)
	--itAd;
	SEQAN_TASSERT(getValue(itAd)==4)
	SEQAN_TASSERT(atEnd(itAd)==false)
	SEQAN_TASSERT(atBegin(itAd)==true)
	goEnd(itAd);
	itAd--;
	SEQAN_TASSERT(getValue(itAd)==2)
	goBegin(itAd);
	TAdjacencyIterator itAd2(itAd);
	TAdjacencyIterator itAd3;
	itAd3 = itAd;
	SEQAN_TASSERT(itAd == itAd2)
	SEQAN_TASSERT(itAd2 == itAd3)
	goEnd(itAd);
	SEQAN_TASSERT(itAd2 != itAd)
	goEnd(itAd2);
	SEQAN_TASSERT(itAd2 == itAd)
	goBegin(itAd2);
	SEQAN_TASSERT(itAd2 != itAd)



//____________________________________________________________________________
// Automaton - Different alphabet
	typedef VertexDescriptor<Graph<Automaton<char> > >::Type VertexDescriptorType;
	typedef EdgeDescriptor<Graph<Automaton<char> > >::Type EdgeDescriptorType;
	Graph<Automaton<char> > automaton;
	VertexDescriptorType rootVertex = addVertex(automaton); // A = 0
	addVertex(automaton); // B = 1
	addVertex(automaton); // C = 2
	addVertex(automaton); // D = 3
	addVertex(automaton); // E = 4
	addVertex(automaton); // F = 5
	addEdge(automaton,0,1,'2');
	addEdge(automaton,1,0,'1');
	addEdge(automaton,4,0,'6');
	addEdge(automaton,0,3,'7');
	addEdge(automaton,1,1,'3');
	addEdge(automaton,1,2,'4');
	addEdge(automaton,5,1,'8');
	addEdge(automaton,2,5,'5');
	addEdge(automaton,3,4,'2');
	addEdge(automaton,5,3,'7');

	VertexDescriptorType succ;
	succ = getSuccessor(automaton,rootVertex,'7');
	SEQAN_TASSERT(succ == 3)
	// Throws an error in debug mode because edge does not exist
	//succ = getSuccessor(automaton,rootVertex,'6');
	succ = getSuccessor(automaton,succ,'2');
	SEQAN_TASSERT(succ == 4)
	succ = getSuccessor(automaton,succ,'6');
	SEQAN_TASSERT(succ == 0)
	// If no map is specified it is assumed that an edge cargo exists!!!
	succ = getSuccessor(automaton,succ,'2');
	SEQAN_TASSERT(succ == 1)

	// Now using shortcuts
	succ = parseString(automaton,rootVertex,"7262");
	SEQAN_TASSERT(succ == 1)
	std::string str = "7262";
	succ = parseString(automaton,rootVertex, str.begin(), str.end());
	SEQAN_TASSERT(succ == 1)
	String<char> str2("7262");
	succ = parseString(automaton,rootVertex, begin(str2), end(str2));
	SEQAN_TASSERT(succ == 1)
	String<char> input("7262");
	succ = parseString(automaton,rootVertex, input);
	SEQAN_TASSERT(succ == 1)

	// Additional cargo
	typedef Graph<Automaton<Dna, short> > TGraph9;
	typedef VertexDescriptor<TGraph9>::Type TVertexDescriptor9;
	typedef EdgeDescriptor<TGraph9>::Type TEdgeDescriptor9;
	typedef Size<TGraph9>::Type TSize9;

	TGraph9 g9;
	addVertex(g9);
	addVertex(g9);
	Dna aDna('a');
	Dna gDna('g');
	TEdgeDescriptor9 edg1 = addEdge(g9,0,1,aDna,12);
	TEdgeDescriptor9 edg2 = addEdge(g9,1,0,gDna,21);
	TGraph9 g10;
	transpose(g9, g10);
	TEdgeDescriptor9 edg1_10 = findEdge(g10, 1, aDna);
	TEdgeDescriptor9 edg1_11 = findEdge(g10, 0, gDna);
	SEQAN_TASSERT(getCargo(edg1)==12)
	SEQAN_TASSERT(getCargo(edg2)==21)
	SEQAN_TASSERT(getCargo(edg1_10)==12)
	SEQAN_TASSERT(getCargo(edg1_11)==21)
}

//////////////////////////////////////////////////////////////////////////////

void Test_WordGraph() {
//____________________________________________________________________________
// Standard automaton: No edge cargo

	typedef Graph<Automaton<Dna, String<Dna>, WordGraph<> > > TWordGraph;
	typedef VertexDescriptor<TWordGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TWordGraph>::Type TEdgeDescriptor;
	

	TWordGraph g;
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	SEQAN_TASSERT(empty(g) == true)

	// Add vertex
	TVertexDescriptor v0 = addVertex(g);
	SEQAN_TASSERT(v0 == 0)
	SEQAN_TASSERT(outDegree(g, v0) == 0)	
	SEQAN_TASSERT(inDegree(g, 0) == 0)
	SEQAN_TASSERT(degree(g, 0) == 0)
	SEQAN_TASSERT(numVertices(g) == 1)
	SEQAN_TASSERT(empty(g) == false)
	addVertex(g);
	addVertex(g);
	TVertexDescriptor v3 = addVertex(g);
	SEQAN_TASSERT(isRoot(g, 0) == true)
	SEQAN_TASSERT(getRoot(g) == 0)
	assignRoot(g,3);
	SEQAN_TASSERT(getRoot(g) == 3)
	SEQAN_TASSERT(isRoot(g, 0) == false)
	SEQAN_TASSERT(isRoot(g, 3) == true)
	root(g) = 2;
	SEQAN_TASSERT(getRoot(g) == 2)
	SEQAN_TASSERT(isRoot(g, 3) == false)
	SEQAN_TASSERT(isRoot(g, 2) == true)

	// Add edge
	TEdgeDescriptor e1 =addEdge(g,v0,v3,"ag");
	SEQAN_TASSERT(findEdge(g,v0,'a') == e1)
	SEQAN_TASSERT(_getId(e1) == 0)
	// First letter -> edge label, all other letters into the cargo
	SEQAN_TASSERT(getCargo(e1) == "g")
	SEQAN_TASSERT(targetVertex(g, e1) == 3)
	SEQAN_TASSERT(sourceVertex(g, e1) == 0) 
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(outDegree(g, v0) == 1)	
	SEQAN_TASSERT(inDegree(g, v0) == 0)
	SEQAN_TASSERT(degree(g, v0) == 1)

	// Add further edges and vertices
	addVertex(g);
	TVertexDescriptor v5 = addVertex(g);
	TEdgeDescriptor e2 =addEdge(g,0,5,"g");
	SEQAN_TASSERT(_getId(e2) == 1) 
	SEQAN_TASSERT(v5 == 5)
	SEQAN_TASSERT(numVertices(g) == 6)
	SEQAN_TASSERT(targetVertex(g, e2) == 5)
	SEQAN_TASSERT(sourceVertex(g, e2) == 0)
	SEQAN_TASSERT(numEdges(g) == 2)
	removeEdge(g,0,5,String<Dna>("g"));
	SEQAN_TASSERT(numEdges(g) == 1)
	e2 =addEdge(g,0,5,"g");
	SEQAN_TASSERT(outDegree(g, v0) == 2)	
	SEQAN_TASSERT(inDegree(g, 5) == 1)
	SEQAN_TASSERT(degree(g, 0) == 2)
	SEQAN_TASSERT(getSuccessor(g, 0, "g") == 5)
	SEQAN_TASSERT(getSuccessor(g, 0, String<Dna>("ag")) == 3)  // The whole edge label or just the first letter
	SEQAN_TASSERT(getSuccessor(g, 0, "a") == getNil<TVertexDescriptor>())
	addVertex(g);
	addVertex(g);
	addEdge(g,3,1,"aggg");
	addEdge(g,3,4,"gg");
	addEdge(g,5,2,"aggg");
	addEdge(g,5,7,"g");
	addEdge(g,7,6,"g");
	SEQAN_TASSERT(parseString(g, 0, "agaggg") == 1)
	SEQAN_TASSERT(parseString(g, 0, "aga") == 3)  // Does not reach 1
	SEQAN_TASSERT(parseString(g, 0, 'g') == 5)
	SEQAN_TASSERT(parseString(g, 0, "ggg") == 6)
	SEQAN_TASSERT(parseString(g, 0, "gaggg") == 2)
	SEQAN_TASSERT(parseString(g, 0, "gagggg") == 2)
	assignRoot(g,0);

	// Output
	// Raw output
	// std::cout << g << ::std::endl;
	// File output
	fstream strm;
	strm.open(TEST_PATH "my_wordgraph.dot", ios_base::out | ios_base::trunc);
	write(strm,g,DotDrawing());
	strm.close();
	// File read
	TWordGraph gTmp;
	strm.open(TEST_PATH "my_wordgraph.dot", ios_base::in);
	read(strm,gTmp,DotDrawing());
	strm.close();

	assignRoot(g,2);
	TWordGraph g_tmp(g);
	SEQAN_TASSERT(numVertices(g_tmp) == 8)
	SEQAN_TASSERT(parseString(g_tmp, 0, "agaggg") == 1)
	SEQAN_TASSERT(inDegree(g_tmp, 5) == 1)
	SEQAN_TASSERT(degree(g_tmp, 0) == 2)
	SEQAN_TASSERT(isRoot(g_tmp, 2) == true)
	TWordGraph g_assign;
	g_assign = g;
	SEQAN_TASSERT(numVertices(g_assign) == 8)
	SEQAN_TASSERT(parseString(g_assign, 0, "agaggg") == 1)
	SEQAN_TASSERT(inDegree(g_assign, 5) == 1)
	SEQAN_TASSERT(degree(g_assign, 0) == 2)

	// Transpose
	transpose(g, g_tmp);
	SEQAN_TASSERT(numVertices(g_tmp) == 8)
	SEQAN_TASSERT(parseString(g_tmp, 2, "aggg") == 5)
	SEQAN_TASSERT(inDegree(g_tmp, 5) == 2)
	SEQAN_TASSERT(outDegree(g_tmp, 0) == 0)
	SEQAN_TASSERT(isRoot(g_tmp, 2) == true)
}


//////////////////////////////////////////////////////////////////////////////

void Test_Tree() {
//____________________________________________________________________________
// Tree without edge cargo

	typedef Graph<Tree<void> > TTree;
	typedef VertexDescriptor<TTree>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TTree>::Type TEdgeDescriptor;
	
	TTree g;
	SEQAN_TASSERT(empty(g) == true)
	createRoot(g);
	TVertexDescriptor rootV = getRoot(g);
	SEQAN_TASSERT(rootV == 0)
	SEQAN_TASSERT(isRoot(g, rootV) == true)
	SEQAN_TASSERT(root(g) == rootV)
	SEQAN_TASSERT(empty(g) == false)
	TVertexDescriptor childC1 = addChild(g,rootV);
	TEdgeDescriptor childC1e = findEdge(g, rootV, childC1);
	SEQAN_TASSERT(_getVertexString(g)[0] == childC1e)
	SEQAN_TASSERT(getIdUpperBound(_getVertexIdManager(g)) == 2)
	SEQAN_TASSERT(getIdUpperBound(_getEdgeIdManager(g)) == 2)
	SEQAN_TASSERT(targetVertex(g, childC1e) == childC1) // Target in a tree = child
	SEQAN_TASSERT(sourceVertex(g, childC1e) == rootV)  // Source in a tree = parent
	SEQAN_TASSERT(childVertex(g, childC1e) == childC1)  // Shortcuts
	SEQAN_TASSERT(parentVertex(g, childC1e) == rootV)
	SEQAN_TASSERT(empty(g) == false)
	SEQAN_TASSERT(outDegree(g, rootV) == 1)
	TVertexDescriptor childC2 = addChild(g,rootV);
	TVertexDescriptor childC3 = addChild(g,rootV);
	SEQAN_TASSERT(outDegree(g, rootV) == 3)
	SEQAN_TASSERT(childC1 == 1)
	SEQAN_TASSERT(childC2 == 2)
	SEQAN_TASSERT(childC3 == 3)
	TVertexDescriptor childC2C1 = addChild(g,childC2);
	TVertexDescriptor childC2C1C1 = addChild(g,childC2C1);
	TVertexDescriptor childC2C1C1C1 = addChild(g,childC2C1C1);
	TVertexDescriptor childC2C1C1C2 = addChild(g,childC2C1C1);
	TVertexDescriptor childC4 = addChild(g,rootV);
	SEQAN_TASSERT(inDegree(g, childC2C1) == 1) 
	SEQAN_TASSERT(outDegree(g, childC2C1) == 1)
	SEQAN_TASSERT(degree(g, childC2C1) == 2)
	SEQAN_TASSERT(numEdges(g) == 8)
	SEQAN_TASSERT(numVertices(g) == 9)
	TEdgeDescriptor childC2C1C1e = findEdge(g, childC2C1C1, childC2C1);
	
	// Raw output
	// std::cout << g << std::endl;
	// File output
	fstream strm;
	strm.open(TEST_PATH "my_tree.dot", ios_base::out | ios_base::trunc);
	write(strm,g,DotDrawing());
	strm.close();
	// File read
	TTree gTmp;
	strm.open(TEST_PATH "my_tree.dot", ios_base::in);
	read(strm,gTmp,DotDrawing());
	strm.close();

	SEQAN_TASSERT(childVertex(g, childC2C1C1e) == childC2C1C1)  
	SEQAN_TASSERT(parentVertex(g, childC2C1C1e) == childC2C1)
	removeChild(g, rootV, childC2);
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(empty(g) == false)
	SEQAN_TASSERT(inDegree(g, rootV) == 0) 
	SEQAN_TASSERT(outDegree(g, rootV) == 3)
	SEQAN_TASSERT(degree(g, rootV) == 3)
	childC2 = addChild(g,rootV);
	childC2C1 = addChild(g,childC2);
	childC2C1C1 = addChild(g,childC2C1);
	childC2C1C1C1 = addChild(g,childC2C1C1);
	childC2C1C1C2 = addChild(g,childC2C1C1);
	removeAllChildren(g, rootV);
	SEQAN_TASSERT(empty(g) == false)
	SEQAN_TASSERT(numEdges(g) == 0)
	SEQAN_TASSERT(numVertices(g) == 1) // Just the root
	SEQAN_TASSERT(inDegree(g, rootV) == 0) 
	SEQAN_TASSERT(outDegree(g, rootV) == 0)
	SEQAN_TASSERT(degree(g, rootV) == 0)
	addChild(g,rootV);addChild(g,rootV);
	SEQAN_TASSERT(empty(g) == false)
	SEQAN_TASSERT(numEdges(g) == 2)
	clearEdges(g);
	SEQAN_TASSERT(empty(g) == false)
	SEQAN_TASSERT(numEdges(g) == 0)
	addChild(g,rootV);addChild(g,rootV);
	clearVertices(g);
	SEQAN_TASSERT(empty(g) == true)
	SEQAN_TASSERT(numEdges(g) == 0)
	createRoot(g);
	childC1 = addChild(g,rootV);
	SEQAN_TASSERT(empty(g) == false)
	SEQAN_TASSERT(numEdges(g) == 1)
	childC3 = addChild(g,rootV);
	childC2 = addChild(g,rootV);
	childC2C1 = addChild(g,childC2);
	childC2C1C1 = addChild(g,childC2C1);
	childC2C1C1C1 = addChild(g,childC2C1C1);
	childC2C1C1C2 = addChild(g,childC2C1C1);
	childC4 = addChild(g,rootV);
	Matrix<unsigned int> mat; 	// Adjacency matrix
	getAdjacencyMatrix(g, mat);
	unsigned int len = length(mat, 0);
	SEQAN_TASSERT(getValue(mat, 0*len+8) == 1)
	SEQAN_TASSERT(getValue(mat, 8*len+0) == 0)
	SEQAN_TASSERT(getValue(mat, 3*len+0) == 0)
	SEQAN_TASSERT(getValue(mat, 0*len+3) == 1)
	SEQAN_TASSERT(getValue(mat, 0*len+4) == 0)
	SEQAN_TASSERT(numEdges(g) == 8)
	SEQAN_TASSERT(numVertices(g) == 9)
	transpose(g); 
	SEQAN_TASSERT(numEdges(g) == 8)
	SEQAN_TASSERT(numVertices(g) == 9)
	TTree g_copy(g);
	SEQAN_TASSERT(numEdges(g) == 8)
	SEQAN_TASSERT(numVertices(g) == 9)
	clear(g_copy);
	g_copy = g;
	SEQAN_TASSERT(numEdges(g) == 8)
	SEQAN_TASSERT(numVertices(g) == 9)
	transpose(g,g_copy);  
	SEQAN_TASSERT(numEdges(g) == 8)
	SEQAN_TASSERT(numVertices(g) == 9)
	removeOutEdges(g,childC2C1C1);
	SEQAN_TASSERT(numEdges(g) == 6)
	SEQAN_TASSERT(numVertices(g) == 7)
	removeVertex(g,childC2C1);
	SEQAN_TASSERT(numEdges(g) == 4)
	SEQAN_TASSERT(numVertices(g) == 5)
	removeInEdges(g,childC2);
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(numVertices(g) == 4)
	removeOutEdges(g,rootV);
	SEQAN_TASSERT(empty(g) == false) // Root is still present
	addVertex(g);
	TEdgeDescriptor my_edge = addEdge(g,0,1);
	removeEdge(g,my_edge);

//____________________________________________________________________________
// Tree with cargo

	typedef Pair<char, int> TPair;
	typedef Tree<TPair> TEdges;
	typedef Graph<TEdges> TCargoGraph;
	typedef VertexDescriptor<TCargoGraph>::Type TVertexDescriptor2;
	typedef EdgeDescriptor<TCargoGraph>::Type TEdgeDescriptor2;

	TCargoGraph g2;
	createRoot(g2);
	SEQAN_TASSERT(numVertices(g2) == 1)
	SEQAN_TASSERT(numEdges(g2) == 0)
	TVertexDescriptor2 ver1 = addChild(g2, getRoot(g2), TPair('a',3));
	SEQAN_TASSERT(ver1 == 1)
	SEQAN_TASSERT(numVertices(g2) == 2)
	TVertexDescriptor2 ver2 = addChild(g2, getRoot(g2));
	SEQAN_TASSERT(ver2 == 2)
	SEQAN_TASSERT(numVertices(g2) == 3)
	TEdgeDescriptor2 ed1 =findEdge(g2,getRoot(g2),ver1);
	TEdgeDescriptor2 ed2 =findEdge(g2,getRoot(g2),ver2);
	SEQAN_TASSERT((getCargo(ed1)).i1 == 'a')
	SEQAN_TASSERT((getCargo(ed1)).i2 == 3)
	SEQAN_TASSERT(targetVertex(g2, ed1) == ver1)
	SEQAN_TASSERT(sourceVertex(g2, ed1) == getRoot(g2))
	SEQAN_TASSERT(numEdges(g2) == 2)
	assignCargo(ed2, TPair('b',4));
	SEQAN_TASSERT((getCargo(ed1)).i1 == 'a')
	SEQAN_TASSERT((getCargo(ed1)).i2 == 3)
	SEQAN_TASSERT((getCargo(ed2)).i1 == 'b')
	SEQAN_TASSERT((getCargo(ed2)).i2 == 4)
	cargo(ed1) = TPair('c',1);
	SEQAN_TASSERT((getCargo(ed1)).i1 == 'c')
	SEQAN_TASSERT((getCargo(ed1)).i2 == 1)
}

//////////////////////////////////////////////////////////////////////////////

void Test_Alignment() {
//____________________________________________________________________________
// Alignment without edge weights
	typedef String<Dna> TString;
	typedef StringSet<TString, IdHolder<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef	Id<TStringSet>::Type TId;

	
	TStringSet str;
	
	TString str0("acaagtaacataaaaaaaaaaaaaaaacccccccccttttttttaaaaa");
	TId id0 = assignValueById(str, str0);

	TString str1("cccaaagggtttttccccccccccccttttttttttaaaaaaagggggggg");
	TId id1 = assignValueById(str, str1);

	TString str2("cacatgtaatcatgggggggggccccccttttaaaaaaaaaaatttt");
	TId id2 = assignValueById(str, str2);


	TGraph g(str);
	SEQAN_TASSERT(getStringSet(g)[0] == str0)
	SEQAN_TASSERT(getStringSet(g)[1] == str1)
	SEQAN_TASSERT(stringSet(g)[2] == str2)
	assignStringSet(g, str);
	SEQAN_TASSERT(getStringSet(g)[0] == str0)
	SEQAN_TASSERT(getStringSet(g)[1] == str1)
	SEQAN_TASSERT(stringSet(g)[2] == str2)
	SEQAN_TASSERT(numEdges(g) == 0)
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(empty(g) == true)

	TVertexDescriptor v0 = addVertex(g, id1,0,2);
	SEQAN_TASSERT(v0 == 0)
	SEQAN_TASSERT(outDegree(g, v0) == 0)	
	SEQAN_TASSERT(inDegree(g, 0) == 0)
	SEQAN_TASSERT(degree(g, 0) == 0)
	SEQAN_TASSERT(numVertices(g) == 1)
	SEQAN_TASSERT(empty(g) == false)
	SEQAN_TASSERT(sequenceId(g, v0) == id1)
	SEQAN_TASSERT(label(g, v0) == "cc")
	SEQAN_TASSERT(segmentBegin(g, v0) == 0)
	SEQAN_TASSERT(segmentLength(g, v0) == 2)

	TVertexDescriptor v1 = addVertex(g, id2, 0, 5);
	TEdgeDescriptor e =addEdge(g,0,1);
	SEQAN_TASSERT(_getVertexString(g)[0] == e)
	SEQAN_TASSERT(getIdUpperBound(_getVertexIdManager(g)) == 2)
	SEQAN_TASSERT(getIdUpperBound(_getEdgeIdManager(g)) == 1)
	SEQAN_TASSERT(v1 == 1)
	SEQAN_TASSERT(numVertices(g) == 2)
	SEQAN_TASSERT(targetVertex(g, e) == 1)
	SEQAN_TASSERT(sourceVertex(g, e) == 0)
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(outDegree(g, v0) == 1)	
	SEQAN_TASSERT(inDegree(g, 1) == 1)
	SEQAN_TASSERT(inDegree(g, 0) == 1)	
	SEQAN_TASSERT(degree(g, 0) == 1)

	// Add more vertices and edges
	addVertex(g, id1, 10, 20);  //2
	TVertexDescriptor v3 = addVertex(g, id2, 1, 2);  //3
	addVertex(g, id1, 1, 3);  //4
	addEdge(g,3,4);
	TEdgeDescriptor my_edge = addEdge(g,3,1);
	addEdge(g,3,0);
	SEQAN_TASSERT(v3 == 3)
	SEQAN_TASSERT(numVertices(g) == 5)
	SEQAN_TASSERT(targetVertex(g, my_edge) == 3)
	SEQAN_TASSERT(sourceVertex(g, my_edge) == 1)
	SEQAN_TASSERT(numEdges(g) == 4)
	SEQAN_TASSERT(outDegree(g, v3) == 3)
	SEQAN_TASSERT(inDegree(g, v3) == 3)
	SEQAN_TASSERT(degree(g, v3) == 3)

	// Graph drawing
	// Raw output
	// File output
	fstream strm;
	strm.open(TEST_PATH "my_alignment.dot", ios_base::out | ios_base::trunc);
	write(strm,g,DotDrawing());
	strm.close();

	// Remove edges
	removeEdge(g,3,1);
	removeEdge(g,0,1);
	SEQAN_TASSERT(numEdges(g) == 2)

	
	// Remove vertices 
	addVertex(g, id2, 3, 6);  //5
	addEdge(g,5,2);
	addEdge(g,2,3);
	addEdge(g,1,3);
	addEdge(g,1,4);
	SEQAN_TASSERT(outDegree(g, 3) == 4)
	SEQAN_TASSERT(outDegree(g, 4) == 2)
	removeVertex(g, v3);
	SEQAN_TASSERT(outDegree(g, 4) == 1)
	SEQAN_TASSERT(outDegree(g, 0) == 0)
	SEQAN_TASSERT(numVertices(g) == 5)
	SEQAN_TASSERT(numEdges(g) == 2)

	// Clear graph
	clearEdges(g);
	SEQAN_TASSERT(numVertices(g) == 5)
	SEQAN_TASSERT(numEdges(g) == 0)
	addEdge(g,2,0);
	addEdge(g,4,1);
	clearVertices(g);
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	addVertex(g,id1,0,1);addVertex(g,id1,1,1);addVertex(g,id1,2,1);
	addVertex(g,id1,3,1);addVertex(g,id1,4,1);
	addEdge(g,2,0);
	addEdge(g,4,1);
	clear(g);
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	addVertex(g,id1,0,1);addVertex(g,id1,1,1);addVertex(g,id1,2,1);
	addVertex(g,id1,3,1);addVertex(g,id1,4,1);
	addEdge(g,2,0);
	addEdge(g,4,1);
	addEdge(g,4,2);
	removeVertex(g,3);
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(outDegree(g, 4) == 2)
	SEQAN_TASSERT(inDegree(g, 4) == 2)

	// Transpose
	transpose(g); 
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(outDegree(g, 4) == 2)
	SEQAN_TASSERT(inDegree(g, 4) == 2)
	TGraph g_copy(g);
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 2)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)
	addVertex(g_copy, id0, 0, 3);
	addEdge(g_copy,3,0);
	g_copy = g;
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 2)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)
	//Copies the graph and transposes just the copy
	transpose(g,g_copy);  // g does not change!
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 2)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)

	// Adjacency matrix
	Matrix<unsigned int> mat;
	getAdjacencyMatrix(g, mat);
	unsigned int len = getIdUpperBound(g.data_align.data_id_managerV);
	SEQAN_TASSERT(getValue(mat,0*len+2) == 1)
	SEQAN_TASSERT(getValue(mat,3*len+2) == 0)
	SEQAN_TASSERT(getValue(mat,0*len+2) == getValue(mat,2*len+0))
	SEQAN_TASSERT(getValue(mat,1*len+4) == getValue(mat,4*len+1))
	SEQAN_TASSERT(getValue(mat,2*len+4) == getValue(mat,4*len+2))

//____________________________________________________________________________
// Alignments with edge weights
	typedef Graph<Alignment<TStringSet> > TGraph2;
	typedef VertexDescriptor<TGraph2>::Type TVertexDescriptor2;
	typedef EdgeDescriptor<TGraph2>::Type TEdgeDescriptor2;


	TGraph2 g2(str);
	SEQAN_TASSERT(numEdges(g2) == 0)
	SEQAN_TASSERT(numVertices(g2) == 0)
	SEQAN_TASSERT(empty(g2) == true)

	TVertexDescriptor2 v20 = addVertex(g2, id1,0, 2);
	TVertexDescriptor2 v21 = addVertex(g2, id2, 0, 5);
	SEQAN_TASSERT(v20 == 0)
	SEQAN_TASSERT(v21 == 1)
	TEdgeDescriptor2 e2 =addEdge(g2,0,1, 100);
	SEQAN_TASSERT(getCargo(e2) == 100)
	SEQAN_TASSERT(numEdges(g2) == 1)
	removeEdge(g2,e2);
	SEQAN_TASSERT(numEdges(g2) == 0)
	e2 =addEdge(g2,0,1, 1005);
	SEQAN_TASSERT(numEdges(g2) == 1)
	removeOutEdges(g2,0);
	SEQAN_TASSERT(numEdges(g2) == 0)
	e2 =addEdge(g2,0,1, 1005);
	SEQAN_TASSERT(findEdge(g2, 0, 1) == e2)
	SEQAN_TASSERT(numEdges(g2) == 1)
	removeInEdges(g2,0);
	SEQAN_TASSERT(numEdges(g2) == 0)



	// Vertex Iterator
	typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	TVertexIterator itV(g);
	SEQAN_TASSERT(atBegin(itV)==true)
	SEQAN_TASSERT(getValue(itV)==0)
	SEQAN_TASSERT(value(itV)==0)
	SEQAN_TASSERT(getValue(itV)==0)
	goNext(itV);
	SEQAN_TASSERT(atBegin(itV)==false)
	SEQAN_TASSERT(getValue(itV)==1)
	++itV;
	SEQAN_TASSERT(getValue(itV)==2)
	SEQAN_TASSERT(atEnd(itV)==false)
	goPrevious(itV);
	SEQAN_TASSERT((*itV)==1)
	SEQAN_TASSERT(atEnd(itV)==false)
	itV--;
	SEQAN_TASSERT(getValue(itV)==0)
	SEQAN_TASSERT(atBegin(itV)==true)

	// OutEdge Iterator
	typedef Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	TOutEdgeIterator it(g, v0);
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==2)
	SEQAN_TASSERT(sourceVertex(g, value(it))==0)
	SEQAN_TASSERT(targetVertex(g, *it)==2)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
	goNext(it);
	SEQAN_TASSERT(atEnd(it)==true)
	SEQAN_TASSERT(atBegin(it)==false)
	goPrevious(it);
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==2)
	--it;
	it--;
	SEQAN_TASSERT(atBegin(it)==true)
	TOutEdgeIterator it2(g, v0);
	TOutEdgeIterator it3;
	it3 = it;
	SEQAN_TASSERT(it == it2)
	SEQAN_TASSERT(it2 == it3)
	goEnd(it);
	SEQAN_TASSERT(it2 != it)
	goEnd(it2);
	SEQAN_TASSERT(it2 == it)
	goBegin(it2);
	SEQAN_TASSERT(it2 != it)
	SEQAN_TASSERT(&g == &hostGraph(it))

	// EdgeIterator
	typedef Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	TEdgeIterator itEdge(g);
	SEQAN_TASSERT(sourceVertex(g, getValue(itEdge))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(itEdge))==2)
	SEQAN_TASSERT(sourceVertex(g, value(itEdge))==0)
	SEQAN_TASSERT(targetVertex(g, *itEdge)==2)
	SEQAN_TASSERT(atEnd(itEdge)==false)
	SEQAN_TASSERT(atBegin(itEdge)==true)
	goNext(itEdge);
	SEQAN_TASSERT(sourceVertex(g, value(itEdge))==1)
	SEQAN_TASSERT(targetVertex(g, *itEdge)==4)
	++itEdge;
	--itEdge;
	SEQAN_TASSERT(sourceVertex(g, value(itEdge))==1)
	SEQAN_TASSERT(targetVertex(g, *itEdge)==4)
	goEnd(itEdge);
	SEQAN_TASSERT(atEnd(itEdge)==true)
	SEQAN_TASSERT(atBegin(itEdge)==false)
	goBegin(itEdge);
	SEQAN_TASSERT(atEnd(itEdge)==false)
	SEQAN_TASSERT(atBegin(itEdge)==true)

	// Adjacency Iterator
	typedef Iterator<TGraph, AdjacencyIterator>::Type TAdjacencyIterator;
	TAdjacencyIterator itAdj(g,2);
	SEQAN_TASSERT(getValue(itAdj) == 4)
	SEQAN_TASSERT(&hostGraph(itAdj) == &g)
	SEQAN_TASSERT(value(itAdj) == 4)
	SEQAN_TASSERT(*itAdj == 4)
	SEQAN_TASSERT(atEnd(itAdj)==false)
	SEQAN_TASSERT(atBegin(itAdj)==true)
	goNext(itAdj);
	SEQAN_TASSERT(*itAdj == 0)
	SEQAN_TASSERT(atEnd(itAdj)==false)
	SEQAN_TASSERT(atBegin(itAdj)==false)
	++itAdj;
	SEQAN_TASSERT(atEnd(itAdj)==true)
	SEQAN_TASSERT(atBegin(itAdj)==false)
	goPrevious(itAdj);--itAdj;
	SEQAN_TASSERT(*itAdj == 4)
	goBegin(itAdj);
	SEQAN_TASSERT(atBegin(itAdj)==true)
	goEnd(itAdj);
	SEQAN_TASSERT(atEnd(itAdj)==true)

	// Bfs Iterator
	typedef Iterator<TGraph, BfsIterator>::Type TBfsIterator;
	TBfsIterator bfsIt(g,2);
	SEQAN_TASSERT(atEnd(bfsIt)==false)
	SEQAN_TASSERT(atBegin(bfsIt)==true)
	++bfsIt;
	SEQAN_TASSERT(getValue(bfsIt) == 4)
	SEQAN_TASSERT(&hostGraph(bfsIt) == &g)
	SEQAN_TASSERT(value(bfsIt) == 4)
	SEQAN_TASSERT(*bfsIt == 4)
	goNext(bfsIt);
	SEQAN_TASSERT(value(bfsIt) == 0)

	// Dfs Iterator
	typedef Iterator<TGraph, DfsPreorder>::Type TDfsPreorder;
	TDfsPreorder dfsIt(g,2);
	SEQAN_TASSERT(atEnd(dfsIt)==false)
	SEQAN_TASSERT(atBegin(dfsIt)==true)
	SEQAN_TASSERT(*dfsIt == 2)
	++dfsIt;
	SEQAN_TASSERT(getValue(dfsIt) == 0)
	SEQAN_TASSERT(&hostGraph(dfsIt) == &g)
	SEQAN_TASSERT(value(dfsIt) == 0)
	SEQAN_TASSERT(*dfsIt == 0)
	goNext(dfsIt);
}


//////////////////////////////////////////////////////////////////////////////

template <typename TStringSet>
void Test_StringSet() {
	typedef	typename Id<TStringSet>::Type TId;

	TStringSet str;
	String<char> bla("a");
	TId id0 = assignValueById(str, bla);
	SEQAN_TASSERT(id0 == 0)
	SEQAN_TASSERT(length(str) == 1)
	SEQAN_TASSERT(str[0] == "a")
	SEQAN_TASSERT(getValueById(str, id0) == "a")
	String<char> bla1("b");
	TId id1 = assignValueById(str, bla1);
	SEQAN_TASSERT(id1 == 1)
	SEQAN_TASSERT(str[1] == "b")
	SEQAN_TASSERT(length(str) == 2)
	SEQAN_TASSERT(getValueById(str, id1) == "b")
	String<char> bla2("c");
	TId id2 = assignValueById(str, bla2);
	SEQAN_TASSERT(id2 == 2)
	SEQAN_TASSERT(str[2] == "c")
	SEQAN_TASSERT(length(str) == 3)
	SEQAN_TASSERT(getValueById(str, id2) == "c")
	String<char> bla3("d");
	TId id3 = assignValueById(str, bla3);
	SEQAN_TASSERT(id3 == 3)
	SEQAN_TASSERT(str[3] == "d")
	SEQAN_TASSERT(length(str) == 4)
	SEQAN_TASSERT(getValueById(str, id3) == "d")
	removeValueById(str,id1);
	SEQAN_TASSERT(getValueById(str, id0) == "a")
	SEQAN_TASSERT(getValueById(str, id1) == "")
	SEQAN_TASSERT(getValueById(str, id2) == "c")
	SEQAN_TASSERT(getValueById(str, id3) == "d")
	SEQAN_TASSERT(length(str) == 3)
	removeValueById(str,id2);
	SEQAN_TASSERT(getValueById(str, id0) == "a")
	SEQAN_TASSERT(getValueById(str, id1) == "")
	SEQAN_TASSERT(getValueById(str, id2) == "")
	SEQAN_TASSERT(getValueById(str, id3) == "d")
	SEQAN_TASSERT(length(str) == 2)
	removeValueById(str,id3);
	SEQAN_TASSERT(getValueById(str, id0) == "a")
	SEQAN_TASSERT(getValueById(str, id1) == "")
	SEQAN_TASSERT(getValueById(str, id2) == "")
	SEQAN_TASSERT(getValueById(str, id3) == "")
	SEQAN_TASSERT(length(str) == 1)
	id1 = assignValueById(str, bla1);
	id2 = assignValueById(str, bla2);
	id3 = assignValueById(str, bla3);	
	SEQAN_TASSERT(getValueById(str, id0) == "a")
	SEQAN_TASSERT(getValueById(str, id1) == "b")
	SEQAN_TASSERT(getValueById(str, id2) == "c")
	SEQAN_TASSERT(getValueById(str, id3) == "d")
	SEQAN_TASSERT(length(str) == 4)
	TStringSet subSet;
	String<unsigned int> ids;
	appendValue(ids, id1);
	appendValue(ids, id3);
	subset(str, subSet, ids);
	SEQAN_TASSERT(getValueById(subSet, id0) == "")
	SEQAN_TASSERT(getValueById(subSet, id1) == "b")
	SEQAN_TASSERT(getValueById(subSet, id2) == "")
	SEQAN_TASSERT(getValueById(subSet, id3) == "d")
	SEQAN_TASSERT(subSet[0] == "b")
	SEQAN_TASSERT(subSet[1] == "d")
	SEQAN_TASSERT(length(subSet) == 2)
	TStringSet subSet2;
	String<unsigned int> ids2;
	appendValue(ids2, id3);
	subset(subSet, subSet2, ids2);
	SEQAN_TASSERT(getValueById(subSet2, id0) == "")
	SEQAN_TASSERT(getValueById(subSet2, id1) == "")
	SEQAN_TASSERT(getValueById(subSet2, id2) == "")
	SEQAN_TASSERT(getValueById(subSet2, id3) == "d")
	SEQAN_TASSERT(length(subSet2) == 1)
	clear(subSet);
	assignValueById(subSet,subSet2, id3);
	assignValueById(subSet,str, id0);
	SEQAN_TASSERT(valueById(subSet, id0) == "a")
	SEQAN_TASSERT(getValueById(subSet, id1) == "")
	SEQAN_TASSERT(getValueById(subSet, id2) == "")
	SEQAN_TASSERT(valueById(subSet, id3) == "d")
	SEQAN_TASSERT(length(subSet) == 2)
	String<char> bla5("f");
	assignValueById(subSet, bla5, id3);
	SEQAN_TASSERT(getValueById(subSet, id3) == "f")
	SEQAN_TASSERT(length(subSet) == 2)
	assignValueById(subSet, bla3, 20);
	SEQAN_TASSERT(getValueById(subSet, 20) == "d")
	SEQAN_TASSERT(length(subSet) == 3)
}



//////////////////////////////////////////////////////////////////////////////
template <typename TGraphType>
void Test_VertexIterator() {
//____________________________________________________________________________
// Graph InternalVertexIterator
	typedef Graph<TGraphType> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	
	TGraph g;
	TVertexDescriptor v0 = addVertex(g);
	TVertexDescriptor v1 = addVertex(g);
	addEdge(g,v0,v1,'t');
	addVertex(g); //2
	addVertex(g); //3
	addVertex(g); //4

	//Tricky case -> id 0 is released
	removeVertex(g,v0);
	removeVertex(g,3);
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	TVertexIterator it(g);
	SEQAN_TASSERT(atBegin(it)==true)
	SEQAN_TASSERT(getValue(it)==1)
	SEQAN_TASSERT(value(it)==1)
	SEQAN_TASSERT(getValue(it)==1)
	goNext(it);
	SEQAN_TASSERT(atBegin(it)==false)
	SEQAN_TASSERT(getValue(it)==2)
	++it;
	SEQAN_TASSERT(getValue(it)==4)
	SEQAN_TASSERT(atEnd(it)==false)
	it++; 
	SEQAN_TASSERT(atEnd(it)==true)
	++it; 
	SEQAN_TASSERT(atEnd(it)==true)
	goPrevious(it);
	// No assignment to vertex iterators
	// *it = 3;
	SEQAN_TASSERT((*it)==4)
	SEQAN_TASSERT(atEnd(it)==false)
	it--;
	SEQAN_TASSERT(getValue(it)==2)
	SEQAN_TASSERT(atBegin(it)==false)
	--it;
	SEQAN_TASSERT(atBegin(it)==true)
	SEQAN_TASSERT(getValue(it)==1)
	--it;
	SEQAN_TASSERT(atBegin(it)==true)
	SEQAN_TASSERT(getValue(it)==1)
	TVertexIterator it2(g);
	TVertexIterator it3;
	it3 = it;
	SEQAN_TASSERT(it == it2)
	SEQAN_TASSERT(it2 == it3)
	goEnd(it);
	SEQAN_TASSERT(it2 != it)
	goEnd(it2);
	SEQAN_TASSERT(it2 == it)
	goBegin(it2);
	SEQAN_TASSERT(it2 != it)
	SEQAN_TASSERT(&hostGraph(it) == &g)
}

//////////////////////////////////////////////////////////////////////////////
void Test_TreeInternalVertexIterator() {
//____________________________________________________________________________
// Tree InternalVertexIterator
	typedef Graph<Tree<void> > TTree;
	typedef VertexDescriptor<TTree>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TTree>::Type TEdgeDescriptor;
	typedef Size<TTree>::Type TSize;
	
	TTree gV;
	TSize numEdges = 8;
	//Parent, Child, Parent, Child, ...
	//The root must be the first vertex
	TVertexDescriptor edges[] = {0,8, 0,3, 0,2, 0,1, 2,4, 4,5, 5,7, 5,6};
	addEdges(gV,edges, numEdges);

	typedef Iterator<TTree, VertexIterator>::Type TVertexIterator;
	TVertexIterator itV(gV);
	SEQAN_TASSERT(atBegin(itV)==true)
	SEQAN_TASSERT(getValue(itV)==0)
	SEQAN_TASSERT(value(itV)==0)
	SEQAN_TASSERT(getValue(itV)==0)
	goNext(itV);
	SEQAN_TASSERT(atBegin(itV)==false)
	SEQAN_TASSERT(getValue(itV)==1)
	++itV;
	SEQAN_TASSERT(getValue(itV)==2)
	SEQAN_TASSERT(atEnd(itV)==false)
	goPrevious(itV);
	SEQAN_TASSERT((*itV)==1)
	SEQAN_TASSERT(atEnd(itV)==false)
	itV--;
	SEQAN_TASSERT(getValue(itV)==0)
	SEQAN_TASSERT(atBegin(itV)==true)
}

//////////////////////////////////////////////////////////////////////////////

template <typename TGraphType>
void Test_OutEdgeIterator() {
//____________________________________________________________________________
// Graph InternalOutEdgeIterator
	typedef Graph<TGraphType> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	
	TGraph g;
	TVertexDescriptor v0 = addVertex(g);
	TVertexDescriptor v1 = addVertex(g);
	addVertex(g); //2
	addEdge(g,0,2,'t');
	addEdge(g,v0,v1,'a');
	addVertex(g); //3
	addVertex(g); //4

	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	TOutEdgeIterator it(g, v0);
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==1)
	SEQAN_TASSERT(sourceVertex(g, value(it))==0)
	SEQAN_TASSERT(targetVertex(g, *it)==1)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
	goNext(it);
	// Slow
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==2)
	// Fast
	SEQAN_TASSERT(sourceVertex(it)==0)
	SEQAN_TASSERT(targetVertex(it)==2)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==false)
	++it;
	SEQAN_TASSERT(atEnd(it)==true)
	SEQAN_TASSERT(atBegin(it)==false)
	goPrevious(it);
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==2)
	--it;
	// Slow
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==1)
	// Fast
	SEQAN_TASSERT(sourceVertex(it)==0)
	SEQAN_TASSERT(targetVertex(it)==1)
	it++; 
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==2)
	it--;
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==1)
	SEQAN_TASSERT(atBegin(it)==true)
	TOutEdgeIterator it2(g, v0);
	TOutEdgeIterator it3;
	it3 = it;
	SEQAN_TASSERT(it == it2)
	SEQAN_TASSERT(it2 == it3)
	goEnd(it);
	SEQAN_TASSERT(it2 != it)
	goEnd(it2);
	SEQAN_TASSERT(it2 == it)
	goBegin(it2);
	SEQAN_TASSERT(it2 != it)
	SEQAN_TASSERT(&g == &hostGraph(it))
}

//////////////////////////////////////////////////////////////////////////////

	
template <typename TGraphType>
void Test_EdgeIterator() {
//____________________________________________________________________________
// Graph InternalEdgeIterator
	typedef Graph<TGraphType> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	
	TGraph g;
	TVertexDescriptor v0 = addVertex(g);
	TVertexDescriptor v1 = addVertex(g);
	addVertex(g); //2
	addEdge(g,0,2,'t');
	addEdge(g,v0,v1,'a');
	addVertex(g); //3
	addVertex(g); //4
	addVertex(g); //5
	addEdge(g,3,4,'t');
	addEdge(g,2,3,'a');
	addEdge(g,4,5,'a');

	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	TEdgeIterator it(g);
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==1)
	SEQAN_TASSERT(sourceVertex(g, value(it))==0)
	SEQAN_TASSERT(targetVertex(g, *it)==1)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
	goNext(it);
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==2)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==false)
	++it;
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==false)
	// Slow
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==2)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==3)
	// Fast
	SEQAN_TASSERT(sourceVertex(it)==2)
	SEQAN_TASSERT(targetVertex(it)==3)
	it++;
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==3)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==4)
	it++;
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==4)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==5)
	++it;
	SEQAN_TASSERT(atEnd(it)==true)
	SEQAN_TASSERT(atBegin(it)==false)
	++it;
	++it;
	goPrevious(it);
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==4)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==5)
	--it;
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==3)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==4)
	it--;
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==2)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==3)
	TEdgeIterator it2(g);
	TEdgeIterator it3;
	goBegin(it);
	it3 = it;
	SEQAN_TASSERT(it == it2)
	SEQAN_TASSERT(it2 == it3)
	goEnd(it);
	SEQAN_TASSERT(it2 != it)
	goEnd(it2);
	SEQAN_TASSERT(it2 == it)
	goBegin(it2);
	SEQAN_TASSERT(it2 != it)
	SEQAN_TASSERT(&hostGraph(it) == &g)
}

//////////////////////////////////////////////////////////////////////////////

template <typename TGraphType>
void Test_AdjacencyIterator() {
//____________________________________________________________________________
// Graph InternalAdjacencyIterator
	typedef Graph<TGraphType> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	
	TGraph g;
	TVertexDescriptor v0 = addVertex(g);
	TVertexDescriptor v1 = addVertex(g);
	addVertex(g); //2
	addEdge(g,0,2,'t');
	addEdge(g,v0,v1,'a');
	addVertex(g); //3
	addVertex(g); //4
	addVertex(g); //5
	addEdge(g,3,2,'t');
	addEdge(g,3,4,'a');
	addEdge(g,4,5,'t');

	typedef typename Iterator<TGraph, AdjacencyIterator>::Type TAdjacencyIterator;
	TAdjacencyIterator it(g,3);
	SEQAN_TASSERT(getValue(it) == 4)
	SEQAN_TASSERT(&hostGraph(it) == &g)
	SEQAN_TASSERT(value(it) == 4)
	SEQAN_TASSERT(*it == 4)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==2)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==false)
	++it;
	SEQAN_TASSERT(atEnd(it)==true)
	SEQAN_TASSERT(atBegin(it)==false)
	goBegin(it);
	it++;
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==false)
	it++;
	goPrevious(it);
	SEQAN_TASSERT(getValue(it)==2)
	--it;
	SEQAN_TASSERT(getValue(it)==4)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
	goEnd(it);
	it--;
	SEQAN_TASSERT(getValue(it)==2)
	goBegin(it);
	TAdjacencyIterator it2(g, 3);
	TAdjacencyIterator it3;
	it3 = it;
	SEQAN_TASSERT(it == it2)
	SEQAN_TASSERT(it2 == it3)
	goEnd(it);
	SEQAN_TASSERT(it2 != it)
	goEnd(it2);
	SEQAN_TASSERT(it2 == it)
	goBegin(it2);
	SEQAN_TASSERT(it2 != it)
}

//////////////////////////////////////////////////////////////////////////////

template <typename TGraphType>
void Test_ExternalProperty() {
//____________________________________________________________________________
// Graph external property maps
	typedef Graph<TGraphType> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TGraph>::Type TSize;
	
	TGraph g;
	TVertexDescriptor v0 = addVertex(g);
	TVertexDescriptor v1 = addVertex(g);
	addVertex(g);
	TEdgeDescriptor e2 =addEdge(g,0,2,'t');
	TEdgeDescriptor e1 =addEdge(g,v0,v1,'a');

	
	// Test external property maps
	String<int> dMap;
	resizeVertexMap(g,dMap);

	String<char> eMap;
	resizeEdgeMap(g,eMap);

	assignProperty(dMap, v0, 3);
	assignProperty(dMap, v1, 1);
	assignProperty(eMap, e1, 'a');
	assignProperty(eMap, e2, 'b');
	SEQAN_TASSERT(getProperty(dMap, v0) == 3)
	SEQAN_TASSERT(getProperty(dMap, v1) == 1)
	SEQAN_TASSERT(getProperty(eMap, e2) == 'b')
	SEQAN_TASSERT(getProperty(eMap, e1) == 'a')
	property(dMap, v1) = 2;
	property(eMap, e2) = 'c';
	SEQAN_TASSERT(getProperty(dMap, v1) == 2)
	SEQAN_TASSERT(getProperty(eMap, e2) == 'c')

	String<int> const dMap2(dMap);
	SEQAN_TASSERT(getProperty(dMap2, v0) == 3)
	SEQAN_TASSERT(getProperty(dMap2, v1) == 2)
	SEQAN_TASSERT(property(dMap2, v1) == 2)

	clear(g);
	addVertex(g);addVertex(g);addVertex(g);

	char names[] = {'r', 's','t'};
	String<char> nameMap;
	resizeVertexMap(g,nameMap, names);
	SEQAN_TASSERT(getProperty(nameMap, v0) == 'r')
	SEQAN_TASSERT(getProperty(nameMap, v1) == 's')
}


//////////////////////////////////////////////////////////////////////////////

void Test_Property() {
//____________________________________________________________________________
// Graph properties
	typedef Pair<char, int> TPair;
	typedef Directed<TPair> TEdges;
	typedef Graph<TEdges> TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;

	// Create a simple graph
	TGraph g;
	TVertexDescriptor v0 = addVertex(g);
	TVertexDescriptor v1 = addVertex(g);
	addVertex(g);
	TEdgeDescriptor e1 =addEdge(g,0,2);
	TEdgeDescriptor e2 =addEdge(g,v0,v1);

	// First Variant: Explicit internal map with member Ids
	InternalMap<TPair, 1> eMap1; // This property map is used to access the first member
	InternalMap<TPair, 2> eMap2; // This property map is used to access the second member
	resizeEdgeMap(g,eMap1);
	resizeEdgeMap(g,eMap2);
	assignProperty(eMap1, e1, 'a');
	assignProperty(eMap2, e1, 20);
	assignProperty(eMap1, e2, 'b');
	assignProperty(eMap2, e2, 50);
	SEQAN_TASSERT(getProperty(eMap1, e1) == 'a')
	SEQAN_TASSERT(getProperty(eMap2, e1) == 20)
	SEQAN_TASSERT(getProperty(eMap1, e2) == 'b')
	SEQAN_TASSERT(getProperty(eMap2, e2) == 50)
	// Note: That these properties are stored inside the cargo of each edge
	SEQAN_TASSERT(getCargo(e1).i1 == 'a')
	SEQAN_TASSERT(getCargo(e1).i2 == 20)
	SEQAN_TASSERT(getCargo(e2).i1 == 'b')
	SEQAN_TASSERT(getCargo(e2).i2 == 50)
	assignProperty(eMap1, e1, 'c');
	assignProperty(eMap2, e1, 10);
	SEQAN_TASSERT(getProperty(eMap1, e1) == 'c')
	SEQAN_TASSERT(getProperty(eMap2, e1) == 10)
	SEQAN_TASSERT(property(eMap1, e1) == 'c')
	SEQAN_TASSERT(property(eMap2, e1) == 10)
	InternalMap<TPair, 1> const eMap3(eMap1);
	InternalMap<TPair, 2> const eMap31(eMap2);
	SEQAN_TASSERT(getProperty(eMap3, e1) == 'c')
	SEQAN_TASSERT(getProperty(eMap3, e2) == 'b')
	SEQAN_TASSERT(getProperty(eMap31, e1) == 10)
	SEQAN_TASSERT(property(eMap31, e1) == 10)
	SEQAN_TASSERT(property(eMap3, e2) == 'b')
	// Create a simple graph with unsigned int cargo
	typedef EdgeDescriptor<Graph<Directed<unsigned int> > >::Type TEdgeDescriptor2;
	Graph<Directed<unsigned int> > g2;
	addVertex(g2);
	addVertex(g2);
	TEdgeDescriptor2 edge1 =addEdge(g2,v0,v0);
	addEdge(g2,0,1);
	InternalMap<unsigned int> edgeMap;
	resizeEdgeMap(g2,edgeMap);
	assignProperty(edgeMap, edge1 ,3);
	SEQAN_TASSERT(getProperty(edgeMap, edge1) == 3)
	SEQAN_TASSERT(property(edgeMap, edge1) == 3)
	InternalMap<unsigned int> const edgeMap2(edgeMap);
	SEQAN_TASSERT(getProperty(edgeMap2, edge1) == 3)
	SEQAN_TASSERT(property(edgeMap2, edge1) == 3)

	// Second Variant: Pointer to member using a class
	InternalPointerMap<char TPair:: *, &TPair::i1> eMap4;
	InternalPointerMap<int TPair:: *, &TPair::i2> eMap5;
	resizeEdgeMap(g,eMap4);
	resizeEdgeMap(g,eMap5);
	assignProperty(eMap4, e1, 'c');
	assignProperty(eMap5, e1, 10);
	assignProperty(eMap4, e2, 'd');
	assignProperty(eMap5, e2, 30);
	SEQAN_TASSERT(getProperty(eMap4, e1) == 'c')
	SEQAN_TASSERT(getProperty(eMap5, e1) == 10)
	SEQAN_TASSERT(getProperty(eMap4, e2) == 'd')
	SEQAN_TASSERT(getProperty(eMap5, e2) == 30)
	property(eMap4,e1)='z';
	property(eMap5,e1)=100;
	SEQAN_TASSERT(getProperty(eMap4, e1) == 'z')
	SEQAN_TASSERT(getProperty(eMap5, e1) == 100)
	InternalPointerMap<char TPair:: *, &TPair::i1> const eMap6(eMap4);
	SEQAN_TASSERT(getProperty(eMap6, e1) == 'z')
	SEQAN_TASSERT(getProperty(eMap6, e2) == 'd')
	SEQAN_TASSERT(property(eMap6, e2) == 'd')
	
	// Third Variant: Raw pointer to member
	char TPair:: * pseudo_map = &TPair::i1;
	assignProperty(pseudo_map, e1, 'z');
	assignProperty(pseudo_map, e2, 'w');
	SEQAN_TASSERT(getProperty(pseudo_map, e1) == 'z')
	SEQAN_TASSERT(getProperty(pseudo_map, e2) == 'w')
	property(pseudo_map,e1)='k';
	SEQAN_TASSERT(getProperty(pseudo_map, e1) == 'k')
}

//////////////////////////////////////////////////////////////////////////////

template <typename TGraphType>
void Test_BfsIter() {
//____________________________________________________________________________
// Graph BfsIterator
	typedef Graph<TGraphType> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TGraph>::Type TSize;

	//Create the graph
	TGraph g;
	addVertex(g);addVertex(g);addVertex(g);addVertex(g);
	addEdge(g,0,2,'t');
	addEdge(g,0,1,'a');
	addEdge(g,1,3,'a');
	
	typedef typename Iterator<TGraph, BfsIterator>::Type TBfsIterator;
	TBfsIterator it(g,0);
	SEQAN_TASSERT(getValue(it) == 0)
	SEQAN_TASSERT(&hostGraph(it) == &g)
	SEQAN_TASSERT(value(it) == 0)
	SEQAN_TASSERT(*it == 0)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==1)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==false)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==2)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==3)
	goNext(it);
	SEQAN_TASSERT(atEnd(it)==true)
	SEQAN_TASSERT(atBegin(it)==false)
	goBegin(it);
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
}

//////////////////////////////////////////////////////////////////////////////

void Test_BfsIterator() {
//____________________________________________________________________________
// Graph InternalBfsIterator
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Create the graph
	Graph<> g;
	addVertex(g);addVertex(g);addVertex(g);addVertex(g);
	addVertex(g);addVertex(g);addVertex(g);addVertex(g);
	addEdge(g,1,0);
	addEdge(g,1,5);
	addEdge(g,0,4);
	addEdge(g,5,2);
	addEdge(g,5,6);
	addEdge(g,2,6);
	addEdge(g,2,3);
	addEdge(g,6,3);
	addEdge(g,6,7);
	addEdge(g,3,7);
	typedef Iterator<Graph<>, BfsIterator>::Type TBfsIterator;
	TBfsIterator it(g,1);
	SEQAN_TASSERT(getValue(it) == 1)
	SEQAN_TASSERT(&hostGraph(it) == &g)
	SEQAN_TASSERT(value(it) == 1)
	SEQAN_TASSERT(*it == 1)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==5)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==false)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==0)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==6)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==2)
	++it;
	SEQAN_TASSERT(getValue(it)==4)
	it++;
	SEQAN_TASSERT(getValue(it)==7)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==3)
	goNext(it);
	SEQAN_TASSERT(atEnd(it)==true)
	SEQAN_TASSERT(atBegin(it)==false)
	goBegin(it);
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
	TBfsIterator it2(g, 1);
	TBfsIterator it3;
	it3 = it;
	SEQAN_TASSERT(it == it2)
	SEQAN_TASSERT(it2 == it3)
	goEnd(it);
	SEQAN_TASSERT(it2 != it)
	goEnd(it2);
	SEQAN_TASSERT(it2 == it)
	goBegin(it2);
	SEQAN_TASSERT(it2 != it)
}

//////////////////////////////////////////////////////////////////////////////


template <typename TGraphType>
void Test_DfsPreorderIter() {
//____________________________________________________________________________
// Graph DfsIterator
	typedef Graph<TGraphType> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TGraph>::Type TSize;

	//Create the graph
	TGraph g;
	addVertex(g);addVertex(g);addVertex(g);addVertex(g);
	addEdge(g,0,1,'t');
	addEdge(g,0,2,'a');
	addEdge(g,1,3,'a');
	
	typedef typename Iterator<TGraph, DfsPreorder>::Type TDfsIterator;
	TDfsIterator it(g,0);
	SEQAN_TASSERT(getValue(it) == 0)
	SEQAN_TASSERT(&hostGraph(it) == &g)
	SEQAN_TASSERT(value(it) == 0)
	SEQAN_TASSERT(*it == 0)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==1)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==false)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==3)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==2)
	goNext(it);
	SEQAN_TASSERT(atEnd(it)==true)
	SEQAN_TASSERT(atBegin(it)==false)
	goBegin(it);
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
}


//////////////////////////////////////////////////////////////////////////////
void Test_DfsPreorderIterator() {
//____________________________________________________________________________
// Graph DfsIterator
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Create the graph
	Graph<> g;
	addVertex(g);addVertex(g);addVertex(g);addVertex(g);
	addVertex(g);addVertex(g);addVertex(g);addVertex(g);
	addEdge(g,1,0);
	addEdge(g,1,5);
	addEdge(g,0,4);
	addEdge(g,5,2);
	addEdge(g,5,6);
	addEdge(g,2,6);
	addEdge(g,2,3);
	addEdge(g,6,3);
	addEdge(g,6,7);
	addEdge(g,3,7);

	typedef Iterator<Graph<>, DfsPreorder>::Type TDfsPreorder;
	TDfsPreorder it(g,1);
	SEQAN_TASSERT(getValue(it) == 1)
	SEQAN_TASSERT(&hostGraph(it) == &g)
	SEQAN_TASSERT(value(it) == 1)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==0)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==false)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==4)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==5)
	SEQAN_TASSERT(value(it)==5)
	SEQAN_TASSERT(*it == 5)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==2)
	++it;
	SEQAN_TASSERT(getValue(it)==3)
	it++;
	SEQAN_TASSERT(getValue(it)==7)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==6)
	goNext(it);
	SEQAN_TASSERT(atEnd(it)==true)
	SEQAN_TASSERT(atBegin(it)==false)
	goBegin(it);
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
	TDfsPreorder it2(g, 1);
	TDfsPreorder it3;
	it3 = it;
	SEQAN_TASSERT(it == it2)
	SEQAN_TASSERT(it2 == it3)
	goEnd(it);
	SEQAN_TASSERT(it2 != it)
	goEnd(it2);
	SEQAN_TASSERT(it2 == it)
	goBegin(it2);
	SEQAN_TASSERT(it2 != it)
}


//////////////////////////////////////////////////////////////////////////////

void Test_Oracle() {
	Graph<Automaton<char> > g;
	createOracleOnReverse(g,"announce");
	SEQAN_TASSERT(parseString(g, 0, "e") == 1)
	SEQAN_TASSERT(parseString(g, 0, "ec") == 2)
	SEQAN_TASSERT(parseString(g, 0, "n") == 3)
	SEQAN_TASSERT(parseString(g, 0, "a") == 8)
	SEQAN_TASSERT(parseString(g, 0, "nn") == 7)

	Graph<Automaton<Dna> > g2;
	createOracle(g2,"ATATA");
	SEQAN_TASSERT(parseString(g2, 0, "A") == 1)
	SEQAN_TASSERT(parseString(g2, 0, "T") == 2)
	SEQAN_TASSERT(parseString(g2, 0, "AT") == 2)
}

//////////////////////////////////////////////////////////////////////////////

void Test_Trie() {
	Graph<Automaton<char> > g;
	String<String<unsigned int> > pos;
	String<String<char> > keywords;
	appendValue(keywords, String<char>("announce"));
	appendValue(keywords, String<char>("annual"));
	appendValue(keywords, String<char>("annually"));
	createTrie(g,pos,keywords);

	// Output
	// File output
	fstream strm;
	strm.open(TEST_PATH "my_trie.dot", ios_base::out | ios_base::trunc);
	String<String<char> > nodeMap;
	_createTrieNodeAttributes(g, pos, nodeMap);
	String<String<char> > edgeMap;
	_createEdgeAttributes(g,edgeMap);
	write(strm,g,nodeMap,edgeMap,DotDrawing());
	strm.close();

	Graph<Automaton<Dna> > gDna;
	clear(pos);
	String<String<Dna> > keyw;
	appendValue(keyw, String<Dna>("ATATATA"));
	appendValue(keyw, String<Dna>("TATAT"));
	appendValue(keyw, String<Dna>("ACGATAT"));
	createTrie(gDna,pos,keyw);

	// Output
	// File output
	fstream strm2;
	strm2.open(TEST_PATH "my_trie_dna.dot", ios_base::out | ios_base::trunc);
	clear(nodeMap);
	_createTrieNodeAttributes(gDna, pos, nodeMap);
	clear(edgeMap);
	_createEdgeAttributes(gDna,edgeMap);
	write(strm2,gDna,nodeMap,edgeMap,DotDrawing());
	strm2.close();
}


//////////////////////////////////////////////////////////////////////////////
void Test_BreadthFirstSearch() {
//____________________________________________________________________________
// Breadth-First Search
	typedef Graph<Undirected<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,4, 1,5, 2,5, 2,6, 2,3, 3,6, 3,7, 5,6, 6,7};
	
	//Create the graph
	TGraph g;
	addEdges(g,edges, numEdges);

	// Predecessor and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// Bfs
	breadth_first_search(g, 1, predMap, distMap);
	
	SEQAN_TASSERT(getProperty(distMap, 0) == 1)
	SEQAN_TASSERT(getProperty(distMap, 1) == 0)
	SEQAN_TASSERT(getProperty(distMap, 2) == 2)
	SEQAN_TASSERT(getProperty(distMap, 3) == 3)
	SEQAN_TASSERT(getProperty(distMap, 4) == 2)
	SEQAN_TASSERT(getProperty(distMap, 5) == 1)
	SEQAN_TASSERT(getProperty(distMap, 6) == 2)
	SEQAN_TASSERT(getProperty(distMap, 7) == 3)
	SEQAN_TASSERT(getProperty(predMap, 0) == 1)
	SEQAN_TASSERT(getProperty(predMap, 1) == getNilPredecessor(g))
	SEQAN_TASSERT(getProperty(predMap, 2) == 5)
	SEQAN_TASSERT(getProperty(predMap, 3) == 6)
	SEQAN_TASSERT(getProperty(predMap, 4) == 0)
	SEQAN_TASSERT(getProperty(predMap, 5) == 1)
	SEQAN_TASSERT(getProperty(predMap, 6) == 5)
	SEQAN_TASSERT(getProperty(predMap, 7) == 6)
}

//////////////////////////////////////////////////////////////////////////////

void Test_DepthFirstSearch() {
//____________________________________________________________________________
// Depth-First Search
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 8;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,3, 0,1, 1,4, 2,4, 2,5, 3,1, 4,3, 5,5};

	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);

	// Predecessor and distance map
	String<unsigned int> predMap;
	String<unsigned int> discoveryTimeMap;
	String<unsigned int> finishingTimeMap;

	// Dfs
	depth_first_search(g, predMap, discoveryTimeMap, finishingTimeMap);

	SEQAN_TASSERT(getProperty(discoveryTimeMap, 0) == 1)
	SEQAN_TASSERT(getProperty(discoveryTimeMap, 1) == 2)
	SEQAN_TASSERT(getProperty(discoveryTimeMap, 2) == 9)
	SEQAN_TASSERT(getProperty(discoveryTimeMap, 3) == 4)
	SEQAN_TASSERT(getProperty(discoveryTimeMap, 4) == 3)
	SEQAN_TASSERT(getProperty(discoveryTimeMap, 5) == 10)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 0) == 8)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 1) == 7)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 2) == 12)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 3) == 5)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 4) == 6)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 5) == 11)
	SEQAN_TASSERT(getProperty(predMap, 0) == getNilPredecessor(g))
	SEQAN_TASSERT(getProperty(predMap, 1) == 0)
	SEQAN_TASSERT(getProperty(predMap, 2) == getNilPredecessor(g))
	SEQAN_TASSERT(getProperty(predMap, 3) == 4)
	SEQAN_TASSERT(getProperty(predMap, 4) == 1)
	SEQAN_TASSERT(getProperty(predMap, 5) == 2)
}

//////////////////////////////////////////////////////////////////////////////

void Test_TopologicalSort() {
//____________________________________________________________________________
// Topological Sort
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 9;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,3, 0,1, 1,2, 3,2, 5,7, 5,6, 6,7, 6,3, 8,7};
	
	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);

	// Predecessor and distance map
	String<TVertexDescriptor> order;
	
	// Topological sort
	topological_sort(g, order);

	SEQAN_TASSERT(getValue(order, 0) == 8)
	SEQAN_TASSERT(getValue(order, 1) == 5)
	SEQAN_TASSERT(getValue(order, 2) == 6)
	SEQAN_TASSERT(getValue(order, 3) == 7)
	SEQAN_TASSERT(getValue(order, 4) == 4)
	SEQAN_TASSERT(getValue(order, 5) == 0)
	SEQAN_TASSERT(getValue(order, 6) == 3)
	SEQAN_TASSERT(getValue(order, 7) == 1)
	SEQAN_TASSERT(getValue(order, 8) == 2)
}

//////////////////////////////////////////////////////////////////////////////

void Test_StronglyConnectedComponents() {
//____________________________________________________________________________
// Strongly-Connected-Components

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 14;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {1,0, 0,4, 2,1, 4,1, 5,1, 6,2, 3,2, 2,3, 7,3, 5,4, 6,5, 5,6, 7,6, 7,7};

	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);

	// Predecessor and distance map
	String<unsigned int> component;
	
	// Strongly Connected Components
	strongly_connected_components(g, component);

	SEQAN_TASSERT(getValue(component, 0) == 3)
	SEQAN_TASSERT(getValue(component, 1) == 3)
	SEQAN_TASSERT(getValue(component, 2) == 2)
	SEQAN_TASSERT(getValue(component, 3) == 2)
	SEQAN_TASSERT(getValue(component, 4) == 3)
	SEQAN_TASSERT(getValue(component, 5) == 1)
	SEQAN_TASSERT(getValue(component, 6) == 1)
	SEQAN_TASSERT(getValue(component, 7) == 0)
}

//////////////////////////////////////////////////////////////////////////////

void Test_PrimsAlgorithm() {
//____________________________________________________________________________
// Prim's algorithm
	typedef Graph<Undirected<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 14;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,6, 1,2, 1,6, 2,3, 2,4, 2,8, 3,5, 3,8, 4,6, 4,7, 5,8, 6,7, 7,8};
	unsigned int weights[] =    {4,   8,   8,   11,  7,   2,   4,   9,   14,  7,   6,   10,  1,   2  };

	//Create the graph
	TGraph g;
	addEdges(g,edges, numEdges);
	String<unsigned int> weightMap;
	resizeEdgeMap(g, weightMap, weights);

	// Tree and predecessor map 
	String<TVertexDescriptor> predMap;

	prims_algorithm(g, 0, weightMap, predMap);

	SEQAN_TASSERT(getProperty(predMap, 0) == getNilPredecessor(g))
	SEQAN_TASSERT(getProperty(predMap, 1) == 0)
	SEQAN_TASSERT(getProperty(predMap, 2) == 1)
	SEQAN_TASSERT(getProperty(predMap, 3) == 2)
	SEQAN_TASSERT(getProperty(predMap, 4) == 2)
	SEQAN_TASSERT(getProperty(predMap, 5) == 3)
	SEQAN_TASSERT(getProperty(predMap, 6) == 7)
	SEQAN_TASSERT(getProperty(predMap, 7) == 8)
	SEQAN_TASSERT(getProperty(predMap, 8) == 2)
}


//////////////////////////////////////////////////////////////////////////////

void Test_KruskalsAlgorithm() {
//____________________________________________________________________________
// Kruskal's algorithm
	typedef Graph<Undirected<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 14;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,6, 1,2, 1,6, 2,3, 2,4, 2,8, 3,5, 3,8, 4,6, 4,7, 5,8, 6,7, 7,8};
	unsigned int weights[] =    {4,   8,   8,   11,  7,   2,   4,   9,   14,  7,   6,   10,  1,   2  };

	//Create the graph
	TGraph g;
	addEdges(g,edges, numEdges);
	String<unsigned int> weightMap;
	resizeEdgeMap(g, weightMap, weights);

	// Tree edges
	String<TVertexDescriptor> treeEdges;
	kruskals_algorithm(g, 0, weightMap, treeEdges);

	SEQAN_TASSERT(getValue(treeEdges, 0) == 6)
	SEQAN_TASSERT(getValue(treeEdges, 1) == 7)
	SEQAN_TASSERT(getValue(treeEdges, 2) == 2)
	SEQAN_TASSERT(getValue(treeEdges, 3) == 4)
	SEQAN_TASSERT(getValue(treeEdges, 4) == 7)
	SEQAN_TASSERT(getValue(treeEdges, 5) == 8)
	SEQAN_TASSERT(getValue(treeEdges, 6) == 0)
	SEQAN_TASSERT(getValue(treeEdges, 7) == 1)
	SEQAN_TASSERT(getValue(treeEdges, 8) == 2)
	SEQAN_TASSERT(getValue(treeEdges, 9) == 8)
	SEQAN_TASSERT(getValue(treeEdges, 10) == 2)
	SEQAN_TASSERT(getValue(treeEdges, 11) == 3)
	SEQAN_TASSERT(getValue(treeEdges, 12) == 0)
	SEQAN_TASSERT(getValue(treeEdges, 13) == 6)
	SEQAN_TASSERT(getValue(treeEdges, 14) == 3)
	SEQAN_TASSERT(getValue(treeEdges, 15) == 5)
}

//////////////////////////////////////////////////////////////////////////////

void Test_DagShortestPath() {
//____________________________________________________________________________
// DAG-Shortest Paths
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,2, 0,1, 1,3, 1,2, 2,5, 2,4, 2,3, 3,5, 3,4, 4,5};
	int weights[] =             {3,   5,   6,   2,   2,   4,   7,   1,   -1,  -2};
	
	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);

	String<int> weightMap;
	resizeEdgeMap(g, weightMap, weights);

	// Predecessor map and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// DAG-Shortest path(Graph, sourceVertex_vertex, weightMap, predMap, distMap)
	dag_shortest_path(g,1,weightMap,predMap,distMap);
	
	SEQAN_TASSERT(getProperty(distMap, 0) == (unsigned int) getInfinityDistance(weightMap))
	SEQAN_TASSERT(getProperty(distMap, 1) == 0)
	SEQAN_TASSERT(getProperty(distMap, 2) == 2)
	SEQAN_TASSERT(getProperty(distMap, 3) == 6)
	SEQAN_TASSERT(getProperty(distMap, 4) == 5)
	SEQAN_TASSERT(getProperty(distMap, 5) == 3)
	SEQAN_TASSERT(getProperty(predMap, 0) == getNilPredecessor(g))
	SEQAN_TASSERT(getProperty(predMap, 1) == getNilPredecessor(g))
	SEQAN_TASSERT(getProperty(predMap, 2) == 1)
	SEQAN_TASSERT(getProperty(predMap, 3) == 1)
	SEQAN_TASSERT(getProperty(predMap, 4) == 3)
	SEQAN_TASSERT(getProperty(predMap, 5) == 4)
}

//////////////////////////////////////////////////////////////////////////////

void Test_BellmanFord() {
//____________________________________________________________________________
// Bellman-Ford

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,3, 1,2, 1,3, 2,4, 3,1, 3,2, 3,4, 4,0, 4,2};
	unsigned int weights[] =    {10,  5,   1,   2,   4,   3,   9,   2,   7,   6};

	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);
	
	String<unsigned int> weightMap;
	resizeEdgeMap(g,weightMap, weights);

	// Out parameters of Bellman-Ford: Predecessor map and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// Bellman-Ford
	bool noNegativeCycle = bellman_ford_algorithm(g,0,weightMap,predMap,distMap);

	SEQAN_TASSERT(getProperty(distMap, 0) == 0)
	SEQAN_TASSERT(getProperty(distMap, 1) == 8)
	SEQAN_TASSERT(getProperty(distMap, 2) == 9)
	SEQAN_TASSERT(getProperty(distMap, 3) == 5)
	SEQAN_TASSERT(getProperty(distMap, 4) == 7)
	SEQAN_TASSERT(getProperty(predMap, 0) == getNilPredecessor(g))
	SEQAN_TASSERT(getProperty(predMap, 1) == 3)
	SEQAN_TASSERT(getProperty(predMap, 2) == 1)
	SEQAN_TASSERT(getProperty(predMap, 3) == 0)
	SEQAN_TASSERT(getProperty(predMap, 4) == 3)
	SEQAN_TASSERT(noNegativeCycle == true)
}

//////////////////////////////////////////////////////////////////////////////

void Test_Dijkstra() {
//____________________________________________________________________________
// Dijkstra

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,3, 1,2, 1,3, 2,4, 3,1, 3,2, 3,4, 4,0, 4,2};
	unsigned int weights[] =    {10,  5,   1,   2,   4,   3,   9,   2,   7,   6};

	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);

	String<unsigned int> weightMap;
	resizeEdgeMap(g , weightMap, weights);

	// Out parameters of Dijkstra: Predecessor map and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// Dijkstra
	dijkstra(g,0,weightMap,predMap,distMap);

	SEQAN_TASSERT(getProperty(distMap, 0) == 0)
	SEQAN_TASSERT(getProperty(distMap, 1) == 8)
	SEQAN_TASSERT(getProperty(distMap, 2) == 9)
	SEQAN_TASSERT(getProperty(distMap, 3) == 5)
	SEQAN_TASSERT(getProperty(distMap, 4) == 7)
	SEQAN_TASSERT(getProperty(predMap, 0) == getNilPredecessor(g))
	SEQAN_TASSERT(getProperty(predMap, 1) == 3)
	SEQAN_TASSERT(getProperty(predMap, 2) == 1)
	SEQAN_TASSERT(getProperty(predMap, 3) == 0)
	SEQAN_TASSERT(getProperty(predMap, 4) == 3)
}

//////////////////////////////////////////////////////////////////////////////

void Test_AllPairsShortestPath() {
//____________________________________________________________________________
// All-Pairs Shortest Path

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 9;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,2, 0,4, 1,3, 1,4, 2,1, 3,0, 3,2, 4,3};
	int weights[] =    {3,   8,   -4,  1,   7,   4,   2,   -5,  6};

	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);
	
	String<int> weightMap;	
	resizeEdgeMap(g,weightMap, weights);

	// Out parameter
	Matrix<int> distMat;
	Matrix<TVertexDescriptor> predMat;

	// All-Pairs shortest path
	all_pairs_shortest_path(g,weightMap, distMat, predMat);

	unsigned int len = length(distMat, 0);
	SEQAN_TASSERT(getValue(distMat, 0*len + 0) == 0)
	SEQAN_TASSERT(getValue(distMat, 0*len + 1) == 1)
	SEQAN_TASSERT(getValue(distMat, 0*len + 2) == -3)
	SEQAN_TASSERT(getValue(distMat, 0*len + 3) == 2)
	SEQAN_TASSERT(getValue(distMat, 0*len + 4) == -4)
	SEQAN_TASSERT(getValue(distMat, 1*len + 0) == 3)
	SEQAN_TASSERT(getValue(distMat, 1*len + 1) == 0)
	SEQAN_TASSERT(getValue(distMat, 1*len + 2) == -4)
	SEQAN_TASSERT(getValue(distMat, 1*len + 3) == 1)
	SEQAN_TASSERT(getValue(distMat, 1*len + 4) == -1)
	SEQAN_TASSERT(getValue(distMat, 2*len + 0) == 7)
	SEQAN_TASSERT(getValue(distMat, 2*len + 1) == 4)
	SEQAN_TASSERT(getValue(distMat, 2*len + 2) == 0)
	SEQAN_TASSERT(getValue(distMat, 2*len + 3) == 5)
	SEQAN_TASSERT(getValue(distMat, 2*len + 4) == 3)
	SEQAN_TASSERT(getValue(distMat, 3*len + 0) == 2)
	SEQAN_TASSERT(getValue(distMat, 3*len + 1) == -1)
	SEQAN_TASSERT(getValue(distMat, 3*len + 2) == -5)
	SEQAN_TASSERT(getValue(distMat, 3*len + 3) == 0)
	SEQAN_TASSERT(getValue(distMat, 3*len + 4) == -2)
	SEQAN_TASSERT(getValue(distMat, 4*len + 0) == 8)
	SEQAN_TASSERT(getValue(distMat, 4*len + 1) == 5)
	SEQAN_TASSERT(getValue(distMat, 4*len + 2) == 1)
	SEQAN_TASSERT(getValue(distMat, 4*len + 3) == 6)
	SEQAN_TASSERT(getValue(distMat, 4*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 0*len + 0) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 0*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 0*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 0*len + 3) == 4)
	SEQAN_TASSERT(getValue(predMat, 0*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 1*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 1*len + 1) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 1*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 1*len + 3) == 1)
	SEQAN_TASSERT(getValue(predMat, 1*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 2*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 2*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 2*len + 2) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 2*len + 3) == 1)
	SEQAN_TASSERT(getValue(predMat, 2*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 3*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 3*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 3*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 3*len + 3) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 3*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 4*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 4*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 4*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 4*len + 3) == 4)
	SEQAN_TASSERT(getValue(predMat, 4*len + 4) == getNilPredecessor(g))
}

//////////////////////////////////////////////////////////////////////////////

void Test_FloydWarshall() {
//____________________________________________________________________________
// Floyd Warshall

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 9;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,2, 0,4, 1,3, 1,4, 2,1, 3,0, 3,2, 4,3};
	int weights[] =    {3,   8,   -4,  1,   7,   4,   2,   -5,  6};

	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);
	
	String<int> weightMap;	
	resizeEdgeMap(g,weightMap, weights);

	// Out parameter
	Matrix<int> distMat;
	Matrix<TVertexDescriptor> predMat;

	// Floyd-Warshall
	floyd_warshall(g,weightMap, distMat, predMat);

	unsigned int len = length(distMat, 0);
	SEQAN_TASSERT(getValue(distMat, 0*len + 0) == 0)
	SEQAN_TASSERT(getValue(distMat, 0*len + 1) == 1)
	SEQAN_TASSERT(getValue(distMat, 0*len + 2) == -3)
	SEQAN_TASSERT(getValue(distMat, 0*len + 3) == 2)
	SEQAN_TASSERT(getValue(distMat, 0*len + 4) == -4)
	SEQAN_TASSERT(getValue(distMat, 1*len + 0) == 3)
	SEQAN_TASSERT(getValue(distMat, 1*len + 1) == 0)
	SEQAN_TASSERT(getValue(distMat, 1*len + 2) == -4)
	SEQAN_TASSERT(getValue(distMat, 1*len + 3) == 1)
	SEQAN_TASSERT(getValue(distMat, 1*len + 4) == -1)
	SEQAN_TASSERT(getValue(distMat, 2*len + 0) == 7)
	SEQAN_TASSERT(getValue(distMat, 2*len + 1) == 4)
	SEQAN_TASSERT(getValue(distMat, 2*len + 2) == 0)
	SEQAN_TASSERT(getValue(distMat, 2*len + 3) == 5)
	SEQAN_TASSERT(getValue(distMat, 2*len + 4) == 3)
	SEQAN_TASSERT(getValue(distMat, 3*len + 0) == 2)
	SEQAN_TASSERT(getValue(distMat, 3*len + 1) == -1)
	SEQAN_TASSERT(getValue(distMat, 3*len + 2) == -5)
	SEQAN_TASSERT(getValue(distMat, 3*len + 3) == 0)
	SEQAN_TASSERT(getValue(distMat, 3*len + 4) == -2)
	SEQAN_TASSERT(getValue(distMat, 4*len + 0) == 8)
	SEQAN_TASSERT(getValue(distMat, 4*len + 1) == 5)
	SEQAN_TASSERT(getValue(distMat, 4*len + 2) == 1)
	SEQAN_TASSERT(getValue(distMat, 4*len + 3) == 6)
	SEQAN_TASSERT(getValue(distMat, 4*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 0*len + 0) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 0*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 0*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 0*len + 3) == 4)
	SEQAN_TASSERT(getValue(predMat, 0*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 1*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 1*len + 1) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 1*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 1*len + 3) == 1)
	SEQAN_TASSERT(getValue(predMat, 1*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 2*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 2*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 2*len + 2) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 2*len + 3) == 1)
	SEQAN_TASSERT(getValue(predMat, 2*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 3*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 3*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 3*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 3*len + 3) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 3*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 4*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 4*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 4*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 4*len + 3) == 4)
	SEQAN_TASSERT(getValue(predMat, 4*len + 4) == getNilPredecessor(g))
}


//////////////////////////////////////////////////////////////////////////////

void Test_TransitiveClosure() {
//____________________________________________________________________________
// Transitive Closure

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 5;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {3,0, 1,2, 2,1, 1,3, 3,2};
		
	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);

	// Transitive-Closure
	Matrix<bool> closure;
	transitive_closure(g,closure);
	
	unsigned int len = length(closure, 0);
	SEQAN_TASSERT(getValue(closure, 0*len + 0) == 1)
	SEQAN_TASSERT(getValue(closure, 0*len + 1) == 0)
	SEQAN_TASSERT(getValue(closure, 0*len + 2) == 0)
	SEQAN_TASSERT(getValue(closure, 0*len + 3) == 0)
	SEQAN_TASSERT(getValue(closure, 1*len + 0) == 1)
	SEQAN_TASSERT(getValue(closure, 1*len + 1) == 1)
	SEQAN_TASSERT(getValue(closure, 1*len + 2) == 1)
	SEQAN_TASSERT(getValue(closure, 1*len + 3) == 1)
	SEQAN_TASSERT(getValue(closure, 2*len + 0) == 1)
	SEQAN_TASSERT(getValue(closure, 2*len + 1) == 1)
	SEQAN_TASSERT(getValue(closure, 2*len + 2) == 1)
	SEQAN_TASSERT(getValue(closure, 2*len + 3) == 1)
	SEQAN_TASSERT(getValue(closure, 3*len + 0) == 1)
	SEQAN_TASSERT(getValue(closure, 3*len + 1) == 1)
	SEQAN_TASSERT(getValue(closure, 3*len + 2) == 1)
	SEQAN_TASSERT(getValue(closure, 3*len + 3) == 1)
}

//////////////////////////////////////////////////////////////////////////////

void Test_FordFulkerson() {
//____________________________________________________________________________
// Ford-Fulkerson
	typedef Graph<Directed<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,4, 1,2, 1,4, 2,3, 2,4, 4,1, 4,5, 5,2, 5,3};
	unsigned int capacity[] =    {16,  13,  12,  10,  20,  9,   4,   14,  7,   4};

	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);
	String<unsigned int> capMap;	
	resizeEdgeMap(g,capMap, capacity);

	// Out-parameter
	String<unsigned int> flow;	
	unsigned int valF = ford_fulkerson(g, 0, 3, capMap, flow);
	
	SEQAN_TASSERT(valF == 23)
	TEdgeIterator itEdge(g);
	for(;!atEnd(itEdge);goNext(itEdge)) {
		SEQAN_TASSERT(getProperty(flow, getValue(itEdge)) <= getProperty(capMap, getValue(itEdge)))
	}
}

//////////////////////////////////////////////////////////////////////////////

void Test_Algorithms() {
//____________________________________________________________________________
// Graph Algorithms
	// Elementary graph algorithms
	Test_BreadthFirstSearch();
	Test_DepthFirstSearch();
	Test_TopologicalSort();
	Test_StronglyConnectedComponents();

	// Minimum Spanning Trees
	Test_PrimsAlgorithm();
	Test_KruskalsAlgorithm();

	// Single-Source shortest paths
	Test_DagShortestPath();
	Test_BellmanFord();
	Test_Dijkstra();

	// All-Pairs Shortest paths
	Test_AllPairsShortestPath();
	Test_FloydWarshall();
	Test_TransitiveClosure();

	//Maximum Flow
	Test_FordFulkerson();

	//Todo
	//Matching
}



//////////////////////////////////////////////////////////////////////////////

void Test_TCoffee() {
//____________________________________________________________________________
// Graph TCoffee

	// Test aa groups
	AAGroups gr;
	gr = AminoAcid('T');
	SEQAN_TASSERT(gr == 0)
	gr = Byte(3);
	SEQAN_TASSERT(gr == 2)
	gr = char('c');
	SEQAN_TASSERT(gr == 5)
	gr = Unicode('j');
	SEQAN_TASSERT(gr == 6)

	// Read a t-coffee library: AminoAcid Alphabet
	typedef StringSet<String<AminoAcid>, IdHolder<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int, Default> > TGraph;
	TStringSet strSet;
	TGraph g(strSet);

	fstream strm; // Read the library
	strm.open(TEST_PATH "garfield.lib", ios_base::in);
	read(strm,g,TCoffeeLib());
	strm.close();

	fstream strmW; // Write the library
	strmW.open(TEST_PATH "my_garfield.lib", ios_base::out | ios_base::trunc);
	write(strmW,g,TCoffeeLib());
	strmW.close();

	//std::cout << g << std::endl;

	fstream strm2; // Alignment graph as dot
	strm2.open(TEST_PATH "my_tcoffee.dot", ios_base::out | ios_base::trunc);
	write(strm2,g,DotDrawing());
	strm2.close();

	Matrix<double> score;  // Calculate a distance matrix on an amino acid string set
	getScoringMatrix(stringSet(g), score);
	Matrix<double> distances;
	scoreToDistanceMatrix(score, distances, 1000);
	normalizeMatrix(distances, 100000, 100);
	Graph<Tree<double> > njTreeOut;
	slowNjTree(distances, njTreeOut);
	std::cout << njTreeOut << std::endl;

	// ToDo: Owner of strings!!!!!!!!!!!
	for(unsigned int i=0; i<length(value(g.data_sequence)); ++i) {  	// Delete sequences
		delete &value(g.data_sequence)[i];
	}

	// Read a t-coffee library: Dna Alphabet
	typedef StringSet<String<Dna>, IdHolder<> > TStringSetDna;
	typedef Graph<Alignment<TStringSetDna, unsigned int, Default> > TGraphDna;
	TStringSetDna strSetDna;
	TGraphDna gDna(strSetDna);

	fstream strmDna; // Read the library
	strmDna.open(TEST_PATH "dna_seq.lib", ios_base::in);
	read(strmDna,gDna,TCoffeeLib());
	strmDna.close();

	fstream strmWDna; // Write the library
	strmWDna.open(TEST_PATH "my_dna_seq.lib", ios_base::out | ios_base::trunc);
	write(strmWDna,gDna,TCoffeeLib());
	strmWDna.close();

	//std::cout << g << std::endl;

	Matrix<double> scoreDna;  // Calculate a distance matrix on an amino acid string set
	getScoringMatrix(stringSet(gDna), scoreDna);
	Matrix<double> distancesDna;
	scoreToDistanceMatrix(scoreDna, distancesDna, 1000);
	normalizeMatrix(distancesDna, 100000, 100);
	Graph<Tree<double> > njTree;
	slowNjTree(distancesDna, njTree);
	std::cout << njTree << std::endl;

	// ToDo: Owner of strings!!!!!!!!!!!
	for(unsigned int i=0; i<length(value(gDna.data_sequence)); ++i) {  	// Delete sequences
		delete &value(gDna.data_sequence)[i];
	}
}

//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	Test_IdManager();	// Test Id Manager
	Test_EdgeStump();	// Test EdgeStumps
	Test_StringSet<StringSet<String<char>, IdHolder<GenerousStorage<> > > >();
	Test_StringSet<StringSet<String<char>, IdHolder<TightStorage<> > > >();

	// Test Graph types
	Test_Directed();	// Directed graphs
	Test_Undirected();	// Undirected graphs
	Test_Automaton();	// Automatons
	Test_WordGraph();	// Word Graph
	Test_Tree();		// Trees
	Test_Alignment();	// Alignment graph

	// Test iterators
	Test_VertexIterator<Directed<char> >();
	Test_VertexIterator<Undirected<char> >();
	Test_VertexIterator<Automaton<char> >();
	Test_TreeInternalVertexIterator();
	Test_OutEdgeIterator<Directed<char> >();
	Test_OutEdgeIterator<Undirected<char> >();
	Test_OutEdgeIterator<Tree<char> >();
	Test_OutEdgeIterator<Automaton<char> >();
	Test_EdgeIterator<Directed<char> >();
	Test_EdgeIterator<Undirected<char> >();
	Test_EdgeIterator<Tree<char> >();
	Test_EdgeIterator<Automaton<char> >();
	Test_AdjacencyIterator<Directed<char> >();
	Test_AdjacencyIterator<Undirected<char> >();
	Test_AdjacencyIterator<Tree<char> >();
	Test_AdjacencyIterator<Automaton<char> >();

	// Test property maps
	Test_ExternalProperty<Directed<char> >();
	Test_ExternalProperty<Undirected<char> >();
	Test_ExternalProperty<Tree<char> >();
	Test_ExternalProperty<Automaton<char> >();	
	Test_Property();


	// Specialized automatons
	Test_Oracle();
	Test_Trie();

	// Test bfs and dfs iterator
	Test_BfsIter<Directed<char> >();
	Test_BfsIter<Undirected<char> >();
	Test_BfsIter<Tree<char> >();
	Test_BfsIter<Automaton<char> >();
	Test_BfsIterator();
	Test_DfsPreorderIter<Directed<char> >();
	Test_DfsPreorderIter<Undirected<char> >();
	Test_DfsPreorderIter<Tree<char> >();
	Test_DfsPreorderIter<Automaton<char> >();
	Test_DfsPreorderIterator();

	// Graph algorithms
	Test_Algorithms();


	// T-Coffee
	Test_TCoffee();


//____________________________________________________________________________
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_base.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_idmanager.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_interface.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_edgestump.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_impl_directed.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_impl_undirected.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_impl_automaton.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_impl_wordgraph.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_impl_tree.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_impl_align.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_iterator.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_iterator_vertex.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_iterator_outedge.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_iterator_adjacency.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_iterator_edge.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_property.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_impl_oracle.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_impl_trie.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_drawing.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_iterator_bfs.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_iterator_dfs.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_algorithm.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_algorithm_tcoffee.h");


	SEQAN_TREPORT("TEST END")

	return 0;
}
