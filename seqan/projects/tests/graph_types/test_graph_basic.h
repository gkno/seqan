#ifndef SEQAN_HEADER_TEST_GRAPH_BASIC_H
#define SEQAN_HEADER_TEST_GRAPH_BASIC_H

using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

void Test_IdManager() {
//____________________________________________________________________________
// IdManager
	typedef Id<IdManager<> >::Type TIdType;
	IdManager<> idm;
	
	// Obtain Ids
	TIdType id=obtainId(idm);
	SEQAN_ASSERT_TRUE(id == 0);
	SEQAN_ASSERT_TRUE(idInUse(idm,0) == true);
	id=obtainId(idm);
	SEQAN_ASSERT_TRUE(id == 1);
	id=obtainId(idm);
	SEQAN_ASSERT_TRUE(id == 2);
	id=obtainId(idm);
	SEQAN_ASSERT_TRUE(id == 3);
	id=obtainId(idm);
	SEQAN_ASSERT_TRUE(id == 4);
	SEQAN_ASSERT_TRUE(idInUse(idm,4) == true);
	SEQAN_ASSERT_TRUE(getIdUpperBound(idm) == 5);
	SEQAN_ASSERT_TRUE(idCount(idm) == 5);

	// Release Ids
	releaseId(idm,3);
	SEQAN_ASSERT_TRUE(idInUse(idm,3) == false);
	releaseId(idm,1);
	SEQAN_ASSERT_TRUE(idInUse(idm,1) == false);
	releaseId(idm,2);
	SEQAN_ASSERT_TRUE(idInUse(idm,2) == false);
	SEQAN_ASSERT_TRUE(idCount(idm) == 2);
	SEQAN_ASSERT_TRUE(getIdUpperBound(idm) == 5);
	releaseId(idm,4);  // Now we can shrink id range
	SEQAN_ASSERT_TRUE(idCount(idm) == 1);
	SEQAN_ASSERT_TRUE(getIdUpperBound(idm) == 4);

	// Ids are reused
	id=obtainId(idm);
	id=obtainId(idm);
	id=obtainId(idm);
	id=obtainId(idm);
	SEQAN_ASSERT_TRUE(getIdUpperBound(idm) == 5);
	releaseId(idm,3);
	SEQAN_ASSERT_TRUE(getIdLowerBound(idm) == 0);
	releaseId(idm,0);
	SEQAN_ASSERT_TRUE(idInUse(idm,0) == false);
	SEQAN_ASSERT_TRUE(getIdLowerBound(idm) == 1);
	releaseId(idm,2);
	SEQAN_ASSERT_TRUE(idInUse(idm,2) == false);
	id=obtainId(idm);
	SEQAN_ASSERT_TRUE(id == 2);
	SEQAN_ASSERT_TRUE(idInUse(idm,2) == true);
	
	// Check copy constructor and assignment operator
	IdManager<> idm2(idm);
	SEQAN_ASSERT_TRUE(idCount(idm2) == 3);
	releaseAll(idm2);
	SEQAN_ASSERT_TRUE(idCount(idm2) == 0);
	SEQAN_ASSERT_TRUE(getIdUpperBound(idm2) == 0);
	SEQAN_ASSERT_TRUE(getIdLowerBound(idm2) == 0);

	// Check assignment
	idm2 = idm;
	SEQAN_ASSERT_TRUE(idCount(idm2) == 3);


//____________________________________________________________________________
// Dummy IdManager

	// Dummy IdManager
	typedef Id<IdManager<void> >::Type TIdType;
	IdManager<void> id_dummy;
	
	// Obtain Ids
	TIdType idd=obtainId(id_dummy);
	SEQAN_ASSERT_TRUE(idd == 0);
	idd=obtainId(id_dummy); // Always zero
	SEQAN_ASSERT_TRUE(idd == 0);
	SEQAN_ASSERT_TRUE(idInUse(id_dummy, 1) == false); // Always false
	SEQAN_ASSERT_TRUE(idInUse(id_dummy, 2) == false);
	idd=obtainId(id_dummy);
	SEQAN_ASSERT_TRUE(idd == 0);
	idd=obtainId(id_dummy);
	SEQAN_ASSERT_TRUE(idd == 0);
	idd=obtainId(id_dummy);
	SEQAN_ASSERT_TRUE(idd == 0);
	SEQAN_ASSERT_TRUE(getIdUpperBound(id_dummy) == 5); 
	SEQAN_ASSERT_TRUE(getIdLowerBound(id_dummy) == 0);
	SEQAN_ASSERT_TRUE(idCount(id_dummy) == 5); //But: Correct id count

	// Release Ids
	releaseId(id_dummy,3); 
	SEQAN_ASSERT_TRUE(idCount(id_dummy) == 4);
	releaseId(id_dummy,1); 
	SEQAN_ASSERT_TRUE(idCount(id_dummy) == 3);
	IdManager<void> id_dummy2(id_dummy);
	SEQAN_ASSERT_TRUE(idCount(id_dummy2) == 3);
	idd=obtainId(id_dummy2);
	SEQAN_ASSERT_TRUE(idd == 0);
	id_dummy = id_dummy2;
	SEQAN_ASSERT_TRUE(idCount(id_dummy) == 4);
	releaseAll(id_dummy);
	SEQAN_ASSERT_TRUE(idCount(id_dummy) == 0);
	SEQAN_ASSERT_TRUE(getIdUpperBound(id_dummy) == 0);
	SEQAN_ASSERT_TRUE(getIdLowerBound(id_dummy) == 0);
}

//////////////////////////////////////////////////////////////////////////////

void Test_EdgeStump() {
//____________________________________________________________________________
// Test all cargoless EdgeStumps in a list
	// No cargo, list, source, id
	EdgeStump<void, true, true, true> es1;
	_assignId(&es1, 4);
	SEQAN_ASSERT_TRUE(_getId(&es1) == 4);
	assignTarget(&es1, 5);
	SEQAN_ASSERT_TRUE(getTarget(&es1) == 5);
	target(&es1)=7;
	SEQAN_ASSERT_TRUE(getTarget(&es1) == 7);
	assignSource(&es1, 20);
	SEQAN_ASSERT_TRUE(getSource(&es1) == 20);
	source(&es1)=30;
	SEQAN_ASSERT_TRUE(getSource(&es1) == 30);
	assignCargo(&es1, 15);  //Does nothing, no cargo -> pointer to 0
	SEQAN_ASSERT_TRUE(getCargo(&es1) == (void*) 0);
	SEQAN_ASSERT_TRUE(cargo(&es1) == (void*) 0);
	SEQAN_ASSERT_TRUE(_getId(&es1) == 4);
	EdgeStump<void, true, true, true> const es1_const(es1);
	SEQAN_ASSERT_TRUE(getCargo(&es1_const) == (void*) 0);
	SEQAN_ASSERT_TRUE(cargo(&es1_const) == (void*) 0);
	SEQAN_ASSERT_TRUE(getTarget(&es1_const) == 7);
	SEQAN_ASSERT_TRUE(target(&es1_const) == 7);
	SEQAN_ASSERT_TRUE(source(&es1_const) == 30);
	SEQAN_ASSERT_TRUE(getSource(&es1_const) == 30);
	SEQAN_ASSERT_TRUE(_getId(&es1_const) == 4);
	EdgeStump<void, true, true, true> es11(es1);
	EdgeStump<void, true, true, true> es12(es1);
	nextT(&es1) = &es11;
	SEQAN_ASSERT_TRUE(getNextT(&es1) == &es11);
	assignNextT(&es1, &es12);
	SEQAN_ASSERT_TRUE(getNextT(&es1) == &es12);
	nextS(&es1) = &es11;
	SEQAN_ASSERT_TRUE(getNextS(&es1) == &es11);
	SEQAN_ASSERT_TRUE(getNextT(&es1) == &es12);
	assignNextS(&es1, &es12);
	SEQAN_ASSERT_TRUE(getNextS(&es1) == &es12);

	// No cargo, list, source, no id
	EdgeStump<void, true, true, false> es2;
	_assignId(&es2, 4);
	SEQAN_ASSERT_TRUE(_getId(&es2) == 0); // No id, always 0
	assignTarget(&es2, 5);
	SEQAN_ASSERT_TRUE(getTarget(&es2) == 5);
	target(&es2)=7;
	SEQAN_ASSERT_TRUE(getTarget(&es2) == 7);
	assignSource(&es2, 20);
	SEQAN_ASSERT_TRUE(getSource(&es2) == 20);
	source(&es2)=30;
	SEQAN_ASSERT_TRUE(getSource(&es2) == 30);
	assignCargo(&es2, 15);  //Does nothing, no cargo -> pointer to 0
	SEQAN_ASSERT_TRUE(getCargo(&es2) == (void*) 0);
	SEQAN_ASSERT_TRUE(cargo(&es2) == (void*) 0);
	SEQAN_ASSERT_TRUE(_getId(&es2) == 0);
	EdgeStump<void, true, true, false> const es2_const(es2);
	SEQAN_ASSERT_TRUE(getCargo(&es2_const) == (void*) 0);
	SEQAN_ASSERT_TRUE(cargo(&es2_const) == (void*) 0);
	SEQAN_ASSERT_TRUE(getTarget(&es2_const) == 7);
	SEQAN_ASSERT_TRUE(target(&es2_const) == 7);
	SEQAN_ASSERT_TRUE(source(&es2_const) == 30);
	SEQAN_ASSERT_TRUE(getSource(&es2_const) == 30);
	SEQAN_ASSERT_TRUE(_getId(&es2_const) == 0);
	EdgeStump<void, true, true, false> es21(es2);
	EdgeStump<void, true, true, false> es22(es2);
	nextT(&es2) = &es21;
	SEQAN_ASSERT_TRUE(getNextT(&es2) == &es21);
	assignNextT(&es2, &es22);
	SEQAN_ASSERT_TRUE(getNextT(&es2) == &es22);
	nextS(&es2) = &es21;
	SEQAN_ASSERT_TRUE(getNextS(&es2) == &es21);
	SEQAN_ASSERT_TRUE(getNextT(&es2) == &es22);
	assignNextS(&es2, &es22);
	SEQAN_ASSERT_TRUE(getNextS(&es2) == &es22);

	// No cargo, list, no source, id
	EdgeStump<void, true, false, true> es3;
	_assignId(&es3, 4);
	SEQAN_ASSERT_TRUE(_getId(&es3) == 4);
	assignTarget(&es3, 5);
	SEQAN_ASSERT_TRUE(getTarget(&es3) == 5);
	target(&es3)=7;
	SEQAN_ASSERT_TRUE(getTarget(&es3) == 7);
	assignSource(&es3, 20);
	SEQAN_ASSERT_TRUE(getSource(&es3) == 0);  // No source
	SEQAN_ASSERT_TRUE(source(&es3) == 0); // No source
	assignCargo(&es3, 15);  //Does nothing, no cargo -> pointer to 0
	SEQAN_ASSERT_TRUE(getCargo(&es3) == (void*) 0);
	SEQAN_ASSERT_TRUE(cargo(&es3) == (void*) 0);
	SEQAN_ASSERT_TRUE(_getId(&es3) == 4);
	EdgeStump<void, true, false, true> const es3_const(es3);
	SEQAN_ASSERT_TRUE(getCargo(&es3_const) == (void*) 0);
	SEQAN_ASSERT_TRUE(cargo(&es3_const) == (void*) 0);
	SEQAN_ASSERT_TRUE(getTarget(&es3_const) == 7);
	SEQAN_ASSERT_TRUE(target(&es3_const) == 7);
	SEQAN_ASSERT_TRUE(getSource(&es3_const) == 0);
	SEQAN_ASSERT_TRUE(source(&es3_const) == 0);
	SEQAN_ASSERT_TRUE(_getId(&es3_const) == 4);
	EdgeStump<void, true, false, true> es31(es3);
	EdgeStump<void, true, false, true> es32(es3);
	nextT(&es3) = &es31;
	SEQAN_ASSERT_TRUE(getNextT(&es3) == &es31);
	assignNextT(&es3, &es32);
	SEQAN_ASSERT_TRUE(getNextT(&es3) == &es32);
	assignNextS(&es3, &es32); // No source
	SEQAN_ASSERT_TRUE(nextS(&es3) == 0);
	SEQAN_ASSERT_TRUE(getNextS(&es3) == 0);
	SEQAN_ASSERT_TRUE(getNextT(&es3) == &es32);

	// No cargo, list, no source, no id
	EdgeStump<void, true, false, false> es4;
	_assignId(&es4, 4);
	SEQAN_ASSERT_TRUE(_getId(&es4) == 0); // No id
	assignTarget(&es4, 5);
	SEQAN_ASSERT_TRUE(getTarget(&es4) == 5);
	target(&es4)=7;
	SEQAN_ASSERT_TRUE(getTarget(&es4) == 7);
	assignSource(&es4, 20);
	SEQAN_ASSERT_TRUE(getSource(&es4) == 0); // No source
	SEQAN_ASSERT_TRUE(source(&es4) == 0);
	assignCargo(&es4, 15);  //Does nothing, no cargo -> pointer to 0
	SEQAN_ASSERT_TRUE(getCargo(&es4) == (void*) 0);
	SEQAN_ASSERT_TRUE(cargo(&es4) == (void*) 0);
	SEQAN_ASSERT_TRUE(_getId(&es4) == 0);
	EdgeStump<void, true, false, false> const es4_const(es4);
	SEQAN_ASSERT_TRUE(getCargo(&es4_const) == (void*) 0);
	SEQAN_ASSERT_TRUE(cargo(&es4_const) == (void*) 0);
	SEQAN_ASSERT_TRUE(getTarget(&es4_const) == 7);
	SEQAN_ASSERT_TRUE(target(&es4_const) == 7);
	SEQAN_ASSERT_TRUE(getSource(&es4_const) == 0);
	SEQAN_ASSERT_TRUE(source(&es4_const) == 0);
	SEQAN_ASSERT_TRUE(_getId(&es4_const) == 0);
	EdgeStump<void, true, false, false> es41(es4);
	EdgeStump<void, true, false, false> es42(es4);
	nextT(&es4) = &es41;
	SEQAN_ASSERT_TRUE(getNextT(&es4) == &es41);
	assignNextT(&es4, &es42);
	SEQAN_ASSERT_TRUE(getNextT(&es4) == &es42);
	assignNextS(&es4, &es42);
	SEQAN_ASSERT_TRUE(getNextS(&es4) == 0);
	SEQAN_ASSERT_TRUE(nextS(&es4) == 0);

//____________________________________________________________________________
// Test all EdgeStumps with cargo in a list
	// Cargo, list, source, id
	EdgeStump<unsigned int, true, true, true> es5;
	_assignId(&es5, 4);
	SEQAN_ASSERT_TRUE(_getId(&es5) == 4);
	assignTarget(&es5, 5);
	SEQAN_ASSERT_TRUE(getTarget(&es5) == 5);
	target(&es5)=7;
	SEQAN_ASSERT_TRUE(getTarget(&es5) == 7);
	assignSource(&es5, 20);
	SEQAN_ASSERT_TRUE(getSource(&es5) == 20);
	source(&es5)=30;
	SEQAN_ASSERT_TRUE(getSource(&es5) == 30);
	assignCargo(&es5, 15); 
	SEQAN_ASSERT_TRUE(getCargo(&es5) == 15);
	SEQAN_ASSERT_TRUE(cargo(&es5) == 15);
	SEQAN_ASSERT_TRUE(_getId(&es5) == 4);
	EdgeStump<unsigned int, true, true, true> const es5_const(es5);
	SEQAN_ASSERT_TRUE(getCargo(&es5_const) == 15);
	SEQAN_ASSERT_TRUE(cargo(&es5_const) == 15);
	SEQAN_ASSERT_TRUE(getTarget(&es5_const) == 7);
	SEQAN_ASSERT_TRUE(target(&es5_const) == 7);
	EdgeStump<unsigned int, true, true, true> es51(es5);
	EdgeStump<unsigned int, true, true, true> es52(es5);
	nextT(&es5) = &es51;
	SEQAN_ASSERT_TRUE(getNextT(&es5) == &es51);
	assignNextT(&es5, &es52);
	SEQAN_ASSERT_TRUE(getNextT(&es5) == &es52);
	nextS(&es5) = &es51;
	SEQAN_ASSERT_TRUE(getNextS(&es5) == &es51);
	SEQAN_ASSERT_TRUE(getNextT(&es5) == &es52);
	assignNextS(&es5, &es52);
	SEQAN_ASSERT_TRUE(getNextS(&es5) == &es52);

	// Cargo, list, source, no id
	EdgeStump<unsigned int, true, true, false> es6;
	_assignId(&es6, 4);
	SEQAN_ASSERT_TRUE(_getId(&es6) == 0); // No id, always 0
	assignTarget(&es6, 5);
	SEQAN_ASSERT_TRUE(getTarget(&es6) == 5);
	target(&es6)=7;
	SEQAN_ASSERT_TRUE(getTarget(&es6) == 7);
	assignSource(&es6, 20);
	SEQAN_ASSERT_TRUE(getSource(&es6) == 20);
	source(&es6)=30;
	SEQAN_ASSERT_TRUE(getSource(&es6) == 30);
	assignCargo(&es6, 15); 
	SEQAN_ASSERT_TRUE(getCargo(&es6) == 15);
	SEQAN_ASSERT_TRUE(cargo(&es6) == 15);
	SEQAN_ASSERT_TRUE(_getId(&es6) == 0);
	EdgeStump<unsigned int, true, true, false> const es6_const(es6);
	SEQAN_ASSERT_TRUE(getCargo(&es6_const) == 15);
	SEQAN_ASSERT_TRUE(cargo(&es6_const) == 15);
	SEQAN_ASSERT_TRUE(getTarget(&es6_const) == 7);
	SEQAN_ASSERT_TRUE(target(&es6_const) == 7);
	EdgeStump<unsigned int, true, true, false> es61(es6);
	EdgeStump<unsigned int, true, true, false> es62(es6);
	nextT(&es6) = &es61;
	SEQAN_ASSERT_TRUE(getNextT(&es6) == &es61);
	assignNextT(&es6, &es62);
	SEQAN_ASSERT_TRUE(getNextT(&es6) == &es62);
	nextS(&es6) = &es61;
	SEQAN_ASSERT_TRUE(getNextS(&es6) == &es61);
	SEQAN_ASSERT_TRUE(getNextT(&es6) == &es62);
	assignNextS(&es6, &es62);
	SEQAN_ASSERT_TRUE(getNextS(&es6) == &es62);

	// Cargo, list, no source, id
	EdgeStump<unsigned int, true, false, true> es7;
	_assignId(&es7, 4);
	SEQAN_ASSERT_TRUE(_getId(&es7) == 4);
	assignTarget(&es7, 5);
	SEQAN_ASSERT_TRUE(getTarget(&es7) == 5);
	target(&es7)=7;
	SEQAN_ASSERT_TRUE(getTarget(&es7) == 7);
	assignSource(&es7, 20);
	SEQAN_ASSERT_TRUE(getSource(&es7) == 0);  // No source
	SEQAN_ASSERT_TRUE(source(&es7) == 0); // No source
	assignCargo(&es7, 15); 
	SEQAN_ASSERT_TRUE(getCargo(&es7) == 15);
	SEQAN_ASSERT_TRUE(cargo(&es7) == 15);
	SEQAN_ASSERT_TRUE(_getId(&es7) == 4);
	EdgeStump<unsigned int, true, false, true> const es7_const(es7);
	SEQAN_ASSERT_TRUE(getCargo(&es7_const) == 15);
	SEQAN_ASSERT_TRUE(cargo(&es7_const) == 15);
	SEQAN_ASSERT_TRUE(getTarget(&es7_const) == 7);
	SEQAN_ASSERT_TRUE(target(&es7_const) == 7);
	EdgeStump<unsigned int, true, false, true> es71(es7);
	EdgeStump<unsigned int, true, false, true> es72(es7);
	nextT(&es7) = &es71;
	SEQAN_ASSERT_TRUE(getNextT(&es7) == &es71);
	assignNextT(&es7, &es72);
	SEQAN_ASSERT_TRUE(getNextT(&es7) == &es72);
	assignNextS(&es7, &es72); // No source
	SEQAN_ASSERT_TRUE(nextS(&es7) == 0);
	SEQAN_ASSERT_TRUE(getNextS(&es7) == 0);
	SEQAN_ASSERT_TRUE(getNextT(&es7) == &es72);

	// Cargo, list, no source, no id
	EdgeStump<unsigned int, true, false, false> es8;
	_assignId(&es8, 4);
	SEQAN_ASSERT_TRUE(_getId(&es8) == 0); // No id
	assignTarget(&es8, 5);
	SEQAN_ASSERT_TRUE(getTarget(&es8) == 5);
	target(&es8)=7;
	SEQAN_ASSERT_TRUE(getTarget(&es8) == 7);
	assignSource(&es8, 20);
	SEQAN_ASSERT_TRUE(getSource(&es8) == 0); // No source
	SEQAN_ASSERT_TRUE(source(&es8) == 0);
	assignCargo(&es8, 15); 
	SEQAN_ASSERT_TRUE(getCargo(&es8) == 15);
	SEQAN_ASSERT_TRUE(cargo(&es8) == 15);
	SEQAN_ASSERT_TRUE(_getId(&es8) == 0);
	EdgeStump<unsigned int, true, false, false> const es8_const(es8);
	SEQAN_ASSERT_TRUE(getCargo(&es8_const) == 15);
	SEQAN_ASSERT_TRUE(cargo(&es8_const) == 15);
	SEQAN_ASSERT_TRUE(getTarget(&es8_const) == 7);
	SEQAN_ASSERT_TRUE(target(&es8_const) == 7);
	EdgeStump<unsigned int, true, false, false> es81(es8);
	EdgeStump<unsigned int, true, false, false> es82(es8);
	nextT(&es8) = &es81;
	SEQAN_ASSERT_TRUE(getNextT(&es8) == &es81);
	assignNextT(&es8, &es82);
	SEQAN_ASSERT_TRUE(getNextT(&es8) == &es82);
	assignNextS(&es8, &es82);
	SEQAN_ASSERT_TRUE(getNextS(&es8) == 0);
	SEQAN_ASSERT_TRUE(nextS(&es8) == 0);

//____________________________________________________________________________
// Test all cargoless EdgeStumps in an array
	// No cargo, no list, source, id
	EdgeStump<void, false, true, true> es9;
	_assignId(&es9, 4);
	SEQAN_ASSERT_TRUE(_getId(&es9) == 4);
	assignTarget(&es9, 5);
	SEQAN_ASSERT_TRUE(getTarget(&es9) == 5);
	target(&es9)=7;
	SEQAN_ASSERT_TRUE(getTarget(&es9) == 7);
	assignSource(&es9, 20);
	SEQAN_ASSERT_TRUE(getSource(&es9) == 20);
	source(&es9)=30;
	SEQAN_ASSERT_TRUE(getSource(&es9) == 30);
	assignCargo(&es9, 15);  //Does nothing, no cargo -> pointer to 0
	SEQAN_ASSERT_TRUE(getCargo(&es9) == (void*) 0);
	SEQAN_ASSERT_TRUE(cargo(&es9) == (void*) 0);
	SEQAN_ASSERT_TRUE(_getId(&es9) == 4);
	EdgeStump<void, false, true, true> const es9_const(es9);
	SEQAN_ASSERT_TRUE(getCargo(&es9_const) == (void*) 0);
	SEQAN_ASSERT_TRUE(cargo(&es9_const) == (void*) 0);
	SEQAN_ASSERT_TRUE(getTarget(&es9_const) == 7);
	SEQAN_ASSERT_TRUE(target(&es9_const) == 7);

	// No cargo, no list, source, no id
	EdgeStump<void, false, true, false> es10;
	_assignId(&es10, 4);
	SEQAN_ASSERT_TRUE(_getId(&es10) == 0); // No id
	assignTarget(&es10, 5);
	SEQAN_ASSERT_TRUE(getTarget(&es10) == 5);
	target(&es10)=7;
	SEQAN_ASSERT_TRUE(getTarget(&es10) == 7);
	assignSource(&es10, 20);
	SEQAN_ASSERT_TRUE(getSource(&es10) == 20);
	source(&es10)=30;
	SEQAN_ASSERT_TRUE(getSource(&es10) == 30);
	assignCargo(&es10, 15);  //Does nothing, no cargo -> pointer to 0
	SEQAN_ASSERT_TRUE(getCargo(&es10) == (void*) 0);
	SEQAN_ASSERT_TRUE(cargo(&es10) == (void*) 0);
	SEQAN_ASSERT_TRUE(_getId(&es10) == 0);
	EdgeStump<void, false, true, false> const es10_const(es10);
	SEQAN_ASSERT_TRUE(getCargo(&es10_const) == (void*) 0);
	SEQAN_ASSERT_TRUE(cargo(&es10_const) == (void*) 0);
	SEQAN_ASSERT_TRUE(getTarget(&es10_const) == 7);
	SEQAN_ASSERT_TRUE(target(&es10_const) == 7);

	// No cargo, no list, no source, id
	EdgeStump<void, false, false, true> es_11;
	_assignId(&es_11, 4);
	SEQAN_ASSERT_TRUE(_getId(&es_11) == 4);
	assignTarget(&es_11, 5);
	SEQAN_ASSERT_TRUE(getTarget(&es_11) == 5);
	target(&es_11)=7;
	SEQAN_ASSERT_TRUE(getTarget(&es_11) == 7);
	assignSource(&es_11, 20);
	SEQAN_ASSERT_TRUE(getSource(&es_11) == 0); // No source
	SEQAN_ASSERT_TRUE(source(&es_11) == 0);
	assignCargo(&es_11, 15);  //Does nothing, no cargo -> pointer to 0
	SEQAN_ASSERT_TRUE(getCargo(&es_11) == (void*) 0);
	SEQAN_ASSERT_TRUE(cargo(&es_11) == (void*) 0);
	SEQAN_ASSERT_TRUE(_getId(&es_11) == 4);
	EdgeStump<void, false, false, true> const es_11_const(es_11);
	SEQAN_ASSERT_TRUE(getCargo(&es_11_const) == (void*) 0);
	SEQAN_ASSERT_TRUE(cargo(&es_11_const) == (void*) 0);
	SEQAN_ASSERT_TRUE(getTarget(&es_11_const) == 7);
	SEQAN_ASSERT_TRUE(target(&es_11_const) == 7);

	
	// No cargo, no list, no source, no id
	EdgeStump<void, false, false, false> es_12;
	_assignId(&es_12, 4);
	SEQAN_ASSERT_TRUE(_getId(&es_12) == 0);
	assignTarget(&es_12, 5);
	SEQAN_ASSERT_TRUE(getTarget(&es_12) == 5);
	target(&es_12)=7;
	SEQAN_ASSERT_TRUE(getTarget(&es_12) == 7);
	assignSource(&es_12, 20);
	SEQAN_ASSERT_TRUE(getSource(&es_12) == 0); // No source
	SEQAN_ASSERT_TRUE(source(&es_12) == 0);
	assignCargo(&es_12, 15);  //Does nothing, no cargo -> pointer to 0
	SEQAN_ASSERT_TRUE(getCargo(&es_12) == (void*) 0);
	SEQAN_ASSERT_TRUE(cargo(&es_12) == (void*) 0);
	SEQAN_ASSERT_TRUE(_getId(&es_12) == 0);
	EdgeStump<void, false, false, false> const es_12_const(es_12);
	SEQAN_ASSERT_TRUE(getCargo(&es_12_const) == (void*) 0);
	SEQAN_ASSERT_TRUE(cargo(&es_12_const) == (void*) 0);
	SEQAN_ASSERT_TRUE(getTarget(&es_12_const) == 7);
	SEQAN_ASSERT_TRUE(target(&es_12_const) == 7);

//____________________________________________________________________________
// Test all EdgeStumps with cargo in an array
	// cargo, no list, source, id
	EdgeStump<unsigned int, false, true, true> es_13;
	_assignId(&es_13, 4);
	SEQAN_ASSERT_TRUE(_getId(&es_13) == 4);
	assignTarget(&es_13, 5);
	SEQAN_ASSERT_TRUE(getTarget(&es_13) == 5);
	target(&es_13)=7;
	SEQAN_ASSERT_TRUE(getTarget(&es_13) == 7);
	assignSource(&es_13, 20);
	SEQAN_ASSERT_TRUE(getSource(&es_13) == 20);
	source(&es_13)=30;
	SEQAN_ASSERT_TRUE(getSource(&es_13) == 30);
	assignCargo(&es_13, 15);  
	SEQAN_ASSERT_TRUE(getCargo(&es_13) == 15);
	SEQAN_ASSERT_TRUE(cargo(&es_13) == 15);
	SEQAN_ASSERT_TRUE(_getId(&es_13) == 4);
	EdgeStump<unsigned int, false, true, true> const es_13_const(es_13);
	SEQAN_ASSERT_TRUE(getCargo(&es_13_const) == 15);
	SEQAN_ASSERT_TRUE(cargo(&es_13_const) == 15);
	SEQAN_ASSERT_TRUE(getTarget(&es_13_const) == 7);
	SEQAN_ASSERT_TRUE(target(&es_13_const) == 7);

	// cargo, no list, source, no id
	EdgeStump<unsigned int, false, true, false> es_14;
	_assignId(&es_14, 4);
	SEQAN_ASSERT_TRUE(_getId(&es_14) == 0); // No id
	assignTarget(&es_14, 5);
	SEQAN_ASSERT_TRUE(getTarget(&es_14) == 5);
	target(&es_14)=7;
	SEQAN_ASSERT_TRUE(getTarget(&es_14) == 7);
	assignSource(&es_14, 20);
	SEQAN_ASSERT_TRUE(getSource(&es_14) == 20);
	source(&es_14)=30;
	SEQAN_ASSERT_TRUE(getSource(&es_14) == 30);
	assignCargo(&es_14, 15); 
	SEQAN_ASSERT_TRUE(getCargo(&es_14) == 15);
	SEQAN_ASSERT_TRUE(cargo(&es_14) == 15);
	SEQAN_ASSERT_TRUE(_getId(&es_14) == 0);
	EdgeStump<unsigned int, false, true, false> const es_14_const(es_14);
	SEQAN_ASSERT_TRUE(getCargo(&es_14_const) == 15);
	SEQAN_ASSERT_TRUE(cargo(&es_14_const) == 15);
	SEQAN_ASSERT_TRUE(getTarget(&es_14_const) == 7);
	SEQAN_ASSERT_TRUE(target(&es_14_const) == 7);

	// cargo, no list, no source, id
	EdgeStump<unsigned int, false, false, true> es_15;
	_assignId(&es_15, 4);
	SEQAN_ASSERT_TRUE(_getId(&es_15) == 4);
	assignTarget(&es_15, 5);
	SEQAN_ASSERT_TRUE(getTarget(&es_15) == 5);
	target(&es_15)=7;
	SEQAN_ASSERT_TRUE(getTarget(&es_15) == 7);
	assignSource(&es_15, 20);
	SEQAN_ASSERT_TRUE(getSource(&es_15) == 0); // No source
	SEQAN_ASSERT_TRUE(source(&es_15) == 0);
	assignCargo(&es_15, 15); 
	SEQAN_ASSERT_TRUE(getCargo(&es_15) == 15);
	SEQAN_ASSERT_TRUE(cargo(&es_15) == 15);
	SEQAN_ASSERT_TRUE(_getId(&es_15) == 4);
	EdgeStump<unsigned int, false, false, true> const es_15_const(es_15);
	SEQAN_ASSERT_TRUE(getCargo(&es_15_const) == 15);
	SEQAN_ASSERT_TRUE(cargo(&es_15_const) == 15);
	SEQAN_ASSERT_TRUE(getTarget(&es_15_const) == 7);
	SEQAN_ASSERT_TRUE(target(&es_15_const) == 7);

	// cargo, no list, no source, no id
	EdgeStump<unsigned int, false, false, false> es_16;
	_assignId(&es_16, 4);
	SEQAN_ASSERT_TRUE(_getId(&es_16) == 0);
	assignTarget(&es_16, 5);
	SEQAN_ASSERT_TRUE(getTarget(&es_16) == 5);
	target(&es_16)=7;
	SEQAN_ASSERT_TRUE(getTarget(&es_16) == 7);
	assignSource(&es_16, 20);
	SEQAN_ASSERT_TRUE(getSource(&es_16) == 0); // No source
	SEQAN_ASSERT_TRUE(source(&es_16) == 0);
	assignCargo(&es_16, 15);  
	SEQAN_ASSERT_TRUE(getCargo(&es_16) == 15);
	SEQAN_ASSERT_TRUE(cargo(&es_16) == 15);
	SEQAN_ASSERT_TRUE(_getId(&es_16) == 0);
	EdgeStump<unsigned int, false, false, false> const es_16_const(es_16);
	SEQAN_ASSERT_TRUE(getCargo(&es_16_const) == 15);
	SEQAN_ASSERT_TRUE(cargo(&es_16_const) == 15);
	SEQAN_ASSERT_TRUE(getTarget(&es_16_const) == 7);
	SEQAN_ASSERT_TRUE(target(&es_16_const) == 7);

	// Special tree edge stump, in a tree the target is the child id
	EdgeStump<unsigned int, true, true, false, TreeTag> es17;
	assignTarget(&es17, 5);
	SEQAN_ASSERT_TRUE(_getId(&es17) == 5);
	EdgeStump<unsigned int, true, true, false, TreeTag> const es17_const(es17);
	SEQAN_ASSERT_TRUE(_getId(&es17_const) == 5);
}


//////////////////////////////////////////////////////////////////////////////

void Test_GraphBasics() {
	Test_IdManager();	// Test Id Manager
	Test_EdgeStump();	// Test EdgeStumps
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(test_graph_basics)
{
	Test_GraphBasics();
}


#endif

