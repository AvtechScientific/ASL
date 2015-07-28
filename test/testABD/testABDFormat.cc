/*
 * Advanced Simulation Library <http://asl.org.il>
 * 
 * Copyright 2015 Avtech Scientific <http://avtechscientific.com>
 *
 *
 * This file is part of Advanced Simulation Library (ASL).
 *
 * ASL is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, version 3 of the License.
 *
 * ASL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with ASL. If not, see <http://www.gnu.org/licenses/>.
 *
 */


/**
	\example testABDFormat.cc
 */

#include "writers/aslABDFormat.h"
#include "aslUtilities.h"
#include "math/aslVectors.h"
#include "data/aslBlocks.h"

bool testNumbers()
{
	cout << "Test of Numbers..." << flush;
	unsigned int aui(3);
	int ai(-2);	
	float af(5);
	double ad(4);	
	
	asl::ABDFileOut afO("test.abd");
	afO << aui << ai << af << ad;
	afO.close();
		
	asl::ABDFileIn afI("test.abd");
	unsigned int bui(0);
	int bi(0);	
	float bf(0);
	double bd(0);	

	afI >> bui >> bi >> bf >> bd;

	bool status((aui==bui) && (ai==bi) && (af==bf) && (ad==bd));
	asl::errorMessage(status);

	return status;
}

bool testAVec()
{
	cout << "Test of AVec..." << flush;
	asl::Block b(asl::makeAVec (10,15),0.1,asl::makeAVec (.1,1.));

	asl::ABDFileOut afO("test.abd");
	afO << b;
	afO.close();
		
	asl::ABDFileIn afI("test.abd");
	asl::Block bn;

	afI >> bn;

	bool status((b.getSize() == bn.getSize()) &&
	            (b.dx == bn.dx) &&
	            (b.position == bn.position));
	asl::errorMessage(status);

	return status;		
}

bool testString()
{
	cout << "Test of String..." << flush;
	std::string b("Hello!!");

	asl::ABDFileOut afO("test.abd");
	afO << b;
	afO.close();
		
	asl::ABDFileIn afI("test.abd");
	std::string bn;

	afI >> bn;

	bool status(b == bn);
	asl::errorMessage(status);

	return status;		
}


bool testBlock()
{
	cout << "Test of Block..." << flush;
	asl::AVec<int> ai(asl::makeAVec(2,3));	
	asl::AVec<float> af(asl::makeAVec(2.f,3.f));	
	asl::AVec<double> ad(asl::makeAVec(4.,5.));	

	asl::ABDFileOut afO("test.abd");
	afO << ai << af << ad;
	afO.close();
		
	asl::ABDFileIn afI("test.abd");
	asl::AVec<int> bi(1);	
	asl::AVec<float> bf(1);
	asl::AVec<double> bd(1);	

	afI >> bi >> bf >> bd;

	bool status((ai == bi) && (af == bf) && (ad == bd));
	asl::errorMessage(status);

	return status;
}

int main()
{
	bool allTestsPassed(true);

	allTestsPassed &= testNumbers();
	allTestsPassed &= testAVec();
	allTestsPassed &= testString();
	allTestsPassed &= testBlock();

	return allTestsPassed ? EXIT_SUCCESS : EXIT_FAILURE;
}