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


#ifndef ASLMATFORMAT_H
#define ASLMATFORMAT_H


#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <vector>


using  namespace std;

namespace acl{
	class ElementBase;
	typedef std::shared_ptr<ElementBase> Element;
}

namespace asl {
	class AbstractData;
	class Probe;
	typedef shared_ptr<AbstractData> SPAbstractData;
	/// writes \p data in a MatLab file  
	/**
		 \ingroup IO
	*/
	void writeMAT(const string &fileName, const AbstractData & data, const string & name);

	/// writes \p data in a MatLab file  
	/**
		 \ingroup IO
		 \param fileName name of the file;
		 \param data is vector which contains pairs of the corresponding 
		 data and its name
	*/
	void writeMAT(const string &fileName, vector<pair<SPAbstractData,string>> data);

	/// writes \p data in a MatLab file  
	/**
		 \ingroup IO
		 \param fileName name of the file;
		 \param data is vector which contains pairs of the corresponding 
		 data and its name
	*/
	void writeMAT(const string &fileName, vector<pair<acl::Element,string>> data);

	
	/// writes \p probe values in a MatLab file  
	/**
		 \ingroup IO
	*/
	void writeMAT(const string &fileName, Probe & probe, unsigned int component, const string & name);

}// asl

#endif //ASLVTKFORMAT_H

