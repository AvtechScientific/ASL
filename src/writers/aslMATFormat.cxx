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


#include "aslMATFormat.h"

#include <matio.h>

#include <utilities/aslMATLABCasters.h>
#include <data/aslDataWrapper.h>
#include <data/aslProbe.h>

#include "math/aslVectors.h"
#include "acl/acl.h"

using  namespace std;

namespace asl {

	void writeMAT(const string &fileName, const AbstractData & data, const string & name)
	{
		mat_t* matFile(Mat_Create(fileName.c_str(),"This file was created by ASL <http://asl.org.il>"));

		unsigned int nComp(data.getDContainer().size());
		std::vector<string> names(nComp,name);
		for (unsigned int i(0); i < nComp; ++i)
			names[i]+="-"+numToStr(i);

		SPMatVar var(castMATLABCellArray(data,names));
		Mat_VarWrite(matFile,var->var, 0);
        Mat_Close(matFile);
	}


	void writeMAT(const string &fileName, vector<pair<SPAbstractData,string>> data)
	{
		mat_t* matFile(Mat_Create(fileName.c_str(),"This file was created by ASL <http://asl.org.il>"));

		for(unsigned int k(0); k<data.size(); ++k)
		{
			unsigned int nComp(data[k].first->getDContainer().size());
			std::vector<string> names(nComp,data[k].second);
			for (unsigned int i(0); i < nComp; ++i)
				names[i]+="-"+numToStr(i);

			SPMatVar var(castMATLABCellArray(*data[k].first,names));
			Mat_VarWrite(matFile,var->var, 0);
		}
        Mat_Close(matFile);		
	}

	void writeMAT(const string &fileName, vector<pair<acl::Element,string>> data)
	{
		mat_t* matFile(Mat_Create(fileName.c_str(),"This file was created by ASL <http://asl.org.il>"));

		for(unsigned int k(0); k<data.size(); ++k)
		{
			SPMatVar var(castMATLABCellArray(data[k].first,data[k].second));
			Mat_VarWrite(matFile,var->var, 0);
		}
        Mat_Close(matFile);		
	}
	
	
	void writeMAT(const string &fileName,Probe & probe, unsigned int component, const string & name)
	{
		mat_t* matFile(Mat_Create(fileName.c_str(),"This file was created by ASL <http://asl.org.il>"));

		SPMatVar var(castMATLABCellArray(probe,component,name));
		Mat_VarWrite(matFile,var->var, 0);
        Mat_Close(matFile);
	}
		
} //asl


