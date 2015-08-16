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


#include <stdexcept>
#include "aslUtilities.h"


using namespace std;

namespace asl
{

	template <typename T> T strToNum(string s)
	{
		istringstream i(s);
		T x;
		char c;

		if (!(i >> x))
			errorMessage("strToNum() - unable to convert " + s + " to the requested type");

		if (i.get(c))
			warningMessage("strToNum() - " + s + " contains a character");

/*
		ToDo: check whether it is needed at all.
		Commented out since it causes cyclic dependency aslcommon <-> aslacl through "aclUtilites.h"

		int decimal(0);
		if (acl::typeToTypeID<T>() == acl::TYPE_INT)
		{
			i >> decimal;
			if (decimal != 0)
				warningMessage("strToNum() - " + s + " is not an integer");
		}
*/
		return x;
	}


	template cl_int strToNum<cl_int>(string s);
	template cl_uint strToNum<cl_uint>(string s);
	template cl_long strToNum<cl_long>(string s);
	template cl_float strToNum<cl_float>(string s);
	template cl_double strToNum<cl_double>(string s);


	void errorMessage(cl_int status, const char *errorMessage)
	{
		if (status != CL_SUCCESS)
		{
			throw std::logic_error("ASL ERROR: " + string(errorMessage) + " (" + to_string(status) + ")." );
		}
	}


	void errorMessage(cl_int status, const string &errorMessage)
	{
		if (status != CL_SUCCESS)
		{
			throw std::logic_error("ASL ERROR: " + errorMessage + " (" + to_string(status) + ")." );
		}
	}
	

	void errorMessage(const char *errorMessage)
	{
			throw std::logic_error("ASL ERROR: " + string(errorMessage));	
	}


	void errorMessage(const string & errorMessage)
	{
			throw std::logic_error("ASL ERROR: " + errorMessage);	
	}


	void errorMessage(bool status)
	{
		if (status)
			cout << " Ok" << endl;
		else
			cerr << " ERROR" << endl;
	}

	
	void warningMessage(const char *warningMessage)
	{
		cout << "ASL WARNING: " << warningMessage << "." << endl;
	}


	void warningMessage(const string & warningMessage)
	{
		cout << "ASL WARNING: " << warningMessage << "." << endl;
	}

	
	string warningString(const char *warningMessage)
	{
		return string("ASL WARNING: ") + warningMessage + ".";
	}

} // asl namespace