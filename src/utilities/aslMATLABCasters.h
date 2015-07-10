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


#ifndef ASLMATLABCASTERS_H
#define ASLMATLABCASTERS_H

#include <memory>
#include <string>
#include <vector>

/**
	\defgroup MATLABInterfacing MATLAB Interfacing
	\ingroup Interfacing
*/

namespace acl{
	class ElementBase;
	typedef std::shared_ptr<ElementBase> Element;
}

struct matvar_t;
namespace asl{
	template <class T> class AVec;
}

/// the matio c++ wrapper \ingroup MATLABInterfacing
namespace matiocpp{
	class Var{
		private:
			template <class T> matvar_t* init(T *d, const asl::AVec<int> & size, const std::string &name="");
		public:
			matvar_t* var;
			bool freeArray;
			inline Var(matvar_t* v,bool freeArr=false):var(v){}
			template <class T> Var(std::vector<T> & d, const std::string &name="");			
			template <class T> Var(T *d, unsigned int size, const std::string &name="",bool freeArr=true);
			template <class T> Var(T *d, const asl::AVec<int> & size, const std::string &name="",bool freeArr=true);
//			Var(T *matiocpp, const asl::AVec<int> & size, const std::string &name=""); ///!!!
			~Var();
	};
}

namespace asl{
	class Block;
	class AbstractData;
	class Probe;
	
	/// @{
	/// \ingroup MATLABInterfacing
	typedef std::shared_ptr<matiocpp::Var> SPMatVar;
	
	SPMatVar castMATLABCellArray(acl::Element source, const std::string &name="");
	SPMatVar castMATLABCellArray(acl::Element source, const AVec<int> & size, const std::string &name="");
	SPMatVar castMATLABCellArray(const AbstractData & d, const std::vector<std::string> &names=
	                          								std::vector<std::string>(0));	
	SPMatVar castMATLABCellArray(Probe & p, unsigned int component, const std::string &name="");
	/// }@
}  //namespace acl


#endif // ASLMATLABCASTERS_H
