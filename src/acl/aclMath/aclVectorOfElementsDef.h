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


#ifndef ACLVECTORELEMENTSDEF_H
#define ACLVECTORELEMENTSDEF_H

#include <memory>
#include <vector>
#include "../aclTypes.h"

using namespace std;

namespace cl
{
	class CommandQueue;
}

namespace acl
{
	
	class ElementBase;
	class MemBlock;
	
	typedef std::shared_ptr<ElementBase> Element;
	typedef std::shared_ptr<MemBlock> ElementData;

	class VectorOfElements;
	typedef std::shared_ptr<cl::CommandQueue> CommandQueue;
	
	/// The class represents several ::ElementData
	/**
		 \ingroup ComplexDataTypes
		 Basicaly the class is identical to the std::vector class. The only 
		 difference is the assigment operator and mathematical operations.
	 */
	class VectorOfElementsData: public std::vector<ElementData>
	{
			/// checks whether all Elements have sizes compatible
			/// with each other and reside on the same device
			bool checkCompatibility() const;
		public:
			/// checks whether all Elements have sizes compatible with \p n
			bool checkSizesCompatibility(unsigned int n) const;
			VectorOfElementsData();
			explicit VectorOfElementsData(unsigned int n);
			/**
			 \param n number of ElementData
			 \param s size of each ElementData
			 \param T type of ElementData: e.g. float double etc.
			 \param queue defines a device where the DataElements should be placed
			 */
			template <typename T> VectorOfElementsData(unsigned int n,
			                                           unsigned int s,
			                                           T);
			template <typename T> VectorOfElementsData(unsigned int n,
			                                           unsigned int s,
			                                           T,
			                                           CommandQueue queue);

			VectorOfElements operator=(const VectorOfElements & a) const;
			VectorOfElements operator=(const VectorOfElementsData & a) const;
			void resizeElements(unsigned int n);
	};

	/// The class represents several ::Element
	/**
		 \ingroup ComplexDataTypes
		 Basicaly the class is identical to the std::vector class. The only 
		 difference is the assigment operator and mathematical operations.
	 */
	class VectorOfElements: public std::vector<Element>
	{
			/// checks whether all Elements have sizes compatible
			/// with each other and reside on the same device
			bool checkCompatibility() const;
		public:						
			/// checks whether all Elements have sizes compatible with \p n
			bool checkSizesCompatibility(unsigned int n) const;
			VectorOfElements();
			explicit VectorOfElements(unsigned int n);
			VectorOfElements(const VectorOfElementsData & a);
			VectorOfElements operator=(const VectorOfElements & a) const;
	};

	typedef shared_ptr<VectorOfElementsData> SPVectorOfElementsData;
	typedef shared_ptr<VectorOfElements> SPVectorOfElements;

	inline bool compatibleSizes (unsigned int s, const VectorOfElements & a);

	///\related VectorOfElements
	acl::TypeID getElementType(const VectorOfElements & a, unsigned int i = 0);
	
	///\related VectorOfElements
	unsigned int getElementsSize(const VectorOfElements & a);
	
//------------------------- Implementation -----------------
	
	inline bool compatibleSizes(unsigned int s, const VectorOfElements & a)
	{
		return a.checkSizesCompatibility(s);
	}

	
}  // namespace acl

#endif // ACLVECTORELEMENTSDEF_H
