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


#ifndef ASLTEMPLATESEXTRAS_H
#define ASLTEMPLATESEXTRAS_H

#include "aslTemplates.h"
#include "aslMatrices.h"

namespace acl {
	class VectorOfElements;
}

namespace asl {


	/// Defines andditionl features related to a VectorTemplate
	/**
		\ingroup Templates
		contains list of edges for a template. this can be used in level 
		set and iso-surface extraction algorithms
	*/
	class VTObjects
	{
		private:
			void initCellMatrices();
		public:
			const VectorTemplate * vt;
			std::vector<unsigned int> edgePoint1;
			std::vector<unsigned int> edgePoint2;
			std::vector<AVec<int>> elementaryCells;
			std::vector<AMatr<>> cellMatrices;
			
			VTObjects(const VectorTemplate * vt, 
			          const std::vector<unsigned int> & ep1, 
			          const std::vector<unsigned int> & ep2,
			          const std::vector<AVec<int>> & elCells);

			/// computes gradient within the elementary cell \p ic and  values in the corners \p val  
			acl::VectorOfElements cellGradient(const acl::VectorOfElements & val, 
			                                   unsigned int ic) const;
			/// fill \p points by corner coordinates of the cell \p ic 
			void getCellPoints(unsigned int ic, std::vector<AVec<>> & points) const;
	};

	inline unsigned int nD(const VTObjects vto);

	/// Returns vtObjecs corresponding to the given VectorTemplate \ingroup Templates
	const VTObjects* vtObject(const VectorTemplate *);

	/// returns coefficient of the laplace operator corresponding to the direction of a cell edge
	/**
		 \related VTObjects
	 */
	double edgeWeight(const VTObjects & vto, unsigned int iEl, unsigned int i, unsigned int j);
	
	
	///Vector template
	/**	
	 \ingroup TemplatesNN
		 
	 \image html t5.png
	*/	
	const VTObjects & d2q5Objs();	

	///Vector template
	/**	
	 \ingroup TemplatesNN
		 
	 \image html t7.png
	*/
	const VTObjects & d3q7Objs();	

	///Vector template
	/**	
	 \ingroup TemplatesNN
		 
	 \image html t9.png
	*/	
	const VTObjects & d2q9Objs();	

	///Vector template
	/**	
	 \ingroup TemplatesNN
		 
	 \image html t15.png
	*/	
	const VTObjects & d3q15Objs();	
	
	///Vector template
	/**	
	 \ingroup TemplatesNN
		 
	 \image html t19.png
	*/	
	const VTObjects & d3q19Objs();	

// ----------------------------- Implementation -------------------------

	inline unsigned int nD(const VTObjects & vto)
	{
		return nD(*vto.vt);
	}

}// asl

#endif // ASLTEMPLATESEXTRAS_H
