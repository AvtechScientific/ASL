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


#ifndef ASLLEVELSET_H
#define ASLLEVELSET_H

#include "aslNumMethod.h"
#include "acl/aclMath/aclVectorOfElementsDef.h"


namespace acl{
	class Kernel;
	class VectorOfElementsData;
}

namespace asl
{
	class VectorTemplate;
	class VTObjects;
	template <typename V> class DataWithGhostNodes;
	typedef DataWithGhostNodes<acl::VectorOfElementsData> 
		DataWithGhostNodesACLData;
	typedef std::shared_ptr<DataWithGhostNodesACLData> 
		SPDataWithGhostNodesACLData;
	typedef DataWithGhostNodes<acl::VectorOfElements> 
		DataWithGhostNodesACL;
	typedef std::shared_ptr<DataWithGhostNodesACL> 
		SPDataWithGhostNodesACL;
	class TemplateVE;
	
	/// Numerical method which computes evolution of an interface
	/**
		 \ingroup InterfaceTracking
		 \ingroup NumMethods
		 \ingroup LevelSet		 
	 */
	class LevelSet: public NumMethod
	{
		public:
			typedef SPDataWithGhostNodesACLData Data;
			typedef SPDataWithGhostNodesACL DataGen;
			/// generates expression resulting true if there is a boundary within the element
			/**
				 \param iEl element number 
			*/
			acl::VectorOfElements isBoundaryEl(unsigned int iEl);

			/// generates expression resulting true if there is a boundary in this direction
			/**
				 \param iDir direction number 
			*/
			acl::VectorOfElements isBoundaryDir(unsigned int iDir);
			
			/// generates expression resulting the value of the \p field  
			/**
				 \param field vector contains values in each direction
				 \param iDir is the direction
			*/			
			acl::VectorOfElements getValueOnBoundary(acl::VectorOfElements field,
			                                         unsigned int iDir);

			/// generates expression for center of a boundary element
			/**
				 \param iEl the element number
				 computes avarage point of corners of the boundary poligon within the element 

			*/
			acl::VectorOfElements getBoundaryCenter(unsigned int iEl);

			/// generates expression for area of a boundary element
			/**
				 \param iEl the element number
				 computes area of the corresponding boundary within cell iEl
			*/
			vector<acl::Element> gcBoundaryArea(unsigned int iEl, 
			                                    acl::VectorOfElements & center, 
			                                    acl::VectorOfElements & area);

			/// generates expression for area of a boundary element
			/**
				 computes area of the corresponding boundary within cell iEl
			*/
			vector<acl::Element> gcBoundaryArea(acl::VectorOfElements & center, 
			                                    acl::VectorOfElements & area);			
			/// computes coordinates of the surface point on the \p iDir vector
			/**
				 \param iDir the element number

			*/
			acl::VectorOfElements getBoundaryPoint(unsigned int iDir);
			
		protected:
			std::unique_ptr<acl::Kernel> kernel;

			Data distanceField;
			Data distanceFieldInternalData;			
			
			const VectorTemplate* vectorTemplate;
			const VTObjects* vto;

			unique_ptr<TemplateVE> distanceTVE;
			vector<acl::VectorOfElements> lVelocities;
			
			void initKernelPropagation();
			virtual void initVelocityComputation()=0;
			virtual void initDistancesComputation()=0;

		public:			
			LevelSet();
			/**
			 \param df is distance field
			 */
			LevelSet(Data df);
			~LevelSet();

			inline const VectorTemplate* getVectorTemplate() const;

			virtual void init();
			virtual void execute();
	};

	typedef std::shared_ptr<LevelSet> SPLevelSet;

// --------------------------- Implementation --------------------
	
	inline const VectorTemplate* LevelSet::getVectorTemplate() const
	{
		return vectorTemplate;
	}
	
} //asl

#endif //ASLLEVELSET_H
