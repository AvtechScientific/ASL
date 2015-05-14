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


#ifndef ASLBCOND_H
#define ASLBCOND_H

#include "aslNumMethod.h"
#include <data/aslBlocks.h>
#include <acl/aclMath/aclVectorOfElementsDef.h>


namespace acl
{
	class ExpressionContainer;
}

namespace asl
{

	class VectorTemplate;
	template <typename V> class DataWithGhostNodes;
	typedef DataWithGhostNodes<acl::VectorOfElementsData> DataWithGhostNodesACLData;
	typedef std::shared_ptr<DataWithGhostNodesACLData> SPDataWithGhostNodesACLData;
	typedef DataWithGhostNodes<acl::VectorOfElements> DataWithGhostNodesACL;
	typedef std::shared_ptr<DataWithGhostNodesACL> SPDataWithGhostNodesACL;
	class TemplateVE;
	class DistanceFunction;
	typedef std::shared_ptr<DistanceFunction> SPDistanceFunction;

	
	class AbstractDataWithGhostNodes; 
	typedef std::shared_ptr<AbstractDataWithGhostNodes> SPAbstractDataWithGhostNodes;
	
	/// Virtual class describes general interface for boundary conditions
	/**
		  \ingroup BoundaryConditions
	*/
	class BCond: public NumMethod
	{
		protected:
			Block block;
			const VectorTemplate * const templ;

			std::vector<int> indices;
			std::vector<int> directions;
			std::vector<double> fractions;
			std::vector<int> neighbourIndices;
			
			acl::SPVectorOfElementsData indicesACL;
			acl::SPVectorOfElementsData neighbourIndicesACL;
			acl::SPVectorOfElementsData directionsACL;
			acl::SPVectorOfElementsData fractionsACL;
			
			void loadIndicesToACL();
			void loadNeighbourIndicesToACL();
			void loadDirectionsToACL();
			void loadfractionsACL();

		public:
			BCond(const Block & b);
			BCond(const Block & b, const VectorTemplate * const t);
			void addPoint(AVec<int> a,int d=0,double fr=0);			
			const Block & getBlock();
			inline const VectorTemplate * getVT();
			
	};

	typedef std::shared_ptr<BCond> SPBCond;

	/// Virtual class describes general interface for boundary conditions
	/**
		\ingroup BoundaryConditions

		The class is aimed to deal with numerical with computational map. Namely
		it should be a DataWithGhestNodes object with values in each point 
		which can be interpreted as a computation map. The inteprretation law is:
		- computation region: \f$ map(\vec r) > 0 \f$
		- outside region: \f$ map(\vec r) \leq 0 \f$
	*/
	class BCondWithMap:public NumMethod
	{
		protected:
			/// flag whether the point list to be generated or not
			bool pointsListFlag;
			acl::VectorOfElements currentPoint;
			const VectorTemplate * const templ;

			/// block 
			Block bl;
			/// boundary description for the particular BC
			SPAbstractDataWithGhostNodes map;			
			SPDistanceFunction mapDF;
			/// Computational domain which contains all boundaries and the particular boundary as well
			SPAbstractDataWithGhostNodes computationalDomain;
			SPDistanceFunction computationalDomainDF;

			unique_ptr<TemplateVE> mapTVE;
			unique_ptr<TemplateVE> cDomainTVE;

			/// initialize mapTVE and cDomainTVE
			virtual void initMapInfrastructure(acl::ExpressionContainer & ec);
			/// returns expression corresponding to check if the node in i^th direction is ghost one
			/**
				 Requires call of initMapInfrastructure
			*/
			acl::VectorOfElements isGhostNode(unsigned int i);
			/// returns expression corresponding to check if the node in i^th direction is computation one
			/**
				 Requires call of initMapInfrastructure
			*/
			acl::VectorOfElements isComputationNode(unsigned int i);
			/// returns expression corresponding to check if nodes in directions \p ii are computation ones
			/**
				 Requires call of initMapInfrastructure
			*/
			acl::VectorOfElements isComputationNode(const vector<unsigned int> & ii);
			/// returns expression corresponding to check if the current node is ghost one
			acl::VectorOfElements isGhostNode();
			/// returns expression corresponding to check if the current node is computation one
			acl::VectorOfElements isComputationNode();
			
			/**
			 \param m the map
			 \param t the corresponding template
			 */
			BCondWithMap(SPAbstractDataWithGhostNodes m, const VectorTemplate * const vt);
			/**
			 \param m the map
			 \param t the corresponding template
			 */
			BCondWithMap(SPDistanceFunction m, const Block & b, const VectorTemplate * const vt);
			/**
			 \param m the map for the particular boundary
			 \param cd the map for the computational domain
			 \param t the corresponding template
			 */
			BCondWithMap(SPAbstractDataWithGhostNodes m, 
			             SPAbstractDataWithGhostNodes cd, 
			             const VectorTemplate * const vt);			
			/**
			 \param m the map for the particular boundary
			 \param cd the map for the computational domain
			 \param t the corresponding template
			 */
			BCondWithMap(SPAbstractDataWithGhostNodes m, 
			             SPDistanceFunction cd, 
			             const VectorTemplate * const vt);			
			/**
			 \param m the map for the particular boundary
			 \param cd the map for the computational domain
			 \param t the corresponding template
			 \param b the computation block
			 */
			BCondWithMap(SPDistanceFunction m, 
			             SPDistanceFunction cd,
			             const Block & b,
			             const VectorTemplate * const vt);			
		public:
			inline const VectorTemplate * getVT();			
	};

	
	/// Virtual class describes general interface for boundary conditions which connect two datas
	/**
		 \ingroup BoundaryConditions
		 The class contains an explicit list of coneting points. 
		 This gives certain flexibility for definition of the boundary geometry.
	 */
	class BCondConnector:public NumMethod
	{
		protected:
			Block block1;
			Block block2;
			
			std::vector<int> indices1;
			std::vector<int> directions1;
			std::vector<int> indices2;
			std::vector<int> directions2;

			acl::SPVectorOfElementsData indices1ACL;
			acl::SPVectorOfElementsData directions1ACL;
			acl::SPVectorOfElementsData indices2ACL;
			acl::SPVectorOfElementsData directions2ACL;
			
			void loadIndicesToACL();
			void loadDirectionsToACL();
		public:
			BCondConnector(const Block & b1, const Block & b2);
			void addGhostPoint(AVec<int> a1,AVec<int> a2);
			void addGhostPoint(AVec<int> a1,int d1,AVec<int> a2,int d2);
			const Block & getBlock1();
			const Block & getBlock2();
	};


	/// Virtual class describes general interface for boundary conditions
	/**
		 \ingroup BoundaryConditions

		 The class differs from the class BCond by no use of the 
		 explicite connection point lists. The boundary slice defined by first and last 
		 points.
		 
		 \todo Add posibility to define several slices with different 
		 dimensionalities
	 */
	class BCondSlice:public NumMethod
	{
		protected:
			Block block;
			const VectorTemplate * const templ;
			
			int pointB;
			int pointE;
			AVec<int> sliceDimentions;
			AVec<int> sliceIncrements;
			int direction;
						
		public:
			BCondSlice(const Block & b);
			BCondSlice(const Block & b, const VectorTemplate * const t);
			void addGhostSlice(AVec<int> pB,AVec<int> pE, int dir);
			const Block & getBlock();
	};	
	
	/// Virtual class describes general interface for boundary conditions which connect two datas
	/**
		 \ingroup BoundaryConditions

		 The class differs from the class BCondConnector by no use of the 
		 explicite connection point lists. The boundary defined by fist and last 
		 point numbers and by definition of corresponding increments and 
		 dimensionality of the point set. 
		 
		 \todo Add posibility to define several slices with different 
		 dimensionalities
	 */
	class BCondConnectorSlice:public NumMethod
	{
		protected:
			Block block1;
			Block block2;
			const VectorTemplate * const templ;
			
			int point1B;
			int point1E;
			AVec<int> sliceDimentions1;
			AVec<int> sliceIncrements1;
			int direction1;

			int point2B;
			int point2E;
			AVec<int> sliceDimentions2;
			AVec<int> sliceIncrements2;
			int direction2;
						
		public:
			BCondConnectorSlice(const Block & b1, const Block & b2, const VectorTemplate *const t);
			void addGhostSlice1(AVec<int> pB,AVec<int> pE, int dir);
			void addGhostSlice2(AVec<int> pB,AVec<int> pE, int dir);
			const Block & getBlock1();
			const Block & getBlock2();
	};

	/// Virtual class describes general interface for boundary conditions 
	/**
		 The class is sutable for BC with moving boundary. It takes the 
		 list points and corresponding normals as an input. The boundary conditions
		 are computed automaticaly for a given VectorTemplate.
		  \ingroup BoundaryConditions
	*/
	class BCondDynamic:public NumMethod
	{
		protected:
			Block block;
			const VectorTemplate * const templ;
						
			acl::SPVectorOfElementsData pointsACL;
			acl::SPVectorOfElementsData normalsACL;
			
		public:
			BCondDynamic(const Block & b);
			BCondDynamic(const Block & b, const VectorTemplate * const t);
			const Block & getBlock();
			inline const VectorTemplate * getVT();
	};

	/// defines names of walls of a box \ingroup BoundaryConditions
	enum SlicesNames {X0, XE, Y0, YE, Z0, ZE};
	
	///	\ingroup BoundaryConditions 
	///@{
	void addSliceX0(BCond &);
	void addSliceXE(BCond &);
	void addSliceY0(BCond &);
	void addSliceYE(BCond &);
	void addSliceZ0(BCond &);
	void addSliceZE(BCond &);
	///@}

	void addSlices(BCond &, const vector<SlicesNames> &);	
		
	///adds slice only points without directions \ingroup BoundaryConditions 
	///@{
	void addSliceX(BCond &, int x);
	void addSliceY(BCond &, int y);
	void addSliceZ(BCond &, int z);
	///@}
		
	///	\ingroup BoundaryConditions 
	///@{
	void addSliceX0(BCondSlice &);
	void addSliceXE(BCondSlice &);
	void addSliceY0(BCondSlice &);
	void addSliceYE(BCondSlice &);
	void addSliceZ0(BCondSlice &);
	void addSliceZE(BCondSlice &);
	///@}

// --------------------------- Implementation ----------------------
	inline const VectorTemplate * BCond::getVT()
	{
		return templ;
	}
		
	inline const VectorTemplate * BCondWithMap::getVT()
	{
		return templ;
	}
	
}

#endif //ASLBCOND_H
