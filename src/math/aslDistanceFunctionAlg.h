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


#ifndef ASLDISTANCEFUNCTIONALG_H
#define ASLDISTANCEFUNCTIONALG_H

#include <acl/aclMath/aclVectorOfElementsDef.h>
#include <aslUtilities.h>
 
namespace asl {

//	template <typename T> class AVec;
//	class Block;

//	class AbstractDataWithGhostNodes; 
//	typedef std::shared_ptr<AbstractDataWithGhostNodes> SPAbstractDataWithGhostNodes;

	/// returns expression corresponding to check if the node in i^th direction is ghost one
	acl::VectorOfElements isGhostNode(TemplateVE & distanceTVE, unsigned int i);

	/// returns expression corresponding to check if the node in i^th direction is computation one
	acl::VectorOfElements isComputationNode(TemplateVE & distanceTVE, unsigned int i);
	/// returns expression corresponding to number of ghost nodes in a cell \p i
	acl::VectorOfElements nGhostNodesInCell(TemplateVE & distanceTVE, unsigned int i);
	/// returns expression corresponding to check if there is a boundary between nodes i^th and j^th  within cell \p iE
	acl::VectorOfElements isBoundaryBetween(TemplateVE & distanceTVE, 
	                                        unsigned int iE, 
	                                        unsigned int i, 
	                                        unsigned int j);

	/// returns expression corresponding to the relative boundary position in direction \p i
	/**
		\f[ x= v_0/(v_0-v_i) \f]
	*/
	acl::VectorOfElements exBoundaryX(TemplateVE & distanceTVE, unsigned int i);
	
	
	/// generates expresion for center of a boundary element
	/**
		 \param iEl the element number
		 computes avarage point of corners of the boundary poligon within the element 
	*/
	acl::VectorOfElements exBoundaryCenter(TemplateVE & distanceTVE, unsigned int iEl);

	/// generates expresion for area of a boundary element
	/**
		 \param iEl the element number
		 computes area of the corresponding boundary within bulk element \p iEl
	*/
	vector<acl::Element> gcBoundaryArea(TemplateVE & distanceTVE,
	                                    unsigned int iEl, 
	                                    acl::VectorOfElements & center, 
	                                    acl::VectorOfElements & area);

	/// generates expresion for area of a boundary element
	/**
		 \param iEl the element number
		 computes area of the corresponding boundary within cell iEl and divide on number of 
		 ghost points within the cell
	*/
	vector<acl::Element> gcBoundaryAreaPerGhostPoint(TemplateVE & distanceTVE,
	                                                 unsigned int iEl, 
	                                                 acl::VectorOfElements & center, 
	                                                 acl::VectorOfElements & area);
	
	/// generates expresion for area of a boundary element
	/**
		 computes area of the corresponding boundary within all cells
	*/
	vector<acl::Element> gcBoundaryArea(TemplateVE & distanceTVE,
	                                    acl::VectorOfElements & area);			

	/// generates expresion for area of a boundary element
	/**
		 computes area of the corresponding boundary per ghost point within all cells
	*/
	vector<acl::Element> gcBoundaryAreaPerGhostPoint(TemplateVE & distanceTVE,
	                                                 acl::VectorOfElements & area);			
	
} // asl

#endif // ASLDISTANCEFUNCTION
