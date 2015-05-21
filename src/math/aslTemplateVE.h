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


#ifndef TEMPLATEVE_H_INCLUDED
#define TEMPLATEVE_H_INCLUDED

#include <acl/aclMath/aclVectorOfElements.h>
#include "data/aslDataWithGhostNodes.h"
#include "aslTemplatesExtras.h"
#include "aslDistanceFunction.h"

/**
 \defgroup DifferentialOpperators Differential Operrators
 */


namespace asl
{
	/// This class contains VectorOfElements which corresponds to values of a field in littice nodes
	/**
		 This class contains VectorOfElements (VE) with PrivateVariables. 
		 This list should be added to the kernel before use as follows:
			
		 \code 
			 Kernel kernel;
			 TemplateVE a;

			 kernel<<a.initValues;
		 \endcode
	*/
	class TemplateVE
	{
		public:
			const VectorTemplate* vectorTemplate;
			const VTObjects* vto;
			acl::VectorOfElements values;
			acl::VectorOfElements initValues;
			/**
			 \param data is a data field
			 \param vectorT is a tempate vector
			 \param i is the number of component (for a case with multicomponent data)
			 \param bIni incase of true it generates private variable for storage otherwise it uses direct acces  
			 */
			TemplateVE(AbstractDataWithGhostNodes & data, 
			           const VectorTemplate & vectorT, 
			           unsigned int i=0,
			           bool bIni=true);
			/**
			 \param data is a data field
			 \param position contains current position (not index)
			 \param vectorT is a tempate vector
			 \param i is the number of component (for a case with multicomponent data)
			 */
			TemplateVE(DistanceFunction & data,
			           acl::VectorOfElements & position,
			           const VectorTemplate & vectorT, 
			           unsigned int i=0);
			///this is used for matematical operations with TemplateVE
			TemplateVE(const acl::VectorOfElements & val, 
			           const VectorTemplate & vectorT);
			
			TemplateVE();
			acl::VectorOfElements getValue(unsigned int i);
			/** 
			 Note: the init function does not make initialization of values in case
			 whete its length remains the same. this is usefull for reuse of private variables

			 \param data is a data field
			 \param vectorT is a tempate vector
			 \param i is the number of component (for a case with multicomponent data)
			 \param bIni incase of true it generates private variable for storage otherwise it uses direct acces  
			 */
			void init(AbstractDataWithGhostNodes & data, 
			          const VectorTemplate & vectorT, 
			          unsigned int i=0,
			          bool bIni=true);

			/** 
			 Note: the init function does not make initialization of values in case
			 whete its length remains the same. this is usefull for reuse of private variables
			 \todo Fnish!!! not distance functino but Position function
			 */
			void init(DistanceFunction & data,
			          acl::VectorOfElements & position,
			          const VectorTemplate & vectorT, 
			          unsigned int i=0);
	}; // class TemplateVE

	TemplateVE operator+ (const TemplateVE &a, const TemplateVE &b);
	TemplateVE operator- (const TemplateVE &a, const TemplateVE &b);	
	TemplateVE operator* (const TemplateVE &a, const TemplateVE &b);	
	TemplateVE operator/ (const TemplateVE &a, const TemplateVE &b);		
	
	/// differential operator \f$ \Delta a \f$ \ingroup DifferentialOpperators
	acl::VectorOfElements laplas(const TemplateVE & a);

	/// differential operator \f$ \vec\nabla a \f$ \ingroup DifferentialOpperators
	acl::VectorOfElements gradient(const TemplateVE & a);	
	/// differential operator \f$ \partial_x a \f$ \ingroup DifferentialOpperators
	acl::VectorOfElements dx(const TemplateVE & a);
	/// differential operator \f$ \partial_y a \f$ \ingroup DifferentialOpperators
	acl::VectorOfElements dy(const TemplateVE & a);
	/// differential operator \f$ \partial_z a \f$ \ingroup DifferentialOpperators
	acl::VectorOfElements dz(const TemplateVE & a);	
	/// differential operator \f$ \partial_{i} a \f$ \ingroup DifferentialOpperators
	acl::VectorOfElements dxi(const TemplateVE & a, unsigned int i);
	
	/// differential operator \f$ \nabla(a \nabla b) \f$ \ingroup DifferentialOpperators
	acl::VectorOfElements divAgradB(const TemplateVE & a, const TemplateVE & b);	

	/// differential operator \f$ \nabla \cdot \vec a \f$ \ingroup DifferentialOpperators
	acl::VectorOfElements div(const vector<TemplateVE> & a);	
	/// differential operator \f$ \nabla \cdot (\vec a c) \f$ \ingroup DifferentialOpperators
	/**
		 This is a discret approximation discribef by the following expression:
		 \f[ \nabla_\alpha (v_\alpha c) = \frac{1}{2}\sum_i w_i a_{i\alpha} (v_{0\alpha}+v_{i\alpha})(c_0+c_i) ) \f]
		where \f$ w_i \f$ are coefficients corresponding to a gradient ones
	*/
	acl::VectorOfElements divProduct(const vector<TemplateVE> & a, const TemplateVE & c);	
	/// differential operator \f$ \nabla \cdot \vec a \f$ \ingroup DifferentialOpperators
	acl::VectorOfElements div(const TemplateVE & ax, const TemplateVE & ay);	
	/// differential operator \f$ \nabla \cdot \vec a \f$ \ingroup DifferentialOpperators
	acl::VectorOfElements div(const TemplateVE & ax, const TemplateVE & ay, const TemplateVE & az);	

	
	/// differential operator \f$ \nabla_i \nabla_j a \f$ \ingroup DifferentialOpperators
	acl::VectorOfElements dIdJ(unsigned int i, unsigned int j, const TemplateVE & a);

	///  generates expresion for bilinear interpolation the template should be an elementary cell 
	acl::VectorOfElements interpolate(const TemplateVE & a, acl::VectorOfElements e);

	
}// asl

#endif // TEMPLATEVE_H_INCLUDED
