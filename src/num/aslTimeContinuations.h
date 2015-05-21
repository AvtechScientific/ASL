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


#ifndef ASLTIMECONTINUATIONS_H
#define ASLTIMECONTINUATIONS_H

#include "aslNumMethod.h"
#include <acl/aclMath/aclVectorOfElementsDef.h>

namespace acl
{
	class Kernel;
}

namespace asl
{
	class VectorTemplate;
	template <typename V> class DataWithGhostNodes;
	typedef DataWithGhostNodes<acl::VectorOfElementsData> DataWithGhostNodesACLData;
	typedef std::shared_ptr<DataWithGhostNodesACLData> SPDataWithGhostNodesACLData;
	class AbstractDataWithGhostNodes;
	typedef std::shared_ptr<AbstractDataWithGhostNodes> SPAbstractDataWithGhostNodes;

	/// Numerical method that generates temporal extrapolation of the data, Abstract class.
	/**
		 \ingroup NumMethods

	*/
	class TimeContinuations: public NumMethod
	{
		public:
			typedef SPDataWithGhostNodesACLData Data;
		protected:	
			acl::VectorOfElementsData inData;
			double factor;
			unsigned int nStorages;			
			TimeContinuations(Data inD, double factor);
			TimeContinuations(acl::VectorOfElementsData & inD, double factor);
		public:			
			void addData(Data inD);
			void addData(acl::VectorOfElementsData & inD);
			virtual void execute()=0;
			virtual void init()=0;
			/// makes reset of the contiuation (storage) cicle 
			void reset();
	};

	
	/// Numerical method that generates temporal extrapolation of the data with Lagrangian polynoms 
	/**
		 \ingroup NumMethods

		 The method computes Lagrange polinomial extrapolation of order \f$k\f$ in time. \p factor 
		 defines extrapolation length:
		 \f[ u(t+dt*factor) = \sum_{i=0}^n y_i  l_i \f]
		 where 
		 \f[ l_i(x) = \prod_{j=0, j\neq i}^k \frac{x-x_j}{x_i-x_j}\f]
		 The first avalible point has \f$ t_0 = -k$. The last avalible point has \f$ t_k = 0$.
	*/
	class TimeContinPLagrange: public TimeContinuations
	{
		private:	
			acl::VectorOfElements storedData;
			std::vector<std::shared_ptr<acl::Kernel>> kernels;
			unsigned int order;
			std::vector<double> coefs; 

		public:			
			TimeContinPLagrange(Data inD, double f, unsigned int order);
			TimeContinPLagrange(acl::VectorOfElementsData & inD, double f, unsigned int order);
			virtual void execute();
			virtual void init();
	};

	typedef std::shared_ptr<TimeContinPLagrange> SPTimeContinPLagrange;

	/// Numerical method that generates temporal extrapolation of the data with Lagrangian polynoms of fractional argument  
	/**
		 \ingroup NumMethods

		 The method computes Lagrange polinomial extrapolation of order \f$k\f$ in time. 
		 The time is taken in the form \f$ (t+t_s)^{-1}\f$. \p factor  
		 defines extrapolation length. \f$ t_s\f$ is defined as: 
	     \f[ t_s = 2 k+ factor \f]
		 The interpolation polinoms are nothing else but the Lagrange one with
		 \f$ x \f$ defined as \f$ x=(t+t_s)^{-1} \f$
		 \f[ u(t+dt*factor) = \sum_{i=0}^n y_i  l_i \f]
		 where 
		 \f[ l_i(x) = \prod_{j=0, j\neq i}^k \frac{x-x_j}{x_i-x_j}\f]

		 The first avalible point has \f$ t_0 = -k$. The last avalible point has \f$ t_k = 0$. 
		 
	*/	
	class TimeContinPLagrangeFraction: public TimeContinuations
	{
		public:
			typedef SPDataWithGhostNodesACLData Data;
			typedef SPAbstractDataWithGhostNodes Field;
		private:	
			acl::VectorOfElements storedData;
			std::vector<std::shared_ptr<acl::Kernel>> kernels;
			unsigned int order;
			double offset;
			std::vector<double> coefs; 

		public:			
			TimeContinPLagrangeFraction(Data inD, double f, unsigned int order);
			TimeContinPLagrangeFraction(acl::VectorOfElementsData & inD, 
			                            double f, unsigned int order);
			void execute();
			virtual void init();
	};

	typedef std::shared_ptr<TimeContinPLagrangeFraction> SPTimeContinPLagrangeFraction;
	
	
} // asl
#endif // ASLTIMECONTINUATIONS_H
