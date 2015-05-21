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


#ifndef ACLCONSTANT_H
#define ACLCONSTANT_H

#include "../aclElementBase.h"
#include "../aclUtilities.h"

namespace acl
{
	
	template <typename T> class Constant: public ElementBase
	{
		private:
			T value;
			string valueStr;
		public:
			explicit Constant(T v);
			virtual string str(const KernelConfiguration & kernelConfig) const;
			virtual string getName() const;
			virtual string getTypeSignature(const KernelConfiguration & kernelConfig) const;
			virtual string getLocalDeclaration(const KernelConfiguration & kernelConfig) const;
			virtual void addToKernelSource(vector<Element> & arguments,
			                               vector<Element> & localDeclarations) const;
			virtual void setAsArgument(cl::Kernel & kernel, unsigned int argumentIndex) const;
	};




	template <typename T> std::string Constant<T>::getName() const
	{
		return valueStr;
	}


	template <typename T> std::string Constant<T>::str(const KernelConfiguration & kernelConfig) const
	{
		return valueStr;
	}


	template <typename T> void Constant<T>::addToKernelSource(vector<Element> & arguments,
	                                                          vector<Element> & localDeclarations) const
	{
	}
		
		

	template <typename T> void Constant<T>::setAsArgument(cl::Kernel & kernel,
	                                                      unsigned int argumentIndex) const
	{
	}


	template <typename T> string Constant<T>::getTypeSignature(const KernelConfiguration & kernelConfig) const
	{
		return "";
	}


	template <typename T> string Constant<T>::getLocalDeclaration(const KernelConfiguration & kernelConfig) const
	{
		return "";
	}


} // namespace acl

#endif // ACLCONSTANT_H
