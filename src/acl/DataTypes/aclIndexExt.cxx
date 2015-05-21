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


#include "aclIndexExt.h"
#include "../../aslUtilities.h"

namespace acl
{


	IndexExt::IndexExt(unsigned int s):
		ElementBase(false, s, TYPE_UINT)
	{
	}


	std::string IndexExt::getName() const
	{
		return "";
	}


	std::string IndexExt::str(const KernelConfiguration & kernelConfig) const
	{
		if (kernelConfig.vectorWidth == 1)
			return string(INDEX);

		string s("(uint"+asl::numToStr(kernelConfig.vectorWidth)+")( 0u");
		for(unsigned int i(1); i < kernelConfig.vectorWidth; ++i)
			s+= ", " + asl::numToStr(i) + "u";
		s+=" )";

		if (kernelConfig.unaligned)
			return "(" + INDEX + "+" + s + ")"; 
		else
			return "(" + asl::numToStr(kernelConfig.vectorWidth) + "*" + INDEX + "+" + s + ")"; 			
	}


	void IndexExt::addToKernelSource(vector<Element> & arguments, vector<Element> & localDeclarations) const
	{
	}


	void IndexExt::setAsArgument(cl::Kernel & kernel, unsigned int argumentIndex) const
	{
	}


	string IndexExt::getTypeSignature(const KernelConfiguration & kernelConfig) const
	{
		return "";
	}


	string IndexExt::getLocalDeclaration(const KernelConfiguration & kernelConfig) const
	{
		return "";
	}


} // namespace acl
