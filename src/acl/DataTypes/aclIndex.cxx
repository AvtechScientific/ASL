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


#include "aclIndex.h"
#include "../../aslUtilities.h"

namespace acl
{

	Index::Index(unsigned int s):
		ElementBase(false, s, TYPE_UINT)
	{
	}


	std::string Index::getName() const
	{
		return "";
	}


	std::string Index::str(const KernelConfiguration & kernelConfig) const
	{
		return string(INDEX);
	}


	void Index::addToKernelSource(vector<Element> & arguments, vector<Element> & localDeclarations) const
	{
	}


	void Index::setAsArgument(cl::Kernel & kernel, unsigned int argumentIndex) const
	{
	}


	string Index::getTypeSignature(const KernelConfiguration & kernelConfig) const
	{
		return "";
	}


	string Index::getLocalDeclaration(const KernelConfiguration & kernelConfig) const
	{
		return "";
	}


} // namespace acl
