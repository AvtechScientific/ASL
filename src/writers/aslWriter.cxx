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


#include "aslWriter.h"
#include <acl/aclMath/aclVectorOfElementsOperations.h>
#include <acl/Kernels/aclKernel.h>
#include <acl/acl.h>
#include <acl/DataTypes/aclArray.h>
#include <aslGenerators.h>
#include <acl/aclGenerators.h>
#include <data/aslDataWrapper.h>
#include <data/aslDataWithGhostNodes.h>

using namespace acl;

namespace asl
{

	Writer * Writer::current(NULL);

	Writer::Writer(const string & file_):
		file(file_),
		numOfWrites(0)
	{
		enable();
	}


	Writer::~Writer()
	{
		// Deactivates this instance of Writer
		if (current == this)
			current = NULL;
	}


	void Writer::enable()
	{
		// Activates this instance of Writer
		current = this;
	}


	void Writer::addScalars(string name, AbstractData & data)
	{
		if ((scalarFields.size() == 0) && (vectorFields.size() == 0))
			block= make_shared<Block>(data.getBlock());

		// check whether the new block is compatible with the old one
		if (block->getSize() != data.getBlock().getSize())
			errorMessage("Writer::addScalars() - attempt to add AbstractData with incompatible block size");

		for (unsigned int i = 0; i < data.getDContainer().size() ; ++i)
		{
			scalarFields.push_back(make_pair(name + "-" + numToStr(i),
			                                 subVE(data.getDContainer(), i, i)));
		}
	}


	void Writer::addVector(string name, AbstractData & data)
	{
		if ((scalarFields.size() == 0) && (vectorFields.size() == 0))
			block= make_shared<Block>(data.getBlock());

		// check whether the new block is compatible with the old one
		if (block->getSize() != data.getBlock().getSize())
			errorMessage("Writer::addVector() - attempt to add AbstractData with incompatible block size");

		vectorFields.push_back(make_pair(name, data.getDContainer()));
	}


	void Writer::addScalars(string name, VectorOfElementsData & data)
	{
		if ((scalarFields.size() == 0) && (vectorFields.size() == 0))
			errorMessage("Writer::addScalars() - attempt to add VectorOfElementsData before any Block was defined");

		// check whether the new data is compatible with the old block
		if (!compatibleSizes((unsigned int)productOfElements(block->getSize()), data))
			errorMessage("Writer::addScalars() - attempt to add VectorOfElementsData with incompatible block size");

		for (unsigned int i = 0; i < data.size() ; ++i)
		{
			scalarFields.push_back(make_pair(name + "-" + numToStr(i),
			                                 subVE(data, i, i)));
		}
	}


	void Writer::addVector(string name, VectorOfElementsData & data)
	{
		if ((scalarFields.size() == 0) && (vectorFields.size() == 0))
			errorMessage("Writer::addVector() - attempt to add VectorOfElementsData before any Block was defined");

		// check whether the new data is compatible with the old block
		if (!compatibleSizes((unsigned int)productOfElements(block->getSize()), data))
			errorMessage("Writer::addVector() - attempt to add VectorOfElementsData with incompatible block size");

		vectorFields.push_back(make_pair(name, data));
	}


	void Writer::addScalars(string name, 
	                        const VectorOfElements & data, 
	                        Kernel & kernel, 
	                        unsigned int nGhost)
	{
		if ((scalarFields.size() == 0) && (vectorFields.size() == 0))
			errorMessage("Writer::addScalars() - attempt to add VectorOfElements before any Block was defined");

		if (kernel.getQueue().get() == 0)
			errorMessage("Writer::addScalars() - attempt to add VectorOfElements before any Queue was defined in Kernel");

		if (data.size() == 0)
			errorMessage("Writer::addScalars() - attempt to add VectorOfElements with size 0");

		auto queue(kernel.getQueue());
		auto v(generateDataContainerACL_SP(offset(*block,-nGhost), 
			                                   acl::getElementType(data),
	                                           data.size(), 
	                                           nGhost, 
	                                           queue));
		acl::initData(v->getEContainer(), generateVEConstantN(data.size(), 0.));
		       
		for (unsigned int i = 0; i < data.size() ; ++i)
		{
			scalarFields.push_back(make_pair(name + "-" + numToStr(i),
			                                 subVE(v->getDContainer(), i, i)));
		}
		
		kernel << assignmentSafe(v->getSubContainer(), data);
	}


	void Writer::addVector(string name, 
	                       const VectorOfElements & data, 
	                       Kernel & kernel, 
	                       unsigned int nGhost)
	{
		if ((scalarFields.size() == 0) && (vectorFields.size() == 0))
			errorMessage("Writer::addVector() - attempt to add VectorOfElements before any Block was defined");

		if (kernel.getQueue().get() == 0)
			errorMessage("Writer::addVector() - attempt to add VectorOfElements before any Queue was defined in Kernel");

		if (data.size() == 0)
			errorMessage("Writer::addVector() - attempt to add VectorOfElements with size 0");
		
		auto queue(kernel.getQueue());
		auto v(generateDataContainerACL_SP(offset(*block,-nGhost), 
			                                   acl::getElementType(data),
	                                           data.size(), 
	                                           nGhost, 
	                                           queue));
		acl::initData(v->getEContainer(), generateVEConstantN(data.size(), 0.));

		vectorFields.push_back(make_pair(name, v->getDContainer()));

		kernel << assignmentSafe(v->getSubContainer(), data);
	}


} // namespace asl
