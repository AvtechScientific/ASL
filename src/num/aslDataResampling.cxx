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


#include "aslDataResampling.h"
#include <aslGenerators.h>
#include <acl/aclGenerators.h>
#include <acl/acl.h>
#include <math/aslTemplateVE.h>

#include <data/aslDataWithGhostNodes.h>
#include <acl/Kernels/aclKernel.h>
#include <math/aslTemplates.h>
#include <math/aslIndex2Position.h>
#include <acl/aclMath/aclVectorOfElements.h>

#include <acl/Kernels/aclKernelConfigurationTemplates.h>

namespace asl
{

	Block generateCoarsedBlock(const Block & b)
	{
		unsigned int n(nD(b));
		AVec<int> s(b.getSize()/2-Block::DV(n,1));
		double dx(b.dx*2);
		AVec<> pos(b.position+Block::V(n,.75*dx));
		return {s,dx,pos};
	}
	
	DataCoarser::DataCoarser():
		SingleKernelNM(acl::KERNEL_BASIC),
		dataIn(0),
		dataOut(0),
		vectorTemplate(NULL)		
	{
	}


	DataCoarser::DataCoarser(Data dIn):
		SingleKernelNM(acl::KERNEL_BASIC),
		dataIn(dIn),
		dataOut(0),
		vectorTemplate(elementaryCellVT(nD(dIn->getBlock())))
	{
	}

		
	
	void DataCoarser::init0()
	{
		
		acl::TypeID type(getElementType(dataIn->getDContainer()));
		unsigned int gN(dataIn->getGhostBorder());		
		Block bOut(generateCoarsedBlock(dataIn->getInternalBlock()));
		acl::VectorOfElementsData veIn(dataIn->getDContainer());
		unsigned int nc(veIn.size());
		
		dataOut=generateDataContainerACL_SP(bOut,type,nc,gN);
		acl::VectorOfElementsData veOut(dataOut->getDContainer());

		acl::VectorOfElements v(acl::generateVEPrivateVariable(nc, type));

		auto dataIn0(resizeGhostNodes(dataIn,0)); 
			
		TemplateVE tVE(*dataIn0, *vectorTemplate);

		Index2PositionDiscreteACL t2p(dataOut->getBlock());
		(*kernel) << t2p.initPosition;
		acl::VectorOfElements fineIndex(t2p.position * 2 *
		                                dataIn0->getBlock().c2iTransformVector);
		
		(*kernel) << (acl::excerpt(tVE.initValues,
		                           fineIndex));
		/// \todoZeev add support for unsigned long type in Hardware
		/// and then remove (cl_uint) casting
		(*kernel) << (veOut = acl::sumOfElements(tVE.values)/(cl_uint)tVE.values.size());	
	}


	Block generateClippedBlock(const Block & b, 
	                           const AVec<int> & a0, 
	                           const AVec<int> & aE)
	{
		unsigned int nd(nD(a0));
		return {aE-a0+AVec<int>(nd,1), b.dx, b.position+b.dx*Block::V(a0)};
	}
	
	DataClipper::DataClipper():
		SingleKernelNM(acl::KERNEL_BASIC),
		dataIn(0),
		dataOut(0)		
	{
	}


	DataClipper::DataClipper(Data dIn, AVec<int> a0, AVec<int> aE):
		SingleKernelNM(acl::KERNEL_BASIC),
		dataIn(dIn),
		dataOut(0),
		a0(a0),
		aE(aE)
	{
	}
	
	void DataClipper::init0()
	{
		acl::TypeID type(getElementType(dataIn->getDContainer()));
		unsigned int gN(dataIn->getGhostBorder());		
		Block bOut(generateClippedBlock(dataIn->getInternalBlock(), a0, aE));
		acl::VectorOfElementsData veIn(dataIn->getDContainer());
		unsigned int nc(veIn.size());
		
		dataOut=generateDataContainerACL_SP(bOut,type,nc,gN);
		acl::VectorOfElementsData veOut(dataOut->getDContainer());
			
		Index2PositionDiscreteACL t2p(dataOut->getBlock());
		(*kernel) << t2p.initPosition;
		acl::VectorOfElements index((t2p.position + a0)*
		                            dataIn->getBlock().c2iTransformVector);
		
		(*kernel) << (veOut = acl::excerpt(veIn, index));	
	}
		
		
} //asl
