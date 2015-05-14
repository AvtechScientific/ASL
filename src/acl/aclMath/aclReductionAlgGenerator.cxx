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


#include "aclReductionAlgGenerator.h"
#include "aclVectorOfElementsOperations.h"
#include <acl/acl.h>
#include <acl/DataTypes/aclArray.h>
#include <acl/aclGenerators.h>
#include <acl/Kernels/aclKernel.h>

namespace acl
{
	
	inline const unsigned int nLocalGPU(unsigned int length, unsigned int nGroups)
	{
		return std::max(1u,std::min(8u, length/nGroups));
	}
	inline const unsigned int nLocalCPU(unsigned int length, unsigned int nGroups)
	{
		return 1;
	}

	unsigned int getLPerUnit(unsigned int length, unsigned int nUnits)
	{
		return ceil((double)length/nUnits);
	}

	unsigned int getNSaturatedUnits(unsigned int length, unsigned int nUnits)
	{
		return length/getLPerUnit(length,nUnits);
	}
		
	unsigned int getLLastUnit(unsigned int length, unsigned int nUnits)
	{
		return length % getLPerUnit(length,nUnits);
	}
	
	template <typename ResType, enum ReductionOperatorType Op> class GroupReduction
	{
	};

	template <typename ResType> class GroupReduction<ResType, ROT_SUM>
	{
		public:
			static void compute(std::vector<std::vector<ResType>> v, unsigned int l, asl::AVec<ResType> & res)
		    {
				for(unsigned int i(0); i<v.size(); ++i)
				{
					vector<ResType> & vv(v[i]);
					ResType r(vv[0]);
					for(unsigned int j(1); j<l; ++j)
						r+=vv[j];
					res[i]=r;
				}				
			}
			static inline VectorOfElements op(VectorOfElements res, VectorOfElements v, acl::TypeID t)
			{
				return res+=v;
			}
	};
//	template void GroupReduction<double, ROT_SUM>::compute(std::vector<std::vector<double>> v, 
//	                                                      asl::AVec<double> & res); 
//	template void GroupReduction<float, ROT_SUM>::compute(std::vector<std::vector<float>> v, 
//	                                                     asl::AVec<float> & res); 
	
	template <typename ResType> class GroupReduction<ResType, ROT_PRODUCT>
	{
		public:
			static void compute(std::vector<std::vector<ResType>> v, unsigned int l, asl::AVec<ResType> & res)
		    {
				for(unsigned int i(0); i<v.size(); ++i)
				{
					vector<ResType> & vv(v[i]);
					ResType r(vv[0]);
					for(unsigned int j(1); j<l; ++j)
						r*=vv[j];
					res[i]=r;
				}				
			}
			static inline VectorOfElements op(VectorOfElements res, VectorOfElements v, acl::TypeID t)
			{
				return res*=v;
			}
	};
//	template void GroupReduction<double, ROT_PRODUCT>::compute(std::vector<std::vector<double>> v, 
//	                                                          asl::AVec<double> & res); 
//	template void GroupReduction<float, ROT_PRODUCT>::compute(std::vector<std::vector<float>> v, 
//	                                                         asl::AVec<float> & res); 
	
	template <typename ResType> class GroupReduction<ResType, ROT_MINIMUM>
	{
		public:
			static void compute(std::vector<std::vector<ResType>> v, unsigned int l, asl::AVec<ResType> & res)
		    {
				for(unsigned int i(0); i<v.size(); ++i)
				{
					vector<ResType> & vv(v[i]);
					ResType r(vv[0]);
					for(unsigned int j(1); j<l; ++j)
						r=std::min(r,vv[j]);
					res[i]=r;
				}				
			}
			static inline VectorOfElements op(VectorOfElements res, VectorOfElements v, acl::TypeID t)
			{
				return res=min(res,v,t);
			}
	};
//	template void GroupReduction<double, ROT_MINIMUM>::compute(std::vector<std::vector<double>> v, 
//	                                                          asl::AVec<double> & res); 
//	template void GroupReduction<float, ROT_MINIMUM>::compute(std::vector<std::vector<float>> v, 
//	                                                         asl::AVec<float> & res); 

	template <typename ResType> class GroupReduction<ResType, ROT_MAXIMUM>
	{
		public:
			static void compute(std::vector<std::vector<ResType>> v, unsigned int l, asl::AVec<ResType> & res)
		    {
				for(unsigned int i(0); i<v.size(); ++i)
				{
					vector<ResType> & vv(v[i]);
					ResType r(vv[0]);
					for(unsigned int j(1); j<l; ++j)
						r=std::max(r,vv[j]);
					res[i]=r;
				}				
			}
			static inline VectorOfElements op(VectorOfElements res, VectorOfElements v, acl::TypeID t)
			{
				return res=max(res,v,t);
			}			
	};
//	template void GroupReduction<double, ROT_MAXIMUM>::compute(std::vector<std::vector<double>> v, 
//	                                                          asl::AVec<double> & res); 
//	template void GroupReduction<float, ROT_MAXIMUM>::compute(std::vector<std::vector<float>> v, 
//	                                                         asl::AVec<float> & res); 

	template <typename ResType, enum ReductionOperatorType Op> 
		void ReductionAlgGenerator<ResType, Op>::compute()
	{
		if(kernel) 
			kernel->compute();
		unsigned int nC(ve.size());
		unsigned int veLength(ve[0]->getSize());
		unsigned int nSaturatedUnits(getNSaturatedUnits(veLength, nGroups*nLocal));
		unsigned int l(std::min(nSaturatedUnits+1,nGroups*nLocal));
		for (unsigned int i(0); i < nC; ++i)
			acl::copy(groupResACL[i], groupRes[i]);
		GroupReduction<ResType, Op>::compute(groupRes, l, res.v());
	};	
	template void ReductionAlgGenerator<double, ROT_SUM>::compute();
	template void ReductionAlgGenerator<float, ROT_SUM>::compute();
	template void ReductionAlgGenerator<double, ROT_PRODUCT>::compute();
	template void ReductionAlgGenerator<float, ROT_PRODUCT>::compute();
	template void ReductionAlgGenerator<double, ROT_MINIMUM>::compute();
	template void ReductionAlgGenerator<float, ROT_MINIMUM>::compute();
	template void ReductionAlgGenerator<double, ROT_MAXIMUM>::compute();
	template void ReductionAlgGenerator<float, ROT_MAXIMUM>::compute();
	
	template <typename ResType, enum ReductionOperatorType Op> 
		void generateKernelCPU(VectorOfElements ve, VectorOfElements groupResACL, Kernel & k)
	{		
		auto nGroups(k.getGroupsNumber());
		int veLength(ve[0]->getSize());
		unsigned int nLocal(nLocalCPU(veLength,nGroups));	
		unsigned int gSize(std::max(k.getSize(),nLocal));
		auto type(getElementType(ve));
		auto typeI(TYPE_SELECT[type]);
		auto lVal(acl::generateVEPrivateVariable(ve.size(),type));
		auto counter(acl::generateVEPrivateVariable(1,typeI));
		auto counterEnd(acl::generateVEPrivateVariable(1,typeI));
		int lPerGroup(getLPerUnit(veLength,nGroups));
		int lLastGroup(getLLastUnit(veLength, nGroups));
		int nSaturatedGroups(getNSaturatedUnits(veLength, nGroups));

		k<<(counterEnd = select(generateVEConstant(lPerGroup), 
		                        generateVEConstant(lLastGroup), 
		                        generateVEGroupID()==nSaturatedGroups, 
		                        type));
		k<<(lVal = select(excerpt(ve,lPerGroup*generateVEGroupID()), counterEnd>0,type));
		// this line ensuers that the conputations will take place on the 1st item only
		k<<(counterEnd = select(counterEnd, generateVEIndex(gSize)==0, type));
		vector<Element> body;
		body<<(GroupReduction<ResType, Op>::
			               op(lVal,excerpt(ve, 
			                               int(lPerGroup)*generateVEGroupID() + 
			                               counter),type));

		auto loop(elementOperators::forLoop((counter = generateVEConstant(1))[0],
		                                    (counter < counterEnd)[0],
		                                    (counter += generateVEConstant(1))[0],
		                                    body));
		k.addExpression(loop);
		k<<(excerpt(groupResACL,generateVEGroupID()) = lVal);
	}

	template <typename ResType, enum ReductionOperatorType Op> 
		void generateKernelGPU(VectorOfElements ve, VectorOfElements groupResACL, Kernel & k)
	{		
		auto nGroups(k.getGroupsNumber());
		int veLength(ve[0]->getSize());
		unsigned int nLocal(nLocalGPU(veLength,nGroups));
		unsigned int gSize(std::max(k.getSize(),nLocal));
		auto type(getElementType(ve));
		auto typeI(TYPE_SELECT[type]);
		auto lVal(acl::generateVEPrivateVariable(ve.size(),type));
		auto counter(acl::generateVEPrivateVariable(1,typeI));
		auto counterEnd(acl::generateVEPrivateVariable(1,typeI));
		auto nUnits(nGroups*nLocal);
		int lPerUnit(getLPerUnit(veLength,nUnits));
		int lLastUnit(getLLastUnit(veLength, nUnits));
		unsigned int nSaturatedUnits(getNSaturatedUnits(veLength, nUnits));
		auto unitID(generateVEGroupID()*nLocal +generateVEIndex());

		if(nSaturatedUnits<nUnits)
			k<<(counterEnd = select(generateVEConstant(lPerUnit), 
			                        generateVEConstant(lLastUnit), 
			                        unitID==nSaturatedUnits, 
			                        typeI));
		if((nSaturatedUnits+1)<nUnits)
			k<<(counterEnd = select(counterEnd, 
			                        unitID<=nSaturatedUnits, 
			                        typeI));
			
		// this line ensuers that the conputations will take place on first item only
		k<<(counterEnd = select(counterEnd, generateVEIndex(gSize)<nLocal, type));
		k<<(lVal = select(excerpt(ve,lPerUnit*unitID), counterEnd>0,type));
		vector<Element> body;
		body<<(GroupReduction<ResType, Op>::
			               op(lVal,excerpt(ve, lPerUnit*unitID + counter),type));

		auto loop(elementOperators::forLoop((counter = generateVEConstant(1))[0],
		                                    (counter < counterEnd)[0],
		                                    (counter += generateVEConstant(1))[0],
		                                    body));
		k.addExpression(loop);
		k<<(excerpt(groupResACL,unitID) = lVal);
	}
	
	template <typename ResType, enum ReductionOperatorType Op> 
		void ReductionAlgGenerator<ResType, Op>::generateAlg(Kernel & k)
	{
		if(!k.getConfiguration().local)
			asl::errorMessage("ReductionAlgGenerator::generateAlg: The kernel should be local");

		nGroups= k.getGroupsNumber();
		unsigned int veLength(ve[0]->getSize());
		nLocal=(getDeviceType(ve[0]->getQueue())== CL_DEVICE_TYPE_CPU ? 
		        nLocalCPU(veLength,nGroups) : nLocalGPU(veLength,nGroups));
		for(unsigned int i(0); i<ve.size(); ++i)
			groupRes[i]=std::vector<ResType>(nGroups*nLocal);
		copy(generateVEData<ResType>(nGroups*nLocal,ve.size()),groupResACL);

		
		switch (getDeviceType(ve[0]->getQueue()))
		{
			case CL_DEVICE_TYPE_CPU:  
				generateKernelCPU<ResType, Op>(ve,groupResACL,k);
				k.setup();
			break;
			case CL_DEVICE_TYPE_GPU:  
				generateKernelGPU<ResType, Op>(ve,groupResACL,k);
				k.setup();
			break;
			default: errorMessage("ReductionAlgGenerator: device type " + 
			                      numToStr(getDeviceType(k.getQueue())) + " is unknown!");
		}
	}

	template void ReductionAlgGenerator<double, ROT_SUM>::generateAlg(Kernel & k);
	template void ReductionAlgGenerator<float, ROT_SUM>::generateAlg(Kernel & k);
	template void ReductionAlgGenerator<double, ROT_PRODUCT>::generateAlg(Kernel & k);
	template void ReductionAlgGenerator<float, ROT_PRODUCT>::generateAlg(Kernel & k);
	template void ReductionAlgGenerator<double, ROT_MINIMUM>::generateAlg(Kernel & k);
	template void ReductionAlgGenerator<float, ROT_MINIMUM>::generateAlg(Kernel & k);
	template void ReductionAlgGenerator<double, ROT_MAXIMUM>::generateAlg(Kernel & k);
	template void ReductionAlgGenerator<float, ROT_MAXIMUM>::generateAlg(Kernel & k);

	template <typename ResType, enum ReductionOperatorType Op> 
		void ReductionAlgGenerator<ResType, Op>::generateAlg()
	{
		KernelConfiguration kConf(KERNEL_BASIC);
		kConf.local = true;
		
		kernel=std::make_shared<Kernel>(kConf);
		unsigned int groupsNumber = getNComputeUnits(ve[0]->getQueue());
		kernel->setGroupsNumber(groupsNumber);
		
		generateAlg(*kernel);
	};	
	template void ReductionAlgGenerator<double, ROT_SUM>::generateAlg();
	template void ReductionAlgGenerator<float, ROT_SUM>::generateAlg();
	template void ReductionAlgGenerator<double, ROT_PRODUCT>::generateAlg();
	template void ReductionAlgGenerator<float, ROT_PRODUCT>::generateAlg();
	template void ReductionAlgGenerator<double, ROT_MINIMUM>::generateAlg();
	template void ReductionAlgGenerator<float, ROT_MINIMUM>::generateAlg();
	template void ReductionAlgGenerator<double, ROT_MAXIMUM>::generateAlg();
	template void ReductionAlgGenerator<float, ROT_MAXIMUM>::generateAlg();

	template <typename ResType, enum ReductionOperatorType Op> 
		ReductionAlgGenerator<ResType, Op>::ReductionAlgGenerator(VectorOfElements v):
			ve(v),		 
			res(asl::AVec<ResType>(ve.size())),
			groupRes(ve.size())
	{
	}
	
	template <typename ResType> std::shared_ptr<ReductionAlgGenerator<ResType, ROT_SUM>> 
		generateSumAlg(VectorOfElements v)
	{
		return make_shared<ReductionAlgGenerator<ResType, ROT_SUM>>(v);
	}

	template std::shared_ptr<ReductionAlgGenerator<double, ROT_SUM>> 
		generateSumAlg<double>(VectorOfElements v);
	template std::shared_ptr<ReductionAlgGenerator<float, ROT_SUM>> 
		generateSumAlg<float>(VectorOfElements v);

	template <typename ResType> std::shared_ptr<ReductionAlgGenerator<ResType, ROT_MINIMUM>> 
		generateMinAlg(VectorOfElements v)
	{
		return make_shared<ReductionAlgGenerator<ResType, ROT_MINIMUM>>(v);
	}

	template std::shared_ptr<ReductionAlgGenerator<double, ROT_MINIMUM>> 
		generateMinAlg<double>(VectorOfElements v);
	template std::shared_ptr<ReductionAlgGenerator<float, ROT_MINIMUM>> 
		generateMinAlg<float>(VectorOfElements v);
		
	template <typename ResType> std::shared_ptr<ReductionAlgGenerator<ResType, ROT_MAXIMUM>> 
		generateMaxAlg(VectorOfElements v)
	{
		return make_shared<ReductionAlgGenerator<ResType, ROT_MAXIMUM>>(v);
	}

	template std::shared_ptr<ReductionAlgGenerator<double, ROT_MAXIMUM>> 
		generateMaxAlg<double>(VectorOfElements v);
	template std::shared_ptr<ReductionAlgGenerator<float, ROT_MAXIMUM>> 
		generateMaxAlg<float>(VectorOfElements v);

	template <typename ResType> std::shared_ptr<ReductionAlgGenerator<ResType, ROT_PRODUCT>> 
		generateProductAlg(VectorOfElements v)
	{
		return make_shared<ReductionAlgGenerator<ResType, ROT_PRODUCT>>(v);
	}

	template std::shared_ptr<ReductionAlgGenerator<double, ROT_PRODUCT>> 
		generateProductAlg<double>(VectorOfElements v);
	template std::shared_ptr<ReductionAlgGenerator<float, ROT_PRODUCT>> 
		generateProductAlg<float>(VectorOfElements v);
		
} // namespace acl
