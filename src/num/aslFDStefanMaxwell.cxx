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


#include "aslFDStefanMaxwell.h"
#include <aslGenerators.h>
#include <acl/aclGenerators.h>
#include <acl/acl.h>
#include <math/aslTemplateVE.h>
#include <acl/aclMath/aclVectorOfElements.h>

#include <data/aslDataWithGhostNodes.h>
#include <acl/Kernels/aclKernel.h>
#include <math/aslTemplates.h>

#include <acl/Kernels/aclKernelConfigurationTemplates.h>

namespace asl
{
	FDStefanMaxwell::FDStefanMaxwell():
		SingleKernelNM(acl::KERNEL_SIMDUA),
		cData(0),
		cInternalData(0),
		vectorTemplate(NULL)
	{
	}


	FDStefanMaxwell::FDStefanMaxwell(Data c1,
	                                 Data c2,
	                                 const acl::VectorOfElements & dC, 
	                                 const VectorTemplate* vT):
		SingleKernelNM(acl::KERNEL_SIMDUA),
		cData({c1,c2}),
		cInternalData(0u),
		vectorTemplate(vT),
		diffusionCoefficients({vector<acl::VectorOfElements>({acl::generateVEConstant(0.), dC}),
			                  vector<acl::VectorOfElements>({dC,acl::generateVEConstant(0.)})})
	{
	}
		
	void FDStefanMaxwell::init0()
	{
		if (cData.size() != diffusionCoefficients.size())
			errorMessage("FDStefanMaxwell::init - some of compenents are underdefined");
		
		acl::TypeID type(getElementType(cData[0]->getDContainer()));
		unsigned int nC(cData.size());
		unsigned int nd(nD(*vectorTemplate));
		unsigned int nDir(vectorTemplate->vectors.size());
		
		cInternalData.resize(cData.size());
		for (unsigned int i(0); i < nC; ++i)
		{
			cInternalData[i] = clone(cData[i]);
			initData(cInternalData[i]->getEContainer(),cData[i]->getEContainer());
		}
		

		acl::VectorOfElements cnew(acl::generateVEPrivateVariable(nC, type));
		acl::VectorOfElements lambda(acl::generateVEPrivateVariable(nC, type));
		acl::VectorOfElements cTotal(acl::generateVEPrivateVariable(1, type));
		acl::MatrixOfElements aM(nC,nC);
		acl::VectorOfElements flux(acl::generateVEPrivateVariable(nC, type));
		
		vector<TemplateVE> velocityT(nd);
		if(velocity.get() != 0)
			for(unsigned int i(0); i < nd; ++i)
			{
				velocityT[i].init(*velocity, *vectorTemplate, i,false);
			}

		TemplateVE phiT;
		if(efPhi.get() != 0)
			phiT.init(*efPhi, *vectorTemplate, 0, false);
		
		vector<TemplateVE> cT(nC);
		for (unsigned int i(0); i < cData.size(); ++i)
		{
			cT[i].init(*cData[i], *vectorTemplate, 0, false);
			(*kernel)<<(subVE(cnew,i) =cT[i].getValue(0)); 
		}

		for(unsigned int kk(1); kk<nDir; ++kk)
		{
			vector<acl::VectorOfElements> cL;
			acl::VectorOfElements mNablaC(nC);
			
			(*kernel) << (cTotal = acl::generateVEConstant(0));
			//computation of cTotal and Lambda
			for (unsigned int i(0); i < nC; ++i)
			{
				cL.push_back((cT[i].getValue(0)+cT[i].getValue(kk))*.5);
				mNablaC[i]=(cT[i].getValue(0)-cT[i].getValue(kk))[0];
				(*kernel) << (cTotal += cL[i]);
				(*kernel) << (subVE(lambda,i) =acl::generateVEConstant(0.));
				for (unsigned int j(0); j < nC; ++j)
					if(j!=i)
						(*kernel) << (subVE(lambda,i) += cL[i]/diffusionCoefficients[i][j]); 
			}
			(*kernel) << (lambda/=cTotal);
			// eccounting of the electric field
			if(efPhi.get() != 0)
				for (unsigned int i(0); i < nC; ++i)
					mNablaC[i]=(subVE(mNablaC,i)-cL[i]*efCharge[i]*
					            (phiT.getValue(kk)-phiT.getValue(0)))[0];

			// compilation of matrix
			for (unsigned int i(0); i < nC; ++i)
				for (unsigned int j(0); j < nC; ++j)
				{
					if(i==j)
						aM.setElement(i,i, lambda[i]);
					else 
						aM.setElement(i,j, (-cL[i]/cTotal/diffusionCoefficients[i][j])[0]);
				}
			if(dustDiffusionCoefficients.size())
				for (unsigned int i(0); i < nC; ++i)
					aM.setElement(i,i,(1./dustDiffusionCoefficients[i]+aM.getVE(i,i))[0]);
			
			(*kernel) << gcSolveSystem(aM, mNablaC, flux);
			acl::VectorOfElements advFlux(nC);
			if(velocity.get() != 0)
			{
				acl::VectorOfElements velL((velocityT[0].getValue(kk)+velocityT[0].getValue(0)) * 
				                      .5 * vectorTemplate->vectors[kk]);
				for (unsigned int i(1); i < nd; ++i)
					copy(velL+(velocityT[i].getValue(kk)+velocityT[i].getValue(0)) * 
					     .5 * vectorTemplate->vectors[kk], velL);
				for (unsigned int i(0); i < nC; ++i)
					advFlux[i]=(cL[i]*velL)[0];
			}
			else
				copy(acl::generateVEConstantN(nC,0), advFlux);	
			
			(*kernel) << (cnew-=(flux + advFlux)*vectorTemplate->laplasCoefs[kk]);
		}
			
		for (unsigned int i(0); i < nC; ++i)	
			(*kernel) << (assignmentSafe(cInternalData[i]->getSubContainer(), subVE(cnew,i)));				
	}

		
	void FDStefanMaxwell::postProcessing()
	{
		for (unsigned int i(0); i < cData.size(); ++i)
			swapBuffers(cData[i]->getDContainer(), cInternalData[i]->getDContainer());
	}
		

	void FDStefanMaxwell::setDiffusionCoefficient(acl::VectorOfElements dC, 
	                                              unsigned int i,unsigned int j)
	{
		if(i==j)
			errorMessage("FDStefanMaxwell::setDiffusionCoefficient - i and j should have values");
		copy(dC,diffusionCoefficients[i][j]);
		copy(dC,diffusionCoefficients[j][i]);
	}
		
	void FDStefanMaxwell::addComponent(Data c, const Param & dC)
	{
		unsigned int n(diffusionCoefficients.size());
		diffusionCoefficients.resize(n+1);
		for(unsigned int i(0); i<n+1; ++i)
			diffusionCoefficients[i].resize(n+1);
		copy(acl::generateVEConstant(0),diffusionCoefficients[n][n]);
		for(unsigned int i(0); i<n; ++i)
			setDiffusionCoefficient(subVE(dC,i),n,i);
		cData.push_back(c);
	}

	void FDStefanMaxwell::addComponent(Data c, const Param & dC, const Param & q)
	{
		addComponent(c, dC);
		efCharge.push_back(q);
	}

	void FDStefanMaxwell::setVelocity(Field v)
	{
		velocity=v;
	}
		

	void FDStefanMaxwell::setElectricField(Field phi)
	{
		efPhi=phi;
	}

	void FDStefanMaxwell::setDustDiffusionCoefficient(unsigned int i, 
	                                                  const acl::VectorOfElements & dd)
	{
		if (i>=cData.size())
			errorMessage("FDStefanMaxwell::setDustDiffusionCoefficient: i is larger than the number of components");
		
		if (i+1>dustDiffusionCoefficients.size() && i<cData.size())
				dustDiffusionCoefficients.resize(i+1);
		acl::copy(dd,dustDiffusionCoefficients[i]);
	}

	void FDStefanMaxwell::setCharge(unsigned int i, const Param & q)
	{
		if (i>=cData.size())
			errorMessage("FDStefanMaxwell::setCharge: i is larger than the number of components");

		if (i+1>efCharge.size() && i<cData.size())
				efCharge.resize(i+1);
		acl::copy(q,efCharge[i]);
	}
		
	SPFDStefanMaxwell generateFDStefanMaxwell(SPDataWithGhostNodesACLData c1,
	                                          SPDataWithGhostNodesACLData c2,
	                                          double diffusionCoeff,            
	                                          SPAbstractDataWithGhostNodes v, 
	                                          const VectorTemplate* vt)
	{
		auto nm(make_shared<FDStefanMaxwell> (c1, c2, acl::generateVEConstant(diffusionCoeff), vt));
		nm->setVelocity(v);
		return nm;
	}

	SPFDStefanMaxwell generateFDStefanMaxwell(SPDataWithGhostNodesACLData c1,
	                                          SPDataWithGhostNodesACLData c2,
	                                          double diffusionCoeff, 
	                                          const VectorTemplate* vt)
	{
		auto nm(make_shared<FDStefanMaxwell> (c1, c2, acl::generateVEConstant(diffusionCoeff), vt));
		return nm;
	}

	FDStefanMaxwellElectricField::FDStefanMaxwellElectricField(SPFDStefanMaxwell sm, Data phi):
		SingleKernelNM(acl::KERNEL_SIMDUA),
		smSolver(sm),
		phi(phi)			
	{
	}
		
	void FDStefanMaxwellElectricField::init0()
	{

		auto & cData(smSolver->getData());
		auto vt(smSolver->getVectorTemplate());
		acl::TypeID type(getElementType(cData[0]->getDContainer()));
		unsigned int nC(cData.size());
		unsigned int nDir(vt->vectors.size());
		
		phiInternalData = clone(phi);
		initData(phiInternalData->getEContainer(), phi->getEContainer() );

		acl::VectorOfElements phiNew(acl::generateVEPrivateVariable(1, type));
		acl::VectorOfElements lambda(acl::generateVEPrivateVariable(nC, type));
		acl::VectorOfElements cTotal(acl::generateVEPrivateVariable(1, type));
		acl::MatrixOfElements aM(nC,nC);
		acl::VectorOfElements flux(acl::generateVEPrivateVariable(nC, type));
		
		vector<TemplateVE> cT(nC);
		for (unsigned int i(0); i < cData.size(); ++i)
			cT[i].init(*cData[i], *vt, 0, false);

		TemplateVE phiT(*phi, *vt, 0, false);
		(*kernel)<<(phiNew = phiT.getValue(0)+phiS->getSubContainer()*stepFactor); 
		
		for(unsigned int kk(1); kk<nDir; ++kk)
		{
			vector<acl::VectorOfElements> cL;
			acl::VectorOfElements cq(nC);
			acl::VectorOfElements q(nC);
			(*kernel) << (cTotal = acl::generateVEConstant(0));

			//computation of cTotal and Lambda
			for (unsigned int i(0); i < nC; ++i)
			{
				cL.push_back((cT[i].getValue(0)+cT[i].getValue(kk))*.5);
				q[i]=smSolver->getCharge(i)[0];
				cq[i]=(smSolver->getCharge(i)*cL[i])[0];
				(*kernel) << (cTotal += cL[i]);
				(*kernel) << (subVE(lambda,i) =acl:: generateVEConstant(0.));
				for (unsigned int j(0); j < nC; ++j)
					if(j!=i)
						(*kernel) << (subVE(lambda,i) += cL[i]/smSolver->getDiffusionCoefficient(i,j)); 
			}
			(*kernel) << (lambda/=cTotal);

			// compilation of matrix
			for (unsigned int i(0); i < nC; ++i)
				for (unsigned int j(0); j < nC; ++j)
				{
					if(i==j)
						aM.setElement(i,i, lambda[i]);
					else 
						aM.setElement(i,j, (-cL[i]/cTotal/smSolver->getDiffusionCoefficient(i,j))[0]);
				}
//			if(dustDiffusionCoefficients.size())
			for (unsigned int i(0); i < nC; ++i)
				aM.setElement(i,i,(1./smSolver->getDustDiffusionCoefficient(i)+aM.getVE(i,i))[0]);

			(*kernel) << gcSolveSystem(aM, cq, flux);

			(*kernel) << (phiNew+=stepFactor*(q*flux)*
			              (phiT.getValue(kk)-phiT.getValue(0))*
			              vt->laplasCoefs[kk]);//!!.... 
		}
			
		(*kernel) << (assignmentSafe(phiInternalData->getSubContainer(), phiNew));
	}

		
	void FDStefanMaxwellElectricField::postProcessing()
	{
		swapBuffers(phi->getDContainer(), 
		            phiInternalData->getDContainer());
	}

	void FDStefanMaxwellElectricField::setPhiS(Field pS)
	{
		phiS=pS;
	}

		
} //asl
