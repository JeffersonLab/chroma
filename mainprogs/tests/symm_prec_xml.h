/*
 * symm_prec_xml.h
 *
 *  Created on: Oct 18, 2018
 *      Author: bjoo
 */

#ifndef MAINPROGS_TESTS_SYMM_PREC_XML_H_
#define MAINPROGS_TESTS_SYMM_PREC_XML_H_
#include <string>
#include "chroma_config.h"
namespace SymmPrecTesting
{


std::string fermact_xml_asymm =
  "<?xml version='1.0'?>                              \
   <Param>					      \
   <FermionAction>                                    \
     <FermAct>CLOVER</FermAct>                        \
     <Mass>0.1</Mass>				      \
     <clovCoeff>1</clovCoeff>			      \
     <AnisoParam>				      \
       <anisoP>false</anisoP>			      \
       <t_dir>3</t_dir>				      \
       <xi_0>1</xi_0>				      \
       <nu>1</nu>				      \
     </AnisoParam>				      \
    <FermState>					      \
       <Name>STOUT_FERM_STATE</Name>		      \
       <rho>0.14</rho>				      \
       <n_smear>2</n_smear>			      \
       <orthog_dir>3</orthog_dir>		      \
       <FermionBC>				      \
         <FermBC>SIMPLE_FERMBC</FermBC>		      \
         <boundary>1 1 1 -1</boundary>		      \
       </FermionBC>				      \
     </FermState>				      \
   </FermionAction>				      \
  </Param>";

std::string fermact_xml_symm =
  "<?xml version='1.0'?>                              \
   <Param>					      \
   <FermionAction>                                    \
     <FermAct>SEOPREC_CLOVER</FermAct>                        \
     <Mass>0.1</Mass>				      \
     <clovCoeff>1</clovCoeff>			      \
     <AnisoParam>				      \
       <anisoP>false</anisoP>			      \
       <t_dir>3</t_dir>				      \
       <xi_0>1</xi_0>				      \
       <nu>1</nu>				      \
     </AnisoParam>				      \
    <FermState>					      \
       <Name>STOUT_FERM_STATE</Name>		      \
       <rho>0.14</rho>				      \
       <n_smear>2</n_smear>			      \
       <orthog_dir>3</orthog_dir>		      \
       <FermionBC>				      \
         <FermBC>SIMPLE_FERMBC</FermBC>		      \
         <boundary>1 1 1 -1</boundary>		      \
       </FermionBC>				      \
     </FermState>				      \
   </FermionAction>				      \
  </Param>";

std::string fermact_xml_asymm_periodic =
  "<?xml version='1.0'?>                              \
   <Param>					      \
   <FermionAction>                                    \
     <FermAct>CLOVER</FermAct>                        \
     <Mass>0.1</Mass>				      \
     <clovCoeff>1</clovCoeff>			      \
     <AnisoParam>				      \
       <anisoP>false</anisoP>			      \
       <t_dir>3</t_dir>				      \
       <xi_0>1</xi_0>				      \
       <nu>1</nu>				      \
     </AnisoParam>				      \
    <FermState>					      \
       <Name>STOUT_FERM_STATE</Name>		      \
       <rho>0.14</rho>				      \
       <n_smear>2</n_smear>			      \
       <orthog_dir>3</orthog_dir>		      \
       <FermionBC>				      \
         <FermBC>SIMPLE_FERMBC</FermBC>		      \
         <boundary>1 1 1 1</boundary>		      \
       </FermionBC>				      \
     </FermState>				      \
   </FermionAction>				      \
  </Param>";

std::string fermact_xml_symm_periodic =
  "<?xml version='1.0'?>                              \
   <Param>					      \
   <FermionAction>                                    \
     <FermAct>SEOPREC_CLOVER</FermAct>                        \
     <Mass>0.1</Mass>				      \
     <clovCoeff>1</clovCoeff>			      \
     <AnisoParam>				      \
       <anisoP>false</anisoP>			      \
       <t_dir>3</t_dir>				      \
       <xi_0>1</xi_0>				      \
       <nu>1</nu>				      \
     </AnisoParam>				      \
    <FermState>					      \
       <Name>STOUT_FERM_STATE</Name>		      \
       <rho>0.14</rho>				      \
       <n_smear>2</n_smear>			      \
       <orthog_dir>3</orthog_dir>		      \
       <FermionBC>				      \
         <FermBC>SIMPLE_FERMBC</FermBC>		      \
         <boundary>1 1 1 1</boundary>		      \
       </FermionBC>				      \
     </FermState>				      \
   </FermionAction>				      \
  </Param>";

std::string inv_param_syssolver_bicgstab_xml = \
		"<?xml version='1.0'?> \
		<Param> \
		  <InvertParam>\
           	 <invType>BICGSTAB_INVERTER</invType>\
		     <RsdBiCGStab>1.0e-8</RsdBiCGStab> \
		     <MaxBiCGStab>1000</MaxBiCGStab> \
          </InvertParam>\
		</Param>";

#ifdef BUILD_QUDA
std::string inv_param_quda_bicgstab_xml = \
		"<?xml version='1.0'?> \
		<Param> \
		  <InvertParam>\
            <invType>QUDA_CLOVER_INVERTER</invType>\
  	  	  	<CloverParams>\
              <Mass>0.1</Mass>				      \
              <clovCoeff>1</clovCoeff>			      \
              <AnisoParam>				      \
                <anisoP>false</anisoP>			      \
                <t_dir>3</t_dir>				      \
                <xi_0>1</xi_0>				      \
                <nu>1</nu>				      \
              </AnisoParam>\
            </CloverParams>\
            <RsdTarget>1.0e-8</RsdTarget>\
            <Delta>1.0e-1</Delta>\
            <Pipeline>0</Pipeline>\
            <MaxIter>500</MaxIter>\
		    <SolverType>BICGSTAB</SolverType> \
            <RsdToleranceFactor>100.0</RsdToleranceFactor>\
            <AntiPeriodicT>true</AntiPeriodicT>\
            <Verbose>true</Verbose>\
            <AsymmetricLinop>false</AsymmetricLinop>\
            <CudaReconstruct>RECONS_12</CudaReconstruct>\
            <CudaSloppyPrecision>SINGLE</CudaSloppyPrecision>\
            <CudaSloppyReconstruct>RECONS_12</CudaSloppyReconstruct>\
            <AxialGaugeFix>false</AxialGaugeFix>\
            <AutotuneDslash>true</AutotuneDslash>\
            <SolutionCheckP>true</SolutionCheckP>\
          </InvertParam>\
		</Param>";

std::string inv_param_quda_multigrid_xml = \
		"<?xml version='1.0'?> \
		<Param> \
		  <InvertParam>\
            <invType>QUDA_MULTIGRID_CLOVER_INVERTER</invType>\
  	  	  	<CloverParams>\
              <Mass>0.1</Mass>				      \
              <clovCoeff>1</clovCoeff>			      \
              <AnisoParam>				      \
                <anisoP>false</anisoP>			      \
                <t_dir>3</t_dir>				      \
                <xi_0>1</xi_0>				      \
                <nu>1</nu>				      \
              </AnisoParam>\
            </CloverParams>\
               <RsdTarget>1.0e-8</RsdTarget> \
            <Delta>1.0e-1</Delta>\
            <Pipeline>4</Pipeline> \
            <MaxIter>500</MaxIter> \
            <RsdToleranceFactor>100.0</RsdToleranceFactor>\
            <AntiPeriodicT>true</AntiPeriodicT>\
            <SolverType>GCR</SolverType>\
            <Verbose>false</Verbose>\
            <AsymmetricLinop>false</AsymmetricLinop>\
            <CudaReconstruct>RECONS_12</CudaReconstruct>\
            <CudaSloppyPrecision>SINGLE</CudaSloppyPrecision>\
            <CudaSloppyReconstruct>RECONS_12</CudaSloppyReconstruct>\
            <AxialGaugeFix>false</AxialGaugeFix>\
            <AutotuneDslash>true</AutotuneDslash>\
            <MULTIGRIDParams>\
              <Verbosity>true</Verbosity>\
              <Precision>HALF</Precision>\
              <Reconstruct>RECONS_12</Reconstruct>\
              <Blocking>\
                <elem>2 2 2 4</elem>\
              </Blocking>\
              <CoarseSolverType>\
                <elem>CA_GCR</elem>\
              </CoarseSolverType>\
              <CoarseResidual>1.0e-1</CoarseResidual>\
              <MaxCoarseIterations>12</MaxCoarseIterations>\
              <RelaxationOmegaMG>1.0</RelaxationOmegaMG>\
              <SmootherType>\
                <elem>CA_GCR</elem>\
              </SmootherType>\
              <SmootherTol>0.25</SmootherTol>\
              <SmootherSchwarzCycle>1</SmootherSchwarzCycle>\
              <NullVectors>24</NullVectors>\
              <Pre-SmootherApplications>0</Pre-SmootherApplications>\
              <Post-SmootherApplications>8</Post-SmootherApplications>\
              <SubspaceSolver>\
                <elem>CG</elem>\
              </SubspaceSolver>\
              <RsdTargetSubspaceCreate>5e-06</RsdTargetSubspaceCreate>\
              <MaxIterSubspaceCreate>500</MaxIterSubspaceCreate>\
              <MaxIterSubspaceRefresh>500</MaxIterSubspaceRefresh>\
              <OuterGCRNKrylov>20</OuterGCRNKrylov>\
              <PrecondGCRNKrylov>10</PrecondGCRNKrylov>\
              <GenerateNullspace>true</GenerateNullspace>\
              <CheckMultigridSetup>false</CheckMultigridSetup>\
              <GenerateAllLevels>true</GenerateAllLevels>\
              <CycleType>MG_RECURSIVE</CycleType>\
              <SchwarzType>ADDITIVE_SCHWARZ</SchwarzType>\
              <RelaxationOmegaOuter>1.0</RelaxationOmegaOuter>\
              <SetupOnGPU>1</SetupOnGPU>\
            </MULTIGRIDParams>\
            <SubspaceID>mg_subspace</SubspaceID>\
            <SolutionCheckP>true</SolutionCheckP>\
          </InvertParam>\
		</Param>";

std::string inv_param_quda_bicgstab_asymm_xml = \
		"<?xml version='1.0'?> \
		<Param> \
		  <InvertParam>\
            <invType>QUDA_CLOVER_INVERTER</invType>\
  	  	  	<CloverParams>\
              <Mass>0.1</Mass>				      \
              <clovCoeff>1</clovCoeff>			      \
              <AnisoParam>				      \
                <anisoP>false</anisoP>			      \
                <t_dir>3</t_dir>				      \
                <xi_0>1</xi_0>				      \
                <nu>1</nu>				      \
              </AnisoParam>\
            </CloverParams>\
            <RsdTarget>1.0e-8</RsdTarget>\
            <Delta>1.0e-1</Delta>\
            <Pipeline>0</Pipeline>\
            <MaxIter>500</MaxIter>\
		    <SolverType>BICGSTAB</SolverType> \
            <RsdToleranceFactor>100.0</RsdToleranceFactor>\
            <AntiPeriodicT>true</AntiPeriodicT>\
            <Verbose>true</Verbose>\
            <AsymmetricLinop>true</AsymmetricLinop>\
            <CudaReconstruct>RECONS_12</CudaReconstruct>\
            <CudaSloppyPrecision>SINGLE</CudaSloppyPrecision>\
            <CudaSloppyReconstruct>RECONS_12</CudaSloppyReconstruct>\
            <AxialGaugeFix>false</AxialGaugeFix>\
            <AutotuneDslash>true</AutotuneDslash>\
            <SolutionCheckP>true</SolutionCheckP>\
          </InvertParam>\
		</Param>";

std::string inv_param_quda_multigrid_asymm_xml = \
		"<?xml version='1.0'?> \
		<Param> \
		  <InvertParam>\
            <invType>QUDA_MULTIGRID_CLOVER_INVERTER</invType>\
  	  	  	<CloverParams>\
              <Mass>0.1</Mass>				      \
              <clovCoeff>1</clovCoeff>			      \
              <AnisoParam>				      \
                <anisoP>false</anisoP>			      \
                <t_dir>3</t_dir>				      \
                <xi_0>1</xi_0>				      \
                <nu>1</nu>				      \
              </AnisoParam>\
            </CloverParams>\
               <RsdTarget>1.0e-8</RsdTarget> \
            <Delta>1.0e-1</Delta>\
            <Pipeline>4</Pipeline> \
            <MaxIter>500</MaxIter> \
            <RsdToleranceFactor>100.0</RsdToleranceFactor>\
            <AntiPeriodicT>true</AntiPeriodicT>\
            <SolverType>GCR</SolverType>\
            <Verbose>false</Verbose>\
            <AsymmetricLinop>true</AsymmetricLinop>\
            <CudaReconstruct>RECONS_12</CudaReconstruct>\
            <CudaSloppyPrecision>SINGLE</CudaSloppyPrecision>\
            <CudaSloppyReconstruct>RECONS_12</CudaSloppyReconstruct>\
            <AxialGaugeFix>false</AxialGaugeFix>\
            <AutotuneDslash>true</AutotuneDslash>\
            <MULTIGRIDParams>\
              <Verbosity>true</Verbosity>\
              <Precision>HALF</Precision>\
              <Reconstruct>RECONS_12</Reconstruct>\
              <Blocking>\
                <elem>2 2 2 4</elem>\
              </Blocking>\
              <CoarseSolverType>\
                <elem>CA_GCR</elem>\
              </CoarseSolverType>\
              <CoarseResidual>1.0e-1</CoarseResidual>\
              <MaxCoarseIterations>12</MaxCoarseIterations>\
              <RelaxationOmegaMG>1.0</RelaxationOmegaMG>\
              <SmootherType>\
                <elem>CA_GCR</elem>\
              </SmootherType>\
              <SmootherTol>0.25</SmootherTol>\
              <SmootherSchwarzCycle>1</SmootherSchwarzCycle>\
              <NullVectors>24</NullVectors>\
              <Pre-SmootherApplications>0</Pre-SmootherApplications>\
              <Post-SmootherApplications>8</Post-SmootherApplications>\
              <SubspaceSolver>\
                <elem>CG</elem>\
              </SubspaceSolver>\
              <RsdTargetSubspaceCreate>5e-06</RsdTargetSubspaceCreate>\
              <MaxIterSubspaceCreate>500</MaxIterSubspaceCreate>\
              <MaxIterSubspaceRefresh>500</MaxIterSubspaceRefresh>\
              <OuterGCRNKrylov>20</OuterGCRNKrylov>\
              <PrecondGCRNKrylov>10</PrecondGCRNKrylov>\
              <GenerateNullspace>true</GenerateNullspace>\
              <CheckMultigridSetup>false</CheckMultigridSetup>\
              <GenerateAllLevels>true</GenerateAllLevels>\
              <CycleType>MG_RECURSIVE</CycleType>\
              <SchwarzType>ADDITIVE_SCHWARZ</SchwarzType>\
              <RelaxationOmegaOuter>1.0</RelaxationOmegaOuter>\
              <SetupOnGPU>1</SetupOnGPU>\
            </MULTIGRIDParams>\
            <SubspaceID>mg_subspace</SubspaceID>\
            <SolutionCheckP>true</SolutionCheckP>\
          </InvertParam>\
		</Param>";
#endif


std::string inv_param_multi_cg_xml = \
		"<?xml version='1.0'?> \
		<Param> \
		  <InvertParam>\
           	 <invType>CG_INVERTER</invType>\
		     <RsdCG>1.0e-8 1.0e-8 1.0e-8</RsdCG> \
		     <MaxCG>1000</MaxCG> \
          </InvertParam>\
		</Param>";

#ifdef BUILD_QUDA
std::string inv_param_multi_cg_quda_xml = \
		"<?xml version='1.0'?> \
		 <Param>\
		 <InvertParam> \
		   <invType>MULTI_CG_QUDA_CLOVER_INVERTER</invType> \
		   <CloverParams> \
              <Mass>0.1</Mass> \
			  <clovCoeff>1</clovCoeff> \
			  <AnisoParam> \
				 <anisoP>false</anisoP>	  \
				 <t_dir>3</t_dir>	 \
				 <xi_0>1</xi_0>	  \
				 <nu>1</nu>	  \
			  </AnisoParam> \
		    </CloverParams> \
		    <RsdTarget>1e-08 1e-08 1e-08</RsdTarget> \
		    <Delta>1.0e-1</Delta> \
		    <Pipeline>0</Pipeline> \
		    <MaxIter>50000</MaxIter> \
		    <RsdToleranceFactor>100</RsdToleranceFactor> \
            <AntiPeriodicT>true</AntiPeriodicT> \
		    <SolverType>CG</SolverType> \
		    <Verbose>false</Verbose> \
		    <CheckShifts>false</CheckShifts> \
		    <AsymmetricLinop>false</AsymmetricLinop> \
		    <CudaReconstruct>RECONS_12</CudaReconstruct> \
 		    <CudaSloppyPrecision>HALF</CudaSloppyPrecision> \
		    <CudaSloppyReconstruct>RECONS_12</CudaSloppyReconstruct> \
		    <AxialGaugeFix>false</AxialGaugeFix> \
		    <AutotuneDslash>true</AutotuneDslash> \
		  </InvertParam> \
		 </Param>";

std::string inv_param_multi_cg_quda_asymm_xml = \
		"<?xml version='1.0'?> \
		 <Param>\
		 <InvertParam> \
		   <invType>MULTI_CG_QUDA_CLOVER_INVERTER</invType> \
		   <CloverParams> \
              <Mass>0.1</Mass> \
			  <clovCoeff>1</clovCoeff> \
			  <AnisoParam> \
				 <anisoP>false</anisoP>	  \
				 <t_dir>3</t_dir>	 \
				 <xi_0>1</xi_0>	  \
				 <nu>1</nu>	  \
			  </AnisoParam> \
		    </CloverParams> \
		    <RsdTarget>1e-08 1e-08 1e-08</RsdTarget> \
		    <Delta>1.0e-1</Delta> \
		    <Pipeline>0</Pipeline> \
		    <MaxIter>50000</MaxIter> \
		    <RsdToleranceFactor>100</RsdToleranceFactor> \
            <AntiPeriodicT>true</AntiPeriodicT> \
		    <SolverType>CG</SolverType> \
		    <Verbose>false</Verbose> \
		    <CheckShifts>false</CheckShifts> \
		    <AsymmetricLinop>true</AsymmetricLinop> \
		    <CudaReconstruct>RECONS_12</CudaReconstruct> \
 		    <CudaSloppyPrecision>HALF</CudaSloppyPrecision> \
		    <CudaSloppyReconstruct>RECONS_12</CudaSloppyReconstruct> \
		    <AxialGaugeFix>false</AxialGaugeFix> \
		    <AutotuneDslash>true</AutotuneDslash> \
		  </InvertParam> \
		 </Param>";
#endif
}// namespace

#endif /* MAINPROGS_TESTS_SYMM_PREC_XML_H_ */
