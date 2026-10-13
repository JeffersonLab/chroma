#include <string>

namespace MGParamTestsEnv 
{

    std::string eig_xml = "<?xml version='1.0'?> \
    <Param> \
      <EigParams> \
        <EigType>TR_LANCZOS</EigType> \
        <SpectrumType>LR</SpectrumType> \
        <NEv>16</NEv>  \
        <NKr>96</NKr>  \
        <UsePolyAcc>true</UsePolyAcc> \
        <UseNormOp>true</UseNormOp> \
        <PolyDeg>128</PolyDeg> \
        <Amin>0.01</Amin> \
        <Amax>4</Amax> \
      </EigParams> \
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

  std::string inv_param_quda_multigrid_asymm_eig_xml = \
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
              <MgEigDeflation>  \
                <elem>              \
                  <Level>1</Level>  \
                  <MGEigParam>      \
                    <EigType>TR_LANCZOS</EigType> \
                    <UseDagger>false</UseDagger> \
                    <UseNormOp>false</UseNormOp> \
                    <NEv>66</NEv>  \
                    <NConv>64</NConv> \
                    <NKr>198</NKr>  \
                    <Tol>1.0e-4</Tol> \
                    <UsePolyAcc>false</UsePolyAcc> \
                  </MGEigParam> \
                </elem> \
             </MgEigDeflation> \
            </MULTIGRIDParams>\
            <SubspaceID>mg_subspace</SubspaceID>\
            <SolutionCheckP>true</SolutionCheckP>\
          </InvertParam>\
		</Param>";    
};
