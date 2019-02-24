/*
 * quda_cmat_xml.h
 *
 *  Created on: Feb 24, 2019
 *      Author: bjoo
 */

#ifndef MAINPROGS_TESTS_QUDA_CMAT_XML_H_
#define MAINPROGS_TESTS_QUDA_CMAT_XML_H_

#include <string>

namespace SymmPrecTesting
{


std::string inv_param_quda_gcr_twisted_xml = \
		"<?xml version='1.0'?> \
		<Param> \
		  <InvertParam>\
            <invType>QUDA_CLOVER_INVERTER</invType>\
  	  	  	<CloverParams>\
              <Mass>0.1</Mass>				      \
              <clovCoeff>1.0</clovCoeff>			      \
			  <TwistedM>0.2345</TwistedM> \
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
		    <SolverType>GCR</SolverType> \
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
		    <GCRInnerParams> \
              <PreconditionCycle>16</PreconditionCycle> \
              <NKrylov>8</NKrylov> \
              <MaxIterPrecondition>12</MaxIterPrecondition> \
              <RsdPrecondition>1.0e-1</RsdPrecondition> \
              <VerboseP>false</VerboseP> \
              <InvTypePrecondition>MR</InvTypePrecondition> \
              <SchwarzType>ADDITIVE_SCHWARZ</SchwarzType> \
              <PrecPrecondition>HALF</PrecPrecondition> \
              <ReconstructPrecondition>RECONS_12</ReconstructPrecondition> \
            </GCRInnerParams> \
          </InvertParam>\
		</Param>";

std::string inv_param_quda_gcr_twisted_asymm_xml = \
		"<?xml version='1.0'?> \
		<Param> \
		  <InvertParam>\
            <invType>QUDA_CLOVER_INVERTER</invType>\
  	  	  	<CloverParams>\
              <Mass>0.1</Mass>				      \
              <clovCoeff>1.0</clovCoeff>			      \
			  <TwistedM>0.2345</TwistedM> \
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
		    <SolverType>GCR</SolverType> \
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
		    <GCRInnerParams> \
              <PreconditionCycle>16</PreconditionCycle> \
              <NKrylov>8</NKrylov> \
              <MaxIterPrecondition>12</MaxIterPrecondition> \
              <RsdPrecondition>1.0e-1</RsdPrecondition> \
              <VerboseP>false</VerboseP> \
              <InvTypePrecondition>MR</InvTypePrecondition> \
              <SchwarzType>ADDITIVE_SCHWARZ</SchwarzType> \
              <PrecPrecondition>HALF</PrecPrecondition> \
              <ReconstructPrecondition>RECONS_12</ReconstructPrecondition> \
            </GCRInnerParams> \
          </InvertParam>\
		</Param>";

} // namespace

#endif /* MAINPROGS_TESTS_QUDA_CMAT_XML_H_ */
