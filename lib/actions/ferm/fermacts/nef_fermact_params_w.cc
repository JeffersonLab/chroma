// $Id: nef_fermact_params_w.cc,v 3.4 2008-11-10 17:59:07 bjoo Exp $
/*! \file
 *  \brief NEF fermion action parameters
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/nef_fermact_params_w.h"

#include "io/param_io.h"

namespace Chroma
{

  //! Default constructor
  NEFFermActParams::NEFFermActParams() : b5(1), c5(1)
  {	  
    Mass = Real(0);
	OverMass= Real(0);
	N5=1;
	b5[0]=Real(1);
	c5[0]=Real(1);
    u0   = Real(1);
    clovCoeffR = clovCoeffT = Real(0);
    max_norm=0;
    max_norm_usedP=false;

    sub_zero=0;
    sub_zero_usedP=false;
  }

  //! Read parameters
  NEFFermActParams::NEFFermActParams(XMLReader& xml, const string& path)
  {	  
    XMLReader paramtop(xml, path);
	
    // Read the stuff for the action
    if (paramtop.count("Mass") != 0) 
    {
      read(paramtop, "Mass", Mass);
      if (paramtop.count("Kappa") != 0) 
      {
	QDPIO::cerr << "Error: found both a Kappa and a Mass tag" << endl;
	QDP_abort(1);
      }
    }
    else if (paramtop.count("Kappa") != 0)
    {
      Real Kappa;
      read(paramtop, "Kappa", Kappa);
      Mass = kappaToMass(Kappa);    // Convert Kappa to Mass
      QDPIO::cout << "Kappa is " << Kappa << "  Mass is " << Mass << endl << flush;
    }
    else
    {
      QDPIO::cerr << "Error: neither Mass or Kappa found" << endl;
      QDP_abort(1);
    }
	
	//same for Mass5:
    if (paramtop.count("OverMass") != 0) 
    {
      read(paramtop, "OverMass", OverMass);
      if (paramtop.count("Kappa5") != 0) 
      {
	QDPIO::cerr << "Error: found both a Kappa5 and a OverMass tag" << endl;
	QDP_abort(1);
      }
    }
    else if (paramtop.count("Kappa5") != 0)
    {
      Real Kappa5;
      read(paramtop, "Kappa5", Kappa5);
      OverMass = kappaToMass(Kappa5)-Real(1);    // Convert Kappa to Mass
      QDPIO::cout << "Kappa5 is " << Kappa5 << "  OverMass is " << OverMass << endl << flush;
    }
    else
    {
      QDPIO::cerr << "Error: neither OverMass or Kappa5 found" << endl;
      QDP_abort(1);
    }

    // Read optional u0
    if (paramtop.count("u0") != 0) 
      read(paramtop, "u0", u0);
    else {
      u0 = Real(1);
      QDPIO::cout << "u0 is " << u0 << endl << flush;
    }
    // Read optional anisoParam
    if (paramtop.count("AnisoParam") != 0) 
      read(paramtop, "AnisoParam", anisoParam);

    // If aniso, read all clover coeff
    if (anisoParam.anisoP)
    {
      read(paramtop, "clovCoeffR", clovCoeffR);
      read(paramtop, "clovCoeffT", clovCoeffT);
    }
    else
    {
      Real clovCoeff;
      read(paramtop, "clovCoeff", clovCoeff);
      clovCoeffR = clovCoeff;
      clovCoeffT = clovCoeff;
    }
	
	//read Ls and other DWF params:
    if (paramtop.count("N5") != 0) 
    {
		read(paramtop, "N5", N5);
	}
	else{
        QDPIO::cerr << "Error: no N5 found" << endl;
        QDP_abort(1);
	}
	if(paramtop.count("b5") !=0) read(paramtop,"b5",b5);
	else{
		b5.resize(N5);
		for(unsigned int s=0; s<N5; s++) b5[s]=Real(1);
	}
	if(paramtop.count("c5") !=0) read(paramtop,"c5",c5);
	else{
		c5.resize(N5);
		for(unsigned int s=0; s<N5; s++) c5[s]=Real(1);
	}

    if( paramtop.count("MaxNorm") != 0 ){
      read(paramtop, "MaxNorm", max_norm);
      QDPIO::cout << "MaxNorm="<<max_norm<<endl;
      max_norm_usedP=true;
    }
    else { 
      max_norm_usedP=false;
      max_norm=Real(0);
    }

    if( paramtop.count("ZeroEnergy") != 0 ) { 
      read(paramtop, "ZeroEnergy", sub_zero);
      sub_zero_usedP=true;
    }
    else { 
      sub_zero_usedP=false;
      sub_zero=Real(0);
    }

    if( paramtop.count("TwistedM") != 0 ) { 
      twisted_m_usedP = true;
      read(paramtop, "TwistedM", twisted_m);
    }
    else {				       
      twisted_m_usedP = false;
    }

  }

  //! Read parameters
  void read(XMLReader& xml, const string& path, NEFFermActParams& param)
  {	  
    NEFFermActParams tmp(xml, path);	
    param = tmp;
  }


  //! Write parameters
  void write(XMLWriter& xml, const string& path, const NEFFermActParams& param)
  {
    push(xml, path);

    write(xml, "Mass", param.Mass);
	write(xml, "OverMass", param.OverMass);
    write(xml, "u0", param.u0);

    if (param.anisoParam.anisoP)
    {
      write(xml, "clovCoeffR", param.clovCoeffR);
      write(xml, "clovCoeffT", param.clovCoeffT);
    }
    else
    {
      write(xml, "clovCoeff", param.clovCoeffR);
    }
	write(xml, "N5", param.N5);
	write(xml, "b5", param.b5);
	write(xml, "c5", param.c5);

    if (param.max_norm_usedP){
      write(xml, "MaxNorm", param.max_norm);
    }

    if (param.sub_zero_usedP == true ) {
      write(xml, "ZeroEnergy", param.sub_zero);
    }

    if (param.twisted_m_usedP == true ) { 
      write(xml, "TwistedM", param.twisted_m);
    }

    pop(xml);
  }


}

