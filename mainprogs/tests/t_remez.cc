// $Id: t_remez.cc,v 3.3 2008-06-02 16:18:41 bjoo Exp $
/*! \file
 *  \brief Test the Remez code
 */

#include "chroma.h"

using namespace Chroma;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {2,2,2,2};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  Real lower;
  Real upper;
  long prec;
  int degree;
  long power_num;
  long power_den;

  XMLReader xml_in(Chroma::getXMLInputFileName());
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();

  try { 
    XMLReader paramtop(xml_in, "/Remez");
    read(paramtop, "lower", lower);
    read(paramtop, "upper", upper);
    read(paramtop, "powerNum", power_num);
    read(paramtop, "powerDen", power_den);
    read(paramtop, "degree", degree);
    if( paramtop.count("prec") == 1 ) { 
      read(paramtop, "prec", prec);
    }
    else { 
      prec=80;
    }

    push(xml_out, "Remez");
    proginfo(xml_out);  // basic program info
    push(xml_out, "Param");
    xml_out << paramtop;
    pop(xml_out);
  
  }
  catch(const std::string& s) { 
    QDPIO::cout << "Caught Exception reading parameters: " << s << endl;
    QDP_abort(1);
  }
  catch(...) { 
    QDPIO::cout << "Caught unknown exception while processing input " << endl;
    QDP_abort(1);
  }

  bool invertP;

  long multmp = power_num*power_den;
  if( multmp > 0 ) { 
    invertP = false;
    QDPIO::cout << "Either both num and den powers are + or both are -" << endl;
  }
  else { 
    invertP = true;
    QDPIO::cout << "One of num or den powers is -" << endl;
  }

  unsigned long pn = power_num > 0 ? power_num : -power_num ;
  unsigned long pd = power_den > 0 ? power_den : -power_den ;

  
  Remez  remez(lower, upper, prec);
  remez.generateApprox(degree, pn, pd);

  QDPIO::cout << "Start getPFE" << endl;
  RemezCoeff_t pfe;

  if( ! invertP ) {
    pfe = remez.getPFE();
  }
  else { 
    pfe = remez.getIPFE();
  }

  QDPIO::cout << "Finished getPFE" << endl;

  push(xml_out, "PFECoeffs");
  write(xml_out, "norm", pfe.norm);
  write(xml_out, "res", pfe.res);
  write(xml_out, "pole", pfe.pole);
  pop(xml_out);

  QDPIO::cout << "Start getIPFE" << endl;
  RemezCoeff_t ipfe;

  if( !invertP )  {
    ipfe = remez.getIPFE();
  }
  else { 
    ipfe = remez.getPFE();
  }
  QDPIO::cout << "Finished getIPFE" << endl;

  push(xml_out, "IPFECoeffs");
  write(xml_out, "norm", ipfe.norm);
  write(xml_out, "res", ipfe.res);
  write(xml_out, "pole", ipfe.pole);
  pop(xml_out);

  pop(xml_out);

  int N = 20;
  for(int n=0; n < N; ++n)
  {
    Double x = exp(log(lower) + n*(log(upper)-log(lower))/Double(N));
    Double f = remez.evalPFE(x,ipfe);
    Double fn = 1/sqrt(x);
    QDPIO::cout << "x=" << x 
		<< " f(x)=" << f
		<< " fn=" << fn
		<< "   diff=" << Double(f-fn) << endl;
  }

  xml_out.close();

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
