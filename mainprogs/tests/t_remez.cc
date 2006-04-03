// $Id: t_remez.cc,v 3.0 2006-04-03 04:59:16 edwards Exp $
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

  XMLFileWriter xml_out("t_remez.xml");
  push(xml_out, "t_remez");

  proginfo(xml_out);  // basic program info

  Real lower = Real(1.0e-4);
  Real upper = 1.6;
  long prec = 80;
  int degree = 14;
  unsigned long power_num = 1;
  unsigned long power_den = 4;

  Remez  remez(lower, upper, prec);

  remez.generateApprox(degree, power_num, power_den);

  QDPIO::cout << "Start getPFE" << endl;
  RemezCoeff_t pfe = remez.getPFE();
  QDPIO::cout << "Finished getPFE" << endl;

  push(xml_out, "Remez_pfe");
  write(xml_out, "norm", pfe.norm);
  write(xml_out, "res", pfe.res);
  write(xml_out, "pole", pfe.pole);
  pop(xml_out);

  QDPIO::cout << "Start getIPFE" << endl;
  RemezCoeff_t ipfe = remez.getIPFE();
  QDPIO::cout << "Finished getIPFE" << endl;

  push(xml_out, "Remez_ipfe");
  write(xml_out, "norm", ipfe.norm);
  write(xml_out, "res", ipfe.res);
  write(xml_out, "pole", ipfe.pole);
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
