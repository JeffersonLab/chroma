// $Id: t_site_ops_scalar.cc,v 1.1 2008-08-26 13:24:29 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "qdp.h"
#include "unittest.h"

#include "testDslashFull.h"
#include "testSiteDslash.h"

using namespace QDP;

int main(int argc, char **argv)
{
  // Initialize QDP++ with argc, and argv. Set Lattice Dimensions
  const int latdims[] = {4,2,6,12};

  // Initialize UnitTest jig
  TestRunner  tests(&argc, &argv, latdims);

  tests.addTest(new testSiteDslash0PlusForward(), "testSiteDslash0PlusForward" );
  tests.addTest(new testSiteDslash0PlusBackwardAdd(), "testSiteDslash0PlusBackward" );
  tests.addTest(new testSiteDslash1PlusForwardAdd(), "testSiteDslash1PlusForwardAdd" );
  tests.addTest(new testSiteDslash1PlusBackwardAdd(), "testSiteDslash1PlusBackwardAdd" );
  tests.addTest(new testSiteDslash2PlusForwardAdd(), "testSiteDslash2PlusForwardAdd" );
  tests.addTest(new testSiteDslash2PlusBackwardAdd(), "testSiteDslash2PlusBackwardAdd" );
  tests.addTest(new testSiteDslash3PlusForwardAdd(), "testSiteDslash3PlusForwardAdd" );
  tests.addTest(new testSiteDslash3PlusBackwardAddStore(), "testSiteDslash3PlusBackwardAddStore" );

  tests.addTest(new testSiteDslash0MinusForward(), "testSiteDslash0MinusForward" );
  tests.addTest(new testSiteDslash0MinusBackwardAdd(), "testSiteDslash0MinusBackward" );
  tests.addTest(new testSiteDslash1MinusForwardAdd(), "testSiteDslash1MinusForwardAdd" );
  tests.addTest(new testSiteDslash1MinusBackwardAdd(), "testSiteDslash1MinusBackwardAdd" );
  tests.addTest(new testSiteDslash2MinusForwardAdd(), "testSiteDslash2MinusForwardAdd" );
  tests.addTest(new testSiteDslash2MinusBackwardAdd(), "testSiteDslash2MinusBackwardAdd" );
  tests.addTest(new testSiteDslash3MinusForwardAdd(), "testSiteDslash3MinusForwardAdd" );
  tests.addTest(new testSiteDslash3MinusBackwardAddStore(), "testSiteDslash3MinusBackwardAddStore" );

  // Run all tests
  tests.run();

  // Testjig is destroyed
  tests.summary();

}

