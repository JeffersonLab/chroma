// $Id: t_site_ops_parscalar.cc,v 1.1 2008-08-26 13:24:29 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "qdp.h"
#include "unittest.h"

#include "testDslashFull.h"
#include "testDecomp.h"
#include "testDecompHvv.h"
#include "testMvvRecons.h"
#include "testRecons.h"

using namespace QDP;

int main(int argc, char **argv)
{
  // Initialize QDP++ with argc, and argv. Set Lattice Dimensions
  const int latdims[] = {2,4,6,8};

  // Initialize UnitTest jig
  TestRunner  tests(&argc, &argv, latdims);
  tests.addTest(new testDecomp0Minus(), "testDecomp0Minus" );
  tests.addTest(new testDecomp1Minus(), "testDecomp1Minus" );
  tests.addTest(new testDecomp2Minus(), "testDecomp2Minus" );
  tests.addTest(new testDecomp3Minus(), "testDecomp3Minus" );
  tests.addTest(new testDecomp0Plus(), "testDecomp0Plus" );
  tests.addTest(new testDecomp1Plus(), "testDecomp1Plus" );  
  tests.addTest(new testDecomp2Plus(), "testDecomp2Plus" );
  tests.addTest(new testDecomp3Plus(), "testDecomp3Plus" );
  tests.addTest(new testDecompHvv0Plus(), "testDecompHvv0Plus");

  tests.addTest(new testDecompHvv1Plus(), "testDecompHvv1Plus");
  tests.addTest(new testDecompHvv2Plus(), "testDecompHvv2Plus");
  tests.addTest(new testDecompHvv3Plus(), "testDecompHvv3Plus");
  tests.addTest(new testDecompHvv0Minus(), "testDecompHvv0Minus");
  tests.addTest(new testDecompHvv1Minus(), "testDecompHvv1Minus");
  tests.addTest(new testDecompHvv2Minus(), "testDecompHvv2Minus");
  tests.addTest(new testDecompHvv3Minus(), "testDecompHvv3Minus");


  tests.addTest(new testMvvRecons0Plus(), "testMvvRecons0Plus");
  tests.addTest(new testMvvRecons1PlusAdd(), "testMvvRecons1PlusAdd");
  tests.addTest(new testMvvRecons2PlusAdd(), "testMvvRecons2PlusAdd");
  tests.addTest(new testMvvRecons2PlusAddStore(), "testMvvRecons2PlusAddStore");
  tests.addTest(new testMvvRecons3PlusAddStore(), "testMvvRecons3PlusAddStore");


  tests.addTest(new testMvvRecons0Minus(), "testMvvRecons0Minus");
  tests.addTest(new testMvvRecons1MinusAdd(), "testMvvRecons1MinusAdd");
  tests.addTest(new testMvvRecons2MinusAdd(), "testMvvRecons2MinusAdd");
  tests.addTest(new testMvvRecons2MinusAddStore(), "testMvvRecons2MinusAddStore");
  tests.addTest(new testMvvRecons3MinusAddStore(), "testMvvRecons3MinusAddStore");

  tests.addTest(new testRecons4DirPlus(), "testReconsDir4Plus");
  tests.addTest(new testRecons4DirMinus(), "testReconsDir4Minus");

  // Run all tests
  tests.run();

  // Testjig is destroyed
  tests.summary();

}

