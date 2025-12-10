// $Id: time_decomp.cc,v 1.1 2008-03-12 00:48:21 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "qdp.h"
#include "unittest.h"

#include "testvol.h"
#include "timeDecomp.h"

using namespace QDP;

int main(int argc, char **argv)
{
  // Initialize UnitTest jig
  TestRunner  tests(&argc, &argv, nrow_in);
  tests.addTest(new timeDecomp(), "timeDslash" );

  // Run all tests
  tests.run();

  // Testjig is destroyed
  tests.summary();

}

