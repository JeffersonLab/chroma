// $Id: t_dslash.cc,v 1.9 2008-03-04 21:50:19 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "qdp.h"
#include "unittest.h"

#include "testvol.h"
#include "testDslashFull.h"
#include "testSiteDslash.h"

using namespace QDP;

int main(int argc, char **argv)
{
  // Initialize UnitTest jig
  TestRunner  tests(&argc, &argv, nrow_in);
  tests.addTest(new testDslashFull(), "testDslashFull" );

  // Run all tests
  tests.run();

  // Testjig is destroyed
  tests.summary();

}

