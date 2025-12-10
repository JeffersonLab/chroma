// $Id: time_clover.cc,v 1.1 2008-09-24 20:01:01 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "qdp.h"
#include "unittest.h"

#include "testvol.h"
#include "timeClover.h"

using namespace QDP;

int main(int argc, char **argv)
{
  // Initialize UnitTest jig
  TestRunner  tests(&argc, &argv, nrow_in);
  tests.addTest(new timeClover(), "timeClover" );

  // Run all tests
  tests.run();

  // Testjig is destroyed
  tests.summary();

}

