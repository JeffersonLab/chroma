// -*- C++ -*-
// $Id: simpleGramSchmidt.h,v 1.1 2007-09-26 02:46:00 kostas Exp $

#ifndef _Simple_GRAMSMIDT_H
#define _Simple_GRAMSMIDT_H

#include "chromabase.h"

namespace Chroma 
{
  //typedef LatticeFermion T ; // save sometyping 

  void SimpleGramSchmidt(multi1d<LatticeFermion>& vec, 
			 const int f,
			 const int t,
			 const OrderedSubset& sub) ;
  
}
#endif


