// -*- C++ -*-
// $Id: mesplq.h,v 1.4 2005-01-14 18:42:35 edwards Exp $

#ifndef __mesplq_h__
#define __mesplq_h__

namespace Chroma {

void MesPlq(const multi1d<LatticeColorMatrix>& u, Double& w_plaq, Double& s_plaq, 
	    Double& t_plaq, Double& link);

}  // end namespace Chroma

#endif
