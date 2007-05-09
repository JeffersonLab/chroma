#ifndef IMPROVEMENT_TERMS_S
#define IMPROVEMENT_TERMS_S


#include "chromabase.h"

namespace Chroma 
{ 
  void Fat7_Links(multi1d<LatticeColorMatrix>& u, multi1d<LatticeColorMatrix>& u_fat, Real u0);

  void Triple_Links(multi1d<LatticeColorMatrix>& u, multi1d<LatticeColorMatrix>& u_triple, Real u0);

/* Pass parameters to the fat link code **/

class  fat7_param
{
public :
  Real c_1l ;
  Real c_3l ;
  Real c_5l ;
  Real c_7l ;
  Real c_Lepage ;

} ;



void Fat7_Links(multi1d<LatticeColorMatrix> & u,
		multi1d<LatticeColorMatrix> & uf,
		fat7_param & pp) ;

} // End Namespace Chroma


#endif

