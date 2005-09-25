// $Id: ischiral_w.h,v 2.0 2005-09-25 21:04:33 edwards Exp $
#ifndef IS_CHIRAL_VECTOR_H
#define IS_CHIRAL_VECTOR_H

#include <chromabase.h>

namespace Chroma {

enum Chirality { CH_MINUS=-1, CH_NONE=0, CH_PLUS=1};

Chirality isChiralVector( const LatticeFermion& chi);

}  // end namespace Chroma


#endif
