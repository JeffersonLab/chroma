// $Id: gramschm.cc,v 1.1 2004-01-04 21:56:04 edwards Exp $
/*! \file
 *  \brief Gramm-Schmidt orthogonolization
 */

//! Gramm-Schmidt orthogonolization
/*!
 * \ingroup eig
 *
 * Arguments:
 *  \param psi         Pseudofermion field     	       (Modify)
 *  \param vec         subspace wrt orthog     	       (Read)
 *  \param Nvec        Number of vectors               (Read)
 *  \param Npsi        Number of source vectors        (Read) 
 */

void GramSchm(multi1d<LatticeFermion>& psi, int Npsi, 
	      multi1d<LatticeFermion>& vec, int Nvec)
{
  START_CODE("GramSchm");

  for(int s = 0; s < Npsi; ++s)
  {
    for(int i = 0; i < Nvec; ++i)
    {
      Complex xp = innerProduct(vec[i], psi[s]);
      psi[s] -= vec[i] * xp;
    }
  }
        
  END_CODE("GramSchm");
}
