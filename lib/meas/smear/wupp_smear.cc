// $Id: wupp_smear.cc,v 3.2 2009-07-22 02:44:04 edwards Exp $
/*! \file
 *  \brief 3d Laplacian solution on color vector
 */

#error "THIS CODE IS NOT YET READY"

#include "chromabase.h"
#include "meas/smear/wupp_smear.h"
#include "util/ft/sftmom.h"

namespace Chroma 
{

  //! Do a covariant Gaussian smearing of a lattice field
  /*!
   * Arguments:
   *
   *  \param u        gauge field ( Read )
   *  \param chi      color vector field ( Modify )
   *  \param width    width of "shell" wave function ( Read )
   *  \param ItrGaus  number of iterations to approximate Gaussian ( Read )
   *  \param j_decay  direction of decay ( Read )
   */

  template<typename T>
  void gausSmear(const multi1d<LatticeColorMatrix>& u, 
		 T& chi, 
		 const Real& width, int ItrGaus, int j_decay)
  {
    T psi;

    Real ftmp = - (width*width) / Real(4*ItrGaus);
    /* The Klein-Gordon operator is (Lapl + mass_sq), where Lapl = -d^2/dx^2.. */
    /* We want (1 + ftmp * Lapl ) = (Lapl + 1/ftmp)*ftmp */
    Real ftmpi = Real(1) / ftmp;
  
    for(int n = 0; n < ItrGaus; ++n)
    {
      psi = chi * ftmp;
      klein_gord(u, psi, chi, ftmpi, j_decay);
    }
  }


  //! Do a covariant Gaussian smearing of a lattice color vector field
  /*! This is a wrapper over the template definition
   *
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param u        gauge field ( Read )
   *  \param chi      color vector field ( Modify )
   *  \param width    width of "shell" wave function ( Read )
   *  \param ItrGaus  number of iterations to approximate Gaussian ( Read )
   *  \param j_decay  direction of decay ( Read )
   */

  void gausSmear(const multi1d<LatticeColorMatrix>& u, 
		 LatticeColorVector& chi, 
		 const Real& width, int ItrGaus, int j_decay)
  {
    gausSmear<LatticeColorVector>(u, chi, width, ItrGaus, j_decay);
  }




  //! Do a covariant Wuppertal smearing of a color vector field
  /*  u -- gauge field ( Read ) */
  /*  chi -- color vector field ( Modify ) */
  /*  mass_sq  -- mass_sq of Wuppertal "shell" wave function ( Read ) */
  /*  ItrMax -- maximal number of iterations to invert ( Read ) */
  /*  j_decay  -- direction of decay ( Read ) */
  /*  RsdCG  -- residue for CG inverter ( Read ) */

  void wuppSmear(const multi1d<LatticeColorMatrix>& u, 
		 LatticeColorVector& chi, 
		 const Real& mass_sq, int ItrMax, int j_decay, const Real& RsdCG)
  {
    LatticeColorVector p;
    LatticeColorVector r;
    LatticeColorVector ap;
    LatticeReal apa;
    LatticeReal tmp;
    LatticeInteger t_coord;
    LatticeBoolean t_mask;
    Real Rsd;
    multi1d<Real> zero(length);
    multi1d<Real> ab(length);
    multi1d<Double> c(length);
    multi1d<Double> d(length);
    multi1d<Double> dd(length);
    multi1d<Double> chi_norm(length);
    multi1d<Double> cp(length);
    int length;
    int t;
    int k;
    int cb;
    int converged;
    multi1d<Boolean> is_zero(length);
    Boolean any_zero;

    START_CODE();

    length = nrow[j_decay];
    Rsd = RsdCG * RsdCG;

    /* The Klein-Gordon operator is (Lapl + mass_sq), where Lapl = -d^2/dx^2.. */
    /* Use CG to compute psi = (Lapl + mass_sq)^(-1) chi */
    /* Note: Lapl is a 3-d laplacian; each time slice is handled independently. */

                              
    zero = 0;

    /* chi_norm = |chi[0]|^2 */
    chi_norm = 0;
    for(cb = 0; cb < Nsubl; cb++)
      SLICE_SUMSQ(chi(cb), j_decay, chi_norm, ADD);

    dd = 0;
    is_zero = chi_norm == dd;
    FILL(any_zero,FALSE);
    for(t = 0; t < length; t++)
      any_zero = any_zero | is_zero[t];

    /* p[1] = r[0] = chi */
    r = chi;
    p = r;

    /* Cp = |r[0]|^2 = chi_norm */
    cp = chi_norm;

    for(t = 0; t < length; t++)
      chi_norm[t] *= Rsd;

    /* We overwrite chi with the solution; initial solution is zero */
    chi = 0;

    t_coord = Layout::latticeCoordinate(j_decay);

    for(k = 1; k <= ItrMax; k++)
    {
      /* c = |r[k-1]|^2 */
      c = cp;

      /* a[k] = |r[k-1]|^2 / < p[k], Ap[k] > ; */
      /* First compute ap = Ap[k] */
      klein_gord (u, p, ap, mass_sq, j_decay);

      apa = real(trace(adj[p[0]] * ap[0]))
	tmp = real(trace(adj[p[1]] * ap[1]))
	apa += tmp;

      d = sumMulti(apa, timeslice)

	dd = c / d;
      ab = FLOAT(dd);
      if( any_zero )
	copymask(ab, is_zero, zero, REPLACE);

      /* Chi[k] += a[k] p[k] */
      /* SLICE_FILL(apa, ab); */
      apa = 0;
      for(t = 0; t < length; t++)
	if( ! is_zero[t] )
	{
	  FILL(tmp, ab(t));
	  t_mask = t_coord == t;
	  copymask(apa, t_mask, tmp, REPLACE);
	}
      for(cb = 0; cb < Nsubl; cb++)
	chi[cb] += p[cb] * apa;
      /*  r[k] -= a[k] Ap[k] */
      for(cb = 0; cb < Nsubl; cb++)
	r[cb] -= ap[cb] * apa;

      /* cp = |r[k]|^2 */
      cp = 0;
      for(cb = 0; cb < Nsubl; cb++)
	SLICE_SUMSQ(r(cb), j_decay, cp, ADD);

      converged = 0;
      for(t = 0; t < length; t++)
	if( cp[t] <= chi_norm[t] )
	  converged += 1;

      /* push(xml_out,"Residues");
	 write(xml_out, "k", k);
	 write(xml_out, "cp", cp);
	 pop(xml_out); */

      if( converged == length )
      {
	push(xml_out,"Wupp_smear_needed");
	write(xml_out, "k", k);
	pop(xml_out);
	return;
      }

      /* b[k+1] = |r[k]|^2 / |r[k-1]|^2 */
      dd = cp / c;
      ab = FLOAT(dd);
      if( any_zero )
	copymask(ab, is_zero, zero, REPLACE);

      /* p[k+1] = r[k] + b[k+1] p[k] */
      /* SLICE_FILL(apa, ab); */
      apa = 0;
      for(t = 0; t < length; t++)
	if( ! is_zero[t] )
	{
	  FILL(tmp, ab(t));
	  t_mask = t_coord == t;
	  copymask(apa, t_mask, tmp, REPLACE);
	}
      for(cb = 0; cb < Nsubl; cb++)
      {
	p[cb] = p[cb] * apa;
	p[cb] += r[cb];
      }
    }

    push(xml_out,"Wupp_smear_not_converged");
    write(xml_out, "k", k);
    pop(xml_out);
                              
    END_CODE();
  }

} // namespace Chroma
