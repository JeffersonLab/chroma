// $Id: multipole_w.cc,v 1.4 2005-04-04 03:22:40 edwards Exp $
/*! \file
 *  \brief Multipole moments
 *
 *  Compute multipole moment within two hadron states using spherical
 *  harmonic expansion
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/multipole_w.h"

namespace Chroma 
{

  // Write a Multipole_t
  void write(XMLWriter& xml, const string& path, const Multipole_t::ElecMag_t& pole)
  {
    push(xml, path);

    write(xml, "L", pole.L);
    write(xml, "M", pole.M);
    write(xml, "electric", pole.electric);
    write(xml, "magnetic", pole.magnetic);
    
    pop(xml);
  }


  // Write a Multipole_t
  void write(XMLWriter& xml, const string& path, const Multipole_t& pole)
  {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);
    write(xml, "j_decay", pole.j_decay);
    write(xml, "Harmonic", pole.corr);

    pop(xml);
  }

  // Anonymous namespace
  namespace 
  {

    //! Dumb factorial function
    int factorial(int m)
    {
      return (m <= 0) ? 1 : m*factorial(m-1);
    }

    //! Compute the P_{lm} Legendre polynomial
    LatticeReal plgndr(int l, int m, const LatticeReal& x)
    {
//      if (m < 0 || m > l || fabs(x) > 1.0)
      if (m < 0 || m > l)
	QDP_error_exit("Bad arguments in plgndr");

      LatticeReal pmm=1.0;
      if (m > 0) 
      {
	LatticeReal somx2 = sqrt((1.0-x)*(1.0+x));
	LatticeReal fact = 1.0;
	for (int i=1; i <= m; i++) 
	{
	  pmm *= -fact*somx2;
	  fact += 2.0;
	}
      }
      if (l == m)
	return pmm;
      else 
      {
	LatticeReal pmmp1 = x*(2*m+1)*pmm;
	if (l == (m+1))
	  return pmmp1;
	else 
	{
	  LatticeReal pll;
	  for (int ll=m+2; ll <= l; ll++) 
	  {
	    pll = (x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
	    pmm = pmmp1;
	    pmmp1 = pll;
	  }
	  return pll;
	}
      }
    }



    //! Compute  r^L*Y_{LM}
    LatticeComplex rYLM(const multi1d<LatticeReal>& x, int l, int m)
    {
      LatticeComplex ylm;

      if (m < -l || m > l)
	QDP_error_exit("Bad arguments in rYLM");

      if (m < 0)
      {
	// Y_{l,-m} = (-1)^m * conj(Y_{l,m})
	int ss = ((-m) & 1 == 1) ? -1 : 1;
	ylm = ss * conj(rYLM(x,l,-m));
      }
      else
      {
	// Build radius squared
	LatticeReal r_sq = zero;
	for(int j=0; j < x.size(); ++j)
	  r_sq += x[j]*x[j];

	LatticeReal r = where(r_sq > 0.0, sqrt(r_sq), LatticeReal(1));

	// z = r*cos(theta)
	LatticeReal cos_theta = x[2] / r;
	
	// y = r*sin(phi),  x = r*cos(phi),  y/x = tan(phi)
	LatticeReal phi = atan2(LatticeReal(x[1]),LatticeReal(x[0]));

	// Normalization
	Real cnst = sqrt((2*l+1)*factorial(l-m)/(2*twopi*factorial(l+m)));

	// Ylm = const*P_{LM}(cos(theta))*exp(i*m*phi)
	ylm = cnst * plgndr(l,m,cos_theta) * cmplx(cos(m*phi),sin(m*phi));
      }

      return ylm;
    }
    

    //! Compute the deriv(r^L*Y_{LM},x_k)
    LatticeComplex deriv_rYLM(const multi1d<LatticeReal>& x, int l, int m, int k)
    {
      LatticeComplex ylm = 1.0;   // HACK

      return ylm;
    }


    //! Compute the electric density
    LatticeSpinMatrix elec_dens(int L, int M, int j_decay)
    {
      SpinMatrix g_one = 1.0;

      multi1d<LatticeReal> x(Nd-1);
      for(int mu=0, j=0; mu < Nd; ++mu)
      {
	if (mu == j_decay) continue;
	
	x[j++] = Layout::latticeCoordinate(mu);
      }

      // Compute r^L*Y_{LM}(x)*J_4(x) 
      return LatticeSpinMatrix((Gamma(1 << j_decay) * g_one) * rYLM(x,L,M));
    }


    //! Compute the magnetic density
    LatticeSpinMatrix mag_dens(int L, int M, int j_decay)
    {
      LatticeSpinMatrix dens;
      SpinMatrix g_one = 1.0;

      multi1d<LatticeReal> x(Nd-1);
      multi1d<int>         g(Nd-1);
      for(int mu=0, j=0; mu < Nd; ++mu)
      {
	if (mu == j_decay) continue;
	
	x[j] = Layout::latticeCoordinate(mu);
	g[j] = 1 << mu;
	j++;
      }

      // The density is   (\vec{r} \cross \vec{J}(x)) \dot div(r^L*Y_{LM}(x))
      // Here, use local current for \vec{J} so  insertion is  gamma_k
      dens  = (x[1]*(Gamma(g[2])*g_one) - x[2]*(Gamma(g[1])*g_one)) * deriv_rYLM(x,L,M,0);
      dens += (x[2]*(Gamma(g[0])*g_one) - x[0]*(Gamma(g[2])*g_one)) * deriv_rYLM(x,L,M,1);
      dens += (x[0]*(Gamma(g[1])*g_one) - x[1]*(Gamma(g[0])*g_one)) * deriv_rYLM(x,L,M,2);

      return dens;
    }
  } // end anonymous namespace


  //! Compute contractions for multipole moments
  /*!
   * \ingroup hadron
   *
   * \param quark_propagator   quark propagator ( Read )
   * \param seq_quark_prop     sequential quark propagator ( Read )
   * \param GammaInsertion     extra gamma matrix insertion ( Read )
   * \param max_power          max value of L ( Read )
   * \param j_decay            direction of decay ( Read )
   * \param t0                 cartesian coordinates of the source ( Read )
   * \param xml                xml file object ( Write )
   * \param xml_group          string used for writing xml data ( Read )
   */

  void multipole(const LatticePropagator& quark_propagator,
		 const LatticePropagator& seq_quark_prop, 
		 int   GammaInsertion,
		 int   max_power,
		 int   j_decay,
		 int   t0,
		 XMLWriter& xml,
		 const string& xml_group)
  {
    START_CODE();

    if (Nd != 4)
      QDP_error_exit("%s only designed and meaningful for Nd=4",__func__);

    if (j_decay < 0 || j_decay >= Nd)
      QDP_error_exit("Invalid args to multipole");

    // Big struct of data
    Multipole_t pole;
    pole.j_decay = j_decay;

    // Initialize the slow Fourier transform phases
    SftMom phases(0, true, j_decay);

    // Length of lattice in j_decay direction and 3pt correlations fcns
    int length = phases.numSubsets();

    int G5 = Ns*Ns-1;
  
    // Construct the anti-quark propagator from the seq. quark prop.
    LatticePropagator anti_quark_prop = Gamma(G5) * seq_quark_prop * Gamma(G5);

    pole.corr.resize(max_power+1);

    // Loop over possible "L" values
    for(int L = 0; L < pole.corr.size(); ++L)
    {
      pole.corr[L].resize(2*L+1);

      // Loop over possible "M" values
      for(int M = -L, MM = 0; MM < pole.corr[L].size(); ++M, ++MM)
      {
	// The electric multipole moment
	// The density is   r^L*Y_{LM}(x)*J_4(x) 
	// Here, use local current for J_4 so  insertion is  gamma_{j_decay}
	multi1d<DComplex> hsum_elec;
	{
	  LatticeComplex corr_fn = 
	    trace(adj(anti_quark_prop) * elec_dens(L,M,j_decay) * 
		  quark_propagator * Gamma(GammaInsertion));
	  
	  hsum_elec = sumMulti(corr_fn, phases.getSet());
	}

	// The magnetic multipole moment
	// The density is   (\vec{r} \times \vec{J}(x)) \dot div(r^L*Y_{LM}(x))
	// Here, use local current for \vec{J} so  insertion is  gamma_k
	multi1d<DComplex> hsum_mag;
	{
	  LatticeComplex corr_fn =
	    trace(adj(anti_quark_prop) * mag_dens(L,M,j_decay) * 
		  quark_propagator * Gamma(GammaInsertion));

	  hsum_mag = sumMulti(corr_fn, phases.getSet());
	}

	multi1d<Complex> elec_corr(length);
	multi1d<Complex> mag_corr(length);

	for (int t=0; t < length; ++t) 
	{
	  int t_eff = (t - t0 + length) % length;
	  elec_corr[t_eff] = Complex(hsum_elec[t]);
	  mag_corr[t_eff]  = Complex(hsum_mag[t]);
	} // end for(t)

	pole.corr[L][MM].L        = L;
	pole.corr[L][MM].M        = M;
	pole.corr[L][MM].electric = elec_corr;
	pole.corr[L][MM].magnetic = mag_corr;

      } // end for(M)
    }   // end for(L)
   
    // Write the pole struct
    write(xml, xml_group, pole);

    END_CODE();
  }

}  // end namespace Chroma
