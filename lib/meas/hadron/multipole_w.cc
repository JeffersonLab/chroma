// $Id: multipole_w.cc,v 1.1 2005-03-24 04:32:03 edwards Exp $
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
  void write(BinaryWriter& bin, const Multipole_t& pole)
  {
    // Loop over possible "L" values
    for(int L = 0; L < pole.corr.size(); ++L)
    {
      if ((2*L+1) != pole.corr[L].size())
	QDP_error_exit("Multipole_t: internal error");

      // Loop over possible "M" values
      for(int MM = 0; MM < pole.corr[L].size(); ++MM)
      {
	write(bin, pole.corr[L][MM].electric);
	write(bin, pole.corr[L][MM].magnetic);
      }
    }
  }

  // Anonymous namespace
  namespace 
  {
    //! Compute the r^L*Y_{LM}
    LatticeComplex YLM(int L, int M, int j_decay)
    {
      LatticeComplex ylm = 1.0;   // HACK

      return ylm;
    }
    

    //! Compute the deriv(r^L*Y_{LM},x_k)
    LatticeComplex deriv_YLM(int L, int M, int j_decay, int k)
    {
      LatticeComplex ylm = 1.0;   // HACK

      return ylm;
    }


    //! Compute the electric density
    LatticeSpinMatrix elec_dens(int L, int M, int j_decay)
    {
      SpinMatrix g_one = 1.0;

      // Compute r^L*Y_{LM}(x)*J_4(x) 
      return LatticeSpinMatrix((Gamma(1 << j_decay) * g_one) * YLM(L,M,j_decay));
    }


    //! Compute the magnetic density
    LatticeSpinMatrix mag_dens(int L, int M, int j_decay)
    {
      LatticeSpinMatrix dens;
      SpinMatrix g_one = 1.0;

      multi1d<LatticeReal> x(Nd);
      for(int mu=0; mu < Nd; ++mu)
	x[mu] = Layout::latticeCoordinate(mu);

      // The density is   (\vec{r} \cross \vec{J}(x)) \dot div(r^L*Y_{LM}(x))
      // Here, use local current for \vec{J} so  insertion is  gamma_k
      if (j_decay != Nd-1 || Nd != 4)
      {
	QDPIO::cerr << "Silly implementation limit: expect j_decay=Nd-1 and Nd=4" << endl;
	QDP_abort(1);
      }

      dens  = (x[1]*(Gamma(1 << 2)*g_one) - x[2]*(Gamma(1 << 1)*g_one)) * deriv_YLM(L,M,j_decay,0);
      dens += (x[2]*(Gamma(1 << 0)*g_one) - x[0]*(Gamma(1 << 2)*g_one)) * deriv_YLM(L,M,j_decay,1);
      dens += (x[0]*(Gamma(1 << 1)*g_one) - x[1]*(Gamma(1 << 0)*g_one)) * deriv_YLM(L,M,j_decay,2);

      return dens;
    }
  } // end anonymous namespace


  //! Compute contractions for multipole moments
  /*!
   * \ingroup hadron
   *
   * \param pole               structures holding formfactors ( Write )
   * \param quark_propagator   quark propagator ( Read )
   * \param seq_quark_prop     sequential quark propagator ( Read )
   * \param GammaInsertion     extra gamma matrix insertion ( Read )
   * \param j_decay            direction of decay ( Read )
   * \param t0                 cartesian coordinates of the source ( Read )
   */

  void multipole(Multipole_t& pole,
		 const LatticePropagator& quark_propagator,
		 const LatticePropagator& seq_quark_prop, 
		 int   GammaInsertion,
		 int   max_power,
		 int   j_decay,
		 int   t0)
  {
    START_CODE();

    // Initialize the slow Fourier transform phases
    SftMom phases(0, true, j_decay);

    // Length of lattice in j_decay direction and 3pt correlations fcns
    int length = phases.numSubsets();

    int G5 = Ns*Ns-1;
  
    // Construct the anti-quark propagator from the seq. quark prop.
    LatticePropagator anti_quark_prop = Gamma(G5) * seq_quark_prop * Gamma(G5);

    pole.corr.resize(max_power);

    // Loop over possible "L" values
    for(int L = 0; L < max_power; ++L)
    {
      pole.corr[L].resize(2*L+1);

      // Loop over possible "M" values
      for(int M = -L, MM = 0; M <= L; ++M, ++MM)
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

	pole.corr[L][MM].electric = elec_corr;
	pole.corr[L][MM].magnetic = mag_corr;

      } // end for(M)
    }   // end for(L)
    
    END_CODE();
  }

}  // end namespace Chroma
