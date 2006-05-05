// $Id: formfac_w.cc,v 3.1 2006-05-05 03:07:20 edwards Exp $
/*! \file
 *  \brief Form-factors 
 *
 *  Form factors constructed from a quark and a sequential quark propagator
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/formfac_w.h"

namespace Chroma 
{

  /*!
   * Structures for hadron parts
   *
   * \ingroup hadron
   *
   * @{
   */

  // Read a momenta struct
  void read(BinaryReader& bin, FormFac_momenta_t& mom)
  {
    read(bin, mom.magic);
    if (mom.magic != 20301)
    {
      QDPIO::cerr << "read(FormFac_momenta_t): magic number invalid" << endl;
      QDP_abort(1);
    }
    read(bin, mom.inser_mom);
    read(bin, mom.local_current);
    read(bin, mom.nonlocal_current);
  }

  // 
  void read(BinaryReader& bin, FormFac_insertion_t& mes)
  {
    read(bin, mes.gamma_value);
    read(bin, mes.momenta);
  }

  // 
  void read(BinaryReader& bin, FormFac_insertions_t& form)
  {
    read(bin, form.output_version);
    read(bin, form.formFac);
  }

  // Write a momenta struct
  void write(BinaryWriter& bin, const FormFac_momenta_t& mom)
  {
    int magic = 20301;
    write(bin, magic);
    write(bin, mom.inser_mom);
    write(bin, mom.local_current);
    write(bin, mom.nonlocal_current);
  }

  // 
  void write(BinaryWriter& bin, const FormFac_insertion_t& mes)
  {
    write(bin, mes.gamma_value);
    write(bin, mes.momenta);
  }

  // 
  void write(BinaryWriter& bin, const FormFac_insertions_t& form)
  {
    write(bin, form.output_version);
    write(bin, form.formFac);
  }

  /*! @} */  // end of group hadron


  //! Compute contractions for current insertion 3-point functions.
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions!
   *
   * \param form               structures holding formfactors ( Write )
   * \param u                  gauge fields (used for non-local currents) ( Read )
   * \param quark_propagator   quark propagator ( Read )
   * \param seq_quark_prop     sequential quark propagator ( Read )
   * \param gamma_insertion    extra gamma insertion at source ( Read )
   * \param phases             fourier transform phase factors ( Read )
   * \param t0                 cartesian coordinates of the source ( Read )
   */

  void FormFac(FormFac_insertions_t& form,
	       const multi1d<LatticeColorMatrix>& u, 
	       const LatticePropagator& quark_propagator,
	       const LatticePropagator& seq_quark_prop, 
	       int gamma_insertion,
	       const SftMom& phases,
	       int t0)
  {
    START_CODE();

    // Length of lattice in j_decay direction and 3pt correlations fcns
    int length = phases.numSubsets();

    int G5 = Ns*Ns-1;
  
    // Construct the anti-quark propagator from the seq. quark prop.
    LatticePropagator anti_quark_prop = Gamma(G5) * seq_quark_prop * Gamma(G5);

    // Rough timings (arbitrary units):
    //   Variant 1: 120
    //   Variant 2: 140
    // See previous cvs versions (before 1.10) for Variant 2 - only keeping Variant 1

    form.formFac.resize(Nd*Nd);

    // Loop over gamma matrices of the insertion current of insertion current
    for(int gamma_value = 0; gamma_value < Nd*Nd; ++gamma_value)
    {
      //  For the case where the gamma value indicates we are evaluating either
      //  the vector or axial vector currents, we will also evaluate
      //  the non-local currents.  The non-local vector current is the conserved
      //  current.  The non-local axial vector current would be partially
      //  conserved but for the Wilson term.  In these cases we will set
      //  mu = corresponding direction.  In all other cases, we will set mu = -1.

      bool compute_nonlocal;
      int mu;

      switch(gamma_value){
      case  1:
      case 14:
	mu = 0;
	compute_nonlocal = true;
	break;
      case  2:
      case 13:
	mu = 1;
	compute_nonlocal = true;
	break;
      case  4:
      case 11:
	mu = 2;
	compute_nonlocal = true;
	break;
      case  8:
      case  7:
	mu = 3;
	compute_nonlocal = true;
	break;
      default:
	mu = -1;
	compute_nonlocal = false;
      }

      // The local non-conserved vector-current matrix element 
      LatticeComplex corr_local_fn =
	trace(adj(anti_quark_prop) * Gamma(gamma_value) * quark_propagator * Gamma(gamma_insertion));

      multi2d<DComplex> hsum, hsum_nonlocal;
      hsum = phases.sft(corr_local_fn);

      // Construct the non-local current matrix element 
      //
      // The form of J_mu = (1/2)*[psibar(x+mu)*U^dag_mu*(1+gamma_mu)*psi(x) -
      //                           psibar(x)*U_mu*(1-gamma_mu)*psi(x+mu)]
      // NOTE: the 1/2  is included down below in the sumMulti stuff
      LatticeComplex corr_nonlocal_fn;
      if(compute_nonlocal){
	corr_nonlocal_fn =
	  trace(adj(u[mu] * shift(anti_quark_prop, FORWARD, mu)) *
		(quark_propagator + Gamma(gamma_value) * quark_propagator) * Gamma(gamma_insertion));
	LatticePropagator tmp_prop1 = u[mu] *
	  shift(quark_propagator, FORWARD, mu);
	corr_nonlocal_fn -= trace(adj(anti_quark_prop) *
				  (tmp_prop1 - Gamma(gamma_value) * tmp_prop1) * Gamma(gamma_insertion));

	hsum_nonlocal = phases.sft(corr_nonlocal_fn);
      }

  
      form.formFac[gamma_value].gamma_value = gamma_value;
      form.formFac[gamma_value].momenta.resize(phases.numMom());  // hold momenta output

      // Loop over insertion momenta and print out results
      for(int inser_mom_num=0; inser_mom_num<phases.numMom(); ++inser_mom_num) 
      {
	form.formFac[gamma_value].momenta[inser_mom_num].inser_mom = phases.numToMom(inser_mom_num);

	multi1d<Complex> local_cur3ptfn(length); // always compute
	multi1d<Complex> nonlocal_cur3ptfn;
	if (compute_nonlocal)
	  nonlocal_cur3ptfn.resize(length);      // possibly compute

	for (int t=0; t < length; ++t) 
	{
	  int t_eff = (t - t0 + length) % length;

	  local_cur3ptfn[t_eff] = Complex(hsum[inser_mom_num][t]);
	  if (compute_nonlocal)
	    nonlocal_cur3ptfn[t_eff] = 0.5 * Complex(hsum_nonlocal[inser_mom_num][t]);

	} // end for(t)

	form.formFac[gamma_value].momenta[inser_mom_num].local_current    = local_cur3ptfn;
	form.formFac[gamma_value].momenta[inser_mom_num].nonlocal_current = nonlocal_cur3ptfn;

      } // end for(inser_mom_num)
    } // end for(gamma_value)
                            
    END_CODE();
  }

}  // end namespace Chroma
