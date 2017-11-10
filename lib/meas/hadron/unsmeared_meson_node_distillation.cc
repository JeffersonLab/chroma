/*! \file
 * \brief Meson nodes
 */

#include "meas/hadron/unsmeared_meson_node_distillation.h"
#include "meas/smear/displace.h"
#include "util/ft/single_phase.h"

#include "util/ferm/diractodr.h"
#include "util/ferm/twoquark_contract_ops.h"

namespace Chroma
{
 
  //----------------------------------------------------------------------------------
  // Utility functions
  namespace
  {
    //! Error output
    StandardOutputStream& operator<<(StandardOutputStream& os, const multi1d<int>& d)
    {
      if (d.size() > 0)
      {
	os << d[0];

	for(int i=1; i < d.size(); ++i)
	  os << " " << d[i];
      }

      return os;
    }

    //! Error output
    StandardOutputStream& operator<<(StandardOutputStream& os, const std::vector<int>& d)
    {
      if (d.size() > 0)
      {
	os << d[0];

	for(int i=1; i < d.size(); ++i)
	  os << " " << d[i];
      }

      return os;
    }

    //----------------------------------------------------------------------------
    // Lazy conversion
    ComplexD quickConvert(const std::complex<double>& src)
    {
      return cmplx(Real(src.real()),Real(src.imag()));
    }


    // Lazy conversion
    std::complex<double> quickConvert(const ComplexD& src)
    {
      return std::complex<double>(toDouble(real(src)),toDouble(imag(src)));
    }
  }


  //---------------------------------------------------------------------------------- 
  //! Constructor
  UnsmearedMesonNodeDistillation::UnsmearedMesonNodeDistillation(const HadronVertex_t& key_, 
								 const MapSingleHadronQuarkFlavorColorSpin_t& deriv_coeffs,
								 const std::map<char,std::string>& flav_map_,
								 int num_vecs_)
    : key(key_), flav_map(flav_map_), num_vecs(num_vecs_)
  {
    // Whether this is mesons (derivatives) or amplitudes (displacements)
    // we set the derivative type
    if (key.type == 'M')
      use_derivP = true;
    else if (key.type == 'A')
      use_derivP = false;
    else
    {
      QDPIO::cerr << __func__ << ": unsupported node type = " << key.type << "\n";
      exit(1);
    }

    //
    // Rotation from DR to DP
    //
    SpinMatrix diracToDRMat(DiracToDRMat());
    diracToDrMatPlus = convertTwoQuarkSpin(diracToDRMat);

    // Finish the gamma_5 transformation of the left solution vectors with a left multiplication of
    // gamma_5 here. Note, the right multiplication of the gamma_5 is done in the insertion.
    // The adjoint in done in the lattice inner product. These gamma_5's are in the DR basis
    diracToDrMatMinus = convertTwoQuarkSpin(adj(diracToDRMat) * Gamma(15));

    //
    // Convert the spin basis from Dirac-Pauli (DP) to Degrand-Rossi (DR) since the props are in that basis
    //
    // Loop over each unique derivative
    std::vector<KeySingleHadronQuarkFlavorColor_t> k1;
    deriv_coeffs.keys(k1);

    // Sanity check - do not expect more than one (generic) flavor key
    if (k1.size() != 1)
    {
      QDPIO::cerr << __func__ << ": internal error - expected 1 flavor key\n";
      QDP_abort(1);
    }

    // Now go to the next level - derivs
    std::vector<KeySingleHadronQuarkDeriv_t>  k2;
    std::vector<MapSingleHadronQuarkSpin_t>   v2;
    deriv_coeffs[k1[0]].keysAndValues(k2,v2);

    for(int d = 0; d < k2.size(); ++d)
    {
      // Get the spin terms
      std::vector<KeySingleHadronQuarkSpin_t>  k3;
      std::vector<ComplexD>                    v3;
      v2[d].keysAndValues(k3,v3);

      // Re-populate a spin matrix
      SpinMatrix spin_dp = zero;

      for(int j = 0; j < k3.size(); ++j)
      {
	// Spin elements
	int ls = k3[j].diracs[1];
	int rs = k3[j].diracs[2];
	Complex weight = v3[j];
	
	pokeSpin(spin_dp, weight, ls, rs);
      }

      // Now rotate from DP to DR
      SpinMatrix spin_dr = diracToDRMat * spin_dp * adj(diracToDRMat);

      // We know the left vectors will be gamma_5 transformed. 
      // Do the right multiplication with gamma_5 onto the left solution vectors by instead multiplying
      // the gamma_5 into the spin components of the insertion
      SpinMatrix g5_spin_dr = Gamma(Ns*Ns-1) * spin_dr;

      // Convert to a sparse version in DR
      std::vector<MatrixSpinRep_t> sparse_g5_spin_dr = convertTwoQuarkSpin(g5_spin_dr);

      // Lowest level - dirac spins
      MapSingleHadronQuarkSpin_t dirac_map;

      // Loop over each spin index. The derivative is unaffected.
      for(std::vector<MatrixSpinRep_t>::const_iterator mm = sparse_g5_spin_dr.begin();
	  mm != sparse_g5_spin_dr.end();
	  ++mm)
      {
	dirac_map.insert(KeySingleHadronQuarkSpin_t(mm->left, mm->right), mm->op);
      }

      // Second level - derivs
      MapSingleHadronQuarkDeriv_t deriv_map;
      deriv_map.insert(k2[d], dirac_map);

      // Third level - flavor
      MapSingleHadronQuarkFlavorColorSpin_t tmp_flavor_map;
      tmp_flavor_map.insert(k1[0], deriv_map);

      // Now, add this into the final map
      coeffs = coeffs + tmp_flavor_map;

    } // for d

  }


  //----------------------------------------------------------------------------------
  //! Contract to form a hadron node
  void UnsmearedMesonNodeDistillation::constructNode(Hadron::HadronDistOperatorRep& obj, 
						     const KeyHadronNode_t& hadron_node, 
						     DistillationSolnCache& soln_cache,
						     DispSolnCache& disp_cache) const
  {
    QDPIO::cout << __func__ << ": entering" << std::endl;

    StopWatch swatch;
    swatch.reset(); swatch.start();

    // For the moment, this node is fixed to be square
    const int left_num_vecs  = num_vecs;
    const int right_num_vecs = num_vecs;

    // Basic initialization
    obj.setNumDils(left_num_vecs, right_num_vecs);

    // Sanity checks
    if (hadron_node.quarks.size() != 2)
    {
      QDPIO::cerr << __func__ << ": expected key size for meson op to be 2\n";
      QDP_abort(1);
    }

    // Sanity checks
    if (hadron_node.quarks[1].t_slice != hadron_node.quarks[2].t_slice)
    {
      QDPIO::cerr << __func__ << ": expect time-slice of quarks to be the same\n";
      QDP_abort(1);
    }

    // Flip momentum in case of a creation op
    multi1d<int> mom;
    if (hadron_node.vertex.creation_op)
      mom = -1 * hadron_node.vertex.mom;
    else
      mom =  hadron_node.vertex.mom;


    // Flip momentum to agree with FT convention \int dx exp(-ipx)
    LatticeComplex phase = singlePhase(-1 * mom);

    // Fix this.
    multi1d<int> mom_to_use = -1 * mom;

    //
    // Read in all the left vectors and hold them in memory
    //
    multi2d< multi1d<LatticeColorVector> > left_solns(Ns,Ns);

    for(int spin_l=0; spin_l < Ns; ++spin_l)
    {
      for(int ls=0; ls < Ns; ++ls)
      {
	left_solns(ls,spin_l).resize(left_num_vecs);

	for(int dl = 0; dl < left_num_vecs; ++dl)
	{
	  // No displacements on left vectors
	  left_solns(ls,spin_l)[dl] = soln_cache.getVector(hadron_node.quarks[1], dl, ls, spin_l);
	}
      }
    }

    //
    // Loop over each unique derivative
    //
    std::vector<KeySingleHadronQuarkFlavorColor_t> k1;
    coeffs.keys(k1);

    // Sanity check - do not expect more than one (generic) flavor key
    if (k1.size() != 1)
    {
      QDPIO::cerr << __func__ << ": internal error - did not expect more than 1 flavor key\n";
      QDP_abort(1);
    }

    // Now go to the next level - derivs
    std::vector<KeySingleHadronQuarkDeriv_t>  k2;
    std::vector<MapSingleHadronQuarkSpin_t>   v2;
    coeffs[k1[0]].keysAndValues(k2,v2);

    //
    // Big loop over derivatives
    //
    for(int d = 0; d < k2.size(); ++d)
    {
      // Get the spin terms
      std::vector<KeySingleHadronQuarkSpin_t>  k3;
      std::vector<ComplexD>                    v3;
      v2[d].keysAndValues(k3,v3);

      // Scan the right spin coefficients of the insertion and see which are active
      multi1d<bool> active(Ns);
      active = false;
      for(int j = 0; j < k3.size(); ++j)
      {
	// Spin elements. These are in the DR basis - same as the props.
	int rs = k3[j].diracs[2];
	active[rs] = true;
      }

      
      //
      // Take a derivative of the right vector
      // Fold in phase which is actually independent of the color components, but still saves work.
      //
      multi2d< multi1d<LatticeColorVector> > right_solns(Ns,Ns);

      for(int spin_r=0; spin_r < Ns; ++spin_r)
      {
	for(int rs=0; rs < Ns; ++rs)
	{
	  // Only need to do all this work for the right spin if we know we are going to use them
	  if (! active[rs]) {continue;}

	  // Do the deed
	  right_solns(rs,spin_r).resize(right_num_vecs);

	  for(int dr = 0; dr < right_num_vecs; ++dr)
	  {
	    // Left-right derivative is turned into a displacement function acting on the right vector
	    right_solns(rs,spin_r)[dr] = phase * disp_cache.getDispVector(use_derivP, mom_to_use,
									  hadron_node.quarks[2], dr, rs, spin_r,
									  k2[d].derivs[2].deriv);
	  } // dr
	} // rs
      } // spin_r


      //
      // Loop over spins of the insertion
      //
      for(int j = 0; j < k3.size(); ++j)
      {
	// Spin elements. These are in the DR basis - same as the props.
	int ls = k3[j].diracs[1];
	int rs = k3[j].diracs[2];
	ComplexD weight = v3[j];

	// Big loop of outer spin indices. We know this is dense.
	for(int spin_r=0; spin_r < Ns; ++spin_r)
	{
	  for(int spin_l=0; spin_l < Ns; ++spin_l)
	  {
	    Hadron::HadronDistOperatorRep::Spin_t& obj_spin = obj.spin(spin_l, spin_r);

	    // Loop over right dilution
	    for(int dr = 0; dr < right_num_vecs; ++dr)
	    {
	      // Loop over left dilution
	      for(int dl = 0; dl < left_num_vecs; ++dl)
	      {
		// Inner product over color indices. 
		// Accumulate with a phase and weight.
		obj_spin.dilution.at(dl,dr) += quickConvert(weight * innerProduct(left_solns(ls,spin_l)[dl], right_solns(rs,spin_r)[dr]));
	      } // dl
	    } // dr
	  } // j
	} // spin_l
      } // spin_r
    } // d


    //
    // Rotate from DR to DP basis
    //
    Hadron::HadronDistOperatorRep obj_tmp;
    multiplyRep(obj_tmp, diracToDrMatMinus, obj);
    obj.clear();
    multiplyRep(obj, obj_tmp, diracToDrMatPlus);


#if 1
    if (1)
    {
      // Big loop of outer spin indices. We know this is dense.
      for(int spin_r=0; spin_r < Ns; ++spin_r)
      {
	for(int spin_l=0; spin_l < Ns; ++spin_l)
	{
	  Hadron::HadronDistOperatorRep::Spin_t& obj_spin = obj.spin(spin_l, spin_r);

	  // Loop over right dilution
	  for(int dr = 0; dr < right_num_vecs; ++dr)
	  {
	    // Loop over left dilution
	    for(int dl = 0; dl < left_num_vecs; ++dl)
	    {
	      QDPIO::cout << __func__ << ": deriv= " << k2[0].derivs[2].deriv
			  << " dl= " << dl << " dr= " << dr << " sl= " << spin_l << " sr= " << spin_r
			  << "  v= " << quickConvert(obj_spin.dilution.at(dl,dr)) << std::endl;
	    } // dl
	  } // dr
	} // spin_l
      } // spin_r
    }
#endif


    // Had enough
    swatch.stop(); 
    QDPIO::cout << __func__ << ": time = " << swatch.getTimeInSeconds() << " secs " << std::endl;

  } // constructNode
    
} // namespace Chroma
