// $Id: dilute_zN_eigvec_source_const.cc,v 3.2 2009-07-12 00:45:36 jbulava Exp $
/*! \file
 *  \brief Random ZN wall source construction
 */

#include "chromabase.h"
#include "handle.h"

#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"

#include "meas/smear/quark_smearing_aggregate.h"
#include "meas/smear/quark_smearing_factory.h"

#include "meas/smear/quark_displacement_aggregate.h"
#include "meas/smear/quark_displacement_factory.h"

#include "meas/smear/simple_quark_displacement.h"
#include "meas/smear/no_quark_displacement.h"

#include "meas/sources/source_const_factory.h"
#include "meas/sources/dilute_zN_eigvec_source_const.h"
#include "meas/sources/zN_src.h"

#include "meas/inline/io/named_objmap.h"

#include "util/ft/sftmom.h"

#include "util/ferm/subset_vectors.h"
namespace Chroma
{
  // Read parameters
  void read(XMLReader& xml, const string& path, DiluteZNEigVecQuarkSourceConstEnv::Params& param)
  {
    DiluteZNEigVecQuarkSourceConstEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const DiluteZNEigVecQuarkSourceConstEnv::Params& param)
  {
    param.writeXML(xml, path);
  }



  // Hooks to register the class
  namespace DiluteZNEigVecQuarkSourceConstEnv
  {
    // Anonymous namespace
    namespace
    {
      //! Callback function
      QuarkSourceConstruction<LatticeFermion>* createFerm(XMLReader& xml_in,
							  const std::string& path)
      {
	return new SourceConst<LatticeFermion>(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;

      //! Name to be used
      const std::string name("RAND_DILUTE_EIGVEC_ZN_SOURCE");
    }  // end namespace

    //! Return the name
    std::string getName() {return name;}

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheFermSourceConstructionFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }


    //! Initialize
    Params::Params()
    {
    }


    //! Read parameters
    Params::Params(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	break;

      default:
	QDPIO::cerr << __func__ << ": parameter version " << version 
		    << " unsupported." << endl;
	QDP_abort(1);
      }

      read(paramtop, "ran_seed", ran_seed);
      read(paramtop, "N", N);
      read(paramtop, "j_decay", j_decay);
      read(paramtop, "t_sources", t_sources);

      read(paramtop, "eigen_vec_id", eigen_vec_id);
      read(paramtop, "eigen_vectors", eigen_vectors);
      read(paramtop, "spin_mask", spin_mask);
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;

      write(xml, "version", version);
      write(xml, "ran_seed", ran_seed);
      write(xml, "N", N);
      write(xml, "j_decay", j_decay);
      write(xml, "t_sources", t_sources);

      write(xml, "eigen_vec_id", eigen_vec_id);
      write(xml, "eigen_vectors", eigen_vectors);
      write(xml, "spin_mask", spin_mask);

      pop(xml);
    }

	
    void fill_laph_subspace_zN( LatticeLAPHSubSpace_t& laph_in, 
				const Seed& rng_seed, const int& N)
    {

      //Obtain the current seed
      Seed curr_seed;
      QDP::RNG::savern(curr_seed);

      //Seed the random number generator
      QDP::RNG::setrn(rng_seed);
	

      //Fill the struct
      for (int t = 0 ; t < laph_in.time_slices.size(); ++t)
      {
	for (int s = 0 ; s < laph_in.time_slices[t].spins.size() ; ++s)
	  for (int v = 0 ; v < laph_in.time_slices[t].spins[s].lap_eigs.size() ; ++v)
	  {
	    laph_in.time_slices[t].spins[s].lap_eigs[v].val = zN_rng(N);
	  }

	//Return the seed to its previous value
	QDP::RNG::setrn(curr_seed);

      }
    }

    //! Construct the source
    template<>
    LatticeFermion
    SourceConst<LatticeFermion>::operator()(const multi1d<LatticeColorMatrix>& u) const
    {
      QDPIO::cout << "Eigenvector-Diluted random complex ZN source" << endl;

      int Nt = QDP::Layout::lattSize()[params.j_decay];
      
      //
      // Sanity checks
      //
      if (params.spin_mask.size() > Ns)
      {
	QDPIO::cerr << name << ": spin mask size incorrect 1" << endl;
	QDP_abort(1);
      }

      if (params.spin_mask.size() == 0)
      {
	QDPIO::cerr << name << ": spin mask size incorrect 2" << endl;
	QDP_abort(1);
      }

      if (params.t_sources.size() > Nt)
      {
	QDPIO::cerr << name << ": time sources size incorrect 1" << endl;
	QDP_abort(1);
      }

      if (params.t_sources.size() == 0)
      {
	QDPIO::cerr << name << ": time sources size incorrect 2" << endl;
	QDP_abort(1);
      }

      //check that there are no repeats, and that each is < Nt
      for (int t = 0 ; t < params.t_sources.size() ; ++t) {

	if ( (t > 0) && (params.t_sources[t] == params.t_sources[0]) ) {
	  QDPIO::cerr << "ERROR: repeat in t_sources" << endl;
	  QDP_abort(1);
	}
	
	if (params.t_sources[t] >= Nt) {

	  QDPIO::cerr << "ERROR: invalid component in t_sources "  
		      << params.t_sources[t] <<  " >= " << Nt << endl;
	  QDP_abort(1);
	}
      }

      XMLBufferWriter eig_vecs_xml;
      
      //Attempt to get eigenvectors from the named object map
      try {
	
	TheNamedObjMap::Instance().getData< Handle< MapObject<int,EVPair<LatticeColorVector> > > >(params.eigen_vec_id);	
	
	TheNamedObjMap::Instance().get(params.eigen_vec_id).getRecordXML(eig_vecs_xml);
	
      }
      catch( std::bad_cast )  {
	QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	QDP_abort(1);
      }
      catch (const string& e)  {
	QDPIO::cerr << name << ": map call failed: " << e << endl;
	QDP_abort(1);
      }
      
      const MapObject<int,EVPair<LatticeColorVector> >& eigen_vecs = 
	*(TheNamedObjMap::Instance().getData< Handle< MapObject<int,EVPair<LatticeColorVector> > > >(params.eigen_vec_id));

      int n_ev = eigen_vecs.size();
      
      //Sanity checks on the eigenvector dilutions
      int n_ev_dil = params.eigen_vectors.size();
			
      //First, ensure there are not more than n_ev elements
      if ( n_ev_dil > n_ev) {

	QDPIO::cerr << "ERROR: n_ev_dil > n_ev" << endl;
	QDP_abort(1);
      }

      //Check that there are no repeats in the vectors, also that each 
      //element is not larger than n_ev
      for (int v = 0 ; v < n_ev_dil ; ++v) {
	if ( (v > 0) && (params.eigen_vectors[v] == params.eigen_vectors[0]) ) {
	  QDPIO::cerr << "ERROR: repeat in eigen_vectors" << endl;
	  QDP_abort(1);
	}

	if (params.eigen_vectors[v] >= n_ev) {
	  QDPIO::cerr << "ERROR: invalid component in eigen_vectors" << endl;
	  QDP_abort(1);
	}
      }

      //params.writeXML(xml_out, "Input");
      
      //params.writeXML(xml_out, "EigVecsXML");
      
      //
      // Finally, do something useful
      //
      
      /*
      // Save current seed
      Seed ran_seed;
      QDP::RNG::savern(ran_seed);
      
      // Set the seed to desired value
      QDP::RNG::setrn(params.ran_seed);
      */


      SftMom phases(0, true, params.j_decay);
      
      LatticeLAPHSubSpace_t laph_noise(n_ev, Nt);
      fill_laph_subspace_zN(laph_noise, params.ran_seed, params.N);

      QDPIO::cout << "Created LapH Noise " << endl;
      
      LatticeFermion dil_source = zero;
      
      for (int t0 = 0 ; t0 < params.t_sources.size() ; ++t0) {
	int curr_t = params.t_sources[t0];
	
	for (int s = 0 ; s < params.spin_mask.size() ; ++s) {
	  int curr_s = params.spin_mask[s];
	  
	  for (int v = 0 ; v < params.eigen_vectors.size() ; ++v) {
	    int curr_v = params.eigen_vectors[v];
	    
	    const Complex& curr_n = laph_noise.time_slices[curr_t].spins[curr_s].lap_eigs[curr_v].val;
	    
	    LatticeFermion temp = zero;
	    //Should only be poking on a single timeslice, Robert Help!!!
	   
	    EVPair<LatticeColorVector> tvec; eigen_vecs.get(curr_v, tvec);
	    pokeSpin(temp, curr_n * tvec.eigenVector, curr_s);
	    
	    dil_source[phases.getSet()[curr_t]] += temp;
	  }//v
	  
	}//s
	
      }//t0
      
      // Reset the seed
      //QDP::RNG::setrn(ran_seed);
      
      return dil_source;
    }
    
  } // end namespace
  
} // end namespace Chroma
