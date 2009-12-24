// -*- C++ -*-
// $Id: dilute_zN_eigvec_source_const.h,v 3.2 2009-07-12 00:45:36 jbulava Exp $
/*! \file
 *  \brief Random Z(N) source construction using dilution in eigenvector 
 *  space
 *
 */

#ifndef __dilute_zN_eigvec_source_const_h__
#define __dilute_zN_eigvec_source_const_h__

#include "meas/sources/source_construction.h"
#include "io/xml_group_reader.h"

namespace Chroma
{

  //! Dilute Z(N) quark source namespace, parameters, and classes
  /*! @ingroup sources */
  namespace DiluteZNEigVecQuarkSourceConstEnv
  {
    bool registerAll();

    //! Return the name
    std::string getName();
  
    //! Random complex Z(N) sources using dilution
    /*! @ingroup sources */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    

      Seed            ran_seed;             /*!< Set the seed to this value */
      int             N;                    /*!< Z(N) */
      
      std::string     eigen_vec_id;         /*! The ID of the eigenvectors in the named object map*/
      multi1d<int>    eigen_vectors;        /*!< Eigenvectors which have support on this dilution projector */
      multi1d<int>    spin_mask;            /*!< Spins which have support on this dilution projector */
      int             j_decay;              /*!< decay direction */
     
      multi1d<int>    t_sources;            /*!< timeslices on which this source has non-zero support */
    };



    struct LatticeLAPHSubSpace_t
    {
			
      LatticeLAPHSubSpace_t(int nev, int nt)
	{
				
	  time_slices.resize(nt);
	  for (int t = 0 ; t < nt ; ++t)
	  {
	    time_slices[t].spins.resize(Ns);
	    for (int s = 0 ; s < Ns ; ++s)
	    {
	      time_slices[t].spins[s].lap_eigs.resize(nev);
	    }
	  }
	}

			
      struct Timeslice_t
      {
	struct Spin_t
	{
	  struct LapEig_t
	  {
	    Complex val;
	  };

	  multi1d<LapEig_t> lap_eigs;

	};

	multi1d<Spin_t> spins;

      };

      multi1d<Timeslice_t> time_slices;

    };


    void fill_laph_subspace_zN( const LatticeLAPHSubSpace_t& laph_in, 
				const Seed& rng_seed, const int& N);


    //! Random complex Z(N) sources using dilution
    /*! @ingroup sources
     *
     * Create a random Z(N) using dilution
     */
    template<typename T>
    class SourceConst : public QuarkSourceConstruction<T>
    {
    public:
      //! Full constructor
      SourceConst(const Params& p) : params(p) {}

      //! Construct the source
      T operator()(const multi1d<LatticeColorMatrix>& u) const;

    private:
      //! Hide partial constructor
      SourceConst() {}

    private:
      Params  params;   /*!< source params */
    };

  }  // end namespace DiluteZNEigVecQuarkSourceConstEnv


  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, DiluteZNEigVecQuarkSourceConstEnv::Params& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const DiluteZNEigVecQuarkSourceConstEnv::Params& param);

}  // end namespace Chroma


#endif
