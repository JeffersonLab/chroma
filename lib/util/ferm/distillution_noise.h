// -*- C++ -*-
/*! \file
 * \brief Support for distillution - random time-slices and quark line noises
 */

#ifndef __distillution_noise_h__
#define __distillution_noise_h__

#include "chromabase.h"

namespace Chroma
{
  //---------------------------------------------------------------------
  //! Keys for lattice noise
  /*! 
   * \ingroup ferm 
   *
   *  Key is used to generate a random number sequence
   */
  struct DistQuarkLines_t
  {
    int          num_vecs;       /*!< Number of eigenvectors */
    int          quark_line;     /*!< Quark line */
    bool         annih;          /*!< Annihilation quark line */
    std::string  mass;           /*!< Mass label */
  };

 
  //---------------------------------------------------------------------
  //! Lattice origin
  /*! \ingroup ferm */
  class DistillutionNoise
  {
  public:
    //! Default constructor
    DistillutionNoise(const string& ensemble, const string& sequence_label, int decay_dir);

    //! Destructor
    ~DistillutionNoise() {} 

    //! Return the ensemble
    virtual string getEnsemble() const {return ensemble;}

    //! Return the sequence
    virtual string getSequence() const {return seqno;}

    //! Return the decay direction
    virtual int getDecayDir() const {return decay_dir;}

    //! Return the actual time origin
    virtual int getOrigin() const {return t_origin;}

    //! Convenience - get shifted time
    /*!< Returns the shifted time taking account of periodicity */
    virtual int getTime(int t_slice) const;

    //! Return the RN-s needed for this line
    /*! Indexing is  (t_slice,vector_num)  */
    virtual multi2d<Complex> getRNG(const DistQuarkLines_t& info) const;

  private:
    string  ensemble;          /*!< Ensemble used for seed of RNG */
    string  seqno;             /*!< Sequence label used for seed of RNG */
    int     decay_dir;         /*!< Decay direction */
    int     t_origin;          /*!< Actual time origin */
  };

} // namespace Chroma

#endif
