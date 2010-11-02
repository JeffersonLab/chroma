/*! \file
 * \brief Support for distillution - random time-slices and quark line noises
 */

#include "util/ferm/distillution_noise.h"
#include "qdp_rannyu.h"

namespace Chroma 
{
  namespace
  {
    //---------------------------------------------------------------------
    // NOTE: snarfed this from  filedb/filehash/ffdb_file_hash.c
    /*
     * Fowler/Noll/Vo hash
     *
     * The basis of the hash algorithm was taken from an idea sent by email to the
     * IEEE Posix P1003.2 mailing list from Phong Vo (kpv@research.att.com) and
     * Glenn Fowler (gsf@research.att.com).  Landon Curt Noll (chongo@toad.com)
     * later improved on their algorithm.
     *
     * The magic is in the interesting relationship between the special prime
     * 16777619 (2^24 + 403) and 2^32 and 2^8.
     *
     * This hash produces the fewest collisions of any function that we've seen so
     * far, and works well on both numbers and strings.
     *
     * PUBLIC: unsigned int __ham_func5 __P((DB *, const void *, unsigned int));
     */
    unsigned int
    __ham_func5(const void* key, unsigned int len)
    {
      const unsigned char* k = (const unsigned char*)key;
      const unsigned char* e = k + len;

      unsigned int h;
      for (h = 0; k < e; ++k) {
	h *= 16777619;
	h ^= *k;
      }
      return (h);
    }



    //---------------------------------------------------------------------
    //! Convert a hash to a rannyu seed
    multi1d<int> hashToSeed(unsigned int hash)
    {
      multi1d<int> seed(4);

      // Hash is only 32 bits but the seed is 48 bits, so fake a bit
      seed[0] = hash & 4095;
      seed[1] = 17;
      hash >>= 12;
      seed[2] = hash & 4095;
      hash >>= 12;
      seed[3] = hash & 4095;

      return seed;
    }


    //---------------------------------------------------------------------
    //! Convert a string to a rannyu seed
    multi1d<int> stringToSeed(const string& output)
    {
      unsigned int hash = __ham_func5(output.c_str(), output.length());
      return hashToSeed(hash);
    }

  } // end namespace



  //---------------------------------------------------------------------
  // Lattice origin
  DistillutionNoise::DistillutionNoise(const string& ensemble_, const string& sequence_label_, int decay_dir_) :
    ensemble(ensemble_), seqno(sequence_label_), decay_dir(decay_dir_)
  {
    if (decay_dir != Nd-1)
    {
      QDPIO::cerr << __func__ << ": only support decay_dir = Nd-1\n";
      QDP_abort(1);
    }
    
    // Use the ensemble and seqno to make a unique seed. This will be used to form the origin.
    BinaryBufferWriter bin;
    write(bin, ensemble);
    write(bin, seqno);

    // Use the no-side-effect version of the RNG with input state, and output random num.
    RANNYU::RNGState_t rng;
    rng.seed = stringToSeed(bin.str());
    RANNYU::random(rng);

    // Make a new origin
    int Lt = Layout::lattSize()[decay_dir];
    t_origin = int(Lt * rng.ran) % Lt;  // The mod makes sure the RNG cannot be 1.000
  }

  //---------------------------------------------------------------------
  //! Convenience - get shifted time
  /*!< Returns the shifted time taking account of periodicity */
  int DistillutionNoise::getTime(int t_slice) const
  {
    return (t_slice + t_origin) % Layout::lattSize()[decay_dir];
  }


  //---------------------------------------------------------------------
  //! Return the RN-s needed for this line
  multi2d<Complex> DistillutionNoise::getRNG(const DistQuarkLines_t& info) const
  {
    // Use tags to make a unique
    BinaryBufferWriter bin;
    write(bin, ensemble);
    write(bin, seqno);
    writeDesc(bin, info.quark_line);
    write(bin, (info.annih) ? 1 : 0);

    // Use the no-side-effect version of the RNG with input state, and output random num.
    RANNYU::RNGState_t rng;
    rng.seed = stringToSeed(bin.str());

    // Loop over the time and vectors to produce the full random numbers
    int Lt = Layout::lattSize()[decay_dir];

    multi2d<Complex> eta(Lt, info.num_vecs);

    for(int t=0; t < Lt; ++t)
    {
      for(int i=0; i < info.num_vecs; ++i)
      {
	RANNYU::random(rng);
	eta(t,i) = rng.ran;
      }
    }

    return eta;
  }


}  // end namespace Chroma
