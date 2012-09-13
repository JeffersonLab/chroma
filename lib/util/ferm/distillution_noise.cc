/*! \file
 * \brief Support for distillution - random time-slices and quark line noises
 */

#include "util/ferm/distillution_noise.h"
#include "util/ferm/crc48.h"
#include "qdp_rannyu.h"

namespace Chroma 
{
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


    //---------------------------------------------------------------------
    //! Convert a hash to a rannyu seed
    multi1d<int> hashToSeed(const CRC48::CRC48_t& hash)
    {
      multi1d<int> seed(4);

      // The CRC48 is six units of 8-bits. Pack this into four 12-bit units
      // First pack two-sets of 8-bits into 16-bit units
      unsigned int h01 = hash.crc[0] | (hash.crc[1] << 8); 
      unsigned int h23 = hash.crc[2] | (hash.crc[3] << 8);
      unsigned int h45 = hash.crc[4] | (hash.crc[5] << 8);

      seed[0]  = h01 & 0xfff;   // use first 12
      h01    >>= 12;            // now only have 4 bits
      h23    <<= 4;             // bump up by 4 bits
      h23     |= h01;           // now have 20 bits
      seed[1]  = h23 & 0xfff;   // use first 12
      h23    >>= 12;            // now only have 8 bits
      h45    <<= 8;             // bump up by 8 bits
      h45     |= h23;           // now have 24 bits
      seed[2]  = h45 & 0xfff;   // used first 12 bits
      seed[3]  = (h45 >> 12);   // use remaining 12 bits

      return seed;
    }


    //---------------------------------------------------------------------
    //! Convert a string to a rannyu seed
    multi1d<int> stringToSeed(const string& output)
    {
      CRC48::CRC48_t crc;

      CRC48::initCRC48(crc);
      CRC48::calcCRC48(crc, output.c_str(), output.length());

      return hashToSeed(crc);
    }


    //---------------------------------------------------------------------
    // Z(N)-rng
    Complex z4rng(RANNYU::RNGState_t& rng)
    {
      RANNYU::random(rng);

      Real twopiN = 0.25 * Chroma::twopi;
      Real theta = twopiN * floor(4*rng.ran);

      return cmplx(cos(theta),sin(theta));
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

    // Throw out the first few RNG's
    for(int i=0; i < 10; ++i)
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
//    return t_slice;
  }


  //---------------------------------------------------------------------
  //! Return the RN-s needed for this line
  multi2d<Complex> DistillutionNoise::getRNG(const DistQuarkLines_t& info) const
  {
    // Use tags to make a unique
    BinaryBufferWriter bin;
    write(bin, ensemble);
    write(bin, seqno);
    write(bin, info.quark_line);
    write(bin, (info.annih) ? 1 : 0);
    write(bin, info.mass);

    // Use the no-side-effect version of the RNG with input state, and output random num.
    RANNYU::RNGState_t rng;
    rng.seed = stringToSeed(bin.str());

    // Throw out the first few RNG's
    for(int i=0; i < 20; ++i)
      RANNYU::random(rng);

    // Loop over the time and vectors to produce the full random numbers
    int Lt = Layout::lattSize()[decay_dir];

    multi2d<Complex> eta(Lt, info.num_vecs);

    for(int t=0; t < Lt; ++t)
    {
      for(int i=0; i < info.num_vecs; ++i)
      {
	eta(t,i) = z4rng(rng);
      }
    }

    return eta;
  }


}  // end namespace Chroma
