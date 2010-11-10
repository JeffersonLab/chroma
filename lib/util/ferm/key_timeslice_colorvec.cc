/*! \file
 * \brief Key for time-sliced color eigenvectors
 */

#include "util/ferm/key_timeslice_colorvec.h"

namespace Chroma 
{ 
  //----------------------------------------------------------------------------
  // Support for the keys of prop colo
  bool operator<(const KeyTimeSliceColorVec_t& a, const KeyTimeSliceColorVec_t& b)
  {
    multi1d<int> lgaa(2);
    lgaa[0] = a.t_slice;
    lgaa[1] = a.colorvec;

    multi1d<int> lgbb(2);
    lgbb[0] = b.t_slice;
    lgbb[1] = b.colorvec;

    return (lgaa < lgbb);
  }



  //----------------------------------------------------------------------------
  // KeyTimeSliceDist read
  void read(BinaryReader& bin, KeyTimeSliceColorVec_t& param)
  {
    read(bin, param.t_slice);
    read(bin, param.colorvec);
  }

  // KeyTimeSliceDist write
  void write(BinaryWriter& bin, const KeyTimeSliceColorVec_t& param)
  {
    write(bin, param.t_slice);
    write(bin, param.colorvec);
  }

  //! KeyTimeSliceDist reader
  void read(XMLReader& xml, const std::string& path, KeyTimeSliceColorVec_t& param)
  {
    XMLReader paramtop(xml, path);
    
    read(paramtop, "t_slice", param.t_slice);
    read(paramtop, "colorvec", param.colorvec);
  }

  // KeyTimeSliceDist writer
  void write(XMLWriter& xml, const std::string& path, const KeyTimeSliceColorVec_t& param)
  {
    push(xml, path);

    write(xml, "t_slice", param.t_slice);
    write(xml, "colorvec", param.colorvec);

    pop(xml);
  }

} // namespace Chroma
