/*! \file
 * \brief Key for time-sliced color eigenvectors
 */

#include "util/ferm/key_timeslice_colorvec.h"

namespace Chroma 
{ 
  //----------------------------------------------------------------------------
  //! Diagnostics
  StandardOutputStream& operator<<(StandardOutputStream& os, const KeyTimeSliceColorVec_t& param)
  {
    os << "KeyTimeSliceColorVec_t:";
    os << " t_slice= " << param.t_slice;
    os << " colorvec= " << param.colorvec;
    os << endl;

    return os;
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
