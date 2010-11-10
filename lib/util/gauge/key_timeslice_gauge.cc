/*! \file
 * \brief Key for time-sliced gauge fields
 */

#include "util/gauge/key_timeslice_gauge.h"

namespace Chroma 
{ 
  //----------------------------------------------------------------------------
  // Support for the keys of prop colo
  bool operator<(const KeyTimeSliceGauge_t& a, const KeyTimeSliceGauge_t& b)
  {
    multi1d<int> lgaa(2);
    lgaa[0] = a.t_slice;
    lgaa[1] = a.dir;

    multi1d<int> lgbb(2);
    lgbb[0] = b.t_slice;
    lgbb[1] = b.dir;

    return (lgaa < lgbb);
  }



  //----------------------------------------------------------------------------
  // KeyTimeSliceDist read
  void read(BinaryReader& bin, KeyTimeSliceGauge_t& param)
  {
    read(bin, param.t_slice);
    read(bin, param.dir);
  }

  // KeyTimeSliceDist write
  void write(BinaryWriter& bin, const KeyTimeSliceGauge_t& param)
  {
    write(bin, param.t_slice);
    write(bin, param.dir);
  }

  //! KeyTimeSliceDist reader
  void read(XMLReader& xml, const std::string& path, KeyTimeSliceGauge_t& param)
  {
    XMLReader paramtop(xml, path);
    
    read(paramtop, "t_slice", param.t_slice);
    read(paramtop, "dir", param.dir);
  }

  // KeyTimeSliceDist writer
  void write(XMLWriter& xml, const std::string& path, const KeyTimeSliceGauge_t& param)
  {
    push(xml, path);

    write(xml, "t_slice", param.t_slice);
    write(xml, "dir", param.dir);

    pop(xml);
  }

} // namespace Chroma
