#ifndef CHROMA_COLORVEC_MOD_READER_H
#define CHROMA_COLORVEC_MOD_READER_H
#include "chromabase.h"
#include "qdp_map_obj_disk_multiple.h"
#include "qdp_disk_map_slice.h"
#include "util/ferm/key_timeslice_colorvec.h"
#include <memory>
#include <sstream>
#include <stdexcept>

namespace Chroma {
// Historical SUPERB MOD writers omit ProgramInfo and store floats. Native
// WRITE_TIMESLICE_MAP_OBJECT_DISK records its actual precision in ProgramInfo.
inline int colorvecMODPrecision(const std::string& metadata) {
 std::istringstream stream(metadata);XMLReader xml(stream);
 std::string id;read(xml,"/MODMetaData/id",id);
 if(id!="eigenVecsTimeSlice")throw std::runtime_error("MOD does not contain time-slice color vectors");
 int bits=32;
 if(xml.count("/MODMetaData/storage_precision")>0)read(xml,"/MODMetaData/storage_precision",bits);
 else if(xml.count("/MODMetaData/ProgramInfo/code_version/basePrecision")>0)
  read(xml,"/MODMetaData/ProgramInfo/code_version/basePrecision",bits);
 if(bits!=32 && bits!=64)throw std::runtime_error("Unsupported color-vector MOD precision");
 return bits;
}

// Native-precision view of either legacy float or double MOD color vectors.
// Never use a float cache as an intermediate when the input stores doubles.
class ColorvecMODReader {
 using F=QDP::MapObjectDiskMultiple<KeyTimeSliceColorVec_t,TimeSliceIO<LatticeColorVectorF>>;
 using D=QDP::MapObjectDiskMultiple<KeyTimeSliceColorVec_t,TimeSliceIO<LatticeColorVectorD>>;
 std::unique_ptr<F> f_;
 std::unique_ptr<D> d_;
 int debug_=0;
public:
 ColorvecMODReader()=default;
 ColorvecMODReader(const ColorvecMODReader&)=delete;
 ~ColorvecMODReader(){close();}
 void open(const std::vector<std::string>& files) {
  if(files.empty())throw std::runtime_error("Empty color-vector MOD list");
  close();int precision=0;
  for(const auto& file:files) {
   F probe;probe.open({file});std::string metadata;probe.getUserdata(metadata);probe.close();
   const int current=colorvecMODPrecision(metadata);
   if(precision && precision!=current)throw std::runtime_error("Mixed color-vector MOD precisions");
   precision=current;
  }
  if(precision==64) {
   if(sizeof(REAL)<8)throw std::runtime_error("Double MOD requires a double-precision Chroma build");
   d_=std::make_unique<D>();d_->open(files);d_->setDebug(debug_);
  } else {f_=std::make_unique<F>();f_->open(files);f_->setDebug(debug_);}
 }
 void close(){if(f_){f_->close();f_.reset();}if(d_){d_->close();d_.reset();}}
 void setDebug(int level){debug_=level;if(f_)f_->setDebug(level);if(d_)d_->setDebug(level);}
 int getDebug()const{return debug_;}
 int getUserdata(std::string& metadata)const{return d_?d_->getUserdata(metadata):f_->getUserdata(metadata);}
 void keys(std::vector<KeyTimeSliceColorVec_t>& keys)const{if(d_)d_->keys(keys);else f_->keys(keys);}
 int size()const{std::vector<KeyTimeSliceColorVec_t> all;keys(all);return all.size();}
 bool exist(const KeyTimeSliceColorVec_t& key)const{return d_?d_->exist(key):f_->exist(key);}
 int get(const KeyTimeSliceColorVec_t& key,TimeSliceIO<LatticeColorVector>& value)const {
  if(d_){LatticeColorVectorD field=zero;TimeSliceIO<LatticeColorVectorD> slice(field,value.getTimeSlice());int rc=d_->get(key,slice);if(!rc)value.getObject()=field;return rc;}
  LatticeColorVectorF field=zero;TimeSliceIO<LatticeColorVectorF> slice(field,value.getTimeSlice());int rc=f_->get(key,slice);if(!rc)value.getObject()=field;return rc;
 }
};
}
#endif
