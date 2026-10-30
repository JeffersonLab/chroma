// Standalone Superb test: complex-double construction -> float storage -> double readback.
#include "superbblas.h"
#include <complex>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
int main(int argc,char** argv) {
  using namespace superbblas;
  using D=std::complex<double>;using F=std::complex<float>;
  if(argc!=2) throw std::runtime_error("Usage: roundtrip temporary-output.s3t");
  Coor<1> dim{{5}},zero{{0}};
  std::array<Coor<1>,2> part{{zero,dim}};
  std::vector<D> a{{1./3.,-1./7.},{1e-20,1e20},{0,0},{-2.1,3.7},{1.0000001,-.9999999}},b(5);
  const char* meta="<DBMetaData><elemental_basis>YLM</elemental_basis><storage_precision>32</storage_precision></DBMetaData>";
  Storage_handle h;
  auto ctx=createCpuContext();
  create_storage<1,F>(dim,SlowToFast,argv[1],meta,std::strlen(meta),checksum_type::BlockChecksum,&h);
  append_blocks<1,F>(&part,1,dim,h,SlowToFast);
  const D* in=a.data();
  save<1,1,D,F>(D(1),&part,1,"i",zero,dim,dim,&in,&ctx,"i",zero,h,SlowToFast);
  close_storage<1,F>(h);
  open_storage<1,F>(argv[1],false,&h);
  check_storage<1,F>(h);
  D* out=b.data();
  load<1,1,F,D>(F(1),h,"i",zero,dim,&part,1,"i",zero,dim,&out,&ctx,SlowToFast,Copy);
  close_storage<1,F>(h);
  for(int i=0;i<5;++i) if(b[i]!=D(F(a[i]))) throw std::runtime_error("roundtrip mismatch");
  values_datatype dtype;std::vector<char> metadata;std::vector<IndexType> dims;
  read_storage_header(argv[1],SlowToFast,dtype,metadata,dims);
  if(dtype!=CFLOAT || dims!=std::vector<IndexType>{5}) throw std::runtime_error("wrong header");
  std::cout<<"PASS: double-to-float S3T write, header dtype, checksums, and double readback\n";
}
