#ifndef QDPJIT_MEMORY_WRAPPER_H
#define QDPJIT_MEMORY_WRAPPER_H

#include "chromabase.h"
#include "qdp_config.h"

#ifdef QDP_IS_QDPJIT
#ifdef QDPJIT_IS_QDPJITNVVM
#define  GetMemoryPtr(a)  QDP_get_global_cache().get_kernel_args( a , false )
#define  GetMemoryPtr2(o1,o2,i1,i2) {std::vector<QDPCache::ArgKey> _i_n_ = {i1,i2}; auto _o_u_t_ = QDP_get_global_cache().get_kernel_args( _i_n_ , false ); o1 = _o_u_t_[0]; o2 = _o_u_t_[1];}
#define  GetMemoryPtrGauge(out,in) {std::vector<QDPCache::ArgKey> _i_n_ = {in[0].getId(),in[1].getId(),in[2].getId(),in[3].getId()}; auto _o_u_t_ = QDP_get_global_cache().get_kernel_args( _i_n_ , false ); out[0]=_o_u_t_[0];out[1]=_o_u_t_[1];out[2]=_o_u_t_[2];out[3]=_o_u_t_[3];}
#define  GetMemoryPtrClover(i1,i2,i3,i4) {std::vector<QDPCache::ArgKey> _i_n_ = {i1,i2,i3,i4}; auto _o_u_t_ = QDP_get_global_cache().get_kernel_args( _i_n_ , false ); clover[0]=_o_u_t_[0]; clover[1]=_o_u_t_[1]; cloverInv[0]=_o_u_t_[2]; cloverInv[1]=_o_u_t_[3];}
#else
// OLD JIT/PTX wrapper
#define  GetMemoryPtr(a)  QDPCache::Instance().getDevicePtr( a )
#endif

#endif /* QDP_IS_QDPJIT */

#endif /* QDPJIT_MEMORY_WRAPPER */ 
