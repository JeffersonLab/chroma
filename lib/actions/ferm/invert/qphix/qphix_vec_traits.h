#ifndef QPHIX_VEC_TRAITS_H
#define QPHIX_VEC_TRAITS_H

#include "chroma_config.h"

#ifdef BUILD_QPHIX
namespace Chroma {
  namespace QPhiXVecTraits {

    template<typename T>
    struct VecTraits { 
      static const int Vec=1;
      static const int Soa=1;
      static const bool compress12=false;
    };

    template<typename TOuter,typename TInner>
    struct MixedVecTraits { 
      static const int Vec=1;
      static const int Soa=1;
      static const bool compress12=false;
      static const int VecInner=1;
      static const int SoaInner=1;

    };

#if defined(CHROMA_QPHIX_ARCH_SSE)
#warning QPHIX for SSE
    // AVX Traits:
    template<>
    struct VecTraits<float> { 
      static const int Vec=4;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
    };

    template<>
    struct VecTraits<double> { 
      static const int Vec=2;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 

    };

    // SSE Traits:
    template<>
    struct MixedVecTraits<double,double> { 
      static const int Vec=2;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
      static const int VecInner=2;
      static const int SoaInner=CHROMA_QPHIX_INNER_SOALEN;
    };

    template<>
    struct MixedVecTraits<double,float> { 
      static const int Vec=2;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
      static const int VecInner=4;
      static const int SoaInner=CHROMA_QPHIX_INNER_SOALEN;
    };

   template<>
    struct MixedVecTraits<float,float> { 
      static const int Vec=4;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
      static const int VecInner=4;
      static const int SoaInner=CHROMA_QPHIX_INNER_SOALEN;
    };

#endif

    // Templates
#if defined(CHROMA_QPHIX_ARCH_AVX) || defined(CHROMA_QPHIX_ARCH_AVX2)
#warning QPHIX for AVX and AVX2
    // AVX Traits:
    template<>
    struct VecTraits<float> { 
      static const int Vec=8;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
    };

    template<>
    struct VecTraits<double> { 
      static const int Vec=4;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 

    };

    // AVX Traits:
    template<>
    struct MixedVecTraits<double,double> { 
      static const int Vec=4;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
      static const int VecInner=4;
      static const int SoaInner=CHROMA_QPHIX_INNER_SOALEN;
    };

    template<>
    struct MixedVecTraits<double,float> { 
      static const int Vec=4;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
      static const int VecInner=8;
      static const int SoaInner=CHROMA_QPHIX_INNER_SOALEN;
    };

   template<>
    struct MixedVecTraits<float,float> { 
      static const int Vec=8;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
      static const int VecInner=8;
      static const int SoaInner=CHROMA_QPHIX_INNER_SOALEN;
    };


#endif

#if defined(CHROMA_QPHIX_ARCH_MIC) || defined(CHROMA_QPHIX_ARCH_AVX512)
#warning QPhiX for MIC or AVX512
    // MIC Traits
    template<>
    struct VecTraits<float> { 
      static const int Vec=16;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
    };

    template<>
    struct VecTraits<double> { 
      static const int Vec=8;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
    };

    // MIC Traits
    template<>
    struct MixedVecTraits<double,double> { 
      static const int Vec=8;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
      static const int VecInner=8;
      static const int SoaInner=CHROMA_QPHIX_INNER_SOALEN;

    };
    template<>
    struct MixedVecTraits<double,float> { 
      static const int Vec=8;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
      static const int VecInner=16;
      static const int SoaInner=CHROMA_QPHIX_INNER_SOALEN;
    };
    template<>
    struct MixedVecTraits<double,QPhiX::half> { 
      static const int Vec=8;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
      static const int VecInner=16;
      static const int SoaInner=CHROMA_QPHIX_INNER_SOALEN;
    };
    template<>
    struct MixedVecTraits<float,float> { 
      static const int Vec=16;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
      static const int VecInner=16;
      static const int SoaInner=CHROMA_QPHIX_INNER_SOALEN;

    };
    template<>
    struct MixedVecTraits<float,QPhiX::half> { 
      static const int Vec=16;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
      static const int VecInner=16;
      static const int SoaInner=CHROMA_QPHIX_INNER_SOALEN;
    };

#endif

#if defined CHROMA_QPHIX_ARCH_QPX
#warning QPhiX for QPX
     // QPX Traits
     template<>
     struct VecTraits<double> {
        static const int Vec=4;
        static const int Soa=CHROMA_QPHIX_SOALEN;
        static const bool compress12=CHROMA_QPHIX_COMPRESS12;
     };

    template<>
    struct MixedVecTraits<double,double> { 
      static const int Vec=4;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
      static const int VecInner=4;
      static const int SoaInner=CHROMA_QPHIX_INNER_SOALEN;

    };
    template<>
    struct MixedVecTraits<double,float> { 
      static const int Vec=4;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
      static const int VecInner=4;
      static const int SoaInner=CHROMA_QPHIX_INNER_SOALEN;

    };

    template<>
      struct MixedVecTraits<float,float> { 
      static const int Vec=4;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
      static const int VecInner=4;
      static const int SoaInner=CHROMA_QPHIX_INNER_SOALEN;

    };

#endif


  }; // namespace QPhiXVecTraits
}; // namespace Chroma


#endif // BUILD_QPHIX


#endif


