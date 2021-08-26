#include <chromabase.h>
#if defined(ARCH_PARSCALAR) && defined(MPI_VERSION)
#include <mpi.h>
#endif

namespace Chroma {

#if defined(ARCH_PARSCALAR) && defined(MPI_VERSION)
  int localRankCommSplit()
  {
    char hostname[256];

    int np_global=0;
    int np_local=0;
    int rank_global=0;
    int rank_local=0;

    MPI_Comm_size(MPI_COMM_WORLD, &np_global);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_global);

    MPI_Comm nodeComm;
    MPI_Comm_split_type( MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank_global,
                     MPI_INFO_NULL, &nodeComm );

    MPI_Comm_size(nodeComm, &np_local);
    MPI_Comm_rank(nodeComm, &rank_local);

    MPI_Comm_free(&nodeComm);
    return rank_local;
  }

  int localRank() 
  {
    int ret_val; 
    // Try
    char *rank = getenv( "PMI_RANK" );

    // Try MVapich
    if( ! rank ) {
        rank = getenv( "MV2_COMM_WORLD_LOCAL_RANK"  );
    }

    // Try OpenMPI
    if( ! rank ) {
       rank = getenv( "OMPI_COMM_WORLD_LOCAL_RANK" );
    }

    // Try SLURM
    if( ! rank ) {
       rank = getenv( "SLURM_LOCALID" );
    }

    if ( rank )  {
       ret_val = atoi( rank );
    }
    else { 
        ret_val = localRankCommSplit();
    } 
    
    QDPIO::cout << "local rank =  " << ret_val << std::endl;  
    return ret_val;
  }

#elif defined ARCH_SCALAR
  int localRank()
  {
    return 0;
  }
#endif
}; // Namespace
