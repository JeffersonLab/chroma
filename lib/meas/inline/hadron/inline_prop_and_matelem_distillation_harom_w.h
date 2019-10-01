// -*- C++ -*-
/*! \file
 * \brief Compute the propagator from distillation
 *
 * Propagator calculation in distillation
 */

#ifndef __inline_prop_and_matelem_distillation_harom_w_h__
#define __inline_prop_and_matelem_distillation_harom_w_h__

#ifndef QDP_IS_QDPJIT_NO_NVPTX

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
#include "io/xml_group_reader.h"
#include <vector>
#include <map>
#include <utility>

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlinePropAndMatElemDistillationHaromEnv 
  {
    bool registerAll();

    //! Parameter structure
    /*! \ingroup inlinehadron */ 
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);

      unsigned long     frequency;

      struct Param_t
      {
	struct Contract_t
	{
	  int           num_vecs;       /*!< Number of color vectors to use */
	  int           decay_dir;      /*!< Decay direction */
	  multi1d<int>  t_sources;      /*!< Array of time slice sources for props */
	  int           Nt_forward;     /*!< Time-slices in the forward direction */

	  // not supported, always set to zero !!!
	  int           Nt_backward;    /*!< Time-slices in the backward direction */
	  
	  std::string   mass_label;     /*!< Some kind of mass label */

	  int           num_tries;      /*!< In case of bad things happening in the solution vectors, do retries */

	  bool          do_inversions;
	  bool          check_results;
	  bool          cache_eigs;

	  std::vector<std::string> fifo;
	  int           nodes_per_cn;   // (QDP++) nodes per compute node
	};

	ChromaProp_t    prop;
	Contract_t      contract;
      };

      struct NamedObject_t
      {
	std::string                 gauge_id;        /*!< Gauge field */
	std::vector<std::string>    colorvec_files;  /*!< Eigenvectors in mod format */
	std::string                 prop_op_file;    /*!< File name for propagator matrix elements */
      };

      Param_t           param;
      NamedObject_t     named_obj;
      std::string       xml_file;       /*!< Alternate XML file pattern */
    };


    //! Inline task for the propagator from distillation
    /*! \ingroup inlinehadron */
    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const Params& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the measurement
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    protected:
      //! Do the measurement
      void func(const unsigned long update_no,
		XMLWriter& xml_out); 

    private:
      Params params;
    };



    class Comms
    {
    public:
      ~Comms();

      void add_receive_from( int node , int size );
      void add_send_to( int node , int size );
      void finishSetup();

      void* get_sendbuf( int node );
      void* get_recvbuf( int node );
      int   get_recvbuf_size( int node );
      bool  exists_recvbuf( int node );

      void cleanup();
      void qmp_wait();
      void send_receive();

      void fake_comms();

    private:
      bool setup_started = false;
      bool setup_finished = false;

      std::vector< std::pair< void* , int > > buffers;
      std::vector< QMP_msgmem_t > msgmem;
      std::vector< QMP_msghandle_t > msghandle;
      QMP_msghandle_t mh;

      std::map< int , void* > map_sendbuf;
      std::map< int , std::pair< void* , int > > map_recvbuf;
    }; // Comms




    template<class T> class Buffer;



    template<class T>
    class TSCollect
    {
      template<class T1> struct ContainedType{};
      template<class T1> struct ContainedType<OLattice<T1> > { typedef T1 Type_t; };
      
    public:
      typedef typename ContainedType<T>::Type_t TypeSend_t;


      template < class S , class R >
      void writeout( S& sendbuf , R& recvsize )
      {
	for ( int node = 0 ; node < Layout::numNodes() ; ++node )
	  {
	    if ( Layout::nodeNumber() == node )
	      {
		printf("messages from node %d\n",node);
		  
		for( auto it = sendbuf.begin() ; it != sendbuf.end() ; ++it )
		  {
		    printf("to node %d\n",it->first);

		    auto& tslices = it->second;
		    for( auto mi = tslices.begin() ; mi != tslices.end() ; ++mi )
		      {
			assert( mi->second.size() > 0 );
			printf("ts = %d, n = %d, size = %d bytes\n", mi->first , mi->second.size() , sizeof(mi->second.at(0)) * mi->second.size() );
			// std::cout << "ts " << mi->first
			// 	  << " n = " << mi->second.size() 
			// 	  << " size = " << sizeof(mi->second.at(0)) * mi->second.size() << " bytes\n";
#if 0
			for ( int n = 0 ; n < mi->second.size() ; ++n )
			  {
			    std::cout << mi->second[n].elem().elem(0).real() << " "
				      << mi->second[n].elem().elem(1).real() << " "
				      << mi->second[n].elem().elem(2).real() << " "
				      << mi->second[n].elem().elem(0).imag() << "\n";
			  }
#endif
		      }
		  }

		for( auto it = recvsize.begin() ; it != recvsize.end() ; ++it )
		  {
		    int size = 0; // sizeof(int);  // number of records

		    auto& tslices = it->second;
		    for( auto ts = tslices.begin() ; ts != tslices.end() ; ++ts )
		      {
			size += ts->second;
		      }
		    printf("recv from node %d, %d bytes\n" , it->first , size );
		    //std::cout << "recv from node " << it->first << ", " << size << " bytes\n";
		  }

	      }
	    QMP_barrier();
	  }
      }


      // multi1d< multi3d< TypeSend_t > >& getRecv() { return _recv; }
      // const multi1d< multi3d< TypeSend_t > >& getRecv() const { return _recv; }


      struct xyz 
      {
	int node_from;
	int linear;
	int x;
	int y;
	int z;
      };


  void exec_prepare( int ts_per_node , int t_source , int Nt_forward , int nodes_per_cn )
  {
    _ts_per_node  = ts_per_node;
    _t_source     = t_source;
    _nodes_per_cn = nodes_per_cn;
    _Nt_forward = Nt_forward;

    const int Nt = Layout::lattSize()[3];
    const int t_end = ( t_source + Nt_forward ) % Nt;

    for(int site=0; site < Layout::vol(); ++site)
      {
	multi1d<int> coord = crtesn(site, Layout::lattSize());

	if (t_source < t_end)
	  {
	    if ( coord[3] < t_source  ||   coord[3] >= t_end )
	      continue;
	  }
	else
	  {
	    if ( t_end <= coord[3]     &&   coord[3] < t_source )
	      continue;
	  }
      
	int node	 = Layout::nodeNumber(coord);
	int linear = Layout::linearSiteIndex(coord);

	int tcorr = ( coord[3] - t_source + Nt ) % Nt;
	int to_node = ( tcorr / ts_per_node ) * nodes_per_cn ;

	if (Layout::nodeNumber() == node)
	  {
	    multi1d<int> size_ts(3);
	    size_ts[0] = Layout::subgridLattSize()[0];
	    size_ts[1] = Layout::subgridLattSize()[1];
	    size_ts[2] = Layout::subgridLattSize()[2];
	    int voln_ts = size_ts[0] * size_ts[1] * size_ts[2];

	    multi1d<int> coord_ts(3);
	    coord_ts[0] = coord[0] % Layout::subgridLattSize()[0];
	    coord_ts[1] = coord[1] % Layout::subgridLattSize()[1];
	    coord_ts[2] = coord[2] % Layout::subgridLattSize()[2];

	    int linear_ts = local_site( coord_ts , size_ts );
	  
	    assert( linear_ts >= 0 );
	    assert( linear_ts < voln_ts );

	    _sendbuf[ to_node ][ tcorr ].resize( voln_ts );

	    _sendbuf_linear[ to_node ][ tcorr ].resize( voln_ts );
	    _sendbuf_linear[ to_node ][ tcorr ][ linear_ts ] = linear;
	    //sendbuf[ to_node ][ tcorr ][ linear_ts ] = lattice.elem(linear);
	  }
	if (Layout::nodeNumber() == to_node)
	  {
	    _recvsize[ node ][ tcorr ] += sizeof( TypeSend_t );
	    _receiving_node = true;
	  }
      }


    if (_receiving_node)
      {
	// First, resize the receive vector
	//
	//_recv.resize(ts_per_node);

	// for ( int i = 0 ; i < ts_per_node ; ++i )
	//   _recv[i].resize( Layout::lattSize()[0],Layout::lattSize()[1],Layout::lattSize()[2] );


	// Next, prepare the cache for the receive vector
	//
	assert( Layout::nodeNumber() % _nodes_per_cn == 0 );
	    
	_recv_vector_xyz.resize( _ts_per_node );

	for ( int h = 0 ; h < _ts_per_node ; ++h )
	  {
	    int ts_vol = Layout::lattSize()[0] * Layout::lattSize()[1] * Layout::lattSize()[2];

	    _recv_vector_xyz[ h ].resize( ts_vol );

	    // int ts    = (Layout::nodeNumber() * ts_per_node + h + t_source ) % Nt;
	    // int tcorr = (Layout::nodeNumber() * ts_per_node + h ) % Nt;

	    int ts    = ( Layout::nodeNumber() / _nodes_per_cn * _ts_per_node + h + _t_source ) % Nt;
	    int tcorr = ( Layout::nodeNumber() / _nodes_per_cn * _ts_per_node + h ) % Nt;
	  
	    // if (Layout::primaryNode())
	    //   std::cout << "node " << Layout::nodeNumber() << " receives ts " << ts << "\n";

	    for( int site_ts = 0 ; site_ts < ts_vol ; ++site_ts )
	      {
		multi1d<int> size_ts(3);
		size_ts[0] = Layout::lattSize()[0];
		size_ts[1] = Layout::lattSize()[1];
		size_ts[2] = Layout::lattSize()[2];

		multi1d<int> coord_ts = crtesn( site_ts , size_ts );

		multi1d<int> coord(4);
		coord[0] = coord_ts[0];
		coord[1] = coord_ts[1];
		coord[2] = coord_ts[2];
		coord[3] = ts;

		int node_from = Layout::nodeNumber(coord);
		//int linear    = Layout::linearSiteIndex(coord);

		//std::cout << linear << " from node " << node_from << "\n";

		// if (site_ts == 0)
		//   if (Layout::primaryNode())
		// 	printf("node %d recvs  ts %d  from node %d\n",Layout::nodeNumber(), ts ,  node_from);

		//int ts_data = ts_vol * sizeof(TypeSend_t);

		//Buffer<TypeSend_t> buf( comms.get_recvbuf( node_from ) );
		//TypeSend_t* buf = static_cast<TypeSend_t*>( comms.get_recvbuf( node_from ) );

		multi1d<int> coord_ts_subgrid(3);
		coord_ts_subgrid[0] = coord_ts[0] % Layout::subgridLattSize()[0];
		coord_ts_subgrid[1] = coord_ts[1] % Layout::subgridLattSize()[1];
		coord_ts_subgrid[2] = coord_ts[2] % Layout::subgridLattSize()[2];

		multi1d<int> size_ts_subgrid(3);
		size_ts_subgrid[0] = Layout::subgridLattSize()[0];
		size_ts_subgrid[1] = Layout::subgridLattSize()[1];
		size_ts_subgrid[2] = Layout::subgridLattSize()[2];

		int linear = local_site( coord_ts_subgrid , size_ts_subgrid );

		//_recv(h)( coord_ts[0] , coord_ts[1] , coord_ts[2] ) = buf.get_buf(tcorr)[ linear ];

		_recv_vector_xyz[ h ][ site_ts ].node_from = node_from;
		_recv_vector_xyz[ h ][ site_ts ].linear = linear;
		_recv_vector_xyz[ h ][ site_ts ].x = coord_ts[0];
		_recv_vector_xyz[ h ][ site_ts ].y = coord_ts[1];
		_recv_vector_xyz[ h ][ site_ts ].z = coord_ts[2];

	      }
	  }



      }

    _prep_done = true;
  }



      void exec( const T& lattice , std::vector< TypeSend_t * >& buf_shm )
      {
	//printf("exec1 enter %d\n",Layout::nodeNumber());

	//const int ts_per_node = 4; // ts_per_node
	const int Nt = Layout::lattSize()[3];
	// const int t_end = ( t_source + Nt_forward ) % Nt;

	//std::map< int , std::map< int , std::vector< TypeSend_t > > > sendbuf;  // send[to_node][tcorr][elems]
	//std::map< int , std::map< int , int > > recvsize;                       // recv[from_node][tcorr] size in bytes
	//bool receiving_node = false;

	StopWatch sniss2;
	sniss2.reset();
	sniss2.start();


#if 1
	//
	// First fill sendbuf and recvsize
	//
	for( auto it = _sendbuf_linear.begin() ; it != _sendbuf_linear.end() ; ++it )
	  {
	    auto& tslices = it->second;
	    for( auto mi = tslices.begin() ; mi != tslices.end() ; ++mi )
	      {
		int to_node = it->first;
		int tcorr = mi->first;

		for ( int n = 0 ; n < mi->second.size() ; ++n )
		  _sendbuf[ to_node ][ tcorr ][ n ] = lattice.elem( _sendbuf_linear[ to_node ][ tcorr ][ n ] );

	      }
	  }

	//printf("exec2 enter %d\n",Layout::nodeNumber());

	QMP_barrier();
	
	sniss2.stop();

	QDPIO::cout << "Time to fill sendbuf = " << sniss2.getTimeInSeconds() << std::endl;

#endif

#if 1
	if (!comms_setup) {
	  QDPIO::cout << "Setting up comms for distributing timeslices to distinct compute nodes.\n";
	  prepare_comms(_sendbuf,_recvsize);
	  comms_setup = true;
	}
#endif
  
	do_comms(_sendbuf);


#if 0
	if (1)
	  writeout(_sendbuf,_recvsize);
#endif



#if 1
	// int vvv = Layout::lattSize()[0] * Layout::lattSize()[1] * Layout::lattSize()[2];
	// std::vector< TypeSend_t > tmp;
	// tmp.resize( vvv );
	// for (int i=0;i<vvv;++i)
	//   zero_rep(tmp[i]);


	StopWatch sniss1;
	sniss1.reset();
	sniss1.start();

#if 1
	if ( _receiving_node )
	  {
	    assert( Layout::nodeNumber() % _nodes_per_cn == 0 );
	    
	    std::vector< Buffer<TypeSend_t> > buffers(Layout::numNodes());
	    for ( int i = 0 ; i < Layout::numNodes() ; ++i )
	      {
		if ( comms.exists_recvbuf( i ) )
		  {
		    buffers[ i ].set_buf( comms.get_recvbuf( i ) );
		    buffers[ i ].set_maxts( _Nt_forward );
		  }
	      }

	    for ( int h = 0 ; h < _ts_per_node ; ++h )
	      {
		int ts_vol = Layout::lattSize()[0] * Layout::lattSize()[1] * Layout::lattSize()[2];

		int tcorr = ( Layout::nodeNumber() / _nodes_per_cn * _ts_per_node + h ) % Nt;

		for( int site_ts = 0 ; site_ts < ts_vol ; ++site_ts )
		  {
		    auto& access = _recv_vector_xyz[h][site_ts];

		    //_recv(h)( access.x , access.y , access.z ) = tmp[ access.linear ];

		    // Buffer<TypeSend_t> buf( comms.get_recvbuf( access.node_from ) );
		    // _recv(h)( access.x , access.y , access.z ) = buf.get_buf(tcorr)[ access.linear ];

		    //_recv(h)( access.x , access.y , access.z ) = 

		    buf_shm[ h ][ site_ts ] = buffers[access.node_from].get_buf(tcorr)[ access.linear ];


		  }
	      }
	  }
#else
	if ( _receiving_node )
	  {
	    assert( Layout::nodeNumber() % _nodes_per_cn == 0 );
	    
	    for ( int h = 0 ; h < _ts_per_node ; ++h )
	      {
		int ts_vol = Layout::lattSize()[0] * Layout::lattSize()[1] * Layout::lattSize()[2];

		// int ts    = (Layout::nodeNumber() * ts_per_node + h + t_source ) % Nt;
		// int tcorr = (Layout::nodeNumber() * ts_per_node + h ) % Nt;

		int ts    = ( Layout::nodeNumber() / _nodes_per_cn * _ts_per_node + h + _t_source ) % Nt;
		int tcorr = ( Layout::nodeNumber() / _nodes_per_cn * _ts_per_node + h ) % Nt;
	  
		// if (Layout::primaryNode())
		//   std::cout << "node " << Layout::nodeNumber() << " receives ts " << ts << "\n";

		for( int site_ts = 0 ; site_ts < ts_vol ; ++site_ts )
		  {
		    multi1d<int> size_ts(3);
		    size_ts[0] = Layout::lattSize()[0];
		    size_ts[1] = Layout::lattSize()[1];
		    size_ts[2] = Layout::lattSize()[2];

		    multi1d<int> coord_ts = crtesn( site_ts , size_ts );

		    multi1d<int> coord(4);
		    coord[0] = coord_ts[0];
		    coord[1] = coord_ts[1];
		    coord[2] = coord_ts[2];
		    coord[3] = ts;

		    int node_from = Layout::nodeNumber(coord);
		    //int linear    = Layout::linearSiteIndex(coord);

		    //std::cout << linear << " from node " << node_from << "\n";

		    // if (site_ts == 0)
		    //   if (Layout::primaryNode())
		    // 	printf("node %d recvs  ts %d  from node %d\n",Layout::nodeNumber(), ts ,  node_from);

		    //int ts_data = ts_vol * sizeof(TypeSend_t);

		    Buffer<TypeSend_t> buf( comms.get_recvbuf( node_from ) );
		    //TypeSend_t* buf = static_cast<TypeSend_t*>( comms.get_recvbuf( node_from ) );

		    multi1d<int> coord_ts_subgrid(3);
		    coord_ts_subgrid[0] = coord_ts[0] % Layout::subgridLattSize()[0];
		    coord_ts_subgrid[1] = coord_ts[1] % Layout::subgridLattSize()[1];
		    coord_ts_subgrid[2] = coord_ts[2] % Layout::subgridLattSize()[2];

		    multi1d<int> size_ts_subgrid(3);
		    size_ts_subgrid[0] = Layout::subgridLattSize()[0];
		    size_ts_subgrid[1] = Layout::subgridLattSize()[1];
		    size_ts_subgrid[2] = Layout::subgridLattSize()[2];

		    int linear = local_site( coord_ts_subgrid , size_ts_subgrid );

		    _recv(h)( coord_ts[0] , coord_ts[1] , coord_ts[2] ) = buf.get_buf(tcorr)[ linear ];
		  }
	      }
	  }
#endif


	QMP_barrier();
	
	sniss1.stop();

	QDPIO::cout << "Time constructing receive vector = " << sniss1.getTimeInSeconds() << std::endl;

#endif


      } // exec
  
    private:
  
      int local_site(const multi1d<int>& coord, const multi1d<int>& latt_size)
      {
	int order = 0;

	for(int mmu=latt_size.size()-1; mmu >= 1; --mmu)
	  order = latt_size[mmu-1]*(coord[mmu] + order);

	order += coord[0];

	return order;
      }


      template < class S , class R >
      void prepare_comms( S& sendbuf , R& recvsize )
      {
	for( auto it = sendbuf.begin() ; it != sendbuf.end() ; ++it )
	  {
	    int size = sizeof(int); // total number of records
      
	    auto& tslices = it->second;
	    for( auto ts = tslices.begin() ; ts != tslices.end() ; ++ts )
	      {
		assert( ts->second.size() > 0 );
		size += sizeof(ts->second.at(0)) * ts->second.size() + sizeof(int);  // The timeslice information (0 .. Nt_forward)
	      }
	  
	    comms.add_send_to( it->first , size );
	    //std::cout << "send added : " <<  it->first << " " <<  size << "\n";
	  }
  
	for( auto it = recvsize.begin() ; it != recvsize.end() ; ++it )
	  {
	    int size = sizeof(int); // total number of records
      
	    auto& tslices = it->second;
	    for( auto ts = tslices.begin() ; ts != tslices.end() ; ++ts )
	      {
		assert( ts->second > 0 );
		size += ts->second + sizeof(int);  // The timeslice information (0 .. Nt_forward)
	      }

	    comms.add_receive_from( it->first , size );    // total number of records

	    //std::cout << "recv added : " <<  it->first << " " <<  size << "\n";

	    // int size = it->second + sizeof(int);
	    // std::cout << "recv added : " <<  it->first << " " <<  size << "\n";
	  }

	comms.finishSetup();
      }




      template < class S >
      void do_comms( S& sendbuf )
      {
  QMP_barrier();
  StopWatch sniss1;
  sniss1.reset();
  sniss1.start();

	for( auto it = sendbuf.begin() ; it != sendbuf.end() ; ++it )
	  {
	    char* sb = static_cast<char*>(comms.get_sendbuf( it->first ));
	    int pos = 0;

	    int recnum = it->second.size();
	    memcpy( sb + pos , &recnum , sizeof(int) );
	    pos += sizeof(int);
      
	    auto& tslices = it->second;
	    for( auto ts = tslices.begin() ; ts != tslices.end() ; ++ts )
	      {
		int tsnum = ts->first;
		memcpy( sb + pos , &tsnum , sizeof(int) );
		pos += sizeof(int);
	  
		memcpy( sb + pos , ts->second.data() , sizeof(ts->second.at(0)) * ts->second.size() );
		pos += sizeof(ts->second.at(0)) * ts->second.size();
	      }
	  }
	//QDPIO::cout << "MPI start..\n";
	// //comms.fake_comms();

  sniss1.stop();

  QMP_barrier();

  QDPIO::cout << "Time for comms copy = " << sniss1.getTimeInSeconds() << std::endl;

  StopWatch sniss2;
  sniss2.reset();
  sniss2.start();

  
	comms.send_receive();
	comms.qmp_wait();

  sniss2.stop();

  QMP_barrier();

  QDPIO::cout << "Time for comms MPI = " << sniss2.getTimeInSeconds() << std::endl;

	//QDPIO::cout << "MPI done.\n";
      }

  
    private:

      bool _receiving_node = false;
      bool _prep_done      = false;
      int _ts_per_node = -1;
      int _t_source = -1;
      int _nodes_per_cn = -1;
      int _Nt_forward = -1;

      std::map< int , std::map< int , std::vector< TypeSend_t > > > _sendbuf;  // send[to_node][tcorr][elems]

      std::map< int , std::map< int , std::vector< int > > > _sendbuf_linear;  // send[to_node][tcorr][linear]
      std::map< int , std::map< int , int > > _recvsize;                       // recv[from_node][tcorr] size in bytes

      //multi1d< multi3d< TypeSend_t > > _recv;

      std::vector< std::vector< struct xyz > > _recv_vector_xyz;                // [ts_per_node][linear]

      bool comms_setup = false;
      Comms comms;
    }; // TSCollect





    //
    // This class processes data buffers with the following format
    // int(no of records)  int(timeslice)  data  int(timeslice)  data ...
    //
    template<class T>
    class Buffer
    {
    public:
      Buffer( void* buf ): buf(static_cast<char*>(buf)) {}
      Buffer(): buf(NULL) {}
      ~Buffer() 
      {
	if (maxtsset)
	  delete[] cache;
      }

      void set_buf( void* _buf )
      { 
	buf = static_cast<char*>(_buf);
      }

      void set_maxts( int maxts )
      {
	cache = new T* [maxts];
	for ( int i = 0 ; i < maxts ; ++i )
	  cache[ i ] = NULL;
	maxtsset=true;
      }

      T* get_buf( int ts )
      {
	assert( maxtsset );

	if (cache[ts] != NULL)
	  return cache[ts];

	assert( buf != NULL );

	multi1d<int> size_ts(3);
	size_ts[0] = Layout::subgridLattSize()[0];
	size_ts[1] = Layout::subgridLattSize()[1];
	size_ts[2] = Layout::subgridLattSize()[2];
	int voln_ts = size_ts[0] * size_ts[1] * size_ts[2];

	int rec_bytes = voln_ts * sizeof(T) + sizeof(int);

	char* cur = buf;

	int recs = *(int*)cur;
	cur += sizeof(int);

	// if (Layout::primaryNode())
	//   printf("buffer has %d records, looking for ts = %d\n",recs,ts);

	bool found=false;
	char* iter = cur;
	for ( ; iter < cur + recs * rec_bytes ; iter += rec_bytes )
	  {
	    int iter_ts = *(int*)iter;

	    // if (Layout::primaryNode())
	    //   printf("ts %d\n",iter_ts);
	
	    if (iter_ts == ts) {
	      found=true;
	      break;
	    }
	  }
	if (!found)
	  {
	    // if (Layout::primaryNode())
	    printf("node %d, ts %d not found!\n",Layout::nodeNumber(),ts);
	    
	    printf("avail ts:\n");
	    char* iter = cur;
	    for ( ; iter < cur + recs * rec_bytes ; iter += rec_bytes )
	      {
		int iter_ts = *(int*)iter;
		printf("%d\n",iter_ts);
	      }
	    printf("end\n");

	    QDP_error_exit("Buffer::get_buf(ts) error: ts not found");
	  }
	iter += sizeof(int);

	cache[ ts ] = static_cast<T*>( static_cast<void*>(iter) );
	
	return static_cast<T*>( static_cast<void*>(iter) );
      }
  
    private:
      char* buf;
      T** cache;
      bool maxtsset = false;
    }; // class Buffer

    

  } // namespace 


}

#endif // QDPJIT

#endif // Chroma
