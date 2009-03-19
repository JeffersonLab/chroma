// $Id: group_baryon_operator_w.cc,v 1.27 2009-03-19 17:17:20 mcneile Exp $
/*! \file
 *  \brief Construct group baryon operators
 */


#include "chromabase.h"

#include "meas/hadron/group_baryon_operator_w.h"
#include "meas/hadron/baryon_operator_factory_w.h"
#include "meas/hadron/baryon_operator_aggregate_w.h"

#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/quark_smearing_aggregate.h"

#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"

#include "meas/sources/source_smearing_aggregate.h"
#include "meas/sources/source_smearing_factory.h"
#include "meas/sinks/sink_smearing_aggregate.h"
#include "meas/sinks/sink_smearing_factory.h"

#include "meas/sources/dilutezN_source_const.h"
#include "meas/sources/zN_src.h"
#include "meas/smear/quark_source_sink.h"
#include "meas/smear/displacement.h"

#include "util/ferm/diractodr.h"
#include "util/ft/sftmom.h"

#include <map>

using std::map;

namespace Chroma
{ 

  void displacementSub(const multi1d<LatticeColorMatrix>& u, 
		       LatticeFermion& chi, 
		       int length, int dir)
  {
    LatticeFermion tmp;
    if (length > 0)
      for(int n = 0; n < length; ++n)
	{
	  tmp = shift(chi, FORWARD, dir);
	  chi = u[dir] * tmp;
	}
    else // If length = or < 0.  If length == 0, does nothing.
      for(int n = 0; n > length; --n)
	{
	  tmp = shift(adj(u[dir])*chi, BACKWARD, dir);
	  chi = tmp;
	}
  }
  


	  //! Baryon sequential sources
	  /*! \ingroup hadron */
	  namespace GroupBaryonOperatorEnv
	  {			
			//	=======
	  	//	Readers
			//	=======

			void read( XMLReader& xml, const string& path, 
			           GroupBaryonOperatorEnv::Params::Qprop_t::Solutions_t& input )
	  	{ 
				XMLReader inputtop(xml, path);
	  	  read(inputtop, "Soln_file_names", input.soln_file_names);
	  	}
	  	void read( XMLReader& xml, const string& path, 
			           GroupBaryonOperatorEnv::Params::Qprop_t& input )
	  	{ 
				XMLReader inputtop(xml, path);
	  	  read(inputtop, "Quarks", input.solns);
	  	}
	  	void read( XMLReader& xml, const string& path, 
			           GroupBaryonOperatorEnv::Params::dilution_t& input )
	  	{ 
				XMLReader inputtop(xml, path);
				input.spatial_mask_size.resize(3);
	  	  read(inputtop, "spatial_mask_size", input.spatial_mask_size);
				input.spatial_mask.resize(3);
	  	  read(inputtop, "spatial_mask", input.spatial_mask);
				input.spin_mask.resize(4);
	  	  read(inputtop, "spin_mask", input.spin_mask);
				input.color_mask.resize(3);
	  	  read(inputtop, "color_mask", input.color_mask);
	  	}
	  	void read( XMLReader& xml, const string& path, 
			           GroupBaryonOperatorEnv::Params& input )
	  	{ 
				XMLReader inputtop(xml, path);
	  	  read(inputtop, "Qprops", input.qprop);
	  	  read(inputtop, "Cfg", input.gaugestuff.cfg);
	  	  read(inputtop, "DilutionScheme", input.dilution);
	  	}
	  	// Writer
	  	void write( XMLWriter& xml, const string& path, 
			            const GroupBaryonOperatorEnv::Params& param )
	  	{
	  	  param.writeXML( xml, path );
	  	}

		//	============
		//	Params stuff
		//	============
		
    //! Initialize
    Params::Params()
    {}

    //! Read parameters
    Params::Params( XMLReader& xml, const std::string& path )
    {
      XMLReader paramtop( xml, path );
			// ============================
			// number of quarks fixed to 3
			// j_decay dir is 3 (time)
			// space-time is 4-dimensions
			// ============================
			int NumberofQuarks = 3;
			j_decay = 3;
			nrow.resize( 4 );
      int version;
      read( paramtop, "Param/Name", name );
      read( paramtop, "Param/version", version );
			QDPIO::cout<< name << " Version " << version <<endl;
      switch ( version )
      {
        case 2:
          break;
        default:
          QDPIO::cerr << name << ": parameter version " << version << " unsupported." << endl;
          QDP_abort( 1 );
      }
			read( paramtop, "Param/LattSize", nrow );
      read( paramtop, "Param/Baryon_type", param.baryon_operator );

      read( paramtop, "Param/SourceSmearing", source_smearing );
      read( paramtop, "Param/SinkSmearing", sink_smearing );

      gaugestuff.link_smearing = readXMLGroup( paramtop, "Param/LinkSmearing", "LinkSmearingType" );
			read( paramtop, "Param/InputFileName", InputFileName );
			QDPIO::cout<< "Main input file is " << InputFileName <<endl;
			read( paramtop, "Cfg", gaugestuff.cfg );

			dilution.resize( NumberofQuarks );
			for(int i=0; i < NumberofQuarks; ++i) dilution[ i ].N = 4; // Z(4) noise
			for(int i=0; i < NumberofQuarks; ++i) dilution[ i ].j_decay = 3; // time-direction

		} // end Params::Params


    // Writer
    void Params::writeXML( XMLWriter& xml, const string& path ) const
    {
      push( xml, path );
      int version = 1;
      write( xml, "version", version );
      write( xml, "nrow", nrow );
      //write( xml, "BaryonOperatorType", GroupBaryonOperatorEnv::name );
			//xml << source_smearing.xml;
      //xml << sink_smearing.xml;
      //xml << link_smearing.xml;
      write( xml, "InputFileName", InputFileName );
      write( xml, "mom2_max", mom2_max );
      write( xml, "j_decay",  j_decay );
      pop( xml );
    } // end void Params::writeXML


		// Used for hybrid list indices, NH and noise indices 
		//  	for indices:  ( ord, i,   j,   k,   which )			
		//  	for NHlimit:  ( ord, NH0, NH1, NH2, which )			
		//  	for NoiseInd: ( ord, 0,   1,   2,   which )			
		int DilSwap( int ord, int i, int j, int k, int which )   
		{  
		  multi1d<int> ijk( 3 );															   
		  switch ( ord )																			   
		  { 																									   
		    case 0:  // [012] ijk 														   
		      ijk[ 0 ] = i; 																	   
		      ijk[ 1 ] = j; 																	   
		      ijk[ 2 ] = k; 																	   
		      break;																					   
		    case 1:  // [210] kji 														   
		      ijk[ 0 ] = k; 																	   
		      ijk[ 1 ] = j; 																	   
		      ijk[ 2 ] = i; 																	   
		      break;																					   
		    case 2:  // [021] ikj 														   
		      ijk[ 0 ] = i; 																	   
		      ijk[ 1 ] = k; 																	   
		      ijk[ 2 ] = j; 																	   
		      break;																					   
		    case 3:  // [102] jik 														   
		      ijk[ 0 ] = j; 																	   
		      ijk[ 1 ] = i; 																	   
		      ijk[ 2 ] = k; 																	   
		      break;																					   
		    case 4:  // [120] jki 														   
		      ijk[ 0 ] = j; 																	   
		      ijk[ 1 ] = k; 																	   
		      ijk[ 2 ] = i;
		      break;																					   
		    case 5:  // [201] kij 														   
		      ijk[ 0 ] = k; 																	   
		      ijk[ 1 ] = i; 																	   
		      ijk[ 2 ] = j; 																	   
		      break;
				default:
					cerr<<"say what? "<<ord<<" "<<i<<" "<<j<<" "<<k<<" "<<which<<endl;
					exit(1);
		  } 																									   
		  return ( ijk[ which ] );														   
		}
		//
		// index things back to [012] ordering
		//
		int DilSwapInv( int ord, int i, int j, int k, int which )   
		{  
		  multi1d<int> ijk( 3 );															   
		  switch ( ord )																			   
		  { 																									   
		    case 0:  // [012] ijk 														   
		      ijk[ 0 ] = i; 																	   
		      ijk[ 1 ] = j; 																	   
		      ijk[ 2 ] = k; 																	   
		      break;																					   
		    case 1:  // [210] kji 														   
		      ijk[ 0 ] = k; 																	   
		      ijk[ 1 ] = j; 																	   
		      ijk[ 2 ] = i; 																	   
		      break;																					   
		    case 2:  // [021] ikj 														   
		      ijk[ 0 ] = i; 																	   
		      ijk[ 1 ] = k; 																	   
		      ijk[ 2 ] = j; 																	   
		      break;																					   
		    case 3:  // [102] jik 														   
		      ijk[ 0 ] = j; 																	   
		      ijk[ 1 ] = i; 																	   
		      ijk[ 2 ] = k; 																	   
		      break;																					   
		    case 4:  // [120] jki 														   
		      ijk[ 0 ] = k; 																	   
		      ijk[ 1 ] = i; 																	   
		      ijk[ 2 ] = j; 
		      break;																					   
		    case 5:  // [201] kij 														   
		      ijk[ 0 ] = j; 																	   
		      ijk[ 1 ] = k; 																	   
		      ijk[ 2 ] = i; 																	   
		      break;
				default:
					cerr<<"Say what? "<<ord<<" "<<i<<" "<<j<<" "<<k<<" "<<which<<endl;
					exit(1);
		  }
		  return ( ijk[ which ] );														   
		}
		
		//	====================
		//	GroupBaryonQQQ stuff
		//	====================
		
    //! Constructor
    GroupBaryonQQQ::GroupBaryonQQQ()
    {
      // The spin basis matrix to goto Dirac
      //spin_rotate_mat = adj( DiracToDRMat() );
		}
    void GroupBaryonQQQ::init( const Params& p_ )
    {
			myparams = p_;
			quark.resize(3);
      // The spin basis matrix to goto Dirac
      spin_rotate_mat = adj( DiracToDRMat() );
    } // GroupBaryonQQQ::init

    //! Full constructor
    GroupBaryonQQQ::GroupBaryonQQQ( const Params& p ) :
        myparams( p )
		{
      // The spin basis matrix to goto Dirac
      spin_rotate_mat = adj( DiracToDRMat() );
      // Factory constructions
      try
      {
        // Create the source quark smearing object
        {
          std::istringstream xml_s( myparams.source_smearing.source.xml );
          XMLReader smeartop( xml_s );
          const string smear_path = "/SourceSmearing";
          QDPIO::cout << "Source quark smearing type = " 
					            << myparams.source_smearing.source.id << endl;
		      Handle< QuarkSourceSink<LatticeFermion> > sourceSmearing(
							TheFermSourceSmearingFactory::Instance().createObject(
				               myparams.source_smearing.source.id, 
											 smeartop,
											 myparams.source_smearing.source.path, 
											 myparams.gaugestuff.u )
							);
        }
        // Create the sink quark smearing object
        {
          std::istringstream xml_s( myparams.sink_smearing.sink.xml );
          XMLReader smeartop( xml_s );
          const string smear_path = "/SinkSmearing";
          QDPIO::cout << "Sink quark smearing type = " 
					            << myparams.sink_smearing.sink.id << endl;
					Handle< QuarkSourceSink<LatticeFermion> > sinkSmearing(
							TheFermSinkSmearingFactory::Instance().createObject(
											myparams.sink_smearing.sink.id, 
											smeartop,
											myparams.sink_smearing.sink.path, 
											myparams.gaugestuff.u )
							);
        }
      }
      catch ( const std::string & e )
      {
        QDPIO::cerr << name << ": Caught Exception smearing: " << e << endl;
        QDP_abort( 1 );
      }

    } // GroupBaryonQQQ::GroupBaryonQQQ

    //! Reader
    void ReadTextInput( Params& params,
                        multi1d<GroupBaryonOp>& AB, multi1d<GroupBaryonOp>& CB,
                        multi1d<GroupBaryonQQQ>& AQQQ, multi1d<GroupBaryonQQQ>& CQQQ )
    {
      QDPIO::cout << "Reading input from TEXT file : " << params.InputFileName <<endl;
      TextFileReader reader( params.InputFileName );
// for testing
#ifdef REDUCETOTIMEDILUTION
			params.NdilReduce = 12;
#else
			params.NdilReduce = 1;
#endif
			// Lattice sizes
      // ------------------------------------------------------------------------------------
      reader >> params.nrow[ 0 ] >> params.nrow[ 1 ] >> params.nrow[ 2 ] >> params.nrow[ 3 ];
			// ------------------------------------------------------------------------------------

			params.Nmomenta = 1; // need to change this in the future when we do multiple momenta
			params.mom2_max = 0; // this obviously needs to be changed and read in ...
			params.j_decay  = 3; // this is OK

      // Orderings are more or less hard wired for now
			//  ... to be continued later when I have the time
			
			//read( xml, "NumSourcePerm", params.NsrcOrderings );
			// ------------------------------------------------------------------------------------
      reader >> params.NsrcOrderings;
			// ------------------------------------------------------------------------------------
			// 012    factor in front is  2  2
      // 210    factor in front is -2  2
      // 021    factor in front is  1 -1
      // 102    factor in front is  1 -1
      // 120    factor in front is -1 -1
      // 201    factor in front is -1 -1

      //read( xml, "SrcQuarkIndices", params.SrcOrderings );
      params.SrcOrderings.resize( params.NsrcOrderings );
      for(int i=0; i < params.NsrcOrderings; ++i)
      {
        params.SrcOrderings[ i ].resize( 3 );
				// ------------------------------------------------------------------------------------
        reader >> params.SrcOrderings[ i ][ 0 ] >> params.SrcOrderings[ i ][ 1 ] >> params.SrcOrderings[ i ][ 2 ];
				// ------------------------------------------------------------------------------------
      }
      //
			// This again, is hard-wired to 012 case for now ...
			//
			//read( xml, "NumSinkPerm", params.NsnkOrderings );
			// ------------------------------------------------------------------------------------
      reader >> params.NsnkOrderings;
			// ------------------------------------------------------------------------------------
      params.SnkOrderings.resize( params.NsnkOrderings );
      //read( xml, "SnkQuarkIndices", params.SnkOrderings );
      for(int i=0; i < params.NsnkOrderings; ++i)
      {
        params.SnkOrderings[ i ].resize( 3 );
				// ------------------------------------------------------------------------------------
        reader >> params.SnkOrderings[ i ][ 0 ] >> params.SnkOrderings[ i ][ 1 ] >> params.SnkOrderings[ i ][ 2 ];
				// ------------------------------------------------------------------------------------
      }
      // Hybrid list sizes
      params.NH.resize( params.NsrcOrderings ); 
      for(int i=0; i < params.NH.size(); ++i) params.NH[ i ].resize( 3 ); // L,M,R quarks
      	// Set the hybrid list sizes
      	//   { something like this would be better                  }
      	//   { Ndil = params.named_obj.prop.op[n].soln_files.size() }			
      //XMLReader& xml;
      //read( xml, "Ndilutions", params.NH[0] );
			// This is the default 012 ordering
			// ------------------------------------------------------------------------------------
      reader >> params.NH[ 0 ][ 0 ] >> params.NH[ 0 ][ 1 ] >> params.NH[ 0 ][ 2 ];
			// ------------------------------------------------------------------------------------
      for(int i=1; i < params.NsrcOrderings; ++i)
			{
      	for(int j=0; j < 3; ++j)
					params.NH[ i ][ j ] = params.NH[ 0 ][ (params.SrcOrderings[ i ][ j ]) ];
			}
      //read( xml, "NumOperators", params.Noperators );
      // ------------------------------------------------------------------------------------
      reader >> params.Noperators;
			// ------------------------------------------------------------------------------------
      // The Baryon operator names ; G1g_L3_TDT_25 ...
      //read( xml, "OperatorNames", params.Names );
      params.Names.resize( params.Noperators );
			int nameindex;
      for(int i=0; i < params.Names.size(); ++i)
      {
				// ------------------------------------------------------------------------------------
        reader >> nameindex >> params.Names[ i ];
				// ------------------------------------------------------------------------------------
      }
      //read( xml, "NumDistinctQQQ", params.NQQQs );
			// ------------------------------------------------------------------------------------
      reader >> params.NQQQs;
			// ------------------------------------------------------------------------------------
#ifdef MAKE_SINK_OPERATORS
      AQQQ.resize( params.NQQQs );
      for(int i=0; i < params.NQQQs; ++i)
      {
				AQQQ[ i ].init( params );
			}
#endif
#ifdef MAKE_SOURCE_OPERATORS
      CQQQ.resize( params.NQQQs );
      for(int i=0; i < params.NQQQs; ++i)
      {
				CQQQ[ i ].init( params );
			}
#endif
			//
			// GroupBaryonOperator structures
			//
#ifdef MAKE_SINK_OPERATORS
      AB.resize( params.Noperators );
#endif
#ifdef MAKE_SOURCE_OPERATORS
      CB.resize( params.Noperators );
#endif
			//
			// Using the BaryonOperator_t structure
			//
      for(int n=0; n < params.Noperators; ++n)
      {
#ifdef MAKE_SINK_OPERATORS
        // This is 1 now because all the 
				// contractions are summed over
				// AB[ n ].baryonoperator.orderings.resize( params.NsnkOrderings );
        // for(int ord=0; ord < params.NsnkOrderings; ++ord)
				AB[ n ].baryonoperator.orderings.resize( 1 );
        for(int ord=0; ord < 1; ++ord)
				{
        	AB[ n ].baryonoperator.orderings[ ord ].op.resize( params.NH[ 0 ][ 0 ]/params.NdilReduce, params.NH[ 0 ][ 1 ]/params.NdilReduce, params.NH[ 0 ][ 2 ]/params.NdilReduce );
        	for(int i=0; i < (params.NH[ 0 ][ 0 ]/params.NdilReduce); ++i)
        	  for(int j=0; j < (params.NH[ 0 ][ 1 ]/params.NdilReduce); ++j)
        	    for(int k=0; k < (params.NH[ 0 ][ 2 ]/params.NdilReduce); ++k)
        	    {
        	      AB[ n ].baryonoperator.orderings[ ord ].op( i, j, k ).ind.resize(1);
        	      AB[ n ].baryonoperator.orderings[ ord ].op( i, j, k ).ind[ 0 ].elem.resize( params.Nmomenta, params.nrow[3] );
        	      for(int t=0; t<params.nrow[3]; ++t)
									AB[ n ].baryonoperator.orderings[ ord ].op( i, j, k ).ind[ 0 ].elem( 0, t ) = cmplx( Real(0), Real(0) );
							}
        }
				AB[ n ].baryonoperator.version = 1;
				AB[ n ].baryonoperator.mom2_max = params.mom2_max;
				AB[ n ].baryonoperator.j_decay = params.j_decay;
				//AB[ n ].baryonoperator.quarkOrderings = params.SnkOrderings;
				AB[ n ].baryonoperator.quarkOrderings.resize(1);
				AB[ n ].baryonoperator.quarkOrderings[0].resize(3);
				AB[ n ].baryonoperator.quarkOrderings[0][0] = 0;
				AB[ n ].baryonoperator.quarkOrderings[0][1] = 1;
				AB[ n ].baryonoperator.quarkOrderings[0][2] = 2;
#endif
#ifdef MAKE_SOURCE_OPERATORS
				//CB[ n ].baryonoperator.orderings.resize( params.NsrcOrderings );
        //for(int ord=0; ord < params.NsrcOrderings; ++ord)
				CB[ n ].baryonoperator.orderings.resize( 1 );
        for(int ord=0; ord < 1; ++ord)
				{
        CB[ n ].baryonoperator.orderings[ ord ].op.resize( params.NH[ 0 ][ 0 ]/params.NdilReduce, params.NH[ 0 ][ 1 ]/params.NdilReduce, params.NH[ 0 ][ 2 ]/params.NdilReduce );
        for(int i=0; i < (params.NH[ 0 ][ 0 ]/params.NdilReduce); ++i)
          for(int j=0; j < (params.NH[ 0 ][ 1 ]/params.NdilReduce); ++j)
            for(int k=0; k < (params.NH[ 0 ][ 2 ]/params.NdilReduce); ++k)
            {
              CB[ n ].baryonoperator.orderings[ ord ].op( i, j, k ).ind.resize(1);
              CB[ n ].baryonoperator.orderings[ ord ].op( i, j, k ).ind[ 0 ].elem.resize( params.Nmomenta, params.nrow[3] );
              for(int t=0; t<params.nrow[3]; ++t) CB[ n ].baryonoperator.orderings[ ord ].op( i, j, k ).ind[ 0 ].elem( 0, t ) = cmplx( Real(0), Real(0) );
            }
				}
				CB[ n ].baryonoperator.version = 1;
				CB[ n ].baryonoperator.mom2_max = params.mom2_max;
				CB[ n ].baryonoperator.j_decay = params.j_decay;
				//CB[ n ].baryonoperator.quarkOrderings = params.SrcOrderings;
				CB[ n ].baryonoperator.quarkOrderings.resize(1);
				CB[ n ].baryonoperator.quarkOrderings[0].resize(3);
				CB[ n ].baryonoperator.quarkOrderings[0][0] = 0;
				CB[ n ].baryonoperator.quarkOrderings[0][1] = 1;
				CB[ n ].baryonoperator.quarkOrderings[0][2] = 2;
#endif
      } // end n (operators)
      for(int i=0; i < params.NQQQs; ++i)
      {
				int dL, hash;
				int q0spin,q1spin,q2spin,q0disp,q1disp,q2disp,nbops;
				reader >> hash >> q0spin >> q1spin >> q2spin >> q0disp >> q1disp >> q2disp >> dL >> nbops;
#ifdef MAKE_SINK_OPERATORS
				AQQQ[ i ].quark[ 0 ].spin = q0spin -1;
				AQQQ[ i ].quark[ 1 ].spin = q1spin -1;
				AQQQ[ i ].quark[ 2 ].spin = q2spin -1;
				AQQQ[ i ].quark[ 0 ].displacement = q0disp;
				AQQQ[ i ].quark[ 1 ].displacement = q1disp;
				AQQQ[ i ].quark[ 2 ].displacement = q2disp;
	      AQQQ[ i ].NBaryonOps = nbops;
        for(int j=0; j < 3; ++j)
        {
          if ( AQQQ[ i ].quark[ j ].displacement >= 0 )
					{
						AQQQ[ i ].quark[ j ].disp_len = dL;
					} 
					else 
					{
						AQQQ[ i ].quark[ j ].disp_len = -dL;
					}
				}
#endif
#ifdef MAKE_SOURCE_OPERATORS
				CQQQ[ i ].quark[ 0 ].spin = q0spin -1;
				CQQQ[ i ].quark[ 1 ].spin = q1spin -1;
				CQQQ[ i ].quark[ 2 ].spin = q2spin -1;
				CQQQ[ i ].quark[ 0 ].displacement = q0disp;
				CQQQ[ i ].quark[ 1 ].displacement = q1disp;
				CQQQ[ i ].quark[ 2 ].displacement = q2disp;
	      CQQQ[ i ].NBaryonOps = nbops;
        for(int j=0; j < 3; ++j)
        {
          if ( CQQQ[ i ].quark[ j ].displacement >= 0 )
					{
						CQQQ[ i ].quark[ j ].disp_len = dL;
					} 
					else 
					{
						CQQQ[ i ].quark[ j ].disp_len = -dL;
					}
				}
#endif
#ifdef MAKE_SINK_OPERATORS
        for(int j=0; j < 3; ++j)
        {
          if ( AQQQ[ i ].quark[ j ].displacement == 0 )
          {
            AQQQ[ i ].quark[ j ].disp_dir = 0;
            AQQQ[ i ].quark[ j ].disp_len = 0;
            AQQQ[ i ].quark[ j ].disp_ind = 0;
          }
          else if ( AQQQ[ i ].quark[ j ].displacement > 0 )
          {
            AQQQ[ i ].quark[ j ].disp_dir = AQQQ[ i ].quark[ j ].displacement - 1;
            AQQQ[ i ].quark[ j ].disp_ind = AQQQ[ i ].quark[ j ].displacement;
            AQQQ[ i ].quark[ j ].disp_len = dL;
          }
          else
          {
            AQQQ[ i ].quark[ j ].disp_dir = -AQQQ[ i ].quark[ j ].displacement - 1;
            AQQQ[ i ].quark[ j ].disp_ind = -AQQQ[ i ].quark[ j ].displacement + 3;
            AQQQ[ i ].quark[ j ].disp_len = -dL;
          }
        } // end j
        AQQQ[ i ].coef.resize( AQQQ[ i ].NBaryonOps );
        AQQQ[ i ].whichBaryonOps.resize( AQQQ[ i ].NBaryonOps );
        AQQQ[ i ].baryon.resize( AQQQ[ i ].NBaryonOps );
        //             name                                    coefficient
        // AQQQ[ <mu,nu,tau> ].whichBaryonOps[ n ]    ,   AQQQ[ <mu,nu,tau> ].coef[ n ]
        //
        for(int n=0; n < AQQQ[ i ].NBaryonOps; ++n)
        {
          Real re, im;
					// ------------------------------------------------------------------------------------
          reader >> AQQQ[ i ].whichBaryonOps[ n ];
          reader >> re >> im;
					// ------------------------------------------------------------------------------------
          AQQQ[ i ].coef[ n ] = cmplx( re, im );
					//
					// Storage will be in the GroupBaryonOp class
					//
					AQQQ[ i ].baryon[ n ] = &AB[ (AQQQ[ i ].whichBaryonOps[ n ]) ];
        }
#endif
#ifdef MAKE_SOURCE_OPERATORS
        for(int j=0; j < 3; ++j)
        {
          if ( CQQQ[ i ].quark[ j ].displacement == 0 )
          {
            CQQQ[ i ].quark[ j ].disp_dir = 0;
            CQQQ[ i ].quark[ j ].disp_len = 0;
            CQQQ[ i ].quark[ j ].disp_ind = 0;
          }
          else if ( CQQQ[ i ].quark[ j ].displacement > 0 )
          {
            CQQQ[ i ].quark[ j ].disp_dir = CQQQ[ i ].quark[ j ].displacement - 1;
            CQQQ[ i ].quark[ j ].disp_ind = CQQQ[ i ].quark[ j ].displacement;
            CQQQ[ i ].quark[ j ].disp_len = dL;
          }
          else
          {
            CQQQ[ i ].quark[ j ].disp_dir = -CQQQ[ i ].quark[ j ].displacement - 1;
            CQQQ[ i ].quark[ j ].disp_ind = -CQQQ[ i ].quark[ j ].displacement + 3;
            CQQQ[ i ].quark[ j ].disp_len = -dL;
          }
        } // end j
        CQQQ[ i ].coef.resize( CQQQ[ i ].NBaryonOps );
        CQQQ[ i ].whichBaryonOps.resize( CQQQ[ i ].NBaryonOps );
        CQQQ[ i ].baryon.resize( CQQQ[ i ].NBaryonOps );
        //             name                                    coefficient
        // CQQQ[ <mu,nu,tau> ].whichBaryonOps[ n ]    ,   CQQQ[ <mu,nu,tau> ].coef[ n ]
        //
        for(int n=0; n < CQQQ[ i ].NBaryonOps; ++n)
        {
          Real re, im;
					// ------------------------------------------------------------------------------------
          reader >> CQQQ[ i ].whichBaryonOps[ n ];
          reader >> re >> im;
					// ------------------------------------------------------------------------------------
          CQQQ[ i ].coef[ n ] = cmplx( re, im );
					//
					// Storage will be in the GroupBaryonOp class
					//
					CQQQ[ i ].baryon[ n ] = &CB[ (CQQQ[ i ].whichBaryonOps[ n ]) ];
        }
#endif
      } // i loop
      //
      // Now for the creation operators
      //
      // rotate by gamma4's to get the barred coefficients
      //       index     0  1   2   3
      //      gamma4 = ( 1, 1, -1, -1 )
      // and take the complex conjugate
#ifdef MAKE_SOURCE_OPERATORS
      for(int i=0; i < params.NQQQs; ++i)
      {
        int spinCount, flipsign;
        spinCount = 0;
        flipsign = 0;
        for(int j=0; j < 3; ++j)
        {
          if ( CQQQ[ i ].quark[ j ].spin > 1 ) spinCount++;
        }
        if ( spinCount % 2 ) flipsign = 1;
        //
				//             name                           coefficient
        // CQQQ[ <mu,nu,tau> ].whichBaryonOps[ n ] , CQQQ[ i ].coef[ n ]
        //
        for(int n=0; n < CQQQ[ i ].NBaryonOps; ++n)
        {
          if ( flipsign )
          {
            // cc and -sign
            // CQQQ[ i ].coef[ n ] = cmplx( -re, im );
            // flip sign (cc at the correlator level)
						CQQQ[ i ].coef[ n ] *= cmplx( Real(-1), Real(0) );
          }
					//
					// Storage will be in the GroupBaryonOp class
					//
					CQQQ[ i ].baryon[ n ] = &CB[ CQQQ[ i ].whichBaryonOps[ n ] ];
        } // n : CQQQ[ i ].baryon[ n ] loop
      } // i : CQQQ[] loop
#endif			
			// solution file names
			params.qprop.solns.resize( 3 );
      for(int n=0; n < 3; ++n)
			{
				params.qprop.solns[ n ].soln_file_names.resize( params.NH[ 0 ][ n ] );
      	for(int i=0; i < params.NH[ 0 ][ n ]; ++i)
				{
					// ------------------------------------------------------------------------------------
        	reader >> params.qprop.solns[ n ].soln_file_names[ i ];
					// ------------------------------------------------------------------------------------
				}
			}
      reader.close();

      QDPIO::cout << "Reading input from text file DONE " << endl;
    } // void ReadTextInput
		
		
    //! Construct the source
    LatticeFermion
    GroupBaryonQQQ::operator() ( Params::dilution_t diln ) const
    {
      //QDPIO::cout << "Diluted random complex ZN source" << endl;
      if ( diln.spatial_mask_size.size() != Nd - 1 )
      {
        QDPIO::cerr << name << ": spatial mask size incorrect 1" << endl;
        QDP_abort( 1 );
      }
      if ( diln.spatial_mask.size() == 0 )
      {
        QDPIO::cerr << name << ": spatial mask incorrect 2" << endl;
        QDP_abort( 1 );
      }
      multi1d<int> lookup_dir( Nd - 1 );
      int mu = 0;
      for(int j=0; j < diln.spatial_mask_size.size(); ++j, ++mu)
      {
        if ( j == diln.j_decay ) ++mu;  // bump up to next dir
        lookup_dir[ j ] = mu;
      }
      for(int j=0; j < diln.spatial_mask.size(); ++j)
      {
        if ( diln.spatial_mask[ j ].size() != Nd - 1 )
        {
          QDPIO::cerr << name << ": spatial mask incorrect 3" << endl;
          QDP_abort( 1 );
        }
      }
      for(int c=0; c < diln.color_mask.size(); ++c)
      {
        if ( diln.color_mask[ c ] < 0 || diln.color_mask[ c ] >= Nc )
        {
          QDPIO::cerr << name << ": color mask incorrect 6" << endl;
          QDP_abort( 1 );
        }
      }
      for(int s=0; s < diln.spin_mask.size(); ++s)
      {
        if ( diln.spin_mask[ s ] < 0 || diln.spin_mask[ s ] >= Ns )
        {
          QDPIO::cerr << name << ": spin mask incorrect 7" << endl;
          QDP_abort( 1 );
        }
      }
      // Save current seed
      Seed ran_seed;
      QDP::RNG::savern( ran_seed );
      // Set the seed to desired value
      QDP::RNG::setrn( diln.ran_seed );
      // Create the noisy quark source on the entire lattice
      LatticeFermion quark_noise;
      zN_src( quark_noise, diln.N );
      // This is the filtered noise source to return
      LatticeFermion quark_source = zero;
      // Filter over the color and spin indices first
      for(int s=0; s < diln.spin_mask.size(); ++s)
      {
        int spin_source = diln.spin_mask[ s ];
        LatticeColorVector colvec = peekSpin( quark_noise, spin_source );
        LatticeColorVector dest = zero;
        for(int c=0; c < diln.color_mask.size(); ++c)
        {
          int color_source = diln.color_mask[ c ];
          LatticeComplex comp = peekColor( colvec, color_source );
          pokeColor( dest, comp, color_source );
        }
        pokeSpin( quark_source, dest, spin_source );
      }
      quark_noise = quark_source;  // reset
      // Filter over the spatial sites
      LatticeBoolean mask = false;  // this is the starting mask
      for(int n=0; n < diln.spatial_mask.size(); ++n)
      {
        LatticeBoolean btmp = true;
        for(int j=0; j < diln.spatial_mask[ n ].size(); ++j)
          btmp &= ( Layout::latticeCoordinate( lookup_dir[ j ] ) % diln.spatial_mask_size[ j ] ) == diln.spatial_mask[ n ][ j ];
        mask |= btmp;
      }
      // Filter over the time slices
      mask &= Layout::latticeCoordinate( diln.j_decay ) == diln.t_source;
      // Zap the unused sites
      quark_source = where( mask, quark_noise, Fermion( zero ) );
      // Reset the seed
      QDP::RNG::setrn( ran_seed );
      return quark_source;
    };
		
		
    //! Construct array of maps of displacements
    void
    GroupBaryonQQQ::displaceQuarks( multi1d< map<int, LatticeFermion> >& disp_quarks,
                                    const multi1d<LatticeFermion>& q,
																		int* qindices                                   
                                  ) const
    {
      START_CODE();
      //QDPIO::cout << __PRETTY_FUNCTION__ << ": entering" << endl;
      disp_quarks.resize( 3 );
      for(int i=0; i < disp_quarks.size(); ++i)
      {
        // Make some shorthands to ease my brain
        map<int, LatticeFermion>& disp_q = disp_quarks[ i ];
				const QuarkTerm_t& term_q = quark[ qindices[ i ] ];
        // If no entry, then create a displaced version of the quark
        if ( disp_q.find( term_q.displacement ) == disp_q.end() )
        {
          //QDPIO::cout << __func__
          //		   << ": i=" << i
          //		   << " disp=" << quark[i].displacement << " disp=" << term_q.displacement
          //		   << " len =" << quark[i].disp_len			<< " len=" << term_q.disp_len
          //		   << " dir =" << quark[i].disp_dir			<< " dir=" << term_q.disp_dir
          //		   << " spin=" << quark[i].spin		    	<< " spin=" << term_q.spin
          //		   << endl;
          LatticeFermion qq = q[ i ];
          displacementSub( myparams.gaugestuff.u, qq, term_q.disp_len, term_q.disp_dir );
          disp_q.insert( std::make_pair( term_q.displacement, qq ) );
        }
      } // for i
      //QDPIO::cout << __PRETTY_FUNCTION__ << ": exiting" << endl;
      END_CODE();
    } // void GroupBaryonQQQ::displaceQuarks
    

    //! First smear then displace the quarks
    void
    GroupBaryonQQQ::smearDisplaceQuarks( multi1d< map<int, LatticeFermion> >& disp_quarks,
                                         const LatticeFermion& q1,
                                         const LatticeFermion& q2,
                                         const LatticeFermion& q3,
																				 int* qindices
                                       ) const
    { 
			START_CODE();
      //QDPIO::cout << __PRETTY_FUNCTION__ << ": entering" << endl;
      multi1d<LatticeFermion> q( 3 );
      // Rotate the quarks up-front.
      // NOTE: this step assumes the spin rotation commutes with the
      // quark smearing and the displacement
      // However, we are careful about the ordering of the smearing and
      // the displacement
//      q[ 0 ] = rotateMat() * q1;
//      q[ 1 ] = rotateMat() * q2;
//      q[ 2 ] = rotateMat() * q3;
      q[ 0 ] = q1;
      q[ 1 ] = q2;
      q[ 2 ] = q3;
      // Displace after the smearing
      displaceQuarks( disp_quarks, q, qindices );
      //QDPIO::cout << __PRETTY_FUNCTION__ << ": exiting" << endl;
      END_CODE();
    } // void GroupBaryonQQQ::smearDisplaceQuarks

		
		//! hack to avoid compiler errors
		LatticeComplex GroupBaryonQQQ::operator() ( const LatticeFermion& q1,
		                                            const LatticeFermion& q2,
																								const LatticeFermion& q3
																							) const
		{ 
			START_CODE();
			LatticeComplex d;
			return d;
			END_CODE();
		}		
		multi1d<LatticeComplex> GroupBaryonQQQ::operator() ( const LatticeFermion& q1,
		                                                     const LatticeFermion& q2,
																												 const LatticeFermion& q3,
																												 enum PlusMinus isign ) const
		{ 
			START_CODE();
			multi1d<LatticeComplex> d(1);
			return d;
			END_CODE();
		}
    
    
    // for the creation operator		
    //! Compute the operator
    LatticeComplex 
    GroupBaryonQQQ::operator() ( const LatticeFermion& q1,
                                 const LatticeFermion& q2,
                                 const LatticeFermion& q3,
																 int* qindices 
                               ) const
    { 
			START_CODE();
			if ( Nc != 3 ){    /* Code is specific to Ns=4 and Nc=3. */
			  QDPIO::cerr<<" code only works for Nc=3 and Ns=4\n";
			  QDP_abort(111) ;
			}
#if QDP_NC == 3


      // The result of displace and smearing (in some unspecified order here)
      multi1d< map<int, LatticeFermion> > disp_quarks;
      // Depending on whether this is the sink or source, do the 
      // appropriate combination of smearing and displacing
      smearDisplaceQuarks( disp_quarks, q1, q2, q3, qindices );
      // The return
      LatticeComplex d;

      LatticeColorVector c0 = peekSpin( disp_quarks[ 0 ].find( quark[ qindices[ 0 ] ].displacement ) ->second,
					                              quark[ qindices[ 0 ] ].spin );
      LatticeColorVector c1 = peekSpin( disp_quarks[ 1 ].find( quark[ qindices[ 1 ] ].displacement ) ->second,
					                              quark[ qindices[ 1 ] ].spin );
      LatticeColorVector c2 = peekSpin( disp_quarks[ 2 ].find( quark[ qindices[ 2 ] ].displacement ) ->second,
					                              quark[ qindices[ 2 ] ].spin );
/*
      multi1d<LatticeFermion> q( 3 );
      // Rotate the quarks up-front.
      // NOTE: this step assumes the spin rotation commutes with the
      // quark smearing and the displacement
      // However, we are careful about the ordering of the smearing and
      // the displacement
      q[ 0 ] = rotateMat() * q1;
      q[ 1 ] = rotateMat() * q2;
      q[ 2 ] = rotateMat() * q3;
*/
      // Contract over color indices with antisym tensors
      d = colorContract( c0, c1, c2 );

      //QDPIO::cout << __PRETTY_FUNCTION__ << ": exiting" << endl;
      END_CODE();
      return d;

#else
      LatticeComplex d;
      d = zero ; 
      return d;
#endif

    } // LatticeComplex GroupBaryonQQQ::operator()

    // for the annihilation operator
    multi1d<LatticeComplex>
    GroupBaryonQQQ::operator() ( const LatticeFermion& q1,
                                 const LatticeFermion& q2,
                                 const LatticeFermion& q3,
																 int* qindices,
                                 enum PlusMinus isign ) const
    { 
      START_CODE();
      if ( Nc != 3 ){    /* Code is specific to Ns=4 and Nc=3. */
	QDPIO::cerr<<" code only works for Nc=3 and Ns=4\n";
	QDP_abort(111) ;
      }
#if QDP_NC == 3


      //QDPIO::cout << __PRETTY_FUNCTION__ << ": entering" << endl;
      // The result of displace and smearing (in some unspecified order here)
      multi1d< map<int, LatticeFermion> > disp_quarks;
      // Depending on whether this is the sink or source, do the 
      // appropriate combination of smearing and displacing
      smearDisplaceQuarks( disp_quarks, q1, q2, q3, qindices );
      // The return
      multi1d<LatticeComplex> d( 1 ); //( 6 )
      // Contract over color indices with antisym tensors

      LatticeColorVector c0 = peekSpin( disp_quarks[ 0 ].find( quark[ qindices[ 0 ] ].displacement ) ->second,
					                              quark[ qindices[ 0 ] ].spin );
      LatticeColorVector c1 = peekSpin( disp_quarks[ 1 ].find( quark[ qindices[ 1 ] ].displacement ) ->second,
					                              quark[ qindices[ 1 ] ].spin );
      LatticeColorVector c2 = peekSpin( disp_quarks[ 2 ].find( quark[ qindices[ 2 ] ].displacement ) ->second,
					                              quark[ qindices[ 2 ] ].spin );
/*
      multi1d<LatticeFermion> q( 3 );
      // Rotate the quarks up-front.
      // NOTE: this step assumes the spin rotation commutes with the
      // quark smearing and the displacement
      // However, we are careful about the ordering of the smearing and
      // the displacement
      q[ 0 ] = rotateMat() * q1;
      q[ 1 ] = rotateMat() * q2;
      q[ 2 ] = rotateMat() * q3;
*/
      // Contract over color indices with antisym tensors
      d[ 0 ] = colorContract( c0, c1, c2 );

      //QDPIO::cout << __PRETTY_FUNCTION__ << ": exiting" << endl;
      END_CODE();
      return d;
#else

      multi1d<LatticeComplex> d( 1 ); 
      return d;
#endif

    } // multi1d<LatticeComplex> GroupBaryonQQQ::operator()


    //	===========================
    //	GroupBaryonOp stuff
		//	now only used to store data
    //	===========================
  	//! Serialize baryon_operator object
  	multi1d<Complex> 
		BaryonOperator_t::serialize()
  	{
  	  int orderings_size = orderings.size();
  	  int op_size3 = orderings[ 0 ].op.size3();
  	  int op_size2 = orderings[ 0 ].op.size2();
  	  int op_size1 = orderings[ 0 ].op.size1();
  	  int ind_size = orderings[ 0 ].op( 0, 0, 0 ).ind.size();
  	  int elem_size2 = orderings[ 0 ].op( 0, 0, 0 ).ind[ 0 ].elem.size2();
  	  int elem_size1 = orderings[ 0 ].op( 0, 0, 0 ).ind[ 0 ].elem.size1();
  	  Complex ord_sizes, op_sizes1, op_sizes2, elem_sizes;
  	  ord_sizes = cmplx( Real( orderings_size ), Real( zero ) );
  	  op_sizes1 = cmplx( Real( op_size2 ), Real( op_size1 ) );
  	  op_sizes2 = cmplx( Real( op_size3 ), Real( ind_size ) );
  	  elem_sizes = cmplx( Real( elem_size2 ), Real( elem_size1 ) );
      multi1d<Complex> baryop_1d( 4 + orderings_size * op_size3 * op_size2 * op_size1 * ind_size * elem_size2 * elem_size1 );
  	  //QDPIO::cout << "baryop_size=" << baryop_1d.size() << endl; 
  	  int cnt = 0;
  	  baryop_1d[ cnt++ ] = ord_sizes;
  	  baryop_1d[ cnt++ ] = op_sizes1;
  	  baryop_1d[ cnt++ ] = op_sizes2;
  	  baryop_1d[ cnt++ ] = elem_sizes;
  	  for(int s=0; s < orderings.size(); ++s)              // orderings
  	  {
  	    for(int i=0; i < orderings[ s ].op.size3(); ++i)              // op_l
  	      for(int j=0; j < orderings[ s ].op.size2(); ++j)            // op_m
  	        for(int k=0; k < orderings[ s ].op.size1(); ++k)          // op_r
  	          for(int l=0; l < orderings[ s ].op( i, j, k ).ind.size(); ++l)     // ind
  	            for(int a=0; a < orderings[ s ].op( i, j, k ).ind[ l ].elem.size2(); ++a)     // elem_l
  	              for(int b=0; b < orderings[ s ].op( i, j, k ).ind[ l ].elem.size1(); ++b)   // elem_r
  	                baryop_1d[ cnt++ ] = orderings[ s ].op( i, j, k ).ind[ l ].elem( a, b );
  	  }
  	  if ( cnt != baryop_1d.size() )
  	  {
  	    QDPIO::cerr << ": size mismatch in BaryonOperator_t serialization " << cnt << endl;
  	    QDP_abort( 1 );
  	  }
  	  return baryop_1d;
  	}
		
    //! Anonymous namespace
    namespace
    {

      //-------------------- callback functions ---------------------------------------

      //! Call-back function
      BaryonOperator<LatticeFermion>* groupBaryon( XMLReader& xml_in,
                                                   const std::string& path,
                                                   const multi1d<LatticeColorMatrix>& u )
      {
        return new GroupBaryonQQQ( Params( xml_in, path ) );
      }
      //! Local registration flag
      bool registered = false;
			
    }  // end anonymous namespace

    const std::string name = "GROUP_BARYON";

    //! Register all the factories
    bool registerAll()
    {
      bool success = true;
      if ( ! registered )
      {
        // Required stuff
        success &= LinkSmearingEnv::registerAll();
				success &= QuarkSourceSmearingEnv::registerAll();
				success &= QuarkSinkSmearingEnv::registerAll();

        //! Register all the factories
        success &= Chroma::TheWilsonBaryonOperatorFactory::Instance().registerObject(name, groupBaryon );

        registered = true;
      }
      return success;
    } // registerAll()


  } // namespace GroupBaryonOperatorEnv


}  // namespace Chroma

