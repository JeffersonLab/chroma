// $Id: group_baryon_operator_w.cc,v 1.14 2006-11-01 04:42:22 edwards Exp $
/*! \file
 *  \brief Construct group baryon operators
 */
#include "chroma.h"

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
#include "meas/smear/displacement.h"

#include "meas/sources/dilutezN_source_const.h"
#include "meas/sources/zN_src.h"
#include "meas/smear/quark_source_sink.h"

#include "util/ferm/diractodr.h"
#include "util/ft/sftmom.h"

#include <map>

using std::map;

namespace Chroma
{ 
  //! Baryon sequential sources
  /*! \ingroup hadron */
  namespace GroupBaryonOperatorEnv
  {
    // Readers
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

// =======================================
// number of zN, quarks, j_decay are fixed
// to 4, 3, 3
// =======================================
      int NumberofQuarks = 3;
      name = "all-to-all";
			
      int version;
      read( paramtop, "version", version );

      switch ( version )
      {
      case 2:
	break;

      default:
	QDPIO::cerr << name << ": parameter version " << version
		    << " unsupported." << endl;
	QDP_abort( 1 );
      }
			
      nrow.resize(4);
      read( paramtop, "nrow", nrow );
		
      read( paramtop, "Baryon_type", param.baryon_operator );
      //source_quark_smearing = readXMLGroup( paramtop, "SourceQuarkSmearing", "wvf_kind" );
      //sink_quark_smearing = readXMLGroup( paramtop, "SinkQuarkSmearing", "wvf_kind" );
      //link_smearing = readXMLGroup( paramtop, "LinkSmearing", "LinkSmearingType" );
      read( paramtop, "InputFileName", InputFileName );
      read( paramtop, "mom2_max", mom2_max );
      read( paramtop, "j_decay", j_decay );
			
      read( paramtop, "QProps", qprop );
      gaugestuff.u.resize(Nd);
      read( paramtop, "Cfg", gaugestuff.cfg );
      dilution.resize( NumberofQuarks );
      for(int i=0; i < NumberofQuarks; ++i) dilution[ i ].N = 4; // Z(4) noise
      for(int i=0; i < NumberofQuarks; ++i) dilution[ i ].j_decay = 3; // time-direction
      read( paramtop, "DilutionScheme", dilution );
    } // end Params::Params


    // Writer
    void Params::writeXML( XMLWriter& xml, const string& path ) const
    {
      push( xml, path );

      int version = 1;
      write( xml, "version", version );
      write( xml, "nrow", nrow );
      //write( xml, "BaryonOperatorType", GroupBaryonOperatorEnv::name );
      xml << source_quark_smearing.xml;
      xml << sink_quark_smearing.xml;
      xml << link_smearing.xml;
      write( xml, "InputFileName", InputFileName );
      //
      // from inline_stoch_baryon_w.cc
      //
      write( xml, "mom2_max", mom2_max );
      write( xml, "j_decay",  j_decay );
      pop( xml );
    } // end void Params::writeXML
		
    /* XML input file
       <Qprops>                                    struct Prop_t
       <Qsolutions>				                      { 														                       
       <elem>																	  //! Operators
       <file_names>													  struct Operator_t
       <elem>./zN_prop_q1_t0</elem>				  {
       <elem>./zN_prop_q1_t1</elem>				  	multi1d<std::string> soln_files;
       <elem>./zN_prop_q1_t2</elem>				  };													   
       <elem>./zN_prop_q1_t3</elem>				  std::string          op_file;
       </file_names> 												  multi1d<Operator_t>  op;
       </elem>                                 };
       <elem>                                    
       <file_names>                          
       <elem>./zN_prop_q2_t0</elem>        
       <elem>./zN_prop_q2_t1</elem>        
       <elem>./zN_prop_q2_t2</elem>				
       <elem>./zN_prop_q2_t3</elem>				struct QProp_t
       </file_names> 												{
       </elem> 																  struct QFiles_t
       <elem>																	  {
       <file_names>															multi1d<std::string> soln_file_names;
       <elem>./zN_prop_q3_t0</elem>				  };
       <elem>./zN_prop_q3_t1</elem>				  multi1d<QFiles_t> Quark; // quark index (0,1,2)
       <elem>./zN_prop_q3_t2</elem>				};
       <elem>./zN_prop_q3_t3</elem>				QProp_t Qprop;
       </file_names> 												Qprop.Quark[ 0,1,2 ].SolnFileName[ i ]
       </elem> 																
       </Qsolutions>			                        read( QProp_t& input )  input.Quark(3)
       </Qprops> 																	read( QFiles_t& input ) input.soln_file_names(Ndil)
    */

    //	====================
    //	GroupBaryonQQQ stuff
    //	====================
		
    //! Constructor
    GroupBaryonQQQ::GroupBaryonQQQ()
    {
      // The spin basis matrix to goto Dirac
      spin_rotate_mat = adj( DiracToDRMat() );
    }
    //! Full constructor
    GroupBaryonQQQ::GroupBaryonQQQ( const Params& p, const multi1d<LatticeColorMatrix>& u_ ) :
      myparams( p ), u_smr( u_ )
    {
      // obsolete
      //readCoeffs( coeffs );
      // The spin basis matrix to goto Dirac
      spin_rotate_mat = adj( DiracToDRMat() );
      // Factory constructions
      try
      {
        // Smear the gauge field if needed
        {
          std::istringstream xml_l( myparams.link_smearing.xml );
          XMLReader linktop( xml_l );
          const string link_path = "/LinkSmearing";
          QDPIO::cout << "Link smearing type = " << myparams.link_smearing.id << endl;

          Handle< LinkSmearing >
	    linkSmearing( TheLinkSmearingFactory::Instance().createObject( 
			    myparams.link_smearing.id,
			    linktop,
			    link_path ) );
          ( *linkSmearing ) ( u_smr );
        }

        // Create the source quark smearing object
        {
          std::istringstream xml_s( myparams.source_quark_smearing.xml );
          XMLReader smeartop( xml_s );
          const string smear_path = "/SourceQuarkSmearing";
          QDPIO::cout << "Source quark smearing type = " << myparams.source_quark_smearing.id << endl;

          sourceQuarkSmearing =
            TheFermSmearingFactory::Instance().createObject( myparams.source_quark_smearing.id,
							     smeartop,
							     smear_path );
        }

        // Create the sink quark smearing object
        {
          std::istringstream xml_s( myparams.sink_quark_smearing.xml );
          XMLReader smeartop( xml_s );
          const string smear_path = "/SinkQuarkSmearing";
          QDPIO::cout << "Sink quark smearing type = " << myparams.sink_quark_smearing.id << endl;

          sinkQuarkSmearing =
            TheFermSmearingFactory::Instance().createObject( myparams.sink_quark_smearing.id,
							     smeartop,
							     smear_path );
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
      QDPIO::cout << "Reading input from text file : " << params.InputFileName <<endl;
      TextReader reader( params.InputFileName );
      //	
      // The spin basis matrix to goto Dirac
      //      spin_rotate_mat = adj( DiracToDRMat() );
      // ------------------------------------------------------------------------------------
      // Lattice sizes
      reader >> params.nrow[ 0 ] >> params.nrow[ 1 ] >> params.nrow[ 2 ] >> params.nrow[ 3 ];
      // ------------------------------------------------------------------------------------
      //
      // Hybrid list sizes
      params.NH.resize( 1 ); // only one ordering now [012] ... hard-wired
      for ( int i = 0; i < params.NH.size(); ++i )
        params.NH[ i ].resize( 3 ); // L,M,R quarks

      params.Nmomenta = 1; // need to change this in the future when we do multiple momenta

      // Set the hybrid list sizes
      //   { something like this would be better                  }
      //   { Ndil = params.named_obj.prop.op[n].soln_files.size() }
      //XMLReader& xml;
      //read( xml, "Ndilutions", params.NH[0] );
      // ------------------------------------------------------------------------------------
      reader >> params.NH[ 0 ][ 0 ] >> params.NH[ 0 ][ 1 ] >> params.NH[ 0 ][ 2 ];
      // ------------------------------------------------------------------------------------
      //for (int i = 0; i < params.qprop.solns.size(); ++i)
      //	params.qprop.solns[ i ].soln_file_names.resize( params.NH[ 0 ][ i ] );

#ifdef OLDSTUFF
      // 012
      params.NH[ 0 ][ 0 ] = NdilL;
      params.NH[ 0 ][ 1 ] = NdilM;
      params.NH[ 0 ][ 2 ] = NdilR;
      // 201
      params.NH[ 1 ][ 0 ] = params.NH[ 0 ][ 2 ]; //NdilR;
      params.NH[ 1 ][ 1 ] = params.NH[ 0 ][ 0 ]; //NdilL;
      params.NH[ 1 ][ 2 ] = params.NH[ 0 ][ 1 ]; //NdilM;
      // 210
      params.NH[ 2 ][ 0 ] = params.NH[ 0 ][ 2 ]; //NdilR;
      params.NH[ 2 ][ 1 ] = params.NH[ 0 ][ 1 ]; //NdilM;
      params.NH[ 2 ][ 2 ] = params.NH[ 0 ][ 0 ]; //NdilL;
      // 102
      params.NH[ 3 ][ 0 ] = params.NH[ 0 ][ 1 ]; //NdilM;
      params.NH[ 3 ][ 1 ] = params.NH[ 0 ][ 0 ]; //NdilL;
      params.NH[ 3 ][ 2 ] = params.NH[ 0 ][ 2 ]; //NdilR;
#endif

      //read( xml, "NumSourcePerm", params.NsrcOrderings );
      // ------------------------------------------------------------------------------------
      reader >> params.NsrcOrderings;
      // ------------------------------------------------------------------------------------

      params.SrcOrderings.resize( params.NsrcOrderings );
      for ( int i = 0; i < params.NsrcOrderings; ++i )
      {
        params.SrcOrderings[ i ].resize( 3 );
	// ------------------------------------------------------------------------------------
        reader >> params.SrcOrderings[ i ][ 0 ] >> params.SrcOrderings[ i ][ 1 ] >> params.SrcOrderings[ i ][ 2 ];
	// ------------------------------------------------------------------------------------
      }
      //read( xml, "SrcQuarkIndices", params.SrcOrderings );

      //read( xml, "NumSinkPerm", params.NsnkOrderings );
      reader >> params.NsnkOrderings;

      params.SnkOrderings.resize( params.NsnkOrderings );
      for ( int i = 0; i < params.NsnkOrderings; ++i )
      {
        params.SnkOrderings[ i ].resize( 3 );
	// ------------------------------------------------------------------------------------
        reader >> params.SnkOrderings[ i ][ 0 ] >> params.SnkOrderings[ i ][ 1 ] >> params.SnkOrderings[ i ][ 2 ];
	// ------------------------------------------------------------------------------------
      }
      //read( xml, "SnkQuarkIndices", params.SnkOrderings );

      //read( xml, "NumOperators", params.Noperators );
      // ------------------------------------------------------------------------------------
      reader >> params.Noperators;
      // ------------------------------------------------------------------------------------
			
      params.Names.resize( params.Noperators );
      // The Baryon operator names ; G1g_L3_TDT_25 ...
      int nameindex;
      for ( int i = 0; i < params.Names.size(); ++i )
      {
	// ------------------------------------------------------------------------------------
        reader >> nameindex >> params.Names[ i ];
	// ------------------------------------------------------------------------------------
      }
			
      //read( xml, "NumDistinctQQQ", params.NQQQs );
      // ------------------------------------------------------------------------------------
      reader >> params.NQQQs;
      // ------------------------------------------------------------------------------------
			
      //read( xml, "OperatorNames", params.Names );
      AB.resize( params.Noperators );
      CB.resize( params.Noperators );
      for ( int n = 0; n < params.Noperators; ++n )
      {				
	AB[ n ].termInCorr.resize( 1 );
        CB[ n ].termInCorr.resize( 1 );
      }
      for ( int n = 0; n < params.Noperators; ++n )
      {
        AB[ n ].termInCorr[ 0 ].hlist.resize( params.NH[ 0 ][ 0 ], params.NH[ 0 ][ 1 ], params.NH[ 0 ][ 2 ] );
	//AB[ n ].termInCorr[ 1 ].hlist.resize( params.NH[ 0 ][ 0 ], params.NH[ 0 ][ 1 ], params.NH[ 0 ][ 2 ] );
        CB[ n ].termInCorr[ 0 ].hlist.resize( params.NH[ 0 ][ 0 ], params.NH[ 0 ][ 1 ], params.NH[ 0 ][ 2 ] );
	//CB[ n ].termInCorr[ 1 ].hlist.resize( params.NH[ 0 ][ 0 ], params.NH[ 0 ][ 1 ], params.NH[ 0 ][ 2 ] );
      }
      for ( int n = 0; n < params.Noperators; ++n )
      {
        for ( int i = 0; i < params.NH[ 0 ][ 0 ]; ++i )
          for ( int j = 0; j < params.NH[ 0 ][ 1 ]; ++j )
            for ( int k = 0; k < params.NH[ 0 ][ 2 ]; ++k )
            {
	      //AB[ n ].termInCorr[ 0 ].hlist( i, j, k ).mom.resize( params.Nmomenta );
              AB[ n ].termInCorr[ 0 ].hlist( i, j, k ).mom.resize( params.Nmomenta, params.nrow[3] );
	      //AB[ n ].termInCorr[ 1 ].hlist( i, j, k ).mom.resize( params.Nmomenta );
              CB[ n ].termInCorr[ 0 ].hlist( i, j, k ).mom.resize( params.Nmomenta, params.nrow[3] );
	      //CB[ n ].termInCorr[ 1 ].hlist( i, j, k ).mom.resize( params.Nmomenta );
            }
      }

      AQQQ.resize( params.NQQQs );
      for ( int i = 0; i < params.NQQQs; ++i )
      {
        AQQQ[ i ].myparams = params;
        AQQQ[ i ].u_smr = params.gaugestuff.u;
        AQQQ[ i ].quark.resize( 3 );
        AQQQ[ i ].orderings.resize( params.NsnkOrderings );
        for ( int j = 0; j < params.NsnkOrderings; ++j )
        {
          AQQQ[ i ].orderings[ j ].hlist.resize( params.NH[ j ][ 0 ], params.NH[ j ][ 1 ], params.NH[ j ][ 2 ] );
        }
      }
      CQQQ.resize( params.NQQQs );
      for ( int i = 0; i < params.NQQQs; ++i )
      {
        CQQQ[ i ].myparams = params;
        CQQQ[ i ].u_smr = params.gaugestuff.u;
        CQQQ[ i ].quark.resize( 3 );
        CQQQ[ i ].orderings.resize( params.NsrcOrderings );
        for ( int j = 0; j < params.NsrcOrderings; ++j )
        {
          CQQQ[ i ].orderings[ j ].hlist.resize( params.NH[ j ][ 0 ], params.NH[ j ][ 1 ], params.NH[ j ][ 2 ] );
        }
      }

      for ( int q = 0; q < params.NQQQs; ++q )
      {
        for ( int p = 0; p < params.NsrcOrderings; ++p )
        {
          for ( int i = 0; i < params.NH[ p ][ 0 ]; ++i )
            for ( int j = 0; j < params.NH[ p ][ 1 ]; ++j )
              for ( int k = 0; k < params.NH[ p ][ 2 ]; ++k )
              {
                CQQQ[ q ].orderings[ p ].hlist( i, j, k ).mom.resize( params.Nmomenta, params.nrow[3] );
              }
        }
        for ( int p = 0; p < params.NsnkOrderings; ++p )
        {
          for ( int i = 0; i < params.NH[ p ][ 0 ]; ++i )
            for ( int j = 0; j < params.NH[ p ][ 1 ]; ++j )
              for ( int k = 0; k < params.NH[ p ][ 2 ]; ++k )
              {
                AQQQ[ q ].orderings[ p ].hlist( i, j, k ).mom.resize( params.Nmomenta, params.nrow[3] );
              }
        }
      }

      for ( int i = 0; i < params.NQQQs; ++i )
      {
	int dL, hash;
	// ------------------------------------------------------------------------------------
        reader >> hash
	       >> AQQQ[ i ].quark[ 0 ].spin
	       >> AQQQ[ i ].quark[ 1 ].spin
	       >> AQQQ[ i ].quark[ 2 ].spin
	       >> AQQQ[ i ].quark[ 0 ].displacement
	       >> AQQQ[ i ].quark[ 1 ].displacement
	       >> AQQQ[ i ].quark[ 2 ].displacement
	       >> dL
	       >> AQQQ[ i ].NBaryonOps;
	// ------------------------------------------------------------------------------------

        AQQQ[ i ].QQQIndex = i;

        for ( int j = 0; j < 3; ++j )
        {
          // Make spin index 0 based
          AQQQ[ i ].quark[ j ].spin -= 1;
          if ( AQQQ[ i ].quark[ j ].displacement == 0 )
          {
            AQQQ[ i ].quark[ j ].disp_dir = 0;
            AQQQ[ i ].quark[ j ].disp_len = 0;
            AQQQ[ i ].quark[ j ].disp_ind = 0;
          }
          else if ( AQQQ[ i ].quark[ j ].displacement > 0 )
          {
            AQQQ[ i ].quark[ j ].disp_dir = AQQQ[ i ].quark[ j ].displacement - 1;
            AQQQ[ i ].quark[ j ].disp_len = dL;
            AQQQ[ i ].quark[ j ].disp_ind = AQQQ[ i ].quark[ j ].displacement;
          }
          else
          {
            AQQQ[ i ].quark[ j ].disp_dir = -AQQQ[ i ].quark[ j ].displacement - 1;
            AQQQ[ i ].quark[ j ].disp_len = -dL;
            AQQQ[ i ].quark[ j ].disp_ind = -AQQQ[ i ].quark[ j ].displacement + 3;
          }
        }
        AQQQ[ i ].coef.resize( AQQQ[ i ].NBaryonOps );
        AQQQ[ i ].whichBaryonOps.resize( AQQQ[ i ].NBaryonOps );
        AQQQ[ i ].baryon.resize( AQQQ[ i ].NBaryonOps );
        //             name                                    coefficient
        // AQQQ[ <mu,nu,tau> ].whichBaryonOps[ n ]    ,   AQQQ[ <mu,nu,tau> ].coef[ n ]
        //
        for ( int n = 0; n < AQQQ[ i ].NBaryonOps; ++n )
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
	  AQQQ[ i ].baryon[ n ] = &AB[ AQQQ[ i ].whichBaryonOps[ n ] ];
        }
      } // i : AQQQ[] loop

      //
      // Now for the creation operators
      //
      // rotate by gamma4's to get the barred coefficients
      //       index     0  1   2   3
      //      gamma4 = ( 1, 1, -1, -1 )
      // and take the complex conjugate
      //          displacement directions ??????
      for ( int i = 0; i < params.NQQQs; ++i )
      {
        CQQQ[ i ].QQQIndex = i;
        int spinCount, flipsign;
        spinCount = 0;
        flipsign = 0;
        for ( int j = 0; j < 3; ++j )
        {
          CQQQ[ i ].quark[ j ].displacement = AQQQ[ i ].quark[ j ].displacement;
          CQQQ[ i ].quark[ j ].disp_dir = AQQQ[ i ].quark[ j ].disp_dir; // sign
          CQQQ[ i ].quark[ j ].disp_ind = AQQQ[ i ].quark[ j ].disp_ind;
          CQQQ[ i ].quark[ j ].spin = AQQQ[ i ].quark[ j ].spin;
          if ( CQQQ[ i ].quark[ j ].spin > 1 ) spinCount++;
        }
        if ( spinCount % 2 ) flipsign = 1;

        CQQQ[ i ].NBaryonOps = AQQQ[ i ].NBaryonOps;
        CQQQ[ i ].coef.resize( CQQQ[ i ].NBaryonOps );
        CQQQ[ i ].whichBaryonOps.resize( CQQQ[ i ].NBaryonOps );
        CQQQ[ i ].baryon.resize( CQQQ[ i ].NBaryonOps );
        //
	//             name                           coefficient
        // CQQQ[ <mu,nu,tau> ].whichBaryonOps[ n ] , CQQQ[ i ].coef[ n ]
        //
        for ( int n = 0; n < CQQQ[ i ].NBaryonOps; ++n )
        {
          Real re, im;
          CQQQ[ i ].whichBaryonOps[ n ] = AQQQ[ i ].whichBaryonOps[ n ];
          re = real( AQQQ[ i ].coef[ n ] );
          im = imag( AQQQ[ i ].coef[ n ] );
          if ( flipsign )
          {
            // cc and -sign
            CQQQ[ i ].coef[ n ] = cmplx( -re, im );
          }
          else
          {
            // just cc
            CQQQ[ i ].coef[ n ] = cmplx( re, -im );
          }
	  //
	  // Storage will be in the GroupBaryonOp class
	  //
	  CQQQ[ i ].baryon[ n ] = &CB[ CQQQ[ i ].whichBaryonOps[ n ] ];
        }
      } // i : CQQQ[] loop
      reader.close();
    } // void ReadTextInput
		
		
    //! Construct the source
    LatticeFermion
    GroupBaryonQQQ::operator() ( Params::dilution_t diln ) const
    {
      QDPIO::cout << "Diluted random complex ZN source" << endl;
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
      for ( int j = 0; j < diln.spatial_mask_size.size(); ++j, ++mu )
      {
        if ( j == diln.j_decay ) ++mu;  // bump up to next dir
        lookup_dir[ j ] = mu;
      }
      for ( int j = 0; j < diln.spatial_mask.size(); ++j )
      {
        if ( diln.spatial_mask[ j ].size() != Nd - 1 )
        {
          QDPIO::cerr << name << ": spatial mask incorrect 3" << endl;
          QDP_abort( 1 );
        }
      }
      for ( int c = 0; c < diln.color_mask.size(); ++c )
      {
        if ( diln.color_mask[ c ] < 0 || diln.color_mask[ c ] >= Nc )
        {
          QDPIO::cerr << name << ": color mask incorrect 6" << endl;
          QDP_abort( 1 );
        }
      }
      for ( int s = 0; s < diln.spin_mask.size(); ++s )
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
      for ( int s = 0; s < diln.spin_mask.size(); ++s )
      {
        int spin_source = diln.spin_mask[ s ];
        LatticeColorVector colvec = peekSpin( quark_noise, spin_source );
        LatticeColorVector dest = zero;
        for ( int c = 0; c < diln.color_mask.size(); ++c )
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
      for ( int n = 0; n < diln.spatial_mask.size(); ++n )
      {
        LatticeBoolean btmp = true;
        for ( int j = 0; j < diln.spatial_mask[ n ].size(); ++j )
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
                                    enum PlusMinus isign ) const
    {
      START_CODE();
      //QDPIO::cout << __PRETTY_FUNCTION__ << ": entering" << endl;
      /* old 
	 struct CoeffTerm_t
	 {
	 struct QuarkTerm_t
	 {
	 int  displacement;
	 int  spin;        
	 int  disp_dir;    
	 int  disp_len;    
	 };
	 multi1d<QuarkTerm_t>  quark;
	 Complex               coeff;
	 };
	 multi1d< multi1d< CoeffTerm_t > >  coeffs;
	 const CoeffTerm_t& term = coeffs[ n ][ l ];
      */ 	
      // three displaced quarks
      disp_quarks.resize( 3 );
      for ( int i = 0; i < disp_quarks.size(); ++i )
      {
        // Make some shorthands to ease my brain
        map<int, LatticeFermion>& disp_q = disp_quarks[ i ];
        const QuarkTerm_t& term_q = quark[ i ];
        // If no entry, then create a displaced version of the quark
        if ( disp_q.find( term_q.displacement ) == disp_q.end() )
        {
          //	      QDPIO::cout << __func__
          //		   << ": n=" << n
          //		   << " l=" << l
          //		   << " i=" << i
          //		   << " disp=" << quark[i].displacement << " disp=" << term_q.displacement
          //		   << " len=" << quark[i].disp_len			<< " len=" << term_q.disp_len
          //		   << " dir=" << quark[i].disp_dir			<< " dir=" << term_q.disp_dir
          //		   << endl;
          LatticeFermion qq = q[ i ];
          switch ( isign )
          {
	  case PLUS:
	    displacement( u_smr, qq, term_q.disp_len, term_q.disp_dir );
	    break;

	  case MINUS:
	    displacement( u_smr, qq, -term_q.disp_len, term_q.disp_dir );
	    break;
          }
          // Insert
          //disp_q.insert( std::make_pair( term_q.displacement, qq ) );
          disp_q.insert( std::make_pair( term_q.disp_ind, qq ) );
        }
      } // for i

      //QDPIO::cout << __PRETTY_FUNCTION__ << ": exiting" << endl;
      END_CODE();
    } // void GroupBaryonQQQ::displaceQuarks


    //! First displace then smear the quarks
    void
    GroupBaryonQQQ::displaceSmearQuarks( multi1d< map<int, LatticeFermion> >& disp_quarks,
                                         const LatticeFermion& q1,
                                         const LatticeFermion& q2,
                                         const LatticeFermion& q3,
                                         enum PlusMinus isign ) const
    {
      START_CODE();
//      QDPIO::cout << __PRETTY_FUNCTION__ << ": entering" << endl;

      multi1d<LatticeFermion> q( 3 );

      // Rotate the quarks up-front.
      // NOTE: this step assumes the spin rotation commutes with the
      // quark smearing and the displacement
      // However, we are careful about the ordering of the smearing and
      // the displacement
      q[ 0 ] = rotateMat() * q1;
      q[ 1 ] = rotateMat() * q2;
      q[ 2 ] = rotateMat() * q3;

      // Displace
      displaceQuarks( disp_quarks, q, isign );

      // Source smear after the displacements
      for ( int i = 0; i < disp_quarks.size(); ++i )
      {
        // Make some shorthands to ease my brain
        map<int, LatticeFermion>& disp_q = disp_quarks[ i ];

        // Loop over all keys
        for ( std::map<int, LatticeFermion>::const_iterator mm = disp_q.begin();
              mm != disp_q.end();
              ++mm )
        {
          ( *sourceQuarkSmearing ) ( disp_q[ mm->first ], u_smr );
        }
      }

      //QDPIO::cout << __PRETTY_FUNCTION__ << ": exiting" << endl;
      END_CODE();
    } // GroupBaryonQQQ::displaceSmearQuarks


    //! First smear then displace the quarks
    void
    GroupBaryonQQQ::smearDisplaceQuarks( multi1d< map<int, LatticeFermion> >& disp_quarks,
                                         const LatticeFermion& q1,
                                         const LatticeFermion& q2,
                                         const LatticeFermion& q3,
                                         enum PlusMinus isign ) const
    {
      START_CODE();
//      QDPIO::cout << __PRETTY_FUNCTION__ << ": entering" << endl;

      multi1d<LatticeFermion> q( 3 );

      // Rotate the quarks up-front.
      // NOTE: this step assumes the spin rotation commutes with the
      // quark smearing and the displacement
      // However, we are careful about the ordering of the smearing and
      // the displacement
      q[ 0 ] = rotateMat() * q1;
      q[ 1 ] = rotateMat() * q2;
      q[ 2 ] = rotateMat() * q3;

      // Sink smear the quarks
/*
  for ( int i = 0; i < q.size(); ++i ) 
  {
  //QDPIO::cout<<"attempting to smear sink"<<endl;
  ( *sinkQuarkSmearing ) ( q[ i ], u_smr );
  //QDPIO::cout<<"in smearDisplace: smeared"<<endl;
  }
*/
      // Displace after the smearing
      displaceQuarks( disp_quarks, q, isign );

      //QDPIO::cout << __PRETTY_FUNCTION__ << ": exiting" << endl;
      END_CODE();
    } // void GroupBaryonQQQ::smearDisplaceQuarks


    //! Manipulate the quark fields
    void
    GroupBaryonQQQ::quarkManip( multi1d< map<int, LatticeFermion> >& disp_quarks,
                                const LatticeFermion& q1,
                                const LatticeFermion& q2,
                                const LatticeFermion& q3,
                                enum PlusMinus isign ) const
    {
      START_CODE();
//      QDPIO::cout << __PRETTY_FUNCTION__ << ": entering" << endl;

      // Depending on whether this is the sink or source, do the appropriate
      // combination of smearing and displacing
      switch ( isign )
      {
      case PLUS:
	// Sink
	smearDisplaceQuarks( disp_quarks, q1, q2, q3, isign );
	break;

      case MINUS:
	// Source
	displaceSmearQuarks( disp_quarks, q1, q2, q3, isign );
	break;

      default:
	QDPIO::cerr << name << ": illegal isign" << endl;
	QDP_abort( 1 );
      }

      //QDPIO::cout << __PRETTY_FUNCTION__ << ": exiting" << endl;
      END_CODE();
    } // void GroupBaryonQQQ::quarkManip

    
    //! Compute the operator
    multi1d<LatticeComplex>
    GroupBaryonQQQ::operator() ( const LatticeFermion& q1,
                                 const LatticeFermion& q2,
                                 const LatticeFermion& q3,
                                 enum PlusMinus isign ) const
    { 
      START_CODE();
//      QDPIO::cout << __PRETTY_FUNCTION__ << ": entering" << endl;
      // The result of displace and smearing (in some unspecified order here)
      multi1d< map<int, LatticeFermion> > disp_quarks;
      // Depending on whether this is the sink or source, do the appropriate
      // combination of smearing and displacing
      quarkManip( disp_quarks, q1, q2, q3, isign );
      // The return
      multi1d<LatticeComplex> d( 1 ); //( 6 )
      // Contract over color indices with antisym tensors

      LatticeFermion qtemp;
      SftMom phases_nomom(0, true, 3);
      multi1d<Double> mycorr = sumMulti( localNorm2(disp_quarks[0].find( quark[ 0 ].displacement ) ->second), phases_nomom.getSet() );
      QDPIO::cout<<"first_corr="<<mycorr[0]<<mycorr[1]<<mycorr[2]<<mycorr[3]<<endl;
      QDPIO::cout<<"spins "<<quark[ 0 ].spin<<" "<<quark[ 1 ].spin<<" "<<quark[ 2 ].spin<<endl;
      switch ( isign )
      {
      case MINUS:
	// Source			
	//d.resize(6);
	// mu,nu,tau contraction
	d[ 0 ] = colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 0 ].disp_ind ) ->second,
					  quark[ 0 ].spin ),
				peekSpin( disp_quarks[ 1 ].find( quark[ 1 ].disp_ind ) ->second,
					  quark[ 1 ].spin ),
				peekSpin( disp_quarks[ 2 ].find( quark[ 2 ].disp_ind ) ->second,
					  quark[ 2 ].spin ) );
	// tau,nu,mu contraction
	//d[ 3 ] = colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 2 ].disp_ind ) ->second,
	d[ 0 ] -= colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 2 ].disp_ind ) ->second,
					   quark[ 2 ].spin ),
				 peekSpin( disp_quarks[ 1 ].find( quark[ 1 ].disp_ind ) ->second,
					   quark[ 1 ].spin ),
				 peekSpin( disp_quarks[ 2 ].find( quark[ 0 ].disp_ind ) ->second,
					   quark[ 0 ].spin ) );
	//d[ 3 ] *= -2.0;
	d[ 0 ] *= 2.0;
	// nu,mu,tau contraction
	//d[ 1 ] = colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 1 ].disp_ind ) ->second,
	d[ 0 ] += colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 1 ].disp_ind ) ->second,
					   quark[ 1 ].spin ),
				 peekSpin( disp_quarks[ 1 ].find( quark[ 0 ].disp_ind ) ->second,
					   quark[ 0 ].spin ),
				 peekSpin( disp_quarks[ 2 ].find( quark[ 2 ].disp_ind ) ->second,
					   quark[ 2 ].spin ) );
	// mu,tau,nu contraction
	//d[ 2 ] = colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 0 ].disp_ind ) ->second,
	d[ 0 ] += colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 0 ].disp_ind ) ->second,
					   quark[ 0 ].spin ),
				 peekSpin( disp_quarks[ 1 ].find( quark[ 2 ].disp_ind ) ->second,
					   quark[ 2 ].spin ),
				 peekSpin( disp_quarks[ 2 ].find( quark[ 1 ].disp_ind ) ->second,
					   quark[ 1 ].spin ) );
	// nu,tau,mu contraction
	//d[ 4 ]= -colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 1 ].disp_ind ) ->second,
	d[ 0 ] -= colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 1 ].disp_ind ) ->second,
					   quark[ 1 ].spin ),
				 peekSpin( disp_quarks[ 1 ].find( quark[ 2 ].disp_ind ) ->second,
					   quark[ 2 ].spin ),
				 peekSpin( disp_quarks[ 2 ].find( quark[ 0 ].disp_ind ) ->second,
					   quark[ 0 ].spin ) );
	// tau,mu,nu contraction
	//d[ 5 ]= -colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 2 ].disp_ind ) ->second,
	d[ 0 ] -= colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 2 ].disp_ind ) ->second,
					   quark[ 2 ].spin ),
				 peekSpin( disp_quarks[ 1 ].find( quark[ 0 ].disp_ind ) ->second,
					   quark[ 0 ].spin ),
				 peekSpin( disp_quarks[ 2 ].find( quark[ 1 ].disp_ind ) ->second,
					   quark[ 1 ].spin ) );
	break;

      case PLUS:
	// Sink
	//d.resize(1);
	d[ 0 ] = colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 0 ].disp_ind ) ->second,
					  quark[ 0 ].spin ),
				peekSpin( disp_quarks[ 1 ].find( quark[ 1 ].disp_ind ) ->second,
					  quark[ 1 ].spin ),
				peekSpin( disp_quarks[ 2 ].find( quark[ 2 ].disp_ind ) ->second,
					  quark[ 2 ].spin ) );
	mycorr = sumMulti( localNorm2(d[0]), phases_nomom.getSet() );
	QDPIO::cout<<"Source_corr="<<mycorr[0]<<mycorr[1]<<mycorr[2]<<mycorr[3]<<endl;
	break;

      default:
	QDPIO::cerr << name << ": illegal isign" << endl;
	QDP_abort( 1 );
			
      }			
      END_CODE();
      //QDPIO::cout << __PRETTY_FUNCTION__ << ": exiting" << endl;
      return d;
    } // multi1d<LatticeComplex> GroupBaryonQQQ::operator()



    //	====================
    //	GroupBaryonOp stuff
    //	====================
		
    // most of this is obsolete ... see GroupBaryonQQQ now
#if 1
    //! Full constructor
    //GroupBaryonOp::GroupBaryonOp( const Params& p, multi1d<LatticeColorMatrix>& u_ ) :
    GroupBaryonOp::GroupBaryonOp( const Params& p, const multi1d<LatticeColorMatrix>& u_ ) :
      myparams( p ), u_smr( u_ )
    {
      //readCoeffs( coeffs );

      // The spin basis matrix to goto Dirac
      spin_rotate_mat = adj( DiracToDRMat() );

      // Factory constructions
      try
      {
        // Smear the gauge field if needed
        {
          std::istringstream xml_l( myparams.link_smearing.xml );
          XMLReader linktop( xml_l );
          const string link_path = "/LinkSmearing";
          QDPIO::cout << "Link smearing type = " << myparams.link_smearing.id << endl;

          Handle< LinkSmearing >
	    linkSmearing( TheLinkSmearingFactory::Instance().createObject( 
			    myparams.link_smearing.id,
			    linktop,
			    link_path ) );
          ( *linkSmearing ) ( u_smr );
        }

        // Create the source quark smearing object
        {
          std::istringstream xml_s( myparams.source_quark_smearing.xml );
          XMLReader smeartop( xml_s );
          const string smear_path = "/SourceQuarkSmearing";

          sourceQuarkSmearing =
            TheFermSmearingFactory::Instance().createObject( 
	      myparams.source_quark_smearing.id,
	      smeartop,
	      smear_path );
        }

        // Create the sink quark smearing object
        {
          std::istringstream xml_s( myparams.sink_quark_smearing.xml );
          XMLReader smeartop( xml_s );
          const string smear_path = "/SinkQuarkSmearing";

          sinkQuarkSmearing =
            TheFermSmearingFactory::Instance().createObject( 
	      myparams.sink_quark_smearing.id,
	      smeartop,
	      smear_path );
        }
      }
      catch ( const std::string & e )
      {
        QDPIO::cerr << name << ": Caught Exception smearing: " << e << endl;
        QDP_abort( 1 );
      }

    } // GroupBaryonOp full constructor

    //! Construct array of maps of displacements
    void
    GroupBaryonOp::displaceQuarks( multi1d< map<int, LatticeFermion> >& disp_quarks,
                                   const multi1d<LatticeFermion>& q,
                                   enum PlusMinus isign ) const
    {
      START_CODE();

      //      QDPIO::cout << __PRETTY_FUNCTION__ << ": entering" << endl;

      disp_quarks.resize( 3 );
      // three displaced quarks
      for ( int i = 0; i < disp_quarks.size(); ++i )
      {
        // Make some shorthands to ease my brain
        map<int, LatticeFermion>& disp_q = disp_quarks[ i ];
        const QuarkTerm_t& term_q = quark[ i ];

        // If no entry, then create a displaced version of the quark
        if ( disp_q.find( term_q.displacement ) == disp_q.end() )
        {
          LatticeFermion qq = q[ i ];

          switch ( isign )
          {
	  case PLUS:
	    displacement( u_smr, qq, term_q.disp_len, term_q.disp_dir );
	    break;

	  case MINUS:
	    displacement( u_smr, qq, -term_q.disp_len, term_q.disp_dir );
	    break;
          }

          // Insert
          disp_q.insert( std::make_pair( term_q.displacement, qq ) );
        }
      } // for i

      //      QDPIO::cout << __PRETTY_FUNCTION__ << ": exiting" << endl;

      END_CODE();
    } // void GroupBaryonOp::displaceQuarks


    //! First displace then smear the quarks
    void
    GroupBaryonOp::displaceSmearQuarks( multi1d< map<int, LatticeFermion> >& disp_quarks,
                                        const LatticeFermion& q1,
                                        const LatticeFermion& q2,
                                        const LatticeFermion& q3,
                                        enum PlusMinus isign ) const
    {
      START_CODE();

      //      QDPIO::cout << __PRETTY_FUNCTION__ << ": entering" << endl;

      multi1d<LatticeFermion> q( 3 );

      // Rotate the quarks up-front.
      // NOTE: this step assumes the spin rotation commutes with the
      // quark smearing and the displacement
      // However, we are careful about the ordering of the smearing and
      // the displacement
      q[ 0 ] = rotateMat() * q1;
      q[ 1 ] = rotateMat() * q2;
      q[ 2 ] = rotateMat() * q3;

      // Displace
      displaceQuarks( disp_quarks, q, isign );

      // Source smear after the displacements
      for ( int i = 0; i < disp_quarks.size(); ++i )
      {
        // Make some shorthands to ease my brain
        map<int, LatticeFermion>& disp_q = disp_quarks[ i ];

        // Loop over all keys
        for ( std::map<int, LatticeFermion>::const_iterator mm = disp_q.begin();
              mm != disp_q.end();
              ++mm )
        {
          ( *sourceQuarkSmearing ) ( disp_q[ mm->first ], u_smr );
        }
      }

      //      QDPIO::cout << __PRETTY_FUNCTION__ << ": exiting" << endl;

      END_CODE();
    } // end GroupBaryonOp::displaceSmearQuarks


    //! First smear then displace the quarks
    void
    GroupBaryonOp::smearDisplaceQuarks( multi1d< map<int, LatticeFermion> >& disp_quarks,
                                        const LatticeFermion& q1,
                                        const LatticeFermion& q2,
                                        const LatticeFermion& q3,
                                        enum PlusMinus isign ) const
    {
      START_CODE();

      //      QDPIO::cout << __PRETTY_FUNCTION__ << ": entering" << endl;

      multi1d<LatticeFermion> q( 3 );

      // Rotate the quarks up-front.
      // NOTE: this step assumes the spin rotation commutes with the
      // quark smearing and the displacement
      // However, we are careful about the ordering of the smearing and
      // the displacement
      q[ 0 ] = rotateMat() * q1;
      q[ 1 ] = rotateMat() * q2;
      q[ 2 ] = rotateMat() * q3;

      // Sink smear the quarks
      for ( int i = 0; i < q.size(); ++i )
        ( *sinkQuarkSmearing ) ( q[ i ], u_smr );

      // Displace after the smearing
      displaceQuarks( disp_quarks, q, isign );

      //      QDPIO::cout << __PRETTY_FUNCTION__ << ": exiting" << endl;

      END_CODE();
    } // end void GroupBaryonOp::smearDisplaceQuarks
		
		
    //! Manipulate the quark fields
    void
    GroupBaryonOp::quarkManip( multi1d< map<int, LatticeFermion> >& disp_quarks,
                               const LatticeFermion& q1,
                               const LatticeFermion& q2,
                               const LatticeFermion& q3,
                               enum PlusMinus isign ) const
    {
      START_CODE();

      //      QDPIO::cout << __PRETTY_FUNCTION__ << ": entering" << endl;

      // Depending on whether this is the sink or source, do the appropriate
      // combination of smearing and displacing
      switch ( isign )
      {
      case PLUS:
	// Sink
	smearDisplaceQuarks( disp_quarks, q1, q2, q3, isign );
	break;

      case MINUS:
	// Source
	displaceSmearQuarks( disp_quarks, q1, q2, q3, isign );
	break;

      default:
	QDPIO::cerr << name << ": illegal isign" << endl;
	QDP_abort( 1 );
      }

      //      QDPIO::cout << __PRETTY_FUNCTION__ << ": exiting" << endl;

      END_CODE();
    } // end void GroupBaryonOp::quarkManip
		
		
    //! Compute the operator
    multi1d<LatticeComplex>
    GroupBaryonOp::operator() ( const LatticeFermion& q1,
                                const LatticeFermion& q2,
                                const LatticeFermion& q3,
                                enum PlusMinus isign ) const
    { START_CODE();
    //      QDPIO::cout << __PRETTY_FUNCTION__ << ": entering" << endl;
    // The result of displace and smearing (in some unspecified order here)
    multi1d< map<int, LatticeFermion> > disp_quarks;
    // Depending on whether this is the sink or source, do the appropriate
    // combination of smearing and displacing
    quarkManip( disp_quarks, q1, q2, q3, isign );
    // The return
    multi1d<LatticeComplex> d;//[ 6 ];
    // Contract over color indices with antisym tensors
    switch ( isign )
    {
    case MINUS:
      // Source			
      d.resize(6);
      // mu,nu,tau contraction
      d[ 0 ] = colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 0 ].displacement ) ->second,
					quark[ 0 ].spin ),
			      peekSpin( disp_quarks[ 1 ].find( quark[ 1 ].displacement ) ->second,
					quark[ 1 ].spin ),
			      peekSpin( disp_quarks[ 2 ].find( quark[ 2 ].displacement ) ->second,
					quark[ 2 ].spin ) );
      d[ 0 ] *= 2.0;
      // nu,mu,tau contraction
      d[ 1 ] = colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 1 ].displacement ) ->second,
					quark[ 1 ].spin ),
			      peekSpin( disp_quarks[ 1 ].find( quark[ 0 ].displacement ) ->second,
					quark[ 0 ].spin ),
			      peekSpin( disp_quarks[ 2 ].find( quark[ 2 ].displacement ) ->second,
					quark[ 2 ].spin ) );
      // mu,tau,nu contraction
      d[ 2 ] = colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 0 ].displacement ) ->second,
					quark[ 0 ].spin ),
			      peekSpin( disp_quarks[ 1 ].find( quark[ 2 ].displacement ) ->second,
					quark[ 2 ].spin ),
			      peekSpin( disp_quarks[ 2 ].find( quark[ 1 ].displacement ) ->second,
					quark[ 1 ].spin ) );
      // tau,nu,mu contraction
      d[ 3 ] = colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 2 ].displacement ) ->second,
					quark[ 2 ].spin ),
			      peekSpin( disp_quarks[ 1 ].find( quark[ 1 ].displacement ) ->second,
					quark[ 1 ].spin ),
			      peekSpin( disp_quarks[ 2 ].find( quark[ 0 ].displacement ) ->second,
					quark[ 0 ].spin ) );
      d[ 3 ] *= -2.0;
      // nu,tau,mu contraction
      d[ 4 ]= -colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 1 ].displacement ) ->second,
					quark[ 1 ].spin ),
			      peekSpin( disp_quarks[ 1 ].find( quark[ 2 ].displacement ) ->second,
					quark[ 2 ].spin ),
			      peekSpin( disp_quarks[ 2 ].find( quark[ 0 ].displacement ) ->second,
					quark[ 0 ].spin ) );
      // tau,mu,nu contraction
      d[ 5 ]= -colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 2 ].displacement ) ->second,
					quark[ 2 ].spin ),
			      peekSpin( disp_quarks[ 1 ].find( quark[ 0 ].displacement ) ->second,
					quark[ 0 ].spin ),
			      peekSpin( disp_quarks[ 2 ].find( quark[ 1 ].displacement ) ->second,
					quark[ 1 ].spin ) );
      break;

    case PLUS:
      // Sink
      d.resize(1);
      d[ 0 ] = colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 0 ].displacement ) ->second,
					quark[ 0 ].spin ),
			      peekSpin( disp_quarks[ 1 ].find( quark[ 1 ].displacement ) ->second,
					quark[ 1 ].spin ),
			      peekSpin( disp_quarks[ 2 ].find( quark[ 2 ].displacement ) ->second,
					quark[ 2 ].spin ) );
      break;

    default:
      QDPIO::cerr << name << ": illegal isign" << endl;
      QDP_abort( 1 );
			
    }
    END_CODE();
    //      QDPIO::cout << __PRETTY_FUNCTION__ << ": exiting" << endl;
    return d;
    } // multi1d<LatticeComplex> GroupBaryonOp::operator()
#endif
    //! Serialize generalized operator object
    multi1d<Complex> 
    GroupBaryonOp::serialize()
    {
      int termInCorr_size = termInCorr.size();
      int hlist_size3 = termInCorr[ 0 ].hlist.size3();
      int hlist_size2 = termInCorr[ 0 ].hlist.size2();
      int hlist_size1 = termInCorr[ 0 ].hlist.size1();
//  	  int mom_size = termInCorr[ 0 ].hlist( 0, 0, 0 ).mom.size();
      int elem_size2 = termInCorr[ 0 ].hlist( 0, 0, 0 ).mom.size2();
      int elem_size1 = termInCorr[ 0 ].hlist( 0, 0, 0 ).mom.size1();
      // dreadful hack - use a complex to hold an int
      //     if you think that's dreadful ...
      int mom_size = 0;
      Complex termInCorr_sizes, hlist_sizes1, hlist_sizes2, elem_sizes;
      termInCorr_sizes = cmplx( Real( termInCorr_size ), Real( zero ) );
      hlist_sizes1 = cmplx( Real( hlist_size2 ), Real( hlist_size1 ) );
      hlist_sizes2 = cmplx( Real( hlist_size3 ), Real( mom_size ) );
      elem_sizes = cmplx( Real( elem_size2 ), Real( elem_size1 ) );

      multi1d<Complex> baryonOp_1d( 4 + termInCorr_size * hlist_size3 * hlist_size2 * hlist_size1 * elem_size2 * elem_size1 );
      //    QDPIO::cout << "baryonprop_size=" << baryonOp_1d.size() << endl;
      int cnt = 0;
      baryonOp_1d[ cnt++ ] = termInCorr_sizes;
      baryonOp_1d[ cnt++ ] = hlist_sizes1;
      baryonOp_1d[ cnt++ ] = hlist_sizes2;
      baryonOp_1d[ cnt++ ] = elem_sizes;
      for ( int s = 0; s < termInCorr.size(); ++s ) // termInCorr
      {
	for ( int i = 0; i < termInCorr[ s ].hlist.size3(); ++i )      // hlist_l
	  for ( int j = 0; j < termInCorr[ s ].hlist.size2(); ++j )    // hlist_m
	    for ( int k = 0; k < termInCorr[ s ].hlist.size1(); ++k )  // hlist_r
	      //for ( int p = 0; p < termInCorr[ s ].hlist( i, j, k ).mom.size(); ++p )  // mom
	      for ( int a = 0; a < termInCorr[ s ].hlist( i, j, k ).mom.size2(); ++a )    // elem_l
		for ( int b = 0; b < termInCorr[ s ].hlist( i, j, k ).mom.size1(); ++b )  // elem_r
		  baryonOp_1d[ cnt++ ] = termInCorr[ s ].hlist( i, j, k ).mom( a, b );
      }
      if ( cnt != baryonOp_1d.size() )
      {
//  	    QDPIO::cerr << GroupBaryonOperatorEnv::name << ": size mismatch in serialization" << endl;
	QDPIO::cerr << ": size mismatch in serialization" << endl;
	QDP_abort( 1 );
      }
      return baryonOp_1d;
    } // multi1d<Complex> GroupBaryonOp::serialize()
		
		
    //! Anonymous namespace
    namespace
    {

      //-------------------- callback functions ---------------------------------------

      //! Call-back function
      BaryonOperator<LatticeFermion>* groupBaryon( XMLReader& xml_in,
                                                   const std::string& path,
                                                   const multi1d<LatticeColorMatrix>& u )
      {
        return new GroupBaryonQQQ( Params( xml_in, path ), u );
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
        success &= QuarkSmearingEnv::registerAll();
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
