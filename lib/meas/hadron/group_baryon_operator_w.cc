// $Id: group_baryon_operator_w.cc,v 1.10 2006-10-13 15:15:37 juge Exp $
/*! \file
 *  \brief Construct group baryon operators
 */

#include "meas/hadron/group_baryon_operator_w.h"
#include "meas/hadron/baryon_operator_factory_w.h"

#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/quark_smearing_aggregate.h"

#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"

#include "meas/smear/displacement.h"
#include "util/ferm/diractodr.h"

#include <map>

using std::map;

namespace Chroma
{
  // Read parameters
  void read( XMLReader& xml, const string& path, GroupBaryonOperatorEnv::Params& param )
  {
    GroupBaryonOperatorEnv::Params tmp( xml, path );
    param = tmp;
  }
  // Writer
  void write( XMLWriter& xml, const string& path, const GroupBaryonOperatorEnv::Params& param )
  {
    param.writeXML( xml, path );
  }


  //! Baryon sequential sources
  /*! \ingroup hadron */
  namespace GroupBaryonOperatorEnv
  {
		//	============
		//	Params stuff
		//	============
		
    //! Initialize
    Params::Params()
    {}

    //! Read parameters
    Params::Params( XMLReader& xml, const string& path )
    {
      XMLReader paramtop( xml, path );

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
			
      source_quark_smearing = readXMLGroup( paramtop, "SourceQuarkSmearing", "wvf_kind" );
      sink_quark_smearing = readXMLGroup( paramtop, "SinkQuarkSmearing", "wvf_kind" );
      link_smearing = readXMLGroup( paramtop, "LinkSmearing", "LinkSmearingType" );

      read( paramtop, "InputFileName", InputFileName );

      read( paramtop, "mom2_max", mom2_max );
      read( paramtop, "j_decay", j_decay );
      read( paramtop, "Seed_l", seed_l );
      read( paramtop, "Seed_m", seed_m );
      read( paramtop, "Seed_r", seed_r );

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
				 write( xml, "Seed_l",   seed_l );
				 write( xml, "Seed_m",   seed_m );
				 write( xml, "Seed_r",   seed_r );

      pop( xml );
    } // end void Params::writeXML
		
		//	====================
		//	GroupBaryonQQQ stuff
		//	====================
		
    //! Full constructor
    GroupBaryonQQQ::GroupBaryonQQQ( const Params& p, const multi1d<LatticeColorMatrix>& u_ ) :
        params( p ), u_smr( u_ )
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
          std::istringstream xml_l( params.link_smearing.xml );
          XMLReader linktop( xml_l );
          const string link_path = "/LinkSmearing";
          QDPIO::cout << "Link smearing type = " << params.link_smearing.id << endl;

          Handle< LinkSmearing >
          linkSmearing( TheLinkSmearingFactory::Instance().createObject( params.link_smearing.id,
                        linktop,
                        link_path ) );
          ( *linkSmearing ) ( u_smr );
        }

        // Create the source quark smearing object
        {
          std::istringstream xml_s( params.source_quark_smearing.xml );
          XMLReader smeartop( xml_s );
          const string smear_path = "/SourceQuarkSmearing";

          sourceQuarkSmearing =
            TheFermSmearingFactory::Instance().createObject( params.source_quark_smearing.id,
                smeartop,
                smear_path );
        }

        // Create the sink quark smearing object
        {
          std::istringstream xml_s( params.sink_quark_smearing.xml );
          XMLReader smeartop( xml_s );
          const string smear_path = "/SinkQuarkSmearing";

          sinkQuarkSmearing =
            TheFermSmearingFactory::Instance().createObject( params.sink_quark_smearing.id,
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
      QDPIO::cout << "Reading input from text file : " << endl;
      TextReader reader( params.InputFileName );

      //multi1d<GroupBaryonOp> AB; // annihilation
      //multi1d<GroupBaryonOp> CB; // creation
      //multi1d<GroupBaryonQQQ> AQQQ; // annihilation
      //multi1d<GroupBaryonQQQ> CQQQ; // creation
			
      // The spin basis matrix to goto Dirac
      //      spin_rotate_mat = adj( DiracToDRMat() ); do something about this

      // hybrid list sizes
      //params.initialize();
      params.NH.resize( 1 ); // only one ordering now [012] ... hard-wired
      for ( int i = 0; i < params.NH.size(); ++i )
        params.NH[ i ].resize( 3 ); // L,M,R quarks

      params.Nmomenta = 1; // need to change this in the future when we do multiple momenta

      // Set the hybrid list sizes
      //   { something like this would be better                  }
      //   { Ndil = params.named_obj.prop.op[n].soln_files.size() }
      //XMLReader& xml;
      //read( xml, "Ndilutions", params.NH[0] );
      reader >> params.NH[ 0 ][ 0 ] >> params.NH[ 0 ][ 1 ] >> params.NH[ 0 ][ 2 ];

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
      reader >> params.NsrcOrderings;

      params.SrcOrderings.resize( params.NsrcOrderings );
      for ( int i = 0; i < params.NsrcOrderings; ++i )
      {
        params.SrcOrderings[ i ].resize( 3 );
        reader >> params.SrcOrderings[ i ][ 0 ] >> params.SrcOrderings[ i ][ 1 ] >> params.SrcOrderings[ i ][ 2 ];
      }
      //read( xml, "SrcQuarkIndices", params.SrcOrderings );

      //read( xml, "NumSinkPerm", params.NsnkOrderings );
      reader >> params.NsnkOrderings;

      params.SnkOrderings.resize( params.NsnkOrderings );
      for ( int i = 0; i < params.NsnkOrderings; ++i )
      {
        params.SnkOrderings[ i ].resize( 3 );
        reader >> params.SnkOrderings[ i ][ 0 ] >> params.SnkOrderings[ i ][ 1 ] >> params.SnkOrderings[ i ][ 2 ];
      }
      //read( xml, "SnkQuarkIndices", params.SnkOrderings );

      //read( xml, "NumOperators", params.Noperators );
      reader >> params.Noperators;

      params.Names.resize( params.Noperators );
      // The Baryon operator names ; G1g_L3_TDT_25 ...
      for ( int i = 0; i < params.Names.size(); ++i )
      {
        reader >> params.Names[ i ];
      }
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
        //AB[ n ].termInCorr[1].hlist.resize( params.NH[0][0], params.NH[0][1], params.NH[0][2] );
        CB[ n ].termInCorr[ 0 ].hlist.resize( params.NH[ 0 ][ 0 ], params.NH[ 0 ][ 1 ], params.NH[ 0 ][ 2 ] );
        //CB[ n ].termInCorr[1].hlist.resize( params.NH[0][0], params.NH[0][1], params.NH[0][2] );
      }
      for ( int n = 0; n < params.Noperators; ++n )
      {
        for ( int i = 0; i < params.NH[ 0 ][ 0 ]; ++i )
          for ( int j = 0; j < params.NH[ 0 ][ 1 ]; ++j )
            for ( int k = 0; k < params.NH[ 0 ][ 2 ]; ++k )
            {
              AB[ n ].termInCorr[ 0 ].hlist( i, j, k ).mom.resize( params.Nmomenta );
              //AB[ n ].termInCorr[1].hlist( i, j, k ).mom.resize(params.Nmomenta);
              CB[ n ].termInCorr[ 0 ].hlist( i, j, k ).mom.resize( params.Nmomenta );
              //CB[ n ].termInCorr[1].hlist( i, j, k ).mom.resize(params.Nmomenta);
            }
      }

      //read( xml, "NumDistinctQQQ", params.NQQQs );
      reader >> params.NQQQs;

      AQQQ.resize( params.NQQQs );
      for ( int i = 0; i < params.NQQQs; ++i )
      {
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
                CQQQ[ q ].orderings[ p ].hlist( i, j, k ).mom.resize( params.Nmomenta );
              }
        }
        for ( int p = 0; p < params.NsnkOrderings; ++p )
        {
          for ( int i = 0; i < params.NH[ p ][ 0 ]; ++i )
            for ( int j = 0; j < params.NH[ p ][ 1 ]; ++j )
              for ( int k = 0; k < params.NH[ p ][ 2 ]; ++k )
              {
                AQQQ[ q ].orderings[ p ].hlist( i, j, k ).mom.resize( params.Nmomenta );
              }
        }
      }

      for ( int i = 0; i < params.NQQQs; ++i )
      {
        int hash, dL;
        reader >> hash
				       >> AQQQ[ i ].quark[ 0 ].spin
							 >> AQQQ[ i ].quark[ 1 ].spin
							 >> AQQQ[ i ].quark[ 2 ].spin
							 >> AQQQ[ i ].quark[ 0 ].displacement
							 >> AQQQ[ i ].quark[ 1 ].displacement
							 >> AQQQ[ i ].quark[ 2 ].displacement
							 >> dL
							 >> AQQQ[ i ].NBaryonOps;

        AQQQ[ i ].QQQIndex = i;

        for ( int j = 0; j < 3; ++j )
        {
          if ( AQQQ[ i ].quark[ j ].displacement == 0 )
          {
            AQQQ[ i ].quark[ j ].disp_dir = 0;
            AQQQ[ i ].quark[ j ].disp_len = 0;
          }
          else if ( AQQQ[ i ].quark[ j ].displacement > 0 )
          {
            AQQQ[ i ].quark[ j ].disp_dir = AQQQ[ i ].quark[ j ].displacement - 1;
            AQQQ[ i ].quark[ j ].disp_len = dL;
          }
          else
          {
            AQQQ[ i ].quark[ j ].disp_dir = -AQQQ[ i ].quark[ j ].displacement - 1;
            AQQQ[ i ].quark[ j ].disp_len = -dL;
          }
        }
        AQQQ[ i ].coef.resize( AQQQ[ i ].NBaryonOps );
        AQQQ[ i ].whichBaryonOps.resize( AQQQ[ i ].NBaryonOps );
        //             name                                    coefficient
        // AQQQ[ <mu,nu,tau> ].whichBaryonOps[ n ]    ,   AQQQ[ <mu,nu,tau> ].coef[ n ]
        //
        for ( int n = 0; n < AQQQ[ i ].NBaryonOps; ++n )
        {
          Real re, im;
          reader >> AQQQ[ i ].whichBaryonOps[ n ];
          reader >> re >> im;
          AQQQ[ i ].coef[ n ] = cmplx( re, im );
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
          CQQQ[ i ].quark[ j ].disp_dir = AQQQ[ i ].quark[ j ].disp_dir;
          CQQQ[ i ].quark[ j ].spin = AQQQ[ i ].quark[ j ].spin;
          if ( CQQQ[ i ].quark[ j ].spin > 1 )
            spinCount++;
        }
        if ( spinCount % 2 )
          flipsign = 1;

        CQQQ[ i ].NBaryonOps = AQQQ[ i ].NBaryonOps;
        CQQQ[ i ].coef.resize( CQQQ[ i ].NBaryonOps );
        CQQQ[ i ].whichBaryonOps.resize( CQQQ[ i ].NBaryonOps );
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
        }
      } // i : CQQQ[] loop
      reader.close();
    } // void ReadTextInput
		
    //! Construct array of maps of displacements
    void
    GroupBaryonQQQ::displaceQuarks( multi1d< map<int, LatticeFermion> >& disp_quarks,
                                    const multi1d<LatticeFermion>& q,
                                    enum PlusMinus isign ) const
    {
      START_CODE();
      //      cout << __PRETTY_FUNCTION__ << ": entering" << endl;
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
          //	      cout << __func__
          //		   << ": n=" << n
          //		   << " l=" << l
          //		   << " i=" << i
          //		   << " disp=" << term.quark[i].displacement
          //		   << " len=" << term.quark[i].disp_len
          //		   << " dir=" << term.quark[i].disp_dir
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
          disp_q.insert( std::make_pair( term_q.displacement, qq ) );
        }
      } // for i

      //      cout << __PRETTY_FUNCTION__ << ": exiting" << endl;

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

      //      cout << __PRETTY_FUNCTION__ << ": entering" << endl;

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

      //      cout << __PRETTY_FUNCTION__ << ": exiting" << endl;

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

      //      cout << __PRETTY_FUNCTION__ << ": entering" << endl;

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

      //      cout << __PRETTY_FUNCTION__ << ": exiting" << endl;

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

      //      cout << __PRETTY_FUNCTION__ << ": entering" << endl;

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

      //      cout << __PRETTY_FUNCTION__ << ": exiting" << endl;

      END_CODE();
    } // void GroupBaryonQQQ::quarkManip
    
		//! Compute the operator
    multi1d<LatticeComplex>
    GroupBaryonQQQ::operator() ( const LatticeFermion& q1,
                                 const LatticeFermion& q2,
                                 const LatticeFermion& q3,
                                 enum PlusMinus isign ) const
    { START_CODE();
      //      cout << __PRETTY_FUNCTION__ << ": entering" << endl;
      // The result of displace and smearing (in some unspecified order here)
      multi1d< map<int, LatticeFermion> > disp_quarks;
      // Depending on whether this is the sink or source, do the appropriate
      // combination of smearing and displacing
      quarkManip( disp_quarks, q1, q2, q3, isign );
      // The return
      //multi1d<LatticeComplex> d;//( 6 );
      multi1d<LatticeComplex> d(1);
      // Contract over color indices with antisym tensors
      switch ( isign )
      {
        case MINUS:
          // Source			
					//d.resize(6);
					// mu,nu,tau contraction
      		d[ 0 ] = colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 0 ].displacement ) ->second,
      		                                  quark[ 0 ].spin ),
      		                        peekSpin( disp_quarks[ 1 ].find( quark[ 1 ].displacement ) ->second,
      		                                  quark[ 1 ].spin ),
      		                        peekSpin( disp_quarks[ 2 ].find( quark[ 2 ].displacement ) ->second,
      		                                  quark[ 2 ].spin ) );
					// tau,nu,mu contraction
      		//d[ 3 ] = colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 2 ].displacement ) ->second,
      		d[ 0 ] -= colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 2 ].displacement ) ->second,
      		                                   quark[ 2 ].spin ),
      		                         peekSpin( disp_quarks[ 1 ].find( quark[ 1 ].displacement ) ->second,
      		                                   quark[ 1 ].spin ),
      		                         peekSpin( disp_quarks[ 2 ].find( quark[ 0 ].displacement ) ->second,
      		                                   quark[ 0 ].spin ) );
					//d[ 3 ] *= -2.0;
					d[ 0 ] *= 2.0;
					// nu,mu,tau contraction
      		//d[ 1 ] = colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 1 ].displacement ) ->second,
      		d[ 0 ] += colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 1 ].displacement ) ->second,
      		                                   quark[ 1 ].spin ),
      		                         peekSpin( disp_quarks[ 1 ].find( quark[ 0 ].displacement ) ->second,
      		                                   quark[ 0 ].spin ),
      		                         peekSpin( disp_quarks[ 2 ].find( quark[ 2 ].displacement ) ->second,
      		                                   quark[ 2 ].spin ) );
					// mu,tau,nu contraction
      		//d[ 2 ] = colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 0 ].displacement ) ->second,
      		d[ 0 ] += colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 0 ].displacement ) ->second,
      		                                   quark[ 0 ].spin ),
      		                         peekSpin( disp_quarks[ 1 ].find( quark[ 2 ].displacement ) ->second,
      		                                   quark[ 2 ].spin ),
      		                         peekSpin( disp_quarks[ 2 ].find( quark[ 1 ].displacement ) ->second,
      		                                   quark[ 1 ].spin ) );
					// nu,tau,mu contraction
      		//d[ 4 ]= -colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 1 ].displacement ) ->second,
      		d[ 0 ] -= colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 1 ].displacement ) ->second,
      		                                   quark[ 1 ].spin ),
      		                         peekSpin( disp_quarks[ 1 ].find( quark[ 2 ].displacement ) ->second,
      		                                   quark[ 2 ].spin ),
      		                         peekSpin( disp_quarks[ 2 ].find( quark[ 0 ].displacement ) ->second,
      		                                   quark[ 0 ].spin ) );
					// tau,mu,nu contraction
      		//d[ 5 ]= -colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 2 ].displacement ) ->second,
      		d[ 0 ] -= colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 2 ].displacement ) ->second,
      		                                   quark[ 2 ].spin ),
      		                         peekSpin( disp_quarks[ 1 ].find( quark[ 0 ].displacement ) ->second,
      		                                   quark[ 0 ].spin ),
      		                         peekSpin( disp_quarks[ 2 ].find( quark[ 1 ].displacement ) ->second,
      		                                   quark[ 1 ].spin ) );
      		break;

        case PLUS:
          // Sink
					//d.resize(1);
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
      //      cout << __PRETTY_FUNCTION__ << ": exiting" << endl;
      return d;
    } // multi1d<LatticeComplex> GroupBaryonQQQ::operator()



		//	====================
		//	GroupBaryonOp stuff
		//	====================
		
    // most of this is obsolete ... see GroupBaryonQQQ now
#if 1
		//! Full constructor
    GroupBaryonOp::GroupBaryonOp( const Params& p, const multi1d<LatticeColorMatrix>& u_ ) :
        params( p ), u_smr( u_ )
    {
      //readCoeffs( coeffs );

      // The spin basis matrix to goto Dirac
      spin_rotate_mat = adj( DiracToDRMat() );

      // Factory constructions
      try
      {
        // Smear the gauge field if needed
        {
          std::istringstream xml_l( params.link_smearing.xml );
          XMLReader linktop( xml_l );
          const string link_path = "/LinkSmearing";
          QDPIO::cout << "Link smearing type = " << params.link_smearing.id << endl;

          Handle< LinkSmearing >
          linkSmearing( TheLinkSmearingFactory::Instance().createObject( params.link_smearing.id,
                        linktop,
                        link_path ) );
          ( *linkSmearing ) ( u_smr );
        }

        // Create the source quark smearing object
        {
          std::istringstream xml_s( params.source_quark_smearing.xml );
          XMLReader smeartop( xml_s );
          const string smear_path = "/SourceQuarkSmearing";

          sourceQuarkSmearing =
            TheFermSmearingFactory::Instance().createObject( params.source_quark_smearing.id,
                smeartop,
                smear_path );
        }

        // Create the sink quark smearing object
        {
          std::istringstream xml_s( params.sink_quark_smearing.xml );
          XMLReader smeartop( xml_s );
          const string smear_path = "/SinkQuarkSmearing";

          sinkQuarkSmearing =
            TheFermSmearingFactory::Instance().createObject( params.sink_quark_smearing.id,
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

      //      cout << __PRETTY_FUNCTION__ << ": entering" << endl;

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

      //      cout << __PRETTY_FUNCTION__ << ": exiting" << endl;

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

      //      cout << __PRETTY_FUNCTION__ << ": entering" << endl;

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

      //      cout << __PRETTY_FUNCTION__ << ": exiting" << endl;

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

      //      cout << __PRETTY_FUNCTION__ << ": entering" << endl;

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

      //      cout << __PRETTY_FUNCTION__ << ": exiting" << endl;

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

      //      cout << __PRETTY_FUNCTION__ << ": entering" << endl;

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

      //      cout << __PRETTY_FUNCTION__ << ": exiting" << endl;

      END_CODE();
    } // end void GroupBaryonOp::quarkManip
		
		//! Compute the operator
    multi1d<LatticeComplex>
    GroupBaryonOp::operator() ( const LatticeFermion& q1,
                                const LatticeFermion& q2,
                                const LatticeFermion& q3,
                                enum PlusMinus isign ) const
    { START_CODE();
      //      cout << __PRETTY_FUNCTION__ << ": entering" << endl;
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
      //      cout << __PRETTY_FUNCTION__ << ": exiting" << endl;
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
  	  int mom_size = termInCorr[ 0 ].hlist( 0, 0, 0 ).mom.size();
//  	  int elem_size2 = termInCorr[ 0 ].hlist( 0, 0, 0 ).mom[ 0 ].elem.size2();
//  	  int elem_size1 = termInCorr[ 0 ].hlist( 0, 0, 0 ).mom[ 0 ].elem.size1();
  	  // dreadful hack - use a complex to hold an int
  	  Complex termInCorr_sizes, hlist_sizes1, hlist_sizes2;//, elem_sizes;
  	  termInCorr_sizes = cmplx( Real( termInCorr_size ), Real( zero ) );
  	  hlist_sizes1 = cmplx( Real( hlist_size2 ), Real( hlist_size1 ) );
  	  hlist_sizes2 = cmplx( Real( hlist_size3 ), Real( mom_size ) );
//  	  elem_sizes = cmplx( Real( elem_size2 ), Real( elem_size1 ) );

      multi1d<Complex> baryonOp_1d( 4 + termInCorr_size * hlist_size3 * hlist_size2 * hlist_size1 * mom_size );//* elem_size2 * elem_size1 );
  	  //    QDPIO::cout << "baryonprop_size=" << baryonOp_1d.size() << endl;
  	  int cnt = 0;
  	  baryonOp_1d[ cnt++ ] = termInCorr_sizes;
  	  baryonOp_1d[ cnt++ ] = hlist_sizes1;
  	  baryonOp_1d[ cnt++ ] = hlist_sizes2;
//  	  barprop_1d[ cnt++ ] = elem_sizes;
  	  for ( int s = 0; s < termInCorr.size(); ++s ) // termInCorr
  	  {
  	    for ( int i = 0; i < termInCorr[ s ].hlist.size3(); ++i )      // hlist_l
  	      for ( int j = 0; j < termInCorr[ s ].hlist.size2(); ++j )    // hlist_m
  	        for ( int k = 0; k < termInCorr[ s ].hlist.size1(); ++k )  // hlist_r
  	          for ( int p = 0; p < termInCorr[ s ].hlist( i, j, k ).mom.size(); ++p )  // mom
  	            //for ( int a = 0; a < termInCorr[ s ].hlist( i, j, k ).mom[ l ].elem.size2(); ++a )    // elem_l
  	              //for ( int b = 0; b < termInCorr[ s ].hlist( i, j, k ).mom[ l ].elem.size1(); ++b )  // elem_r
  	                baryonOp_1d[ cnt++ ] = termInCorr[ s ].hlist( i, j, k ).mom[ p ];//.elem( a, b );
  	  }
  	  if ( cnt != baryonOp_1d.size() )
  	  {
  	    QDPIO::cerr << GroupBaryonOperatorEnv::name << ": size mismatch in serialization" << endl;
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

        //! Register all the factories
        success &= Chroma::TheWilsonBaryonOperatorFactory::Instance().registerObject( name,
                   groupBaryon );

        registered = true;
      }
      return success;
    } // registerAll()


  } // namespace GroupBaryonOperatorEnv


}  // namespace Chroma
