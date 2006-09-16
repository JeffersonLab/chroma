// $Id: group_baryon_operator_w.cc,v 1.7 2006-09-16 05:20:35 juge Exp $
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
      source_quark_smearing = readXMLGroup( paramtop, "SourceQuarkSmearing", "wvf_kind" );
      sink_quark_smearing = readXMLGroup( paramtop, "SinkQuarkSmearing", "wvf_kind" );
      link_smearing = readXMLGroup( paramtop, "LinkSmearing", "LinkSmearingType" );

      read( paramtop, "main_input_file", mainInputFile );
      //read( paramtop, "operator_names_file", operatorNamesFile );
      
			//read( paramtop, "operator_coeff_file", operator_coeff_file );
      //read( paramtop, "displacement_length", displacement_length );
    }
    // Writer
    void Params::writeXML( XMLWriter& xml, const string& path ) const
    {
      push( xml, path );
      int version = 1;
      write( xml, "version", version );
      write( xml, "BaryonOperatorType", GroupBaryonOperatorEnv::name );
      write( xml, "main_input_file", mainInputFile );
      //write( xml, "operator_names_file", operatorNamesFile );
      //write( xml, "operator_coeff_file", operator_coeff_file );
      //write( xml, "displacement_length", displacement_length );
      xml << source_quark_smearing.xml;
      xml << sink_quark_smearing.xml;
      xml << link_smearing.xml;
      pop( xml );
    }
    //! Full constructor
    GroupBaryon::GroupBaryon( const Params& p, const multi1d<LatticeColorMatrix>& u_ ) :
        params( p ), u_smr( u_ )
    {
      readCoeffs( coeffs );
      // The spin basis matrix to goto Dirac
      rotate_mat = adj( DiracToDRMat() );
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
    }
		// obsolete
    //! Reader
    void GroupBaryon::readCoeffs( multi1d< multi1d< CoeffTerm_t > >& coef )
    {
      TextReader reader( params.operator_coeff_file );
      int num_rows;
      reader >> num_rows;
      coeffs.resize( num_rows );
      for ( int n = 0; n < coeffs.size(); ++n )
      {
        int num_terms;
        reader >> num_terms;
        coeffs[ n ].resize( num_terms );
        for ( int l = 0; l < coeffs[ n ].size(); ++l )
        {
          CoeffTerm_t& term = coeffs[ n ][ l ];
          term.quark.resize( 3 );
          // Make spin index 0 based
          {
            multi1d<int> spin( 3 );
            reader >> spin[ 0 ] >> spin[ 1 ] >> spin[ 2 ];
            for ( int i = 0; i < spin.size(); ++i )
              term.quark[ i ].spin = spin[ i ] - 1;
          }
          // Convert displacements to  disp_dir, disp_len
          {
            multi1d<int> displacement( 3 );
            reader >> displacement[ 0 ] >> displacement[ 1 ] >> displacement[ 2 ];
            for ( int i = 0; i < displacement.size(); ++i )
            {
              term.quark[ i ].displacement = displacement[ i ];
              if ( displacement[ i ] == 0 )
              {
                term.quark[ i ].disp_dir = 0;
                term.quark[ i ].disp_len = 0;
              }
              else if ( displacement[ i ] > 0 )
              {
                term.quark[ i ].disp_dir = term.quark[ i ].displacement - 1;
                term.quark[ i ].disp_len = params.displacement_length;
              }
              else
              {
                term.quark[ i ].disp_dir = -term.quark[ i ].displacement - 1;
                term.quark[ i ].disp_len = -params.displacement_length;
              }
            }
          }
          // Read the garbage around a complex
          {
            Real re, im;
            char lparen, comma, rparen;
            reader >> lparen >> re >> comma >> im >> rparen;
            term.coeff = cmplx( re, im );
          }
        }
      }
      reader.close();
    }
		
		Params::initialize()
		{
			// hybrid list sizes
			NH.resize( 4 ); // [012] [201] [210] [102]
			for(int i=0; i<4; ++i) NH[ i ].resize(3); // 0,1,2
			int term1 = 0;
			int term2 = 1;
			Aorderings.resize(2); // two terms				
			Aorderings[ term1 ].resize(2); // two orderings
			Aorderings[ term2 ].resize(2); // two orderings
			Asign.resize(2); // two terms				
			Asign[ term1 ].resize(2); // two orderings
			Asign[ term2 ].resize(2); // two orderings
			Corderings.resize(2); // two terms				
			Corderings[ term1 ].resize(2); // two orderings
			Corderings[ term2 ].resize(2); // two orderings
			Csign.resize(2); // two terms				
			Csign[ term1 ].resize(2); // two orderings
			Csign[ term2 ].resize(2); // two orderings
			//
			// Signs are hard-wired here ... CHECK THESE
			//
			// Annihilation Operators
			//
			// first term
			Aorderings[ term1 ][ 0 ] = 0; // [012]
			Aorderings[ term1 ][ 1 ] = 1; // [201]
			Asign[ term1 ][ 0 ] =  1.0; // [012] +
			Asign[ term1 ][ 1 ] = -1.0; // [201] -
			// second term
			Aorderings[ term2 ][ 0 ] = 0; // [012]
			Aorderings[ term2 ][ 1 ] = 2; // [210]
			Asign[ term2 ][ 0 ] = -1.0; // [012] -
			Asign[ term2 ][ 1 ] = -1.0; // [210] -
			//
			// Creation Operators
			//
			// first term
			Corderings[ term1 ][ 0 ] = 0; // [012]
			Corderings[ term1 ][ 1 ] = 3; // [102]
			Csign[ term1 ][ 0 ] = 1.0; // [012] +
			Csign[ term1 ][ 1 ] = 1.0; // [102] +
			// second term
			Corderings[ term2 ][ 0 ] = 1; // [201]
			Corderings[ term2 ][ 1 ] = 2; // [210]
			Csign[ term2 ][ 0 ] =  1.0; // [201] +
			Csign[ term2 ][ 1 ] = -1.0; // [210] -
		}

    //! Reader
    void GroupBaryon::newReadInput( Params& params )
    {
			QDPIO::cout<<"Reading input from text file : "<<endl;			
      TextReader reader( params.InputFileName );
			
			params.initialize();			
			multi1d<BaryonOp_t> AB; // annihilation
			multi1d<BaryonOp_t> CB; // creation
			multi1d<QQQOperator_t> AQQQ; // annihilation
			multi1d<QQQOperator_t> CQQQ; // creation
			// Set the hybrid list sizes
			//   { something like this would be better                  }
			//   { Ndil = params.named_obj.prop.op[n].soln_files.size() }
		reader >> params.NdilL >> params.NdilM >> params.NdilR;
			
		reader >> params.NsrcOrderings;

      params.SrcOrderingsresize(3);
      for(int i=0; i < params.NsrcOrderings; ++i)
			{
      	params.SrcOrderings[i].resize(3);
		reader >> params.SrcOrderings[i][0] >> params.SrcOrderings[i][1] >> params.SrcOrderings[i][2];
			} 
			// 012                              201
			params.NH[ 0 ][ 0 ] = params.NdilL; params.NH[ 1 ][ 0 ] = params.NdilR;
			params.NH[ 0 ][ 1 ] = params.NdilM; params.NH[ 1 ][ 1 ] = params.NdilL;
			params.NH[ 0 ][ 2 ] = params.NdilR; params.NH[ 1 ][ 2 ] = params.NdilM;
			// 210                              102
			params.NH[ 2 ][ 0 ] = params.NdilR; params.NH[ 3 ][ 0 ] = params.NdilM;
			params.NH[ 2 ][ 1 ] = params.NdilM; params.NH[ 3 ][ 1 ] = params.NdilL;
			params.NH[ 2 ][ 2 ] = params.NdilL; params.NH[ 3 ][ 2 ] = params.NdilR;

		reader >> params.NsnkOrderings;

      for(int i=0; i < params.NsnkOrderings; ++i)
			{
      	params.SnkOrderings[i].resize(3);
		reader >> params.SnkOrderings[i][0] >> params.SnkOrderings[i][1] >> params.SnkOrderings[i][2];
			}
			
		reader >> params.Noperators;
						
			AB.resize( params.Noperators );
			CB.resize( params.Noperators );
      for(int n=0; n < params.Noperators; ++n)
			{
				AB[ n ].termInCorr.resize( 2 );
				CB[ n ].termInCorr.resize( 2 );
      }
			for(int n=0; n < params.Noperators; ++n)
			{
				AB[ n ].termInCorr[0].hlist.resize( params.NdilL, params.NdilM, params.NdilR );
				AB[ n ].termInCorr[1].hlist.resize( params.NdilL, params.NdilM, params.NdilR );
				CB[ n ].termInCorr[0].hlist.resize( params.NdilL, params.NdilM, params.NdilR );
				CB[ n ].termInCorr[1].hlist.resize( params.NdilL, params.NdilM, params.NdilR );
			}
			for(int n=0; n < params.Noperators; ++n)
			{
				for(int i=0; i < params.NdilL; ++i)
				for(int j=0; j < params.NdilM; ++j)
				for(int k=0; k < params.NdilR; ++k)
				{
					AB[ n ].termInCorr[0].hlist( i, j, k ).resize(params.Nmomenta);
					AB[ n ].termInCorr[1].hlist( i, j, k ).resize(params.Nmomenta);
					CB[ n ].termInCorr[0].hlist( i, j, k ).resize(params.Nmomenta);
					CB[ n ].termInCorr[1].hlist( i, j, k ).resize(params.Nmomenta);
				}
			}

		reader >> params.NQQQs;

			AQQQ.resize( params.NQQQs );
			for(int i=0; i < params.NQQQs; ++i) 
			{
				AQQQ[i].initialize();
				AQQQ[i].quark.resize(3);
				AQQQ[i].orderings.resize(params.NsnkOrderings);
				for(int j=0; j < params.NsnkOrderings; ++j) 
				{
					AQQQ[i].orderings[j].resize(params.Ndil[j][0],params.Ndil[j][1],params.Ndil[j][2]);
				}
			}
			CQQQ.resize( params.NQQQs );
			for(int i=0; i < params.NQQQs; ++i) 
			{
				CQQQ[i].initialize();
				CQQQ[i].quark.resize(3);
				CQQQ[i].orderings.resize(params.NsrcOrderings);
				for(int j=0; j < params.NsrcOrderings; ++j) 
				{
					CQQQ[i].orderings[j].resize(params.Ndil[j][0],params.Ndil[j][1],params.Ndil[j][2]);
				}
			}
			for(int q=0; q < params.NQQQs; ++q) 
			{
				for(int p=0; p < params.NsrcOrderings; ++p)
				{ 
					for(int i=0; i < params.Ndil[p][0]; ++i) 
					for(int j=0; j < params.Ndil[p][1]; ++j) 
					for(int k=0; k < params.Ndil[p][2]; ++k) 
					{
						CQQQ[q].orderings[p].hlist(i,j,k).resize(params.Nmomenta);
					}
				}
				for(int p=0; p < params.NsnkOrderings; ++p)
				{ 
					for(int i=0; i < params.Ndil[p][0]; ++i) 
					for(int j=0; j < params.Ndil[p][1]; ++j) 
					for(int k=0; k < params.Ndil[p][2]; ++k) 
					{
						AQQQ[q].orderings[p].hlist(i,j,k).resize(params.Nmomenta);
					}
				}
			}

      for(int i=0; i < params.NQQQs; ++i)
			{
				int hash, dL;
				reader >> hash 
				       >> AQQQ[i].quark[ 0 ].spin 
							 >> AQQQ[i].quark[ 1 ].spin 
							 >> AQQQ[i].quark[ 2 ].spin
               >> AQQQ[i].quark[ 0 ].displacement
							 >> AQQQ[i].quark[ 1 ].displacement
							 >> AQQQ[i].quark[ 2 ].displacement 
							 >> dL 
							 >> AQQQ[i].NBaryonOps;
							         
				for(int j=0; j < 3; ++j)
        {
          if ( disp[ j ] == 0 )
          {
            AQQQ[i].quark[j].disp_dir = 0;
            AQQQ[i].quark[j].disp_len = 0;
          }
          else if ( disp[ j ] > 0 )
          {
            AQQQ[i].quark[ j ].disp_dir = AQQQ[i].quark[ j ].displacement - 1;
            AQQQ[i].quark[ j ].disp_len = dL;
          }
          else
          {
            AQQQ[i].quark[ j ].disp_dir = -AQQQ[i].quark[ j ].displacement - 1;
            AQQQ[i].quark[ j ].disp_len = -dL;
          }
        }
				AQQQ[i].coef.resize(AQQQ[i].NBaryonOps);
				AQQQ[i].whichBaryonOp.resize(AQQQ[i].NBaryonOps);
				//             name                                    coefficient
				// AQQQ[ <mu,nu,tau> ].whichBaryonOp[ n ]    ,   AQQQ[ <mu,nu,tau> ].coef[ n ]
				//
				for(int n=0; n < AQQQ[i].NBaryonOps; ++n)
				{
					Real re, im;
					reader >> AQQQ[i].whichBaryonOp[ n ];
					reader >> re >> im;
					AQQQ[i].coef[ n ] = cmplx(re,im);
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
      for(int i=0; i < params.NQQQs; ++i)
			{
				int spinCount, flipsign;
				spinCount = 0; flipsign = 0;
				for(int j=0; j < 3; ++j)
        {
          CQQQ[i].quark[ j ].displacement = AQQQ[i].quark[ j ].displacement;
          CQQQ[i].quark[ j ].disp_dir = AQQQ[i].quark[ j ].disp_dir;
          for(int k=0; k < 4; ++k) CQQQ[i].quark[ j ].spin[ k ] = AQQQ[i].quark[ j ].spin[ k ];
					if( CQQQ[i].quark[ j ].spin[ 2 ] > 1 ) spinCount++;
					if( CQQQ[i].quark[ j ].spin[ 3 ] > 1 ) spinCount++;
        }
				if ( spinCount%2 ) flipsign = 1;
				
				CQQQ[i].NBaryonOps = AQQQ[i].NBaryonOps;
				CQQQ[i].coef.resize(CQQQ[i].NBaryonOps);
				CQQQ[i].whichBaryonOp.resize(CQQQ[i].NBaryonOps);
				
				//             name                           coefficient
				// CQQQ[ <mu,nu,tau> ].whichBaryonOp[ n ] , CQQQ[ i ].coef[ n ]
				//
				for(int n=0; n < CQQQ[i].NBaryonOps; ++n)
				{
					Real re, im;
					CQQQ[i].whichBaryonOp[ n ] = AQQQ[i].whichBaryonOp[ n ];
					re = real( AQQQ[i].coef[ n ] );
					im = imag( AQQQ[i].coef[ n ] );
					if( flipsign ) 
					{
						// cc and -sign
						CQQQ[i].coef[ n ] = cmplx(-re,im);
					}
					else
					{
						// just cc
						CQQQ[i].coef[ n ] = cmplx(re,-im);
					}					
				}
			} // i : CQQQ[] loop
			
			// The Baryon operator names ; G1g_L3_TDT_25 ...
			for(int i=0; i < params.Names.size(); ++i)
			{
				reader >> params.Names[ i ];
			}
      reader.close();
    }
		//
		int DilSwap( int ord, int i, int j, int k, int which )   
		{ // Used for hybrid list indices, NH and noise indices  
			//  	for indices:  ( ord, i,   j,   k,   which )				   
			//  	for NHlimit:  ( ord, NH0, NH1, NH2, which )				   
			//  	for NoiseInd: ( ord, 0,   1,   2,   which )				   
		  multi1d<int> ijk( 3 );															   
		  switch ( ord )																			   
		  { 																									   
		    case 0:  // [012] ijk 														   
		      ijk[ 0 ] = i; 																	   
		      ijk[ 1 ] = j; 																	   
		      ijk[ 2 ] = k; 																	   
		      break;																					   
		    case 1:  // [201] kij 														   
		      ijk[ 0 ] = k; 																	   
		      ijk[ 1 ] = i; 																	   
		      ijk[ 2 ] = j; 																	   
		      break;																					   
		    case 2:  // [210] kji 														   
		      ijk[ 0 ] = k; 																	   
		      ijk[ 1 ] = j; 																	   
		      ijk[ 2 ] = i; 																	   
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
		    case 5:  // [021] ikj 														   
		      ijk[ 0 ] = i; 																	   
		      ijk[ 1 ] = k; 																	   
		      ijk[ 2 ] = j; 																	   
		      break;																					   
		  } 																									   
		  return ( ijk[ which ] );														   
		} 																										   

	  //! Serialize generalized operator object
	  multi1d<Complex> BaryonOp_t::serialize()
	  {
	    int op_size3 = hlist.size3();
	    int op_size2 = hlist.size2();
	    int op_size1 = hlist.size1();
	    int mom_size = hlist( 0, 0, 0 ).mom.size();
	    // dreadful hack - use a complex to hold an int
	    Complex op_sizes1, op_sizes2;
	    op_sizes1 = cmplx( Real( op_size1 ), Real( op_size2 ) );
	    op_sizes2 = cmplx( Real( op_size3 ), Real( mom_size ) );

	    multi1d<Complex> baryonOp_1d( 2 + op_size3 * op_size2 * op_size1 * mom_size );

	    //    QDPIO::cout << "baryprop_size=" << baryprop_1d.size() << endl;

	    int cnt = 0;
	    baryonOp_1d[ cnt++ ] = op_sizes1;
	    baryonOp_1d[ cnt++ ] = op_sizes2;

	    for(int i=0; i < op.size3(); ++i)  // op_l
	    {
	      for(int j=0; j < op.size2(); ++j)  // op_m
	      {
	        for(int k=0; k < op.size1(); ++k)  // op_r
	        {
	          for(int p=0; p < hlist( i, j, k ).mom.size1(); ++p)  // elem_r mom number
	          {
	            baryonOp_1d[ cnt++ ] = hlist( i, j, k ).mom( p );
	          }
	        }
	      }
	    }
	    if ( cnt != baryonOp_1d.size() )
	    {
	      QDPIO::cerr << InlineStochGroupBaryonEnv::name << ": size mismatch in serialization" << endl;
	      QDP_abort( 1 );
	    }
	    return baryonOp_1d;
	  }
		/*  how the QQQ contraction part will look like (spin mu = m, spin nu = n, spin tau = t )
        //
				//  (i)   (j)   (k)   Example: ordering=3     (j)*   (i)*   (k)*
				// q   * q   * q                           eta    eta    eta     saved in hlist(j,i,k)  
				//   L    M      R                            [1]    [0]    [2]
				//                                     _ijk
				// combine with [012] to make operator O     (eg one contribution from <mu,nu,tau>)
				//                                      (1)
				//  _ijk                                   
				//  O    += c   * Csign[ term1 ][ ord0 ] * qqq_op[mnt].orderings[ ord0 ].hlist(i1,j1,k1) 
				//   (1)     mnt
				//  _ijk                                   
				//  O    += c   * Csign[ term1 ][ ord3 ] * qqq_op[mnt].orderings[ ord3 ].hlist(i2,j2,k2) 
				//   (1)     mnt
				//  
				// where
				//        i1 = DilSwap( ord0, i, j, k, 0 ) = i     i2 = DilSwap( ord3, i, j, k, 0 ) = j
				//        j1 = DilSwap( ord0, i, j, k, 1 ) = j     j2 = DilSwap( ord3, i, j, k, 1 ) = i
				//        k1 = DilSwap( ord0, i, j, k, 2 ) = k     k2 = DilSwap( ord3, i, j, k, 2 ) = k
				//
        
				for(srcOrd=0; srcOrd < NsrcOrderings=4; ++srcOrd)
				
				for(term=0; term < 2; ++term)
				
				const QuarkSourceSolutions_t& q0 = quarks[ COrderings[ term1 ][ 0 ] ]; 
        const QuarkSourceSolutions_t& q1 = quarks[ COrderings[ term1 ][ 1 ] ];
        const QuarkSourceSolutions_t& q2 = quarks[ COrderings[ term1 ][ 2 ] ];

        qqq_op[mnt].orderings[ srcord ].hlist.resize( q0.dilutions.size(),
                                                      q1.dilutions.size(),
                                                      q2.dilutions.size() );
				
				for ( int i = 0; i < q0.dilutions.size(); ++i ) 
        {
          for ( int j = 0; j < q1.dilutions.size(); ++j )
          {
            for ( int k = 0; k < q2.dilutions.size(); ++k )
            {
              LatticeComplex bar = ( *baryonOperator ) ( q0.dilutions[ i ].source,
                                                         q1.dilutions[ j ].source,
                                                         q2.dilutions[ k ].source,
                                                         MINUS );
              qqq_op[mnt].orderings[ srcord ].hlist( i, j, k ).elem = phases.sft( bar );
            } // end for k
          } // end for j
        } // end for i
		*/
  	void QQQOperator_t::sum( multi1d<multi1d<BaryonOp_t> >& B, 
  	                         int NdilL, int NdilM, int NdilR )
  	{
  	  DComplex dctemp, dctemp2;
			for(int Bops=0; Bops < NBaryonOps; ++Bops)
  	  {
  	    for(int term=0; term < 2; ++term) // hard-wired
  	    {
  	      for(int ord=0; ord < 2; ++ord) // hard-wired
  	      {
  	        // this should probably be done at input level ...
						dctemp.re = ( B[ whichBaryonOps[Bops] ][ term ].sign[ord] > 0 )
						           ?  B[ whichBaryonOps[Bops] ][ term ].coeff[ QQQIndex ].re 
											 : -B[ whichBaryonOps[Bops] ][ term ].coeff[ QQQIndex ].re;
  	        dctemp.im = ( B[ whichBaryonOps[Bops] ][ term ].sign[ord] > 0 )
						           ?  B[ whichBaryonOps[Bops] ][ term ].coeff[ QQQIndex ].im 
											 : -B[ whichBaryonOps[Bops] ][ term ].coeff[ QQQIndex ].im;
  	        for(int i=0; i < NdilL; ++i)
  	        for(int j=0; j < NdilM; ++j)
  	        for(int k=0; k < NdilR; ++k)
  	        {
  	          int ii = DilSwap( ord , i, j, k, 0 );
  	          int jj = DilSwap( ord , i, j, k, 1 );
  	          int kk = DilSwap( ord , i, j, k, 2 );
  	          for(int p=0; p < param.Nmomenta; ++p)
  	          {
  	            dctemp2 = orderings[ord].hlist(ii,jj,kk).elem[ p ];
  	            B[ whichBaryonOps[Bops] ][ term ].termInCorr[0].hlist(i,j,k).mom[p].re 
								 += ( dctemp.re * dctemp2.re - dctemp.im * dctemp2.im );
  	            B[ whichBaryonOps[Bops] ][ term ].termInCorr[0].hlist(i,j,k).mom[p].im 
								 += ( dctemp.re * dctemp2.im + dctemp.im * dctemp2.re );
  	          }
  	        }
  	      }
  	    }
  	  }
  	}
		
    //! Compute the operator
    //multi1d<LatticeComplex>
    LatticeComplex
    QQQOperator_t::operator() ( const LatticeFermion& q1,
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
      LatticeComplex d;
      //for ( int n = 0; n < NBaryonOps; ++n )
      //{
        // Contract over color indices with antisym tensors
        d = colorContract( peekSpin( disp_quarks[ 0 ].find( quark[ 0 ].displacement ) ->second,
                                     quark[ 0 ].spin ),
                           peekSpin( disp_quarks[ 1 ].find( quark[ 1 ].displacement ) ->second,
                                     quark[ 1 ].spin ),
                           peekSpin( disp_quarks[ 2 ].find( quark[ 2 ].displacement ) ->second,
                                     quark[ 2 ].spin ) );
      //}
      END_CODE();
      //      cout << __PRETTY_FUNCTION__ << ": exiting" << endl;
      return d;
    }
    //! Construct array of maps of displacements
    void
    QQQOperator_t::displaceQuarks( multi1d< map<int, LatticeFermion> >& disp_quarks,
                                   const multi1d<LatticeFermion>& q,
                                   enum PlusMinus isign ) const
    {
      START_CODE();
      //      cout << __PRETTY_FUNCTION__ << ": entering" << endl;
      disp_quarks.resize( 3 );
			// quark[ 0, 1, 2 ]
      for ( int i = 0; i < disp_quarks.size(); ++i )
      {
        // Make some shorthands to ease my brain
        map<int, LatticeFermion>& disp_q = disp_quarks[ i ];
        //const CoeffTerm_t::QuarkTerm_t& term_q = term.quark[ i ];
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
    }

    //! Construct array of maps of displacements
    void
    GroupBaryon::displaceQuarks( multi1d< map<int, LatticeFermion> >& disp_quarks,
                                 const multi1d<LatticeFermion>& q,
                                 enum PlusMinus isign ) const
    {
      START_CODE();

      //      cout << __PRETTY_FUNCTION__ << ": entering" << endl;

      disp_quarks.resize( 3 );

      for ( int n = 0; n < coeffs.size(); ++n )
      {
        for ( int l = 0; l < coeffs[ n ].size(); ++l )
        {
          const CoeffTerm_t& term = coeffs[ n ][ l ];

          for ( int i = 0; i < disp_quarks.size(); ++i )
          {
            // Make some shorthands to ease my brain
            map<int, LatticeFermion>& disp_q = disp_quarks[ i ];
            const CoeffTerm_t::QuarkTerm_t& term_q = term.quark[ i ];

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
        } // for l
      }  // for n

      //      cout << __PRETTY_FUNCTION__ << ": exiting" << endl;

      END_CODE();
    }


    //! First displace then smear the quarks
    void
    GroupBaryon::displaceSmearQuarks( multi1d< map<int, LatticeFermion> >& disp_quarks,
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
    }


    //! First smear then displace the quarks
    void
    GroupBaryon::smearDisplaceQuarks( multi1d< map<int, LatticeFermion> >& disp_quarks,
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
    }

    //! Manipulate the quark fields
    void
    GroupBaryon::quarkManip( multi1d< map<int, LatticeFermion> >& disp_quarks,
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
    }

    //! Compute the operator
    multi1d<LatticeComplex>
    GroupBaryon::operator() ( const LatticeFermion& q1,
                              const LatticeFermion& q2,
                              const LatticeFermion& q3,
                              enum PlusMinus isign ) const
    {
      START_CODE();

      //      cout << __PRETTY_FUNCTION__ << ": entering" << endl;

      // The result of displace and smearing (in some unspecified order here)
      multi1d< map<int, LatticeFermion> > disp_quarks;

      // Depending on whether this is the sink or source, do the appropriate
      // combination of smearing and displacing
      quarkManip( disp_quarks, q1, q2, q3, isign );

      // The return
      multi1d<LatticeComplex> d( coeffs.size() );

      for ( int n = 0; n < coeffs.size(); ++n )
      {
        d[ n ] = zero;

        for ( int l = 0; l < coeffs[ n ].size(); ++l )
        {
          const CoeffTerm_t& term = coeffs[ n ][ l ];

          // Contract over color indices with antisym tensors
          LatticeComplex b_oper =
            colorContract( peekSpin( disp_quarks[ 0 ].find( term.quark[ 0 ].displacement ) ->second,
                                     term.quark[ 0 ].spin ),
                           peekSpin( disp_quarks[ 1 ].find( term.quark[ 1 ].displacement ) ->second,
                                     term.quark[ 1 ].spin ),
                           peekSpin( disp_quarks[ 2 ].find( term.quark[ 2 ].displacement ) ->second,
                                     term.quark[ 2 ].spin ) );

          d[ n ] += term.coeff * b_oper;
        }
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
        return new GroupBaryon( Params( xml_in, path ), u );
      }

    }  // end anonymous namespace


    //! Baryon operators
    /*! \ingroup hadron */
    bool registerAll( void )
    {
      bool success = true;

      // Required stuff
      success &= LinkSmearingEnv::registered;
      success &= QuarkSmearingEnv::registered;

      //! Register all the factories
      success &= Chroma::TheWilsonBaryonOperatorFactory::Instance().registerObject( name,
                 groupBaryon );

      return success;
    }

    const std::string name = "GROUP_BARYON";

    const bool registered = registerAll();

  } // namespace BaryonOperatorCallMapEnv


}  // end namespace Chroma
