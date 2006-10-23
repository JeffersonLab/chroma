/*
 *  \brief Group Theory based Baryon Operator Construction with All2All
 */
#include "chroma.h"
#include "chromabase.h"
#include "meas/hadron/group_baryon_operator_w.h"
using namespace Chroma;
using namespace GroupBaryonOperatorEnv;
//
//! Make a group theory based baryon operator with all-to-all quark 
//! propagators using the dilution method of TrinLat (2004)
//
/*! \defgroup GroupBaryonOperator
 *  \ingroup main
 *
 * Main program for making a group baryon operator
 */
int main( int argc, char *argv[] )
{
	// ==================================
  // Put the machine into a known state
	// ==================================
  Chroma::initialize( &argc, &argv );
  START_CODE();
  QDPIO::cout << "GroupBaryonAll2All: Group Baryon Operators with All2All" << endl;
	// ===============
  // Read input data
	// ===============
  // Instantiate xml reader
  XMLReader xml_in( Chroma::getXMLInputFileName() );
	const std::string path = "/GroupBaryonAll2All";
  Params params( xml_in, path );
	// ==================
  // Setup lattice size
	// ==================
  Layout::setLattSize( params.nrow );
  Layout::create();
	// ============
  // XML business
	// ============
	XMLFileWriter& xml_out = Chroma::getXMLOutputInstance(); // XMLDAT by default
  push( xml_out, "GroupBaryonAll2AllRun" );
  	proginfo( xml_out ); // Print out basic program info
  	write( xml_out, "InputCopy", xml_in ); // save a copy of the input
  	xml_out.flush();
	// =============
  // Startup gauge
	// =============
  XMLReader gauge_file_xml, gauge_xml;
  gaugeStartup(gauge_file_xml, gauge_xml, params.gaugestuff.u, params.gaugestuff.cfg);
		// test the config
		Double w_plaq, s_plaq, t_plaq, link;
  	MesPlq(params.gaugestuff.u, w_plaq, s_plaq, t_plaq, link);
		QDPIO::cout<< "gauge meas " << w_plaq <<" "<< s_plaq <<" "<< t_plaq <<" "<< link <<" "<<endl;	
	// ====================
	// qqq's and baryonop's
	// ====================
	multi1d<GroupBaryonQQQ> AQQQ; // annihilation
	multi1d<GroupBaryonQQQ> CQQQ; // creation
	multi1d<GroupBaryonOp>  AOp;
	multi1d<GroupBaryonOp>  COp;
	// ====================================
	// Read in the text file produced by  
	//   Gen_Input_For_All2All_Baryons2.pl
	// called "Chroma_Input.txt"
	// ====================================
	ReadTextInput( params, AOp, COp, AQQQ, CQQQ );
	// ===============
	// Noise Solutions
	// ===============
	/*
	  struct PropSourceConst_t              struct QuarkSourceSolutions_t {
		{ PropSourceConst_t();            		  struct QuarkSolution_t {
		  multi1d<int> getTSrce() const;  		    LatticeFermion source;
		  GroupXML_t       source;        		    LatticeFermion soln;
		  int              j_decay;       		    PropSourceConst_t source_header;
		  int              t_source;      		    ChromaProp_t      prop_header;
		};																		  };
		struct ChromaProp_t {           			  int j_decay;
		  ChromaProp_t();               			  Seed seed;
		  QuarkSpinType   quarkSpinType;			  multi1d<QuarkSolution_t> dilutions;
		  GroupXML_t      fermact;      			};
		  bool            obsvP;        			
		  GroupXML_t      invParam;     			
		};
	  quarks[n] is a QuarkSourceSolutions_t
	  quarks[n].dilutions[i] is a QuarkSolution_t
	  quarks[n].dilutions[i].soln is a LatticeFermion
	  quarks[n].dilutions[i].source_header is a PropSourceConst_t
	  quarks[n].dilutions[i].prop_header is a ChromaProp_t
	*/
  multi1d<Params::QuarkSourceSolutions_t>  quarks( params.qprop.solns.size() );
  QDPIO::cout << "num_quarks= " << params.qprop.solns.size() << endl;
  try
  {
    //QDPIO::cout << "quarks.size= " << quarks.size() << endl;
    for(int n=0; n < quarks.size(); ++n)
    {
			//QDPIO::cout << "Attempt to read solutions for source number=" << n << endl;
			quarks[n].dilutions.resize(params.qprop.solns[n].soln_file_names.size());
			//QDPIO::cout << "dilutions.size= " << quarks[n].dilutions.size() << endl;
			for(int i=0; i < quarks[n].dilutions.size(); ++i)
			{
				XMLReader file_xml, record_xml;
				QDPFileReader from(file_xml, params.qprop.solns[n].soln_file_names[i], QDPIO_SERIAL);
					read(from, record_xml, quarks[n].dilutions[i].soln);
				close(from);
				// need the seed
				read(record_xml, "/Propagator/PropSource/Source/ran_seed", quarks[n].seed);
				params.dilution[ n ].ran_seed = quarks[ n ].seed;
				// read headers ... may not be necessary
	  		read(record_xml, "/Propagator/PropSource",  quarks[n].dilutions[i].source_header);
	  		read(record_xml, "/Propagator/ForwardProp", quarks[n].dilutions[i].prop_header);
				/* some reminders
				myPropSourceConst = quarks[n].dilutions[i].source_header;
				            mystr = myPropSourceConst.source.xml;  // XML
				            mystr = myPropSourceConst.source.id;   // RAND_DILUTE_ZN_SOURCE
				            mystr = myPropSourceConst.source.path; // /Source
				     myChromaProp = quarks[n].dilutions[i].prop_header;
				       myGroupXML = myPropSourceConst.source;
				QDPIO::cout<< "quark number " << n << " : dilution number " << i <<endl;
				QDPIO::cout<< "j_decay="<<myPropSourceConst.j_decay << " t_source=" << myPropSourceConst.t_source <<endl;
				QDPIO::cout<<"seed2:"<<quarks[n].seed<<endl<<endl;
				*/
			}
    }
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error extracting headers: " << e << endl;
    QDP_abort(1);
  }
	// =============
	// Noise Sources ... re-make from the seeds
	// =============
  QDPIO::cout << "re-generating the noise sources" << endl;
	try
	{
	  int N;
		for(int n=0; n < quarks.size(); ++n)
	  {
	    bool first = true;
	    LatticeFermion quark_noise;      // noisy source on entire lattice
	    for(int i=0; i < quarks[ n ].dilutions.size(); ++i)
	    {
	      std::istringstream xml_s( quarks[ n ].dilutions[ i ].source_header.source.xml );
	      XMLReader sourcetop( xml_s );
	      //	QDPIO::cout << "Source = " << quarks[n].dilutions[i].source_header.source.id << endl;
	      if ( quarks[ n ].dilutions[ i ].source_header.source.id != DiluteZNQuarkSourceConstEnv::name )
	      {
	        QDPIO::cerr << "Expected source_type = " << DiluteZNQuarkSourceConstEnv::name << endl;
	        QDP_abort( 1 );
	      }
	      QDPIO::cout << "Quark num= " << n << "  dilution num= " << i << endl;
	      DiluteZNQuarkSourceConstEnv::Params srcParams( sourcetop,
	                                                     quarks[ n ].dilutions[ i ].source_header.source.path );
	      DiluteZNQuarkSourceConstEnv::SourceConst<LatticeFermion> srcConst( srcParams );
	      if ( first )
	      {
	        first = false;
	        quarks[ 0 ].j_decay = srcParams.j_decay;
	        // Grab N
	        N = srcParams.N;
	        // Set the seed to desired value
	        quarks[ n ].seed = srcParams.ran_seed;
	        QDP::RNG::setrn( quarks[ n ].seed );
	        // Create the noisy quark source on the entire lattice
	        zN_src( quark_noise, N );
	      }
	      // The seeds must always agree - here the seed is the unique id of the source
	      if ( toBool( srcParams.ran_seed != quarks[ n ].seed ) )
	      {
	        QDPIO::cerr << "quark source=" << n << "  dilution=" << i << " seed does not match" << endl;
	        QDP_abort( 1 );
	      }
	      // The N's must always agree
	      if ( toBool( srcParams.N != N ) )
	      {
	        QDPIO::cerr << "quark source=" << n << "  dilution=" << i << " N does not match" << endl;
	        QDP_abort( 1 );
	      }
	      // Use a trick here, create the source and subtract it from the global noisy
	      // Check at the end that the global noisy is zero everywhere.
	      // NOTE: the seed will be set every call
	      quarks[ n ].dilutions[ i ].source = srcConst( params.gaugestuff.u );
	      quark_noise -= quarks[ n ].dilutions[ i ].source;
	    } // end for i

	    Double dcnt = norm2( quark_noise );
	    if ( toDouble( dcnt ) != 0.0 )   // problematic - seems to work with unnormalized sources
	    {
	      QDPIO::cerr << "Noise not saturated by all potential solutions: dcnt=" << dcnt << endl;
	      QDP_abort( 1 );
	    }
	  } // end for n
	} // end try
	catch ( const std::string& e )
	{
	  QDPIO::cerr << ": Caught Exception creating source: " << e << endl;
	  QDP_abort( 1 );
	}
	QDPIO::cout << "Source vectors reconstructed from Seeds" << endl;
	
	// ====================================
	//       Main Part of Program
	// Run through the QQQ combinations and
	// multiply by appropriate coefficients
	// to make the various operators
	// ====================================
	//
	// Create the group baryon operator object
	//
#if 0
	  // Diagnostic
	  {
      for(int n=0; n < quarks.size(); ++n)
      {
				for(int i=0; i < quarks[n].dilutions.size(); ++i)
				{
	    		// Keep a copy of the phases with NO momenta
	    		SftMom phases_nomom(0, true, 3);

	    		multi1d<Double> source_corr = sumMulti(localNorm2(quarks[n].dilutions[i].source), 
								   phases_nomom.getSet());

	    		multi1d<Double> soln_corr = sumMulti(localNorm2(quarks[n].dilutions[i].soln), 
								 phases_nomom.getSet());

	    		push(xml_out, "elem");
	    		write(xml_out, "n", n);
	    		write(xml_out, "i", i);
	    		write(xml_out, "source_corr", source_corr);
	    		write(xml_out, "soln_corr", soln_corr);
	    		pop(xml_out);
				}
			}
	  }
#endif
/*	
  std::istringstream  xml_op(params.param.baryon_operator);
  XMLReader  optop(xml_op);
  const string operator_path = "/BaryonOperator";
  Handle< GroupBaryonQQQ >
    baryonOperator(TheWilsonBaryonOperatorFactory::Instance().createObject
		               (
								    params.param.baryon_operator,
								    optop,
								    operator_path,
								    params.gaugestuff.u 
									 )
									);
*/
	//
	//enum PlusMinus {PLUS = 1, MINUS = -1};
	//
	// Sinks first (only one term)
	int ordering = 0;
	//enum PlusMinus Plus;
	multi1d<LatticeComplex> result(1);
	//SftMom phases_nomom(0, true, quarks[n].dilutions[i].source_header.j_decay);
	// dont average over equivalent momenta
  SftMom phases(params.mom2_max, false, params.j_decay);
	for(int qqq=0; qqq < params.NQQQs; ++qqq)
	{
		for(int i=0; i < params.NH[ordering][0]; ++i)
		for(int j=0; j < params.NH[ordering][1]; ++j)
		for(int k=0; k < params.NH[ordering][2]; ++k)
		{
			result = AQQQ[ qqq ]( quarks[ 0 ].dilutions[ i ].soln,
    	                      quarks[ 1 ].dilutions[ j ].soln,
    	                      quarks[ 2 ].dilutions[ k ].soln,
    	                      PLUS );
#if 1
	  // Diagnostic
	  {
      for(int n=0; n < quarks.size()/quarks.size(); ++n)
      {
				for(int w=0; w < quarks[n].dilutions.size()/quarks[n].dilutions.size(); ++w)
				{
	    		// Keep a copy of the phases with NO momenta
	    		SftMom phases_nomom(0, true, params.j_decay);

	    		multi1d<Double> source_corr = sumMulti( localNorm2(result[0]), 
								                                  phases_nomom.getSet() );
	    		multi1d<Double> q0_corr = sumMulti( localNorm2(quarks[ 0 ].dilutions[ i ].soln), 
								                                  phases_nomom.getSet() );
	    		multi1d<Double> q1_corr = sumMulti( localNorm2(quarks[ 1 ].dilutions[ j ].soln), 
								                                  phases_nomom.getSet() );
	    		multi1d<Double> q2_corr = sumMulti( localNorm2(quarks[ 2 ].dilutions[ k ].soln), 
								                                  phases_nomom.getSet() );

	    		push(xml_out, "Elem");
	    		write(xml_out, "N", n);
	    		write(xml_out, "I", w);
	    		write(xml_out, "Source_corr", source_corr);
	    		write(xml_out, "q0_corr", q0_corr);
	    		write(xml_out, "q1_corr", q1_corr);
	    		write(xml_out, "q2_corr", q2_corr);
	    		pop(xml_out);
				}
			}
	  }
#endif
QDPIO::cout<<"done with "<<i<<" "<<j<<" "<<k<<" results size is "<<result.size()<<endl;
QDPIO::cout<<"NB="<<AQQQ[ qqq ].NBaryonOps<<endl;
QDPIO::cout<<"Np="<<params.Nmomenta<<endl;
	    for(int c=0; c < AQQQ[ qqq ].NBaryonOps; ++c)
			{
//				AQQQ[ qqq ].baryon[ c ]->termInCorr[ 0 ].hlist( i, j, k ).mom(0,0)
//				= cmplx(onee,twoo);
				AQQQ[ qqq ].baryon[ c ]->termInCorr[ 0 ].hlist( i, j, k ).mom
				= phases.sft( result[ 0 ] );
			}
QDPIO::cout<<"done with p=0 computation"<<endl;
	    for(int c=0; c < AQQQ[ qqq ].NBaryonOps; ++c)
			{
	    	for(int p=0; p < params.Nmomenta; ++p)
				{
	    		//for(int t=0; t < params.nrow[ 3 ]; ++t)
					{
//						AQQQ[ qqq ].baryon[ c ]->termInCorr[ 0 ].hlist( i, j, k ).mom( p, 0 ) 
//						*= AQQQ[ qqq ].coef[ c ];
					}
				}
QDPIO::cout<<"done with baryon=0 computation"<<endl;
			}
QDPIO::cout<<AQQQ[ 0 ].baryon[ 0 ]->termInCorr[ 0 ].hlist( i,j,k ).mom( 0, 0 )<<endl;
		}
	}
	
/*	
	std::istringstream xml_op( params.param.baryon_operator );
	XMLReader optop( xml_op );
	const string operator_path = "/GroupBaryonOperator";
	
	// baryon_operator_type NUCLEONS
	Handle< BaryonOperator<LatticeFermion> >
	baryonOperator( TheWilsonBaryonOperatorFactory::Instance().createObject( 
	                params.param.baryon_operator,
	                optop,
	                operator_path,
	                params.gaugestuff.u ) 
								);
	//
	// Permutations of quarks within an operator
	//    only doing one of them now ...
	//    but this will change once we start to do decays 
	int num_orderings = 1;   // number of permutations of the numbers  0,1,2
	multi1d< multi1d<int> > perms( num_orderings );
	{
	  multi1d<int> p( 3 );
	  if ( num_orderings >= 1 )
	  {
	    p[ 0 ] = 0;
	    p[ 1 ] = 1;
	    p[ 2 ] = 2;
	    perms[ 0 ] = p;
	  }
	  if ( num_orderings >= 2 )
	  {
	    p[ 0 ] = 0;
	    p[ 1 ] = 2;
	    p[ 2 ] = 1;
	    perms[ 1 ] = p;
	  }
	  if ( num_orderings >= 3 )
	  {
	    p[ 0 ] = 1;
	    p[ 1 ] = 0;
	    p[ 2 ] = 2;
	    perms[ 2 ] = p;
	  }
	  if ( num_orderings >= 4 )
	  {
	    p[ 0 ] = 1;
	    p[ 1 ] = 2;
	    p[ 2 ] = 0;
	    perms[ 3 ] = p;
	  }
	  if ( num_orderings >= 5 )
	  {
	    p[ 0 ] = 2;
	    p[ 1 ] = 1;
	    p[ 2 ] = 0;
	    perms[ 4 ] = p;
	  }
	  if ( num_orderings >= 6 )
	  {
	    p[ 0 ] = 2;
	    p[ 1 ] = 0;
	    p[ 2 ] = 1;
	    perms[ 5 ] = p;
	  }
	}
	//
	// Annihilation Operator A
	//
	/*
  struct BaryonOperator_t                              
  {
    struct BaryonOperatorInsertion_t
    {
      struct BaryonOperatorIndex_t
      {
        struct BaryonOperatorElement_t
        {
          multi2d<DComplex> elem;
        };
        multi1d<BaryonOperatorElement_t> ind;
      };
      multi3d<BaryonOperatorIndex_t> op;
    };
    multi1d<BaryonOperatorInsertion_t> orderings;
    Seed seed_l, seed_m, seed_r;
    multi1d< multi1d<int> > perms;
    int mom2_max, j_decay;
    multi1d<Complex> serialize();
  }; 
	
	class GroupBaryonOp : public BaryonOperator<LatticeFermion>
	{ public:
    GroupBaryonOp( const Params& p, const multi1d<LatticeColorMatrix>& u );
    std::string Name;  // Name of the Operator ... may also be the filename
    //protected:
    // .termInCorr[ whichterm ].hlist(i,j,k).mom[ p ]
    struct termInCorr_t
    {
      struct hlist_t
      {
        multi1d<DComplex> mom;
      };
      multi3d<hlist_t> hlist;
    };
    multi1d<termInCorr_t> termInCorr; // termInCorr(1); used to be size 2
    GroupBaryonOp(){}
    Params params;
    multi1d<LatticeColorMatrix> u_smr;
    SpinMatrix spin_rotate_mat;
		multi1d<Complex> serialize();
	}; // end class GroupBaryonOp
	*/
/*
//	swatch.start();
	BaryonOperator_t baryon_opA;
	baryon_opA.mom2_max = params.mom2_max;
	baryon_opA.j_decay = params.j_decay;
	baryon_opA.seed_l = quarks[ 0 ].seed;
	baryon_opA.seed_m = quarks[ 1 ].seed;
	baryon_opA.seed_r = quarks[ 2 ].seed;
	baryon_opA.orderings.resize( num_orderings );
	baryon_opA.perms.resize( num_orderings );

	push( xml_out, "OperatorA" );

	// Sanity check
	if ( toBool( baryon_opA.seed_l == baryon_opA.seed_m ) )
	{
	  QDPIO::cerr << "baryon op seeds are the same" << endl;
	  QDP_abort( 1 );
	}
	// Sanity check
	if ( toBool( baryon_opA.seed_l == baryon_opA.seed_r ) )
	{
	  QDPIO::cerr << "baryon op seeds are the same" << endl;
	  QDP_abort( 1 );
	}
	// Sanity check
	if ( toBool( baryon_opA.seed_m == baryon_opA.seed_r ) )
	{
	  QDPIO::cerr << "baryon op seeds are the same" << endl;
	  QDP_abort( 1 );
	}
	// Construct annihilation operator A
	try
	{
	  for ( int ord = 0; ord < baryon_opA.orderings.size(); ++ord )
	  {
	    QDPIO::cout << "Operator A: ordering = " << ord << endl;

	    baryon_opA.perms[ ord ] = perms[ ord ];

	    // Operator construction
	    const QuarkSourceSolutions_t& q0 = quarks[ perms[ ord ][ 0 ] ];
	    const QuarkSourceSolutions_t& q1 = quarks[ perms[ ord ][ 1 ] ];
	    const QuarkSourceSolutions_t& q2 = quarks[ perms[ ord ][ 2 ] ];

	    baryon_opA.orderings[ ord ].op.resize( q0.dilutions.size(),
	                                           q1.dilutions.size(),
	                                           q2.dilutions.size() );

	    for ( int i = 0; i < q0.dilutions.size(); ++i )
	    {
	      for ( int j = 0; j < q1.dilutions.size(); ++j )
	      {
	        for ( int k = 0; k < q2.dilutions.size(); ++k )
	        {
	          multi1d<LatticeComplex> bar = ( *baryonOperator ) ( q0.dilutions[ i ].source,
	                                                              q1.dilutions[ j ].source,
	                                                              q2.dilutions[ k ].source,
	                                                              MINUS );

	          baryon_opA.orderings[ ord ].op( i, j, k ).ind.resize( bar.size() );
	          for ( int l = 0; l < bar.size(); ++l )
	            baryon_opA.orderings[ ord ].op( i, j, k ).ind[ l ].elem = phases.sft( bar[ l ] );

	        } // end for k
	      } // end for j
	    } // end for i
	  } // end for ord
	} // end try
	catch ( const std::string& e )
	{
	  QDPIO::cerr << ": Caught Exception creating group baryon sink operator: " << e << endl;
	  QDP_abort( 1 );
	}
	catch ( ... )
	{
	  QDPIO::cerr << ": Caught generic exception creating group baryon sink operator" << endl;
	  QDP_abort( 1 );
	}
	pop( xml_out ); // OperatorA
//	swatch.stop();
//	QDPIO::cout << "Operator A computed: time= " << swatch.getTimeInSeconds() << " secs" << endl;


	//
	// Creation Operator B    ... probably should have called it C
	//
//	swatch.start();
	BaryonOperator_t baryon_opB;
	baryon_opB.mom2_max = params.mom2_max;
	baryon_opB.j_decay = j_decay;
	baryon_opB.seed_l = quarks[ 0 ].seed;
	baryon_opB.seed_m = quarks[ 1 ].seed;
	baryon_opB.seed_r = quarks[ 2 ].seed;
	baryon_opB.orderings.resize( num_orderings );
	baryon_opB.perms.resize( num_orderings );

	push( xml_out, "OperatorB" );

	// Sanity check
	if ( toBool( baryon_opB.seed_l == baryon_opB.seed_m ) )
	{
	  QDPIO::cerr << "baryon op seeds are the same" << endl;
	  QDP_abort( 1 );
	}
	// Sanity check
	if ( toBool( baryon_opB.seed_l == baryon_opB.seed_r ) )
	{
	  QDPIO::cerr << "baryon op seeds are the same" << endl;
	  QDP_abort( 1 );
	}
	// Sanity check
	if ( toBool( baryon_opB.seed_m == baryon_opB.seed_r ) )
	{
	  QDPIO::cerr << "baryon op seeds are the same" << endl;
	  QDP_abort( 1 );
	}
	//
	// Construct creation operator B
	//
	try
	{
	  for ( int ord = 0; ord < baryon_opB.orderings.size(); ++ord )
	  {
	    QDPIO::cout << "Operator B: ordering = " << ord << endl;

	    baryon_opB.perms[ ord ] = perms[ ord ];

	    // Operator construction
	    const QuarkSourceSolutions_t& q0 = quarks[ perms[ ord ][ 0 ] ];
	    const QuarkSourceSolutions_t& q1 = quarks[ perms[ ord ][ 1 ] ];
	    const QuarkSourceSolutions_t& q2 = quarks[ perms[ ord ][ 2 ] ];

	    baryon_opB.orderings[ ord ].op.resize( q0.dilutions.size(),
	                                           q1.dilutions.size(),
	                                           q2.dilutions.size() );

	    for ( int i = 0; i < q0.dilutions.size(); ++i )
	    {
	      for ( int j = 0; j < q1.dilutions.size(); ++j )
	      {
	        for ( int k = 0; k < q2.dilutions.size(); ++k )
	        {
	          multi1d<LatticeComplex> bar = ( *baryonOperator ) ( q0.dilutions[ i ].soln,
	                                                              q1.dilutions[ j ].soln,
	                                                              q2.dilutions[ k ].soln,
	                                                              PLUS );
	          baryon_opB.orderings[ ord ].op( i, j, k ).ind.resize( bar.size() );
	          for ( int l = 0; l < bar.size(); ++l )
	            baryon_opB.orderings[ ord ].op( i, j, k ).ind[ l ].elem = phases.sft( bar[ l ] );
	        } // end for k
	      } // end for j
	    } // end for i
	  } // end for ord
	} // end try
	catch ( const std::string& e )
	{
	  QDPIO::cerr << ": Caught Exception creating group baryon source operator: " << e << endl;
	  QDP_abort( 1 );
	}
	catch ( ... )
	{
	  QDPIO::cerr << ": Caught generic exception creating group baryon source operator" << endl;
	  QDP_abort( 1 );
	}

	pop( xml_out ); // OperatorB

//	swatch.stop();

//	QDPIO::cout << "Operator B computed: time= " << swatch.getTimeInSeconds() << " secs" << endl;
*/

	// Save the operators
	// ONLY SciDAC output format is supported!
//	swatch.start();
/*
	{
	  XMLBufferWriter file_xml;
	  push( file_xml, "baryon_operator" );
	  file_xml << params.param.baryon_operator;
	  write( file_xml, "Config_info", gauge_xml );
	  pop( file_xml );

	  QDPFileWriter to( file_xml, params.named_obj.prop.op_file,      // are there one or two files???
	                    QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN );

	  // Write the scalar data
	  {
	    XMLBufferWriter record_xml;
	    write( record_xml, "SourceBaryonOperator", baryon_opA );
	    write( to, record_xml, baryon_opA.serialize() );
	  }

	  // Write the scalar data
	  {
	    XMLBufferWriter record_xml;
	    write( record_xml, "SinkBaryonOperator", baryon_opB );
	    write( to, record_xml, baryon_opB.serialize() );
	  }

	  close( to );
	}
*/

//	swatch.stop();

//	QDPIO::cout << "Operators written: time= " << swatch.getTimeInSeconds() << " secs" << endl;

#if 0
  /*
   * Read in a Chroma prop
   */
  LatticePropagator prop;
  XMLReader prop_in_xml;
	XMLReader prop_in_file_xml;

  read( prop_in_xml, "/Propagator/ForwardProp", prop_header );
	
  push( xml_out, "SciDAC_propagator" );
  write( xml_out, "prop_in_file", input.prop.prop_in_file );

  readQprop( prop_in_file_xml, 
	           prop_in_xml, 
						 prop,
             input.prop.prop_in_file, 
						 QDPIO_SERIAL 
						);
  write( xml_out, "File_xml", prop_in_file_xml );
  write( xml_out, "Record_xml", prop_in_xml );
  pop( xml_out );
  // Try to invert this record XML into a source struct
  // Also pull out the id of this source
  ChromaProp_t prop_header;
  PropSourceConst_t source_header;
  try
  {
    read( prop_in_xml, "/Propagator/ForwardProp", prop_header );
    read( prop_in_xml, "/Propagator/PropSource", source_header );
  }
  catch ( const string & e )
  {
    QDPIO::cerr << "Error extracting forward_prop header: " << e << endl;
    throw;
  }
  // Derived from input prop
  int j_decay = source_header.j_decay;
  int t_source = source_header.t_source;
  // Initialize the slow Fourier transform phases
  // This is used to get the time-slice subsets in the j_decay dir
  SftMom phases( 0, true, j_decay );

  // Length of lattice in j_decay direction and 3pt correlations fcns
  int length = phases.numSubsets();

  // Sanity check - write out the propagator (pion) correlator in the j_decay direction
  {
    multi1d<Double> prop_corr = sumMulti( localNorm2( prop ),
                                          phases.getSet() );
    push( xml_out, "Prop_correlator" );
    write( xml_out, "prop_corr", prop_corr );
    pop( xml_out );
  }
  xml_out.flush();
  /*
   * Charge-conjugation and time-reverse the beasty
   */
  {
    /* Time-charge reverse the quark propagators */
    /* S_{CT} = gamma_5 gamma_4 = gamma_1 gamma_2 gamma_3 = Gamma(7) */
    LatticePropagator prop_tmp = - ( Gamma( 7 ) * prop * Gamma( 7 ) );
    // This is a really dumb way to implement this. Shift each slice around.
    // Do nothing on the t=0 slice
    prop[ phases.getSet() [ 0 ] ] = prop_tmp;

    LatticePropagator tmp1, tmp2, tmp3;

    for(int t=1; t < length; ++t)
    {
      int tp = ( length - t ) % length;
      if ( t < tp )
      {
        tmp1[ phases.getSet() [ tp ] ] = prop_tmp;
        for(int k=tp - 1; k >= t; --k )
        {
          tmp2[ phases.getSet() [ k ] ] = shift( tmp1, FORWARD, j_decay );
          tmp1[ phases.getSet() [ k ] ] = tmp2;
        }
        prop[ phases.getSet() [ t ] ] = tmp1;
      }
      else if ( t == tp )
      {
        prop[ phases.getSet() [ t ] ] = prop_tmp;
      }
      else if ( t > tp )
      {
        tmp1[ phases.getSet() [ tp ] ] = prop_tmp;
        for(int k=tp + 1; k <= t; ++k)
        {
          tmp2[ phases.getSet() [ k ] ] = shift( tmp1, BACKWARD, j_decay );
          tmp1[ phases.getSet() [ k ] ] = tmp2;
        }
        prop[ phases.getSet() [ t ] ] = tmp1;
      }
    }
  }
  // Sanity check - write out the propagator (pion) correlator in the j_decay direction
  {
    multi1d<Double> trev_prop_corr = sumMulti( localNorm2( prop ),
                                     phases.getSet() );
    push( xml_out, "TRevProp_correlator" );
    write( xml_out, "trev_prop_corr", trev_prop_corr );
    pop( xml_out );
  }
  /*
   * Now write them thangs...
   */
  {
    XMLBufferWriter prop_out_file_xml;
    push( prop_out_file_xml, "propagator" );
    int id = 0;    // NEED TO FIX THIS - SOMETHING NON-TRIVIAL NEEDED
    write( prop_out_file_xml, "id", id );
    pop( prop_out_file_xml );

    source_header.t_source = ( length - source_header.t_source ) % length;

    XMLBufferWriter prop_out_record_xml;
    push( prop_out_record_xml, "Propagator" );
    write( prop_out_record_xml, "ForwardProp", prop_header );
    write( prop_out_record_xml, "PropSource", source_header );
    {
      QDPIO::cout << "Create config info" << endl;
      XMLReader gauge_xml( prop_in_xml, "/Propagator/Config_info" );
      ostringstream gauge_str;
      gauge_xml.print( gauge_str );
      write( prop_out_record_xml, "Config_info", gauge_str );
      QDPIO::cout << "Done config info" << endl;
    }
    pop( prop_out_record_xml );

    // Write the source
    writeQprop( prop_out_file_xml, prop_out_record_xml, prop,
                input.prop.prop_out_file, input.prop.prop_out_volfmt,
                QDPIO_SERIAL );
  }
#endif

  pop( xml_out );   // GroupBaryonAll2All
  END_CODE();
  // Time to bolt
  Chroma::finalize();
  exit( 0 );
}
