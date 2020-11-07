
#include "chroma.h"
#include "chroma_config.h"
#include <string>
#include "actions/gauge/gaugeacts/rect_gaugeact.h"
#include "actions/gauge/gaugeacts/plaq_plus_spatial_two_plaq_gaugeact.h"
#include "update/molecdyn/integrator/lcm_toplevel_integrator.h"

// Specials
#include "update/molecdyn/hamiltonian/exact_hamiltonian.h"

using namespace Chroma;

#ifdef CHROMA_USE_ITTNOTIFY
#include "ittnotify.h"
#endif

//! To insure linking of code, place the registered code flags here
/*! This is the bit of code that dictates what fermacts are in use */
bool linkageHack(void)
{
	bool foo = true;

	// Gauge Monomials
	foo &= GaugeMonomialEnv::registerAll();

	// Ferm Monomials
	foo &= WilsonTypeFermMonomialAggregrateEnv::registerAll();

	// MD Integrators
	foo &= LCMMDComponentIntegratorAggregateEnv::registerAll();

	// Chrono predictor
	foo &= ChronoPredictorAggregrateEnv::registerAll();

	// Inline Measurements
	foo &= InlineAggregateEnv::registerAll();

	return foo;
}

int main(int argc, char *argv[]) 
{
	Chroma::initialize(&argc, &argv);

	START_CODE();

	// Chroma Init stuff -- Open DATA and XMLDAT
	QDPIO::cout << "Linkage = " << linkageHack() << std::endl;

	// Snarf it all
	XMLReader param_in(Chroma::getXMLInputFileName());
	XMLReader paramtop(param_in, "/MonomialTests");

	multi1d<int> nrow(Nd);

	try {
		read(paramtop, "nrow", nrow);
	}
	catch(const std::string& e) {
		QDPIO::cerr << "Unable to read nrow from XML: " << e << std::endl;
		QDP_abort(1);
	}

	Layout::setLattSize(nrow);
	Layout::create();

	QDPIO::cout << " HELLLOOOO!!!!!!!!"  << std::endl;
	// Dump output
	XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
	XMLFileWriter& xml_log = Chroma::getXMLLogInstance();
	push(xml_out, "MonomialTimings");
	push(xml_log, "MonomialTimings");

	// Read Parameters
	multi1d<int> boundary(Nd);           // Ferm BC's
	std::string monomial_name;           // String for Factory


	Cfg_t cfg;
	try {
		read(paramtop, "./GaugeStartup", cfg);
	}
	catch( const std::string& e ) {
		QDPIO::cerr << " Error reading XML " << e << std::endl;
		QDP_abort(1);
	}

	multi1d<LatticeColorMatrix> u(Nd);
	{
		XMLReader file_xml;
		XMLReader config_xml;

		gaugeStartup(file_xml, config_xml, u, cfg);
	}

	// Try and create an array of monomials:
	try {
		readNamedMonomialArray(paramtop, "./Monomials");
	}
	catch(const std::string& e) {
		QDPIO::cout << "Failed to read monomials " << std::endl;
		QDP_abort(1);
	}
	QDPIO::cout << "ALL MONOMIALS READ " << std::endl << std::flush;

	// Read the list of monomials to test.
	multi1d<std::string> monomial_test_ids;
	multi1d<int> n_iters;
	try {
	read(paramtop, "//MonomialTestIds", monomial_test_ids);
	}
	catch( ... ) {
		QDPIO::cout << "Caught exception: reading monomial_test_ids" <<std::endl;
		QDP_abort(1);
	}

	QDPIO::cout << "MONOMIAL TEST IDS READ " << std::endl << std::flush;

	try {
	read(paramtop, "//n_iters", n_iters);
	}
	catch( const std::string e) {
		QDPIO::cout << "Couldnt read n_iters. Caught exception: "<< e << std::endl;
		QDP_abort(1);
	}

	QDPIO::cout << "ALL parameters read " << std::endl;
	if ( n_iters.size() != monomial_test_ids.size())  {
		if( n_iters.size() == 1 ) {
			QDPIO::cout << "Using same value of n_iters for all monomials" << std::endl;
			int n_iters_val = n_iters[0];
			n_iters.resize( monomial_test_ids.size());
			for(int i=0; i < monomial_test_ids.size(); ++i) {
				n_iters[i]=n_iters_val;
			}
		}
		else {
			QDPIO::cout <<" There are " << monomial_test_ids.size()
				  			  <<" monomials but " << n_iters.size()
							  <<" values for n_iters" << std::endl;
			QDP_abort(1);
		}
	}

	QDPIO::cout << "STARTING MOMENTA " << std::endl << std::flush ;
	// Fictitious momenta for now
	multi1d<LatticeColorMatrix> p(Nd);

	// Get some noise into the momenta...
	for(int mu=0; mu<Nd; mu++) {
		gaussian(p[mu]);
		p[mu] *= sqrt(0.5);
		taproj(p[mu]);
	}

	QDPIO::cout << "create state" << std::endl;

	// Create a field state
	GaugeFieldState gauge_state(p,u);

	QDPIO::cout << "There are " << monomial_test_ids.size() << " monomials to test" << std::endl;
	typedef multi1d<LatticeColorMatrix> P;
	typedef multi1d<LatticeColorMatrix> Q;
	QDPIO::cout << "Binding Monomials" << std::endl;

	multi1d< Chroma::IntegratorShared::MonomialPair > mon_pairs(monomial_test_ids.size());
	try {
		Chroma::IntegratorShared::bindMonomials( monomial_test_ids, mon_pairs );
	}
	catch( const std::string& err ) {
		QDPIO::cout << "Caught exception string:" << err << std::endl;
	}
	StopWatch swatch;
	double seconds;
	P F;

	for( int m = 0; m < monomial_test_ids.size(); m++) {
		// Get a monomial
		Monomial<P,Q>& the_mon = *( mon_pairs[m].mon );
		std::string& mon_name = mon_pairs[m].id;
		QDPIO::cout << "Timing monomial force: " << mon_name << std::endl;
		QDPIO::cout << "Refreshing internal fields: " << std::endl;

		the_mon.refreshInternalFields(gauge_state);

		QDPIO::cout << "Single Application in case one is JIT-ing" << std::endl;
		{

			swatch.reset();
			swatch.start();

			the_mon.dsdq(F,gauge_state);

			swatch.stop();
			seconds = swatch.getTimeInSeconds();
			QDPInternal::globalSum(seconds);
			seconds /= (double)Layout::numNodes();
			QDPIO::cout << "First application took: " << seconds <<" (s) "<< std::endl << std::flush;
		}

		{
			int hits = n_iters[m];
			QDPIO::cout << "Timing monomial force: " << mon_name << " with " << hits << " hits " << std::endl;
			swatch.reset();
			swatch.start();
#ifdef CHROMA_USE_ITTNOTIFY
			__itt_resume();
#endif
			for(int i=0; i < hits; i++) {
				the_mon.dsdq(F,gauge_state);
			}
#ifdef CHROMA_USE_ITTNOTIFY
			__itt_pause();
#endif
			swatch.stop();
			seconds = swatch.getTimeInSeconds();
			QDPInternal::globalSum(seconds);
			seconds /= (double)Layout::numNodes();

			QDPIO::cout << "monomial_id = " << mon_name
					<< ", Nhits = " << hits
					<< ", Time = " << seconds << "(s)"
					<< ", Time per hit = " << seconds/(double)hits << "(s)"<< std::endl;

			push(xml_out, "ForceTiming");
			write(xml_out, "monomial_id", mon_name);
			write(xml_out, "nhits", hits);
			write(xml_out, "time", seconds);
			write(xml_out, "time_per_hit", (seconds/(double)hits));
			pop(xml_out);
		}
	}
#if 0  
		// Time taproj
		seconds = 0;
		hits = 1;
		QDPIO::cout << "Timing taproj()" << std::endl;
		QDPIO::cout << "Calibrating";
		while( seconds < 10 ) {
			hits *= 2;
			swatch.reset();
			swatch.start();
			for(int i=0; i < hits; i++) {
				for(int mu=0; mu < Nd; mu++) {
					taproj(gauge_state.getP()[mu]);
				}
			}
			swatch.stop();
			seconds = swatch.getTimeInSeconds();
			QDPInternal::globalSum(seconds);
			seconds /= (double)Layout::numNodes();
			QDPIO::cout << "." << std::flush;
		}

		QDPIO::cout << " " << seconds << "(s)"<< std::endl;


		QDPIO::cout << "Timing taproj with " << hits << " hits " << std::endl;
		swatch.reset();
		swatch.start();
		for(int i=0; i < hits; i++) {
			for(int mu=0; mu < Nd; mu++) {
				taproj(gauge_state.getP()[mu]);
			}
		}
		swatch.stop();
		seconds = swatch.getTimeInSeconds();
		QDPInternal::globalSum(seconds);
		seconds /= (double)Layout::numNodes();

		QDPIO::cout << "taproj"
				<< ", Nhits = " << hits
				<< ", Time = " << seconds << "(s)"
				<< ", Time per hit = " << seconds/(double)hits << "(s)" << std::endl;

		push(xml_out, "TaprojTiming");
		write(xml_out, "nhits", hits);
		write(xml_out, "time", seconds);
		write(xml_out, "time_per_hit", (seconds/(double)hits));
		pop(xml_out);

		// Time expmat
		seconds = 0;
		hits = 1;
		QDPIO::cout << "Timing expmat()" << std::endl;
		QDPIO::cout << "Calibrating";

		while( seconds < 10 ) {
			hits *= 2;
			swatch.reset();
			swatch.start();
			for(int i=0; i < hits; i++) {
				for(int mu=0; mu < Nd; mu++) {
					expmat(gauge_state.getP()[mu], EXP_EXACT);
				}
			}
			swatch.stop();
			seconds = swatch.getTimeInSeconds();
			QDPInternal::globalSum(seconds);
			seconds /= (double)Layout::numNodes();

			QDPIO::cout << "." << std::flush;
		}
		QDPIO::cout << " " << seconds << "(s)"<< std::endl;

		QDPIO::cout << "Timing expmat with " << hits << " hits " << std::endl;
		swatch.reset();
		swatch.start();
		for(int i=0; i < hits; i++) {
			for(int mu=0; mu < Nd; mu++) {
				expmat(gauge_state.getP()[mu], EXP_EXACT);
			}
		}
		swatch.stop();
		seconds = swatch.getTimeInSeconds();
		QDPInternal::globalSum(seconds);
		seconds /= (double)Layout::numNodes();

		QDPIO::cout << "expmat"
				<< ", Nhits = " << hits
				<< ", Time = " << seconds << "(s)"
				<< ", Time per hit = " << seconds/(double)hits << "(s)" << std::endl;

		push(xml_out, "ExpmatTiming");
		write(xml_out, "nhits", hits);
		write(xml_out, "time", seconds);
		write(xml_out, "time_per_hit", (seconds/(double)hits));
		pop(xml_out);

		// Time reunit
		seconds = 0;
		hits = 1;
		int numbad;
		QDPIO::cout << "Timing reunit()" << std::endl;
		QDPIO::cout << "Calibrating";

		while( seconds < 10 ) {
			hits *= 2;
			swatch.reset();
			swatch.start();
			for(int i=0; i < hits; i++) {
				for(int mu=0; mu < Nd; mu++) {
					reunit((gauge_state.getQ())[mu], numbad, REUNITARIZE_ERROR);
				}
			}
			swatch.stop();
			seconds = swatch.getTimeInSeconds();
			QDPInternal::globalSum(seconds);
			seconds /= (double)Layout::numNodes();

			QDPIO::cout << "." << std::flush;
		}
		QDPIO::cout << " " << seconds << "(s)"<< std::endl;

		QDPIO::cout << "Timing reunit with " << hits << " hits " << std::endl;
		swatch.reset();
		swatch.start();
		for(int i=0; i < hits; i++) {
			for(int mu=0; mu < Nd; mu++) {
				reunit((gauge_state.getQ())[mu], numbad, REUNITARIZE_ERROR);
			}
		}
		swatch.stop();
		seconds = swatch.getTimeInSeconds();
		QDPInternal::globalSum(seconds);
		seconds /= (double)Layout::numNodes();

		QDPIO::cout << "reunit"
				<< ", Nhits = " << hits
				<< ", Time = " << seconds << "(s)"
				<< ", Time per hit = " << seconds/(double)hits << "(s)" << std::endl;

		push(xml_out, "ReunitTiming");
		write(xml_out, "nhits", hits);
		write(xml_out, "time", seconds);
		write(xml_out, "time_per_hit", (seconds/(double)hits));
		pop(xml_out);

#endif

		pop(xml_log);   // t_leapfrog
		pop(xml_out);   // t_leapfrog

		END_CODE();

		// Finish
		Chroma::finalize();
		exit(0);
	}
