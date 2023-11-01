/*! \file
 * \brief Compute the disconnected diagrams with 4D probing
 *
 * Propagator calculation on a colorstd::vector
 */

#include "inline_disco_prob_defl_superb_w.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "fermact.h"
#include "meas/glue/mesplq.h"
#include "meas/hadron/greedy_coloring.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/make_xml_file.h"
#include "qdp_stdio.h"
#include "util/ferm/key_prop_colorvec.h"
#include "util/ferm/key_val_db.h"
#include "util/ferm/map_obj/map_obj_aggregate_w.h"
#include "util/ferm/map_obj/map_obj_factory_w.h"
#include "util/ferm/subset_vectors.h"
#include "util/ferm/superb_contractions.h"
#include "util/ferm/transf.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"

#include <cassert>
#include <string>
#include <complex>

#ifdef BUILD_SB

namespace Chroma 
{ 
  namespace InlineDiscoProbDeflSuperb 
  {
    

    //! Propagator input
    void read(XMLReader& xml, const std::string& path, InlineDiscoProbDeflSuperb::Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "sdb_file", input.sdb_file);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const InlineDiscoProbDeflSuperb::Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "sdb_file", input.sdb_file);

      pop(xml);
    }

    //! Propagator input
    void read(XMLReader& xml, const std::string& path, InlineDiscoProbDeflSuperb::Params::Param_t& param)
    {
      XMLReader inputtop(xml, path);
      
      read(inputtop,"max_path_length",param.max_path_length);
      if(inputtop.count("p2_max")!=0){	
      	read(inputtop,"p2_max",param.p2_max);
	param.use_p_list = false;
	QDPIO::cout<<"Using momenta centered at the origin, with a max of "<<param.p2_max<<std::endl;
      }
      else if(inputtop.count("p_file")!=0){
	read(inputtop,"p_num",param.p_num);
        read(inputtop,"p_file",param.p_file);
	param.use_p_list = true;
        int c1[param.p_num],c2[param.p_num],c3[param.p_num];
	int n,m,l,i=1,lines;
	std::ifstream read(param.p_file);
	while(read>>n>>m>>l){
	  c1[i]=n;
	  c2[i]=m;
	  c3[i]=l;
	  i++;
	}
	lines=i-1;//Total number of momenta.
	param.p_list.resize(lines, Nd - 1);
	for(int mom = 0; mom < lines; mom++)
	{
	  param.p_list[mom][0] = c1[mom];
	  param.p_list[mom][1] = c2[mom];
	  param.p_list[mom][2] = c3[mom];
	  QDPIO::cout<<"Momentum number "<<mom<<" is px "<<param.p_list[mom][0]<<" py "<<param.p_list[mom][1]<<" pz "<<param.p_list[mom][2]<<std::endl;
	}
	QDPIO::cout<<"Using momentum list "<<std::endl;
        QDPIO::cout<<"There are "<<lines<<" momenta in file."<<std::endl;
      }
      else
      {
	QDPIO::cout<<"Could not find valid XML momentum input."<<std::endl;
	QDP_abort(1);
      }
      read(inputtop,"mass_label",param.mass_label);

      read(inputtop,"Propagator",param.prop) ;
      
      if(inputtop.count("use_ferm_state_links")!=0){ 
	read(inputtop,"use_ferm_state_links",param.use_ferm_state_links) ;
        QDPIO::cout<<"Ferm state links set equal to "<<param.use_ferm_state_links<<std::endl;
      }
      else
	param.use_ferm_state_links=false ;

      param.projParam = readXMLGroup(inputtop,"Projector","projectorType") ;

      if(inputtop.count("probing_distance")!=0){ 
	read(inputtop,"probing_distance",param.probing_distance) ;
      }
      else
	param.probing_distance = param.max_path_length;

      if(inputtop.count("probing_power")!=0){ 
	read(inputtop,"probing_power",param.probing_power) ;
      }
      else
	param.probing_power = 2;

      if(inputtop.count("probing_file")!=0){ 
	read(inputtop,"probing_file",param.probing_file) ;
      }

      if(inputtop.count("noise_vectors")!=0){ 
	read(inputtop,"noise_vectors",param.noise_vectors) ;
      }
      else
	param.noise_vectors = 1;

      if(inputtop.count("max_rhs")!=0){ 
	read(inputtop,"max_rhs",param.max_rhs) ;
      }
      else
	param.max_rhs = 1;
     }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const InlineDiscoProbDeflSuperb::Params::Param_t& param)
    {
      push(xml, path);

      write(xml,"max_path_length",param.max_path_length);
      write(xml,"p2_max",param.p2_max);
      write(xml,"mass_label",param.mass_label);

      write(xml,"Propagator",param.prop) ;
      write(xml,"use_ferm_state_links",param.use_ferm_state_links) ;
      write(xml,"probing_distance",param.probing_distance) ;
      write(xml,"probing_power",param.probing_power) ;
      write(xml,"noise_vectors",param.noise_vectors) ;
      write(xml,"max_rhs",param.max_rhs) ;
      xml << param.projParam.xml;

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, InlineDiscoProbDeflSuperb::Params& input)
    {
      InlineDiscoProbDeflSuperb::Params tmp(xml, path);
      input = tmp;
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const InlineDiscoProbDeflSuperb::Params& input)
    {
      push(xml, path);
    
      write(xml, "Param", input.param);
      write(xml, "NamedObject", input.named_obj);

      pop(xml);
    }

    //! Meson operator     
    struct KeyOperator_t
    {
      int t_slice                ; /*!< Meson operator time slice */
      multi1d<int>       disp    ; /*!< Displacement dirs of quark (right)*/
      multi1d<int>       mom     ; /*!< D-1 momentum of this operator */
      std::string  mass_label    ; /*!< Mass label */
      
      KeyOperator_t(){
	mom.resize(Nd-1);
      }
    };
    
    bool operator<(const KeyOperator_t& a, const KeyOperator_t& b) {
      return (a.t_slice != b.t_slice ? a.t_slice < b.t_slice : (
               a.mom != b.mom ? a.mom < b.mom : (
               a.disp != b.disp ? a.disp < b.disp : (
               a.mass_label < b.mass_label))));
    }
    
    template <typename T, typename Q>
    T& operator<<(T& os, const multi1d<Q>& d)
    {
      for (int i=0; i<d.size();i++){
        os << d[i] << " " ;
      }
      return os;
    }

    template <typename T>
    T& operator<<(T& os, const KeyOperator_t& d)
    {
      os << "KeyOperator_t:"
         << " t_slice = " << d.t_slice
         << ", disp = ";
      for (int i=0; i<d.disp.size();i++){
        os << d.disp[i] << " " ;
      }
      os << ", mom = ";
      for (int i=0; i<d.mom.size();i++){
        os << d.mom[i] << " " ;
      }
      os << std::endl;

      return os;
    }
    class ValOperator_t{
    public:
      multi1d<ComplexD> op ;  
      ValOperator_t(){op.resize(Ns*Ns);} // Here go the 16 gamma matrices
      ~ValOperator_t(){}
    } ;

    //-------------------------------------------------------------------------
    //! stream IO
    template <typename T>
    T& operator<<(T& os, const ValOperator_t& d)
    {
      os << "ValOperator_t:\n";
      for (int i=0; i<d.op.size();i++){
	os <<"     gamma["<<i<<"] = "<< d.op[i] << std::endl ;
      }
      
      return os;
    }
    
    struct KeyVal_t{
      SerialDBKey <KeyOperator_t> k ;
      SerialDBData<ValOperator_t> v ;
    };

    //! KeyOperator reader    
    void read(BinaryReader& bin, KeyOperator_t& d){
      read(bin,d.t_slice);
      int n ;
      read(bin,n);
      d.disp.resize(n); 
      read(bin,d.disp);
      d.mom.resize(Nd-1) ;
      read(bin,d.mom);
      read(bin, d.mass_label, 32);
    }
    //! KeyOperator writer
    void write(BinaryWriter& bin, const KeyOperator_t& d){
      write(bin,d.t_slice);
      int n ;
      n = d.disp.size();
      write(bin,n);
      write(bin,d.disp);
      write(bin,d.mom);
      write(bin,d.mass_label);
    }

    //! ValOperator reader    
    void read(BinaryReader& bin, ValOperator_t& d){
      d.op.resize(Ns*Ns);
      read(bin,d.op);
    }
    //! ValOperator writer
    void write(BinaryWriter& bin, const ValOperator_t& d){
      write(bin,d.op);
    }


    /** This is does the noise part of the disconnected diagram **/ 
    void do_shift(LatticeFermion& q_mu, 
        const LatticeFermion& q, const multi1d<LatticeColorMatrix>& u, int mu, int sign)
    {
        if(sign>0)
          q_mu = u[mu] * shift(q, FORWARD, mu);
        else
          q_mu = shift(adj(u[mu]) * q, BACKWARD, mu);
    }

    template <typename T>
    struct detox_aux {
      using rtype = T;
      static T get(const T& e)
      {
	return e;
      }
    };

#ifdef QDP_IS_QDPJIT
    template <typename T>
    struct detox_aux<QDP::Word<T>> {
      using rtype = T;
      static T get(const QDP::Word<T>& w)
      {
	return w.elem();
      }
    };
#endif

    // Return the same value but with the most elemental type
    template <typename T>
    typename detox_aux<T>::rtype detox(const T& e)
    {
      return detox_aux<T>::get(e);
    }

    void do_disco(std::map<KeyOperator_t, ValOperator_t>& db,
		  const std::vector<std::shared_ptr<LatticeFermion>>& qbar,
		  const std::vector<std::shared_ptr<LatticeFermion>>& q, const SftMom& p,
		  const multi1d<LatticeColorMatrix>& u, const int& max_path_length)
    {

      const int Nt = Layout::lattSize()[3];
      const int a = qbar.size();
      assert(qbar.size() == q.size());

      // Copy q and qbar to tensor forms
      std::string q_order = "cxyzXnSst*";
      SB::Tensor<Nd + 6, SB::Complex> qt(
	q_order, SB::latticeSize<Nd + 6>(q_order, {{'n', 1}, {'S', Ns}, {'s', 1}, {'*', a}}));
      for (int i = 0; i < a; ++i)
	SB::asTensorView(*q[i])
	  .rename_dims({{'s', 'S'}})
	  .copyTo(qt.kvslice_from_size({{'*', i}}, {{'*', 1}}));
      std::string qbar_order = "cxyzXNQqt*";
      SB::Tensor<Nd + 6, SB::Complex> qbart(
	qbar_order, SB::latticeSize<Nd + 6>(qbar_order, {{'N', 1}, {'Q', Ns}, {'q', 1}, {'*', a}}));
      for (int i = 0; i < a; ++i)
	SB::asTensorView(*qbar[i])
	  .rename_dims({{'s', 'Q'}})
	  .copyTo(qbart.kvslice_from_size({{'*', i}}, {{'*', 1}}));

      // Construct a vector with the desired contractions
      std::vector<std::vector<int>> disps;
      disps.push_back(std::vector<int>()); // no displacement
      for (int i = 1; i <= max_path_length; ++i)
	disps.push_back(std::vector<int>(i, 3)); // displacements on positive z-dir
      for (int i = 1; i <= max_path_length; ++i)
	disps.push_back(std::vector<int>(i, -3)); // displacements on negative z-dir

      // Put all gamma matrices in a single tensor
      std::vector<SB::Tensor<2, SB::Complex>> gamma_mats;
      {
	for (int g = 0; g < Ns * Ns; ++g)
	{
	  SpinMatrix gmat = Gamma(g) * SB::SpinMatrixIdentity();
	  gamma_mats.push_back(SB::asTensorView(gmat).cloneOn<SB::Complex>(SB::OnDefaultDevice));
	}
      }

      // Contract S and Q with all the gammas, and apply the displacements
      std::string order_out = "gmNndsqt*";
      std::pair<SB::Tensor<9, SB::Complex>, std::vector<int>> r =
	SB::doMomGammaDisp_contractions<9>(u, std::move(qbart), std::move(qt), 0, p, 0, SB::none,
					   gamma_mats, disps, false, order_out);

      // Gather all traces at the master node
      SB::Tensor<9, SB::Complex> con =
	r.first.make_sure(SB::none, SB::OnHost, SB::OnMaster).getLocal();

      const std::vector<int>& disps_perm = r.second;

      // Do the update only on the master node
      if (con)
      {
	std::pair<KeyOperator_t, ValOperator_t> kv;
	kv.first.mom.resize(Nd - 1);
	kv.second.op.resize(Ns * Ns);
	for (int i = 0; i < Ns * Ns; ++i)
	  kv.second.op[i] = 0.0;

	for (int d = 0; d < disps_perm.size(); ++d)
	{
	  // Normalize paths
	  int disp_d_len = disps[disps_perm[d]].size();
	  kv.first.disp.resize(std::max(disp_d_len, 1));
	  for (int i = 0; i < disp_d_len; ++i)
	    kv.first.disp[i] = disps[disps_perm[d]][i];
	  if (disp_d_len == 0)
	    kv.first.disp[0] = 0;

	  for (int mom = 0; mom < p.numMom(); ++mom)
	  {
	    for (int i = 0; i < Nd - 1; ++i)
	      kv.first.mom[i] = p.numToMom(mom)[i];

	    for (int t = 0; t < Nt; ++t)
	    {
	      kv.first.t_slice = t;

	      auto it = db.find(kv.first);
	      if (it == db.end())
		it = db.insert(kv).first;

	      for (int ai = 0; ai < a; ++ai)
	      {
		for (int g = 0; g < Ns * Ns; ++g)
		{
		  std::complex<double> a = con.get({g, mom, 1, 1, d, 1, 1, t, ai});
#ifdef QDP_IS_QDPJIT
		  it->second.op[g].elem().elem().elem().real().elem() += a.real();
		  it->second.op[g].elem().elem().elem().imag().elem() += a.imag();
#else
		  it->second.op[g].elem().elem().elem().real() += a.real();
		  it->second.op[g].elem().elem().elem().imag() += a.imag();
#endif
		}
	      }
	    }
	  }
	}
      }
    }

    // Update the mean and var for each observable in db
    void do_update(std::map< KeyOperator_t, ValOperator_t >& dbmean,
	          std::map< KeyOperator_t, ValOperator_t >& dbvar,
	          const std::map< KeyOperator_t, ValOperator_t >& db, bool first_it)
    {
      for(std::map< KeyOperator_t, ValOperator_t >::const_iterator it=db.begin(); it != db.end(); it++) {
        std::pair<std::map< KeyOperator_t, ValOperator_t >::iterator, bool> itbo;
        // Insert mean
        itbo = dbmean.insert(*it);
	assert(itbo.second == first_it);
        if(!itbo.second ){
          // if insert fails, key already exists, so add result
           itbo.first->second.op += it->second.op;
        }

        // Insert variance
        std::pair<KeyOperator_t, ValOperator_t> kv; 
        kv.first = it->first;
        kv.second.op.resize(it->second.op.size());
        for(int i(0); i<it->second.op.size(); i++)
          kv.second.op[i] = it->second.op[i] * conj(it->second.op[i]);
        itbo = dbvar.insert(kv);
	assert(itbo.second == first_it);
        if(!itbo.second ){
          // if insert fails, key already exists, so add result
          itbo.first->second.op += kv.second.op;
        }
      }
    }

    void show_stats(const std::map< KeyOperator_t, ValOperator_t >& dbmean,
	          const std::map< KeyOperator_t, ValOperator_t >& dbvar,
	          const std::map< KeyOperator_t, ValOperator_t >& dbdet,
                  unsigned int hadamard_normalization, unsigned num_noise)
    {
      if (num_noise <= 1) return;

      const int Nt = Layout::lattSize()[3];
      std::map<KeyOperator_t, std::vector<double>> dbmean_avg, dbdet_avg, dbvar_avg;
      std::map<KeyOperator_t, unsigned int> avg_n; // number of averaged values
      for(std::map< KeyOperator_t, ValOperator_t >::const_iterator it=dbvar.cbegin();it != dbvar.cend(); it++){
        // Average over t_slice and forward/backward directions
        std::pair<KeyOperator_t, std::vector<double>> kv;
        kv.first = it->first;
        kv.first.t_slice = 0;
        for(int k=0;k<it->first.disp.size();k++) kv.first.disp[k] = abs(kv.first.disp[k]);

        // Update dbvar_avg
        // Compute the variance as E[x^2] - E[x]^2
        std::map< KeyOperator_t, ValOperator_t >::const_iterator itmean = dbmean.find(it->first);
        assert(itmean != dbmean.cend());
        kv.second.resize(it->second.op.size());
        for(int i(0);i<it->second.op.size();i++) {
          DComplex a = it->second.op[i] / num_noise - itmean->second.op[i] * conj(itmean->second.op[i]) / num_noise / num_noise / hadamard_normalization / hadamard_normalization;
          kv.second[i] = detox(a.elem().elem().elem().real());
        }
        std::pair<std::map< KeyOperator_t, std::vector<double> >::iterator, bool> itbo = dbvar_avg.insert(kv);
        if(itbo.second ){
          avg_n[kv.first] = 1;
        } else {
          // if insert fails, key already exists, so add result
          for(int i(0);i<kv.second.size();i++) itbo.first->second[i] += kv.second[i];
          avg_n[kv.first]++;
        }

        // Update dbmean_avg
        for(int i(0);i<it->second.op.size();i++) kv.second[i] = abs(std::complex<double>(detox(itmean->second.op[i].elem().elem().elem().real()), detox(itmean->second.op[i].elem().elem().elem().imag()))) / hadamard_normalization / num_noise;
        itbo = dbmean_avg.insert(kv);
        if(!itbo.second){
          for(int i(0);i<it->second.op.size();i++) itbo.first->second[i] += kv.second[i];
        }

        // Update dbdet_avg
        itmean = dbdet.find(it->first);
        if (itmean != dbdet.cend()) {
          for(int i(0);i<it->second.op.size();i++) kv.second[i] = abs(std::complex<double>(detox(itmean->second.op[i].elem().elem().elem().real()), detox(itmean->second.op[i].elem().elem().elem().imag())));
        } else {
          for(int i(0);i<it->second.op.size();i++) kv.second[i] = 0.0;
        }
        itbo = dbdet_avg.insert(kv);
        if(!itbo.second){
          for(int i(0);i<it->second.op.size();i++) itbo.first->second[i] += kv.second[i];
        }
      }
      for(std::map< KeyOperator_t, std::vector<double>>::iterator it=dbvar_avg.begin();it != dbvar_avg.end(); it++) {
        const unsigned int n = avg_n[it->first];
        QDPIO::cout << "DISCO VARIANCE with " << num_noise << " noise vectors key: disp = " << it->first.disp << " mom = " << it->first.mom << "   val: " << std::endl;
        for(int i(0);i<it->second.size();i++) QDPIO::cout << "Gamma[" << i << "]: avg_det = " << dbdet_avg[it->first][i]/n << " avg = " << dbmean_avg[it->first][i]/n << "   var = " << it->second[i]/n << std::endl;
      }
    }
 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }
      
    const std::string name = "DISCO_PROBING_DEFLATION_SUPERB";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= WilsonTypeFermActsEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }


    //----------------------------------------------------------------------------
    // Param stuff
    Params::Params() { frequency = 0; }

    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      try 
      {
	XMLReader paramtop(xml_in, path);

	if (paramtop.count("Frequency") == 1)
	  read(paramtop, "Frequency", frequency);
	else
	  frequency = 1;

	// Parameters for source construction
	read(paramtop, "Param", param);

	// Read in the output propagator/source configuration info
	read(paramtop, "NamedObject", named_obj);

	// Possible alternate XML file pattern
	if (paramtop.count("xml_file") != 0) 
	{
	  read(paramtop, "xml_file", xml_file);
	}
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
	QDP_abort(1);
      }
    }


    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	std::string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "HadaDisco4D");
	write(xml_out, "update_no", update_no);
	write(xml_out, "xml_file", xml_file);
	pop(xml_out);

	XMLFileWriter xml(xml_file);
	func(update_no, xml);
      }
      else
      {
	func(update_no, xml_out);
      }
    }


    // Real work done here
    void 
    InlineMeas::func(unsigned long update_no,
		     XMLWriter& xml_out) 
    {
      typedef LatticeFermion               T;
      typedef multi1d<LatticeColorMatrix>  P;
      typedef multi1d<LatticeColorMatrix>  Q;

      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      // Test and grab a reference to the gauge field
      multi1d<LatticeColorMatrix> u;
      XMLBufferWriter gauge_xml;
      try
      {
	u = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
	TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
	QDP_abort(1);
      }
      catch (const std::string& e) 
      {
	QDPIO::cerr << name << ": std::map call failed: " << e << std::endl;
	QDP_abort(1);
      }

      push(xml_out, "HadamardProp");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": disconnected diagram calculation" << std::endl;

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      write(xml_out, "Input", params);

      // Write out the config header
      write(xml_out, "Config_info", gauge_xml);

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      // Calculate some gauge invariant observables just for info.
      MesPlq(xml_out, "Observables", u);

      std::shared_ptr<Coloring> coloring;
      if (!params.param.probing_file.empty()) {
        QDPIO::cout << "Reading colors from file " << params.param.probing_file << std::endl;
        coloring.reset(new Coloring(params.param.probing_file));
      } else {
	QDPIO::cout << "Generating a " << params.param.probing_distance
		    << "-distance coloring with a power " << params.param.probing_power
		    << std::endl;
	coloring.reset(new Coloring(params.param.probing_distance, params.param.probing_power));
      }
    
      //
      // Initialize fermion action
      //
      std::istringstream  xml_s(params.param.prop.fermact.xml);
      XMLReader  fermacttop(xml_s);
      QDPIO::cout << "FermAct = " << params.param.prop.fermact.id << std::endl;

      Handle< FermionAction<T,P,Q> >
	      S_f(TheFermionActionFactory::Instance().createObject(params.param.prop.fermact.id,
							       fermacttop,
							       params.param.prop.fermact.path));

      Handle< FermState<T,P,Q> > state(S_f->createState(u));

      Handle< SystemSolver<LatticeFermion> > PP = S_f->qprop(state,
							       params.param.prop.invParam);
      Handle< Projector<LatticeFermion> > proj = S_f->projector(state, params.param.projParam); 

      // Initialize the slow Fourier transform phases
      int decay_dir           = Nd-1 ; // hadamard needs this for now
      //Initialize ft differently based on momentum list or max value.
      SftMom ft = params.param.use_p_list ? SftMom(params.param.p_list, decay_dir) : SftMom(params.param.p2_max, false, decay_dir);

      // number of colors
      int Nsrc = coloring->numColors();
      QDPIO::cout << "num colors " << Nsrc << std::endl;

      DComplex tr = 0.0 ;
      DComplex trDef = 0.0 ;

      StopWatch swatch;
      swatch.start();

      // Do the projector part of the trace
      // Loop over the U and V vectors
      QDPIO::cout<<"Now computing the projector contribution"<<std::endl;
      StopWatch swatch_det;
      swatch_det.start();
      std::map< KeyOperator_t, ValOperator_t > dbdet;
      for (int k = 0 ; k < proj->rank() ; k++) {
        // collect dk pairs of vectors
        LatticeFermion vi_lambda, // = v[i]/(u[i]'*Dslash*v[i])
                       ui, vi;     // = u[i], v[i]
	proj->V(k,vi);
        DComplex lambda;
        proj->lambda(k, lambda);
        vi_lambda = vi / lambda;
	proj->U(k,ui);

	std::vector<std::shared_ptr<LatticeFermion>> vi_lambda_sh(
	  1, std::shared_ptr<LatticeFermion>(&vi_lambda, [](LatticeFermion*) {}));
	std::vector<std::shared_ptr<LatticeFermion>> ui_sh(
	  1, std::shared_ptr<LatticeFermion>(&ui, [](LatticeFermion*) {}));
	do_disco(dbdet, ui_sh, vi_lambda_sh, ft,
		 params.param.use_ferm_state_links ? state->getLinks() : u,
		 params.param.max_path_length);
      }
      swatch_det.stop();
      QDPIO::cout << "Projector contribution computed in time= " << swatch_det.getTimeInSeconds() << " secs" << std::endl;
 

      // Loop over the source color and spin, creating the source
      // and calling the relevant propagator routines.
      std::map< KeyOperator_t, ValOperator_t > dbmean, dbvar;
      for (int noise = 0 ; noise < params.param.noise_vectors; noise++) {
        std::map< KeyOperator_t, ValOperator_t > db;

        // doing a new noise vector
        QDPIO::cout << " Doing noise vector " << noise  << std::endl; 

	//generate a random std::vector
	LatticeComplex vec ;
	LatticeReal rnd1, theta;
	random(rnd1); 
	Real twopiN = Chroma::twopi / 4; 
        theta = twopiN * floor(4*rnd1);
	vec = cmplx(cos(theta),sin(theta));

        // All the loops
        const int N_rhs = (params.param.max_rhs + Ns * Nc - 1) / Ns / Nc;
        for (int k1 = 0, dk = std::min(Nsrc, N_rhs); k1 < Nsrc ; k1 += dk, dk = std::min(Nsrc - k1, N_rhs)) {
          // collect (Ns*Nc*dk) pairs of vectors
          std::vector<std::shared_ptr<LatticeFermion>> v_chi(Ns * Nc * dk), v_psi(Ns * Nc * dk), v_q(Ns * Nc * dk);
          for (int col=0; col<v_chi.size(); col++) v_chi[col].reset(new LatticeFermion);
          for (int col=0; col<v_psi.size(); col++) v_psi[col].reset(new LatticeFermion);
          for (int col=0; col<v_q.size(); col++) v_q[col].reset(new LatticeFermion);
          for (int i_v = 0 ; i_v < dk ; i_v++) {
            LatticeInteger hh ; 
            coloring->getVec(hh, k1 + i_v);
            LatticeComplex rv = vec*hh;
            for(int color_source(0);color_source<Nc;color_source++){
              LatticeColorVector vec_srce = zero ;
              pokeColor(vec_srce,rv,color_source) ;
              
              for(int spin_source=0; spin_source < Ns; ++spin_source){
                // Insert a ColorVector into spin index spin_source
                // This only overwrites sections, so need to initialize first
                *v_chi[i_v * Ns * Nc + color_source * Ns + spin_source]  = zero;
                CvToFerm(vec_srce, *v_chi[i_v * Ns * Nc + color_source * Ns + spin_source], spin_source);
                *v_psi[i_v * Ns * Nc + color_source * Ns + spin_source]  = zero;
              } 
            }
          }

          (*PP)(v_psi, std::vector<std::shared_ptr<const LatticeFermion>>(v_chi.begin(), v_chi.end()));
          proj->VUAObliqueProjector(v_q, std::vector<std::shared_ptr<const LatticeFermion>>(v_psi.begin(), v_psi.end()));
          for (int i=0; i<v_psi.size(); ++i)
            *v_q[i] = *v_psi[i] - *v_q[i]; // q <= (I - V*inv(U'*AV*)*U'*A)*quark_soln

          // here the recursive call goes to compute the loops
          // result is ADDED to db
          StopWatch swatch_dots;
          swatch_dots.start();
	  do_disco(db, v_chi, v_q, ft, params.param.use_ferm_state_links ? state->getLinks() : u,
		   params.param.max_path_length);
	  swatch_dots.stop();
	  QDPIO::cout << "Computing inner products " << swatch_dots.getTimeInSeconds() << " secs"
		      << std::endl;
        } // for k1

        // Update dbmean, dbvar
        do_update(dbmean, dbvar, db, noise == 0);

        // Show stats
        show_stats(dbmean, dbvar, dbdet, 1, noise+1);
      } // noise

      // Normalize the traces
      for(std::map< KeyOperator_t, ValOperator_t >::iterator it=dbmean.begin();it != dbmean.end(); it++){
        for(int k=0;k<it->second.op.size();k++){
          it->second.op[k] = it->second.op[k]/toDouble(params.param.noise_vectors);
        }
      }

      // Add the deterministic part to the traces
      for(std::map< KeyOperator_t, ValOperator_t >::iterator it=dbmean.begin();it != dbmean.end(); it++)
        it->second.op += dbdet[it->first].op;

      swatch.stop();
      QDPIO::cout << "Traces were  computed: time= " 
		  << swatch.getTimeInSeconds() 
		  << " secs" << std::endl;

      // write out the results

      if (Layout::nodeNumber() == 0)
      {
	// DB storage
	LocalBinaryStoreDB<LocalSerialDBKey<KeyOperator_t>, LocalSerialDBData<ValOperator_t>>
	  qdp_db;

	// Open the file, and write the meta-data and the binary for this operator
	XMLBufferWriter file_xml;

	push(file_xml, "DBMetaData");
	write(file_xml, "id", std::string("DiscoBlocks"));
	write(file_xml, "lattSize", QDP::Layout::lattSize());
	write(file_xml, "decay_dir", decay_dir);
	write(file_xml, "Params", params.param);
	write(file_xml, "Config_info", gauge_xml);
	pop(file_xml);

	std::string file_str(file_xml.str());
	qdp_db.setMaxUserInfoLen(file_str.size());

	//Slightly modify code to account for changes from multifile write.
	//Be consistent with old mode of filename write.
	std::string file_name = params.named_obj.sdb_file;
	qdp_db.open(file_name, O_RDWR | O_CREAT, 0664);

	qdp_db.insertUserdata(file_str);

	LocalSerialDBKey<KeyOperator_t> key;
	LocalSerialDBData<ValOperator_t> val;
	std::map<KeyOperator_t, ValOperator_t>::iterator it;
	// Store all the data
	for (it = dbmean.begin(); it != dbmean.end(); it++)
	{
	  key.key() = it->first;
	  key.key().mass_label = params.param.mass_label;
	  val.data().op.resize(it->second.op.size());
	  for (int i(0); i < it->second.op.size(); i++)
	    val.data().op[i] = it->second.op[i];
	  qdp_db.insert(key, val);
	}

	qdp_db.close();
      }

      pop(xml_out);  // close last tag

      snoop.stop();
      QDPIO::cout << name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << std::endl;
      
      QDPIO::cout << name << ": ran successfully" << std::endl;
      
      END_CODE(); 
    }

  }// namespace

} // namespace Chroma
// vim: sw=2 sts=2

#endif
