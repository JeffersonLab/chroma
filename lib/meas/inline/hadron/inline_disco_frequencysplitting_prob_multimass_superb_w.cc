/*! \file
 * \brief Compute the disconnected diagrams with 4D probing and Frequency Splitting
 * \Travis Whyte -> thwhyte@wm.edu so you know who to blame
 *
 * Propagator calculation on a colorstd::vector
 */

#include "inline_disco_frequencysplitting_prob_multimass_superb_w.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "actions/ferm/linop/eoprec_clover_linop_w.h"
#include "fermact.h"
#include "meas/glue/mesplq.h"
#include "meas/hadron/greedy_coloring.h"
#include "meas/hadron/interpolation.h"
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
#include "util/ferm/mgproton.h"
#include "util/ferm/transf.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"


#include <cassert>
#include <string>
#include <complex>
#include <iomanip>

//using namespace Interp;

namespace Chroma 
{ 
  namespace InlineDiscoFreqSplitProbMMSuperb 
  {

    typedef std::vector<std::shared_ptr<LatticeFermion>> X1dvector;
    typedef std::vector<X1dvector> X2dmatrix;    


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, InlineDiscoFreqSplitProbMMSuperb::Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "sdb_file", input.sdb_file);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const InlineDiscoFreqSplitProbMMSuperb::Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "sdb_file", input.sdb_file);

      pop(xml);
    }

    //! Propagator input
    void read(XMLReader& xml, const std::string& path, InlineDiscoFreqSplitProbMMSuperb::Params::Param_t& param)
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
      } else {
	QDPIO::cout<<"Could not find valid XML momentum input."<<std::endl;
	QDP_abort(1);
      }
      read(inputtop,"mass_label",param.mass_label);

      read(inputtop,"Propagator",param.prop);
      
      if(inputtop.count("use_ferm_state_links")!=0){ 
	read(inputtop,"use_ferm_state_links",param.use_ferm_state_links) ;
        QDPIO::cout<<"Ferm state links set equal to "<<param.use_ferm_state_links<<std::endl;
      } else {
	param.use_ferm_state_links=false ;
      }


      if(inputtop.count("probing_distance")!=0){ 
	read(inputtop,"probing_distance",param.probing_distance) ;
      } else {
	param.probing_distance = param.max_path_length;
      }

      if(inputtop.count("probing_power")!=0){
	read(inputtop,"probing_power",param.probing_power);
      }else{
	param.probing_power = 2;
      }

      if(inputtop.count("probing_file")!=0){ 
	read(inputtop,"probing_file",param.probing_file) ;
      }

      if(inputtop.count("num_samples")!=0){
	read(inputtop,"num_samples",param.num_samples);
	}else{
	//Default value of num_samples for variance estimation for now...
	param.num_samples = 5;
      }
      if(inputtop.count("max_rhs")!=0){ 
	read(inputtop,"max_rhs",param.max_rhs) ;
      } else {
	param.max_rhs = 1;
     }

      if(inputtop.count("max_rhs_interp")!=0){
        read(inputtop,"max_rhs_interp",param.max_rhs_interp) ;
      } else {
        param.max_rhs_interp = 1;
     }


     if(inputtop.count("use_interpolation")!=0){
	read(inputtop,"use_interpolation",param.use_interpolation);
	}else{
	param.use_interpolation = false;
      }

    if(inputtop.count("num_shifts")!=0){
      read(inputtop,"num_shifts",param.num_shifts);
    }else{
      //will be modified if you did not declare shifts
      param.num_shifts = 0;
    }

     if(inputtop.count("Shifts")!=0){
	read(inputtop,"Shifts",param.shifts);
	}else{
	QDPIO::cout << "No shifts supplied. Setting default shifts and performing interpolation." << std::endl;
	//if you don't have shifts, I am assuming you don't have a reasonable guess/know good shifts
	//so I will set some shifts and then do the interpolation. This will overwrite
	//the param.use_interpolation flag from false to true. Default values of shifts
	//will be 0, .003, .03, .3, 1
	param.shifts.resize(5);
	param.shifts[0] = 0.0;
	param.shifts[1] = 0.003;
	param.shifts[2] = 0.03;
	param.shifts[3] = 0.3;
	param.shifts[4] = 1.0;
	param.use_interpolation = true;
	param.num_shifts = param.shifts.size();
     }

    if (inputtop.count("use_mg")!=0){
        read(inputtop,"use_mg",param.use_mg);
    }else{
        param.use_mg = false;
    }


      if(inputtop.count("noise_vectors")!=0){
        read(inputtop,"noise_vectors",param.noise_vectors) ;
      } else {
	param.noise_vectors.resize(param.shifts.size()+1);
	for (int i = 0; i < param.noise_vectors.size(); i++){
        param.noise_vectors[i] = 1;
	}
      }

     if(inputtop.count("del_s")!=0){
       read(inputtop,"del_s",param.del_s);
       }else{
       //default shift discretization for now
       //you get logspace shifts in regions
       //10^-4 -> 10^-3, 10^-3 -> 10^-2, 10^-2 -> 1
       param.del_s.resize(4);
       param.del_s[0] = -4.0; param.del_s[1] = -3.0; param.del_s[2] = -2.0; 
       param.del_s[3] = 0.0;
      }

      if (inputtop.count("num_bcshifts")!=0){
	read(inputtop,"num_bcshifts",param.num_bcshifts);
      }else{
	param.num_bcshifts.resize(param.del_s.size()-1);
	//you get 20 shifts in each interval
	for (int i = 0; i < param.num_bcshifts.size(); i++){
	    param.num_bcshifts[i] = 20;
	}
      }

	
     if(inputtop.count("display_level_stats")!=0){
	read(inputtop,"display_level_stats",param.display_level_stats);
	}else{
	param.display_level_stats = false;
     }

     if(inputtop.count("use_hpe")!=0){
       read(inputtop,"use_hpe",param.use_hpe);
       }else{
       param.use_hpe = false;
     }

     if(inputtop.count("hpe_power")!=0){
       read(inputtop,"hpe_power",param.hpe_power);
       }else{
       //default value
       param.hpe_power = 4;
     }

      if(inputtop.count("hpe_probing_files")!=0){
        std::string tmp;
        std::string delim = " ";
        read(inputtop,"hpe_probing_files",tmp) ;
        param.hpe_probing_files.resize(param.max_path_length+1);
        size_t pos = 0;
        std::string tmp_file_str;
        int i = 0;
        while ( (pos = tmp.find(delim)) != std::string::npos) {
            tmp_file_str = tmp.substr(0, pos);
            param.hpe_probing_files[i] = tmp_file_str;
            QDPIO::cout << "Parsed hpe probing file " << tmp_file_str << std::endl;
            tmp.erase(0, pos + delim.length());
            i+=1;
        }
         QDPIO::cout << "Parsed hpe probing file " << tmp << std::endl;
         param.hpe_probing_files[i] = tmp;
      }else{
	param.hpe_probing_files.resize(param.max_path_length+1);
	for (int i = 0; i < param.hpe_probing_files.size(); i++){
	    param.hpe_probing_files[i] = "nofile";
	}
      }

      if(inputtop.count("disp_gamma")!=0){
	read(inputtop,"disp_gamma",param.disp_gamma);
	//it will be determined if use_interpolation == true
	//and it won't matter if use_interpolation == false
      }
  
      if (inputtop.count("do_trace_est")!=0){
	 read(inputtop,"do_trace_est",param.do_trace_est);
      }
  

  }
    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const InlineDiscoFreqSplitProbMMSuperb::Params::Param_t& param)
    {
      push(xml, path);

      write(xml,"max_path_length",param.max_path_length);
      write(xml,"p2_max",param.p2_max);
      write(xml,"mass_label",param.mass_label);
      write(xml,"Propagator",param.prop);
      write(xml,"use_ferm_state_links",param.use_ferm_state_links) ;
      write(xml,"probing_distance",param.probing_distance) ;
      write(xml,"probing_power",param.probing_power);
      write(xml,"noise_vectors",param.noise_vectors) ;
      write(xml,"max_rhs",param.max_rhs) ;
      write(xml,"use_interpolation",param.use_interpolation);
      write(xml,"Shifts",param.shifts);
      write(xml,"num_samples",param.num_samples);
      write(xml,"use_hpe",param.use_hpe);
      

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, InlineDiscoFreqSplitProbMMSuperb::Params& input)
    {
      InlineDiscoFreqSplitProbMMSuperb::Params tmp(xml, path);
      input = tmp;
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const InlineDiscoFreqSplitProbMMSuperb::Params& input)
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
         << " , disp = ";
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

    struct StatHolder_t {
      std::vector<std::map< KeyOperator_t, ValOperator_t >> levels_var;
      std::vector<std::map< KeyOperator_t, ValOperator_t >> levels_avg;
      std::vector<int> level_cost;
      std::map< KeyOperator_t, ValOperator_t > total_var;
      std::map< KeyOperator_t, ValOperator_t > total_avg;
    };

    typedef std::vector<std::map<KeyOperator_t, ValOperator_t>> VecMap;
    typedef std::vector<std::vector<std::map<KeyOperator_t, ValOperator_t>>> MatMap;

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

    struct Restart_t {

      int level;
      int noise;
      //int prob_int;
      int num_levels;
      //QDP::Seed rng_seed;
      std::vector<double> count = std::vector<double>(2);
      std::vector<double> level_cost; //needs to be resized!!!
      std::vector<double> shifts; //needs to be resized!!!
      std::map<KeyOperator_t, ValOperator_t> db_trace;
      std::map<KeyOperator_t, ValOperator_t> db_var;

    };

    void send_to_nodes(Restart_t restart)
    {
    QDPIO::cin >> restart.level;
    QDPIO::cin >> restart.noise;
    QDPIO::cin >> restart.num_levels;
    QDPIO::cin >> restart.count[0];
    QDPIO::cin >> restart.count[1];
    for (int i = 0; i < restart.level_cost.size(); i++){
	QDPIO::cin >> restart.level_cost[i];
    }


    }

    void restart_seed_out(const std::string& filename, const QDP::Seed& rng_seed)
    {

    QDP::TextFileWriter out(filename);
    out << rng_seed;
    out.close();
    
    }

    QDP::Seed restart_seed_in(const std::string& filename)
    {
    
    QDP::TextFileReader in(filename);
    QDP::Seed tmp;
    in >> tmp;
    in.close();
    return tmp; 
    }

    void assign_KeyOperator(std::string& str, KeyOperator_t& keyop)
    {

      std::vector<std::string> out_string;
      std::string delim = " ";
      size_t pos = 0;
      std::string substring;
      while ((pos = str.find(delim)) != std::string::npos) {
	    substring = str.substr(0, pos);
	    out_string.push_back(substring);
	    str.erase(0, pos+delim.length());
      }
      out_string.push_back(str);
     
      std::string key = "t_slice";
      auto loc = std::find(out_string.begin(), out_string.end(), key);
      keyop.t_slice = std::stoi(out_string[std::distance(out_string.begin(), loc)+2]);
      key = "disp";
      loc = std::find(out_string.begin(), out_string.end(), key);
      key = "mom";
      auto loc2 = std::find(out_string.begin(), out_string.end(), key);
      std::vector<int> tmpdisp;
      for (int i = std::distance(out_string.begin(), loc)+2; i <= std::distance(out_string.begin(), loc2)-2; i++){
	  tmpdisp.push_back(std::stoi(out_string[i]));
      }
      keyop.disp.resize(tmpdisp.size());
      for (int i = 0; i < keyop.disp.size(); i++){
	  keyop.disp[i] = tmpdisp[i];
      }
      std::vector<int> tmpmom;
      for (int i = std::distance(out_string.begin(), loc2)+2; i < out_string.size()-1; i++){
	tmpmom.push_back(std::stoi(out_string[i]));
      }
      keyop.mom.resize(tmpmom.size());
      for (int i = 0; i < keyop.mom.size(); i++){
	keyop.mom[i] = tmpmom[i];
      }
    }

    Restart_t checkpoint_in(const std::string& tracefile, const std::string& varfile, const int& max_path_length)

    {
    
    const int Nt = Layout::lattSize()[3];
    Restart_t restart;
    std::string line_t;
    std::string line_v;
    std::ifstream file_trace;
    file_trace.open(tracefile, std::ios::in);
    std::ifstream file_var;
    file_var.open(varfile, std::ios::in);
    if (!file_trace.is_open()){
       QDPIO::cout << "Could not open trace checkpoint file...aborting" << std::endl;
       QDP_abort(1);
    }
    if (!file_var.is_open()){
       QDPIO::cout << "Could not open variance checkpoint file...aborting" << std::endl;
       QDP_abort(1);
    }
    getline(file_trace,line_t);
    restart.level = stoi(line_t);
    getline(file_trace,line_t);
    restart.noise = stoi(line_t);
    getline(file_trace,line_t);
    restart.num_levels = stoi(line_t);
    getline(file_trace,line_t);
    restart.count[0] = stod(line_t);
    getline(file_trace,line_t);
    restart.count[1] = stod(line_t);
    restart.level_cost.resize(restart.num_levels);
    restart.shifts.resize(restart.num_levels);
    for (int i = 0; i < restart.level_cost.size(); i++){
	getline(file_trace,line_t);
	restart.level_cost[i] = stod(line_t);
    }
    for (int i = 0; i < restart.shifts.size(); i++){
	getline(file_trace,line_t);
	restart.shifts[i] = stod(line_t);
    }
    QDPIO::cout << "Restart params read in successfully" << std::endl;
    for (int i = 0; i < (2*max_path_length+1)*Nt; i++){
	std::pair<KeyOperator_t, ValOperator_t> kv ;
	std::pair<KeyOperator_t, ValOperator_t> kvv ;
	getline(file_trace,line_t);
	getline(file_var,line_v);
	assign_KeyOperator(line_t, kv.first);
	assign_KeyOperator(line_v, kvv.first);
	getline(file_trace,line_t);
	getline(file_var,line_v);
        auto it = restart.db_trace.find(kv.first);
        if (it == restart.db_trace.end()){ it = restart.db_trace.insert(kv).first; }
        auto itv = restart.db_var.find(kvv.first);
        if (itv == restart.db_var.end()) {itv = restart.db_var.insert(kvv).first;}

	for (int g = 0; g < Ns * Ns; g++){
	    getline(file_trace,line_t);
	    getline(file_var,line_v);
#ifdef QDP_IS_QDPJIT
	    it->second.op[g].elem().elem().elem().real().elem() = stod(line_t);
	    itv->second.op[g].elem().elem().elem().real().elem() = stod(line_v);
	    getline(file_trace,line_t);
	    getline(file_var,line_v);
	    it->second.op[g].elem().elem().elem().imag().elem() = stod(line_t);
	    itv->second.op[g].elem().elem().elem().imag().elem() = stod(line_v);
#else
	    it->second.op[g].elem().elem().elem().real() = stod(line_t);
	    itv->second.op[g].elem().elem().elem().real() = stod(line_v);
	    getline(file_trace,line_t);
	    getline(file_var,line_v);
	    it->second.op[g].elem().elem().elem().imag() = stod(line_t);
	    itv->second.op[g].elem().elem().elem().imag() = stod(line_v);
#endif
	}
	
      }
      QDPIO::cout << "Traces and variances read in successfully" << std::endl;

    
    


    //} //nodelist

    return restart;
    } //func

    void checkpoint_out(const int& level, const int& noise, const QDP::Seed& rng_seed,
			const int& num_levels, const std::vector<double> level_cost,
			const std::vector<double>& count, multi1d<Real>& shifts, std::map<KeyOperator_t, ValOperator_t>& db_trace,
      			const std::map<KeyOperator_t, ValOperator_t>& db_var, const std::string& tracefile, const std::string& varfile)
    {

    int nodelist = Layout::nodeNumber();
    //only master node writes to file
    if (nodelist == 0){
    int prec = 17; //for output to file
    std::ofstream file_trace;
    std::ofstream file_var;
    file_trace.open(tracefile, std::ios::out);
    file_var.open(varfile, std::ios::out);
    if (file_trace.is_open() && file_var.is_open()){
        file_trace << level << '\n';
	file_trace << noise << '\n';
	file_trace << num_levels << '\n';
	for (auto i : count){ file_trace << i << '\n'; }
	for (auto i : level_cost){ file_trace << i << '\n'; }
	for (int i = 0; i < shifts.size(); i++) { file_trace << std::setprecision(std::numeric_limits<double>::max_digits10) << shifts[i] << std::endl; }
	std::map< KeyOperator_t, ValOperator_t >::const_iterator itv = db_var.begin();
	for (std::map< KeyOperator_t, ValOperator_t >::const_iterator it=db_trace.begin(); it != db_trace.end(); it++) {
	    file_trace << it->first << '\n';
	    file_var << itv->first << '\n';
	    for (int g = 0; g < Ns * Ns; g++){
		DComplex a;
		DComplex b;
#ifdef QDP_IS_QDPJIT
		a.elem().elem().elem().real() = it->second.op[g].elem().elem().elem().real().elem();
		a.elem().elem().elem().imag() = it->second.op[g].elem().elem().elem().imag().elem();
		//file_trace << it->second.op[g].elem().elem().elem().real().elem() << '\n';
		file_trace << std::fixed << std::setprecision(prec) << detox(a.elem().elem().elem().real())<< '\n';
		//file_trace << it->second.op[g].elem().elem().elem().imag().elem() << '\n';
		file_trace << std::fixed << std::setprecision(prec) << detox(a.elem().elem().elem().imag()) << '\n';
		b.elem().elem().elem().real() = itv->second.op[g].elem().elem().elem().real().elem();
		b.elem().elem().elem().imag() = itv->second.op[g].elem().elem().elem().imag().elem();
		//file_var << itv->second.op[g].elem().elem().elem().real().elem() << '\n';
		file_var << std::fixed << std::setprecision(prec) << detox(b.elem().elem().elem().real()) << '\n';
		//file_var << itv->second.op[g].elem().elem().elem().imag().elem() << '\n';
		file_var << std::fixed << std::setprecision(prec) << detox(b.elem().elem().elem().imag()) << '\n';
		
#else
		a.elem().elem().elem().real() = it->second.op[g].elem().elem().elem().real();
		a.elem().elem().elem().imag() = it->second.op[g].elem().elem().elem().imag();
		file_trace << std::fixed << std::setprecision(prec) << detox(a.elem().elem().elem().real()) << '\n';
		file_trace << std::fixed << std::setprecision(prec) << detox(a.elem().elem().elem().imag()) << '\n';
		//file_trace << it->second.op[g].elem().elem().elem().real() << '\n';
		//file_trace << it->second.op[g].elem().elem().elem().imag() << '\n';
		b.elem().elem().elem().real() = itv->second.op[g].elem().elem().elem().real();
		b.elem().elem().elem().imag() = itv->second.op[g].elem().elem().elem().imag();
		//file_var << itv->second.op[g].elem().elem().elem().real() << '\n';
		//file_var << itv->second.op[g].elem().elem().elem().imag() << '\n';
		file_var << std::fixed << std::setprecision(prec) << detox(b.elem().elem().elem().real()) << '\n';
		file_var << std::fixed << std::setprecision(prec) << detox(b.elem().elem().elem().imag()) << '\n';
#endif
	    } //g
	    itv++;
	} //it

    } //is_open
    file_trace.close();
    file_var.close();
    } //node

    } //func

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

    void do_disco(std::map<KeyOperator_t, ValOperator_t>& db,
                  const std::vector<std::shared_ptr<LatticeFermion>>& qbar,
                  const std::vector<std::shared_ptr<LatticeFermion>>& q, const SftMom& p,
                  const multi1d<LatticeColorMatrix>& u, const int& max_path_length,
				  const DComplex& cshift, const int& dk)
    {

      const int Nt = Layout::lattSize()[3];
      const int a = qbar.size();
      assert(qbar.size() == q.size());
	  std::shared_ptr<LatticeFermion> q_tmp;
	  q_tmp.reset(new LatticeFermion);
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
      //for (int i = 0; i < a; ++i)
      for (int idk = 0; idk < dk; ++idk){
	for (int col = 0; col < Nc * Ns; ++col){ 
	if (col == 2 || col == 3 || col == 6 || col == 7 || col == 10 || col == 11){
		*q_tmp = -1.0 * cshift * (Gamma(15) * *qbar[idk * Nc * Ns + col]);
        SB::asTensorView(*q_tmp)
          .rename_dims({{'s', 'Q'}})
          .copyTo(qbart.kvslice_from_size({{'*', idk * Nc * Ns + col}}, {{'*', 1}}));
	}else{
		*q_tmp = cshift * (Gamma(15) * *qbar[idk * Nc * Ns + col]);
        SB::asTensorView(*q_tmp)
          .rename_dims({{'s', 'Q'}})
          .copyTo(qbart.kvslice_from_size({{'*', idk * Nc * Ns + col}}, {{'*', 1}}));
	}
      }
    }

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

              //for (int ai = 0; ai < a; ++ai)
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


    void do_disco_hop(std::map< KeyOperator_t, ValOperator_t >& db,
      const LatticeFermion& qbar,
      const LatticeFermion& q,
      const SftMom& p,
      const multi1d<LatticeColorMatrix>& u,
      const multi1d<int>& path,
      const int& max_path_length ){
      
      const int Nt = Layout::lattSize()[3];

      ValOperator_t val ;
      KeyOperator_t key ;
      std::pair<KeyOperator_t, ValOperator_t> kv ; 
      if(path.size()==0){
  kv.first.disp.resize(1);
  kv.first.disp[0] = 0 ;
      }
      else
      kv.first.disp = path ;

    multi2d< multi1d<ComplexD> > foo(p.numMom(), Nt) ;
    for (int t(0); t < Nt; t++)
      for (int m(0); m < p.numMom(); m++)
      foo(m,t).resize(Ns*Ns);
      for(int g(0);g<Ns*Ns;g++){
      LatticeComplex cc = localInnerProduct(qbar,Gamma(g)*q);
  for (int m(0); m < p.numMom(); m++){
          LatticeComplex pcc = p[m]*cc;
          for (int t(0); t < Nt; t++) {
        foo(m,t)[g] = sum(pcc,p.getSet()[t]);
    }
  }
      }

      for (int t(0); t < Nt; t++) {
        kv.first.t_slice = t ;
        for (int m(0); m < p.numMom(); m++){
          for(int i(0);i<(Nd-1);i++)
            kv.first.mom[i] = p.numToMom(m)[i] ;
          
          kv.second.op = foo(m,t);
          std::pair<std::map< KeyOperator_t, ValOperator_t >::iterator, bool> itbo;

          itbo = db.insert(kv);
          if(!itbo.second ){

            for(int i(0);i<kv.second.op.size();i++){
              itbo.first->second.op[i] += kv.second.op[i] ;
            }
          }
        }
      }

    }

    void do_disco_hop(std::map<KeyOperator_t, ValOperator_t>& db,
                  const std::vector<std::shared_ptr<LatticeFermion>>& qbar,
                  const std::vector<std::shared_ptr<LatticeFermion>>& q, const SftMom& p,
                  const multi1d<LatticeColorMatrix>& u, const int& disp)
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
      //for hpe, restrict to just one displacement
      if (disp == 0){
	disps.push_back(std::vector<int>()); // no displacement
      }else{
	disps.push_back(std::vector<int>(disp, 3));
	disps.push_back(std::vector<int>(disp, -3));
      }

      //disps.push_back(std::vector<int>()); // no displacement
     // for (int i = 1; i <= max_path_length; ++i)
     //   disps.push_back(std::vector<int>(i, 3)); // displacements on positive z-dir
      //for (int i = 1; i <= max_path_length; ++i)
     //   disps.push_back(std::vector<int>(i, -3)); // displacements on negative z-dir

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

	

    std::vector<std::vector<double>> retrieve_variances(const MatMap& dbmean, 
							const MatMap& dbvar,
							unsigned int hadamard_normalization,
							unsigned num_noise, int num_shifts,
							int dk, int g)
    {

    std::vector<std::vector<double>> out_vars(num_shifts, std::vector<double>(num_shifts));
    //out_vars.resize(num_shifts, std::vector<double>(num_shifts));
    for (int ii = 0; ii < num_shifts; ii++){
	for (int j = ii; j < num_shifts; j++){
  
	 std::map<KeyOperator_t, std::vector<double>> dbvar_avg;
        std::map<KeyOperator_t, unsigned int> avg_n; // number of averaged values
        for(std::map< KeyOperator_t, ValOperator_t >::const_iterator it=dbvar[ii][j].cbegin();it != dbvar[ii][j].cend(); it++){
	    // Average over t_slice and forward/backward directions
              std::pair<KeyOperator_t, std::vector<double>> kv;
	      kv.first = it->first;
	      kv.first.t_slice = 0;
	      for(int k=0;k<it->first.disp.size();k++) kv.first.disp[k] = abs(kv.first.disp[k]);
	      std::map< KeyOperator_t, ValOperator_t >::const_iterator itmean = dbmean[ii][j].find(it->first);
	      assert(itmean != dbmean[ii][j].cend());
	      kv.second.resize(it->second.op.size());
	      for(int i(0);i<it->second.op.size();i++) {
	      DComplex a = it->second.op[i] / num_noise - itmean->second.op[i] * conj(itmean->second.op[i]) / num_noise / num_noise / hadamard_normalization / hadamard_normalization;
	      kv.second[i] = detox(a.elem().elem().elem().real());
	       }
	      std::pair<std::map< KeyOperator_t, std::vector<double> >::iterator, bool> itbo = dbvar_avg.insert(kv);
	      if(itbo.second ){
	      avg_n[kv.first] = 1;
	      } else {
	      for(int i(0);i<kv.second.size();i++) itbo.first->second[i] += kv.second[i];
	      avg_n[kv.first]++;
	      }
	  }
	  for(std::map< KeyOperator_t, std::vector<double>>::iterator it=dbvar_avg.begin(); it != dbvar_avg.end(); it++) {
	      const unsigned int n = avg_n[it->first];
	      //if you want 0  displacement variances
	      if (it->first.disp.size() == 1 && it->first.disp[0] == 0 && dk == 0){
		  out_vars[ii][j] = it->second[g]/n;
		  QDPIO::cout << "Variance for Displacement = " << dk << " and Gamma = " << g << " : " << it->second[g]/n << std::endl;
		  break;
	      //if you want the first displacement
	      } else if (it->first.disp.size() == 1 && it->first.disp[0] == 3 && dk == 1){
		  out_vars[ii][j] = it->second[g]/n;
		  QDPIO::cout << "Variance for Displacement = " << dk << " and Gamma = " << g << " : " << it->second[g]/n << std::endl;
		  break;
	      //any other displacement, disp.size() = the displacement you want
	      } else if (it->first.disp.size() == dk && (it->first.disp.size() != 1)){
		  out_vars[ii][j] = it->second[g]/n;
		  QDPIO::cout << "Variance for Displacement = " << dk << " and Gamma = " << g << " : " << it->second[g]/n << std::endl;
		  break;
	      }
	  }	
		  
	}
    }
    
    return out_vars;
    }

    std::vector<double> retrieve_variances(const VecMap& dbmean,
					   const VecMap& dbvar,
                                           unsigned int hadamard_normalization,
                                           unsigned num_noise, int num_shifts,
					   int dk, int g)
    {

    std::vector<double> out_vars(num_shifts);
    //out_vars.resize(num_shifts);

    for (int ii = 0; ii < num_shifts; ii++){
	std::map<KeyOperator_t, std::vector<double>> dbvar_avg;
        std::map<KeyOperator_t, unsigned int> avg_n; // number of averaged values
        for(std::map< KeyOperator_t, ValOperator_t >::const_iterator it=dbvar[ii].cbegin();it != dbvar[ii].cend(); it++){
	    // Average over t_slice and forward/backward directions
              std::pair<KeyOperator_t, std::vector<double>> kv;
              kv.first = it->first;
              kv.first.t_slice = 0;
              for(int k=0;k<it->first.disp.size();k++) kv.first.disp[k] = abs(kv.first.disp[k]);
              std::map< KeyOperator_t, ValOperator_t >::const_iterator itmean = dbmean[ii].find(it->first);
              assert(itmean != dbmean[ii].cend());
              kv.second.resize(it->second.op.size());
              for(int i(0);i<it->second.op.size();i++) {
              DComplex a = it->second.op[i] / num_noise - itmean->second.op[i] * conj(itmean->second.op[i]) / num_noise / num_noise / hadamard_normalization / hadamard_normalization;
              kv.second[i] = detox(a.elem().elem().elem().real());
               }
              std::pair<std::map< KeyOperator_t, std::vector<double> >::iterator, bool> itbo = dbvar_avg.insert(kv);
              if(itbo.second ){
              avg_n[kv.first] = 1;
              } else {
              for(int i(0);i<kv.second.size();i++) itbo.first->second[i] += kv.second[i];
              avg_n[kv.first]++;
              }
          }

          for(std::map< KeyOperator_t, std::vector<double>>::iterator it=dbvar_avg.begin();it != dbvar_avg.end(); it++) {
              const unsigned int n = avg_n[it->first];
              if (it->first.disp.size() == 1 && it->first.disp[0] == 0 && dk == 0){
                  out_vars[ii] = it->second[g]/n;
		  QDPIO::cout << "Variance for Displacement = " << dk << " and Gamma = " << g << " : " << it->second[g]/n << std::endl;
                  break;
              } else if (it->first.disp.size() == 1 && it->first.disp[0] == 3 && dk == 1){
                  out_vars[ii] = it->second[g]/n;
		  QDPIO::cout << "Variance for Displacement = " << dk << " and Gamma = " << g << " : " << it->second[g]/n << std::endl;
                  break;
              } else if (it->first.disp.size() == dk && (it->first.disp.size() != 1)){
                  out_vars[ii] = it->second[g]/n;
		  QDPIO::cout << "Variance for Displacement = " << dk << " and Gamma = " << g << " : " << it->second[g]/n << std::endl;
                  break;
              }
          }
	  
    }

    return out_vars;

    }


    double retrieve_variances(const std::map< KeyOperator_t, ValOperator_t >& dbmean,
			      const std::map< KeyOperator_t, ValOperator_t >& dbvar,
                              unsigned int hadamard_normalization,
                              unsigned num_noise,
                              int dk, int g)
    {

    double out_vars;
    std::map<KeyOperator_t, std::vector<double>> dbvar_avg;
    std::map<KeyOperator_t, unsigned int> avg_n; // number of averaged values
    for(std::map< KeyOperator_t, ValOperator_t >::const_iterator it=dbvar.cbegin();it != dbvar.cend(); it++){
    // Average over t_slice and forward/backward directions
              std::pair<KeyOperator_t, std::vector<double>> kv;
              kv.first = it->first;
              kv.first.t_slice = 0;
              for(int k=0;k<it->first.disp.size();k++) kv.first.disp[k] = abs(kv.first.disp[k]);
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
              for(int i(0);i<kv.second.size();i++) itbo.first->second[i] += kv.second[i];
              avg_n[kv.first]++;
              }
          }

          for(std::map< KeyOperator_t, std::vector<double>>::iterator it=dbvar_avg.begin();it != dbvar_avg.end(); it++) {
              const unsigned int n = avg_n[it->first];
              if (it->first.disp.size() == 1 && it->first.disp[0] == 0 && dk == 0){
                  out_vars = it->second[g]/n;
                  QDPIO::cout << "Variance for Displacement = " << dk << " and Gamma = " << g << " : " << it->second[g]/n << std::endl;
                  break;
              } else if (it->first.disp.size() == 1 && it->first.disp[0] == 3 && dk == 1){
                  out_vars = it->second[g]/n;
                  QDPIO::cout << "Variance for Displacement = " << dk << " and Gamma = " << g << " : " << it->second[g]/n << std::endl;
                  break;
              } else if (it->first.disp.size() == dk && (it->first.disp.size() != 1)){
                  out_vars = it->second[g]/n;
                  QDPIO::cout << "Variance for Displacement = " << dk << " and Gamma = " << g << " : " << it->second[g]/n << std::endl;
                  break;
              }
          }

    

    return out_vars;

    }


    void show_stats(const std::map< KeyOperator_t, ValOperator_t >& dbmean,
                  const std::map< KeyOperator_t, ValOperator_t >& dbvar,
                  const std::map< KeyOperator_t, ValOperator_t >& dbdet,
                  unsigned int hadamard_normalization, unsigned num_noise,
		  const int& level)
    {
      if (num_noise <= 1) return;

      const int Nt = Layout::lattSize()[3];
      std::map<KeyOperator_t, std::vector<double>> dbmean_avg, dbdet_avg, dbvar_avg;
      std::map<KeyOperator_t, unsigned int> avg_n; 
      for(std::map< KeyOperator_t, ValOperator_t >::const_iterator it=dbvar.cbegin();it != dbvar.cend(); it++){
        
        std::pair<KeyOperator_t, std::vector<double>> kv;
        kv.first = it->first;
        kv.first.t_slice = 0;
        for(int k=0;k<it->first.disp.size();k++) kv.first.disp[k] = abs(kv.first.disp[k]);


        std::map< KeyOperator_t, ValOperator_t >::const_iterator itmean = dbmean.find(it->first);
        assert(itmean != dbmean.cend());
        kv.second.resize(it->second.op.size());
        for(int i(0);i<it->second.op.size();i++) {
          DComplex a = it->second.op[i] / num_noise - itmean->second.op[i] * conj(itmean->second.op[i]) / num_noise / num_noise;
          kv.second[i] = detox(a.elem().elem().elem().real()) / hadamard_normalization / hadamard_normalization;
        }
        std::pair<std::map< KeyOperator_t, std::vector<double> >::iterator, bool> itbo = dbvar_avg.insert(kv);
        if(itbo.second ){
          avg_n[kv.first] = 1;
        } else {

          for(int i(0);i<kv.second.size();i++) itbo.first->second[i] += kv.second[i];
          avg_n[kv.first]++;
        }
    

        for(int i(0);i<it->second.op.size();i++) kv.second[i] = abs(std::complex<double>(detox(itmean->second.op[i].elem().elem().elem().real()), detox(itmean->second.op[i].elem().elem().elem().imag()))) / hadamard_normalization / num_noise;
        itbo = dbmean_avg.insert(kv);
        if(!itbo.second){
          for(int i(0);i<it->second.op.size();i++) itbo.first->second[i] += kv.second[i];
        }


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
        QDPIO::cout << "DISCO FS VARIANCE with " << num_noise << " noise vectors for level = " << level << ", key: disp = " << it->first.disp << " mom = " << it->first.mom << "   val: " << std::endl;
        for(int i(0);i<it->second.size();i++) QDPIO::cout << "Gamma[" << i << "]: avg_det = " << dbdet_avg[it->first][i]/n << " avg = " << dbmean_avg[it->first][i]/n << "   var = " << it->second[i]/n << std::endl;
      }
    }

      void modifyMass(std::string& xml_str, const Real& m_q, const Real& shifts)
      {
      
      //delimiter is the opening/close  mass tag
      std::string open_tag = "<Mass>";
      std::string close_tag = "</Mass>";
      auto o_pos = xml_str.find(open_tag);
      auto c_pos = xml_str.find(close_tag);
      if(o_pos == std::string::npos){
        return;
      }else{
      //std::string mq0_str = std::to_string(toDouble(m_q0));
      std::string new_mq = std::to_string(toDouble(m_q + shifts));
      int m_q0_size = c_pos - (o_pos + open_tag.size());
      xml_str.replace(o_pos+open_tag.size(), m_q0_size, new_mq);
      }

      }

    void modifySubspaceId(std::string& xml_str, std::string& modifier)
    {
      //tags for mgproto/qop_mg
      std::string open_tag = "<SubspaceId>";
      std::string close_tag = "</SubspaceId>";
      auto o_pos = xml_str.find(open_tag);
      auto c_pos = xml_str.find(close_tag);
      if(o_pos == std::string::npos){
	//if met, then either QUDA MG or no MG
	open_tag = "<SubspaceID>";
	close_tag = "</SubspaceID>";
	auto o_pos = xml_str.find(open_tag);
	auto c_pos = xml_str.find(close_tag);
	if(o_pos == std::string::npos){
	    return;
	}else{
	    int id_size = c_pos - (o_pos + open_tag.size());
	    std::string new_id = xml_str.substr(o_pos+open_tag.size(), id_size);
	    new_id.append(modifier);
	    xml_str.replace(o_pos+open_tag.size(), id_size, new_id);
	}
      }else{
      int id_size = c_pos - (o_pos + open_tag.size());
      std::string new_id = xml_str.substr(o_pos+open_tag.size(), id_size);
      new_id.append(modifier);
      xml_str.replace(o_pos+open_tag.size(), id_size, new_id);
      }
    } 

    std::string retrieveSubspaceId(std::string& xml_str)
    {

      std::string id = "foo";
      //tags for mgproto/qop_mg
      std::string open_tag = "<SubspaceId>";
      std::string close_tag = "</SubspaceId>";
      auto o_pos = xml_str.find(open_tag);
      auto c_pos = xml_str.find(close_tag);
      if(o_pos == std::string::npos){
	  //either not using MG at all or using QUDA MG
	  //Idk why but the xml tags are different
	  //Below are tags for QUDA MG subspaces
	  //At this point, I am assuming that there are no other MG that
	  //are currently in production use
	  open_tag = "<SubspaceID>";
	  close_tag = "</SubspaceID>";
	  auto o_pos = xml_str.find(open_tag);
	  auto c_pos = xml_str.find(close_tag);
	  if(o_pos == std::string::npos){
	    return id;
	  }else{
            int id_size = c_pos - (o_pos + open_tag.size());
	    id = xml_str.substr(o_pos+open_tag.size(), id_size);
	  }
        //return id;
      }else{
      int id_size = c_pos - (o_pos + open_tag.size());
      id = xml_str.substr(o_pos+open_tag.size(), id_size);
      }
      return id;
    }


    void exactHPE(std::map< KeyOperator_t, ValOperator_t >& db, 
		  const InlineDiscoFreqSplitProbMMSuperb::Params::Param_t& param,
		  const multi1d<LatticeColorMatrix>& u)
    {

      typedef LatticeFermion               T;
      typedef multi1d<LatticeColorMatrix>  P;
      typedef multi1d<LatticeColorMatrix>  Q;

      //multi1d<int> new_path_d;
      //multi1d<int> path_d;

      std::vector<std::shared_ptr<Coloring>> hop_coloring(param.max_path_length+1);
      for (int i = 0; i < param.max_path_length+1; i++){
      if (param.hpe_probing_files[i] != "nofile"){
        QDPIO::cout << "Reading colors from file " << param.hpe_probing_files[i] << std::endl;
        hop_coloring[i].reset(new Coloring(param.hpe_probing_files[i]));
      } else {
        QDPIO::cout << "Generating a " << param.hpe_power - 1 << "-distance coloring for the Hopping Term" << std::endl;
        hop_coloring[i].reset(new Coloring(i, param.hpe_power - 1));
      }
      }
  
      Real m_q0 = getMass(param.prop.fermact);
      QDPIO::cout << "Base mass m_q0 = " << toDouble(m_q0) << std::endl;
      std::string ferm_xml = param.prop.fermact.xml;
      modifyMass(ferm_xml, m_q0, param.shifts[param.shifts.size()-1]);
      QDPIO::cout << "Performing exact HPE contribution with m_q = " << toDouble(m_q0 + param.shifts[param.shifts.size()-1]) << std::endl;
      std::istringstream xml_t(ferm_xml);
      XMLReader fermacttop(xml_t);
      
      Handle< FermionAction<T,P,Q> >
              S_f(TheFermionActionFactory::Instance().createObject(param.prop.fermact.id,
                                                               fermacttop,
                                                               param.prop.fermact.path));
      Handle< FermState<T,P,Q> > state(S_f->createState(u));
    
      EvenOddPrecCloverLinOp D;
      CloverFermActParams clov_params(fermacttop, param.prop.fermact.path);
      D.create(state, clov_params);

      int decay_dir           = Nd-1 ;
      SftMom ft = param.use_p_list ? SftMom(param.p_list, decay_dir) : SftMom(param.p2_max, false, decay_dir);


      for (int disp = 0; disp < param.max_path_length+1; disp++){

      int Nvecs = hop_coloring[disp]->numColors();
      const int N_v = (param.max_rhs + Ns * Nc - 1) / Ns / Nc;
      QDPIO::cout << Nvecs << " colors for the hopping term trace calculation for Displacement " << disp << std::endl;

      QDPIO::cout<<"Now computing the hopping contribution"<<std::endl;

      for (int k1 = 0, dk = std::min(Nvecs, N_v); k1 < Nvecs ; k1 += dk, dk = std::min(Nvecs - k1, N_v)) {
	  //LatticeFermion source, psi, chi, eta, tmp;
	  std::vector<std::shared_ptr<LatticeFermion>> source(Ns * Nc * dk), psi(Ns * Nc * dk);
	  std::shared_ptr<LatticeFermion>  chi, eta, tmp;
	  for (int col = 0; col < source.size(); col++) {source[col].reset(new LatticeFermion);}
	  for (int col = 0; col < psi.size(); col++) {psi[col].reset(new LatticeFermion);}
	  //psi.reset(new LatticeFermion);
	  chi.reset(new LatticeFermion);
	  eta.reset(new LatticeFermion);
	  tmp.reset(new LatticeFermion);
	  for (int i_v = 0 ; i_v < dk ; i_v++) {
              LatticeInteger hh ;
              hop_coloring[disp]->getVec(hh, k1 + i_v);
              LatticeReal a = 1.0;
              LatticeComplex rv = cmplx(a * hh, 0);
	      for(int color_source(0);color_source<Nc;color_source++){
                  LatticeColorVector vec_srce = zero ;
                  pokeColor(vec_srce,rv,color_source) ;
		  for(int spin_source=0; spin_source < Ns; ++spin_source){
		      *source[i_v * Ns * Nc + color_source * Ns + spin_source] = zero; 
		      *psi[i_v * Ns * Nc + color_source * Ns + spin_source] = zero; 
		      *chi = zero; 
		      *eta  = zero; 
		      *tmp  = zero;
		      CvToFerm(vec_srce, *source[i_v * Ns * Nc + color_source * Ns + spin_source], spin_source);

		      //do D_ee^{-1} * x_e and D_oo^{-1} * x_o
		      D.evenEvenInvLinOp(*eta, *source[i_v * Ns * Nc + color_source * Ns + spin_source], PLUS);
		      D.oddOddInvLinOp(*eta, *source[i_v * Ns * Nc + color_source * Ns + spin_source], PLUS);
		      *psi[i_v * Ns * Nc + color_source * Ns + spin_source] = *eta;
		      *chi = *eta;
		      for (int p = 1; p < param.hpe_power; ++p){
		      
			  // now do D_ee^{-1} * D_eo * above and
			  //  D_oo^{-1} * D_oe * above
			  D.oddHoppingOp(*tmp, *chi, PLUS);
			  D.evenHoppingOp(*tmp, *chi, PLUS);
			  *chi = *tmp;
			  *psi[i_v * Ns * Nc + color_source * Ns + spin_source] += *tmp;
		      } //p

		      } //spin
		    } //color

		  } //iv

		do_disco_hop(db, source, psi, ft, param.use_ferm_state_links ? state->getLinks() : u, disp);
	    } //k1
		 

         //if (new_path_d.size() < param.max_path_length){
         // new_path_d.resize(new_path_d.size()+1);
	  //}

      } //disp
			  

    } //end of exact hpe

    void exactHPE2(std::map< KeyOperator_t, ValOperator_t >& db,
                  const InlineDiscoFreqSplitProbMMSuperb::Params::Param_t& param,
                  const multi1d<LatticeColorMatrix>& u)
    {

      typedef LatticeFermion               T;
      typedef multi1d<LatticeColorMatrix>  P;
      typedef multi1d<LatticeColorMatrix>  Q;

      std::vector<std::shared_ptr<Coloring>> hop_coloring(param.max_path_length+1);
      for (int i = 0; i < param.max_path_length+1; i++){
      if (param.hpe_probing_files[i] != "nofile"){
        QDPIO::cout << "Reading colors from file " << param.hpe_probing_files[i] << std::endl;
        hop_coloring[i].reset(new Coloring(param.hpe_probing_files[i]));
      } else {
        QDPIO::cout << "Generating a " << param.hpe_power - 1 << "-distance coloring for the Hopping Term" << std::endl;
        hop_coloring[i].reset(new Coloring(i, param.hpe_power - 1));
      }
      }

      Real m_q0 = getMass(param.prop.fermact);
      QDPIO::cout << "Base mass m_q0 = " << toDouble(m_q0) << std::endl;
      std::string ferm_xml = param.prop.fermact.xml;
      modifyMass(ferm_xml, m_q0, param.shifts[param.shifts.size()-1]);
      QDPIO::cout << "Performing exact HPE contribution with m_q = " << toDouble(m_q0 + param.shifts[param.shifts.size()-1]) << std::endl;
      std::istringstream xml_t(ferm_xml);
      XMLReader fermacttop(xml_t);

      Handle< FermionAction<T,P,Q> >
              S_f(TheFermionActionFactory::Instance().createObject(param.prop.fermact.id,
                                                               fermacttop,
                                                               param.prop.fermact.path));
      Handle< FermState<T,P,Q> > state(S_f->createState(u));

      EvenOddPrecCloverLinOp D;
      CloverFermActParams clov_params(fermacttop, param.prop.fermact.path);
      D.create(state, clov_params);

      int decay_dir           = Nd-1 ;
      SftMom ft = param.use_p_list ? SftMom(param.p_list, decay_dir) : SftMom(param.p2_max, false, decay_dir);

      multi1d<int> new_path_d;
      multi1d<int> path_d;

      for (int disp = 0; disp < param.max_path_length+1; disp++){
      int Nvecs = hop_coloring[disp]->numColors();
      const int N_v = (param.max_rhs + Ns * Nc - 1) / Ns / Nc;

      QDPIO::cout << Nvecs << " colors for the hopping term trace calculation for Displacement " << disp << std::endl;

      
      QDPIO::cout<<"Now computing the hopping contribution"<<std::endl;
     for (int k1 = 0, dk = std::min(Nvecs, N_v); k1 < Nvecs ; k1 += dk, dk = std::min(Nvecs - k1, N_v)) {

    LatticeFermion source, psi, chi, eta, tmp;

        for (int i_v = 0 ; i_v < dk ; i_v++) {
      LatticeInteger hh ;
      hop_coloring[disp]->getVec(hh, k1 + i_v);
      LatticeReal a = 1.0;
      LatticeComplex rv = cmplx(a * hh, 0);
      for(int color_source(0);color_source<Nc;color_source++){
          LatticeColorVector vec_srce = zero ;
          pokeColor(vec_srce,rv,color_source) ;
          for(int spin_source=0; spin_source < Ns; ++spin_source){

	source = zero; psi = zero; chi = zero; eta = zero; tmp = zero;

	CvToFerm(vec_srce, source, spin_source);


	D.evenEvenInvLinOp(eta, source, PLUS);
	D.oddOddInvLinOp(eta, source, PLUS);
	psi = eta;
	chi = eta;
	for (int p = 1; p < param.hpe_power; ++p){


	    D.oddHoppingOp(tmp, chi, PLUS);
	    D.evenHoppingOp(tmp, chi, PLUS);
	    chi = tmp;
	    psi += tmp;
	}
            LatticeFermion q; 
            if (disp != 0){
	for(int sign(-1);sign<2;sign+=2){
	    for(int mu(2);mu<3;mu++){
	  for( int i = 0; i < new_path_d.size() ; i++){
	      new_path_d[i] = sign*(mu+1);
	  }

	  q = psi;
	    LatticeFermion q_mu;
	    for (int k = 1; k < disp+1; k++){
	        if (param.use_ferm_state_links){
	      do_shift(q_mu, q, state->getLinks(), mu, sign);
	        } else {
	      do_shift(q_mu, q, u, mu, sign);
	        }
	         q = q_mu;
	    }
	

	    if (param.use_ferm_state_links){
	        do_disco_hop(db, source, q, ft, state->getLinks(), new_path_d, param.max_path_length);
	    }else{
	        do_disco_hop(db, source, q, ft, u, new_path_d, param.max_path_length); 
	    }
  

	    } 


            } 
	 
	 } else { 
	    if (param.use_ferm_state_links){
	    do_disco_hop(db, source, psi, ft, state->getLinks(), new_path_d, param.max_path_length);
	     }else{
	    do_disco_hop(db, source, psi, ft, u, new_path_d, param.max_path_length);
	    }
	}

	} 
           } 
        } 
          } 
       if (new_path_d.size() < param.max_path_length){
       new_path_d.resize(new_path_d.size()+1);
      }
  } //disp

  } //end of func

    std::vector<Interp::MinCosts_t> getOptimalShifts(SB::ShiftedChimeraSolver& PP, InlineDiscoFreqSplitProbMMSuperb::Params::Param_t& param, const multi1d<LatticeColorMatrix>& u)
    {

      typedef LatticeFermion               T;
      typedef multi1d<LatticeColorMatrix>  P;
      typedef multi1d<LatticeColorMatrix>  Q;

      Real tmpmass = getMass(param.prop.fermact);

      std::vector<EvenOddPrecCloverLinOp> D;
      std::vector<bool> do_hpe;
      do_hpe.resize(param.num_shifts);
      int num_clovers = 0;
      std::vector<int> which_clov;
      if (param.use_hpe){
      for (int i = 0; i < param.shifts.size(); i++){
          if ( toDouble(tmpmass + param.shifts[i]) >= 0.0){
              do_hpe[i] = true;
	      which_clov.push_back(i);
	      num_clovers++;
          }else{
              do_hpe[i] = false;
          }
	}
      }else{
      for (int i = 0; i <param.shifts.size(); i++){
	  do_hpe[i] = false;
      }
      }

      if (param.use_hpe){
      D.resize(num_clovers);
      for (int i = 0; i < num_clovers; i++){
      modifyMass(param.prop.fermact.xml, tmpmass, param.shifts[which_clov[i]]);
      std::istringstream  xml_s(param.prop.fermact.xml);
      XMLReader  fermacttop(xml_s);
      Handle< FermionAction<T,P,Q> > S_r(TheFermionActionFactory::Instance().createObject(param.prop.fermact.id,
                                                               fermacttop,
                                                               param.prop.fermact.path));
      Handle< FermState<T,P,Q> > state_r(S_r->createState(u));
      CloverFermActParams clov_params(fermacttop, param.prop.fermact.path);
      D[i].create(state_r, clov_params);
      }
      }

      
      std::shared_ptr<Coloring> coloring;
      if (!param.probing_file.empty()) {
        QDPIO::cout << "Reading colors from file " << param.probing_file << std::endl;
        coloring.reset(new Coloring(param.probing_file));
      } else {
        QDPIO::cout << "Generating a " << param.probing_distance << "-distance coloring" << std::endl;
        coloring.reset(new Coloring(param.probing_distance, param.probing_power));
      }

    //declared before the linear equations solve
    int decay_dir           = Nd-1 ;
    SftMom ft = param.use_p_list ? SftMom(param.p_list, decay_dir) : SftMom(param.p2_max, false, decay_dir);
    int Nsrc = coloring->numColors();
    QDPIO::cout << "num colors " << Nsrc << std::endl;
    QDPIO::cout << "Optimizing z Displacement  = " << param.disp_gamma[0] << " and Gamma = " << param.disp_gamma[1] << std::endl;
    //set to the original mass
    modifyMass(param.prop.fermact.xml, tmpmass, param.shifts[0]);
    modifyMass(param.prop.invParam.xml, tmpmass, param.shifts[0]);
    std::istringstream xml_s(param.prop.fermact.xml);
    XMLReader fermacttop(xml_s);
    
      Handle< FermionAction<T,P,Q> >
              S_f(TheFermionActionFactory::Instance().createObject(param.prop.fermact.id,
                                                               fermacttop,
                                                               param.prop.fermact.path));

      Handle< FermState<T,P,Q> > state(S_f->createState(u));

     //Handle< SystemSolver<LatticeFermion> > PP = S_f->qprop(state,
//							      param.prop.invParam);	
      //SB::ShiftedChimeraSolver PP{param.prop.fermact, param.prop.invParam, u};
      //have to calculate the test shifts or have them provided...
    

      Complex NOne;
      Real a;
      a = -1.0;
      NOne = cmplx(a,0.0);

      Interp::CostHolder_t costs;
      costs.r_variances.resize(param.num_shifts);
      costs.rs_variances.resize(param.num_shifts, std::vector<double>(param.num_shifts));
      costs.level_costs.resize(param.num_shifts); 
      costs.shifts.resize(param.num_shifts);
      QDPIO::cout << "Shifts used to estimate the variance are : " << std::endl;
      for (int i = 0; i < param.num_shifts; i++){
      costs.shifts[i] = toDouble(param.shifts[i]);
      QDPIO::cout << costs.shifts[i] << std::endl;
      }   

      MatMap dbrs_mean(param.num_shifts, VecMap(param.num_shifts));;
      MatMap dbrs_var(param.num_shifts, VecMap(param.num_shifts));;
      VecMap dbr_mean(param.num_shifts);
      VecMap dbr_var(param.num_shifts);

      std::vector<Interp::MinCosts_t> mincosts;
      int pos;
      std::vector<double> level_sum;
      level_sum.resize(param.num_shifts);
      std::vector<double> count(2); //for checkpoint
      //int level = 0; //for checkpoint

      //Restart_t restart;
      //QDP::Seed rng_seed;
      //QDP::RNG::savern(rng_seed);
      //std::string seed_file = "current_seed_" + param.mass_label + "_" + param.probing_file;
      /*if (param.interp_restart && param.use_interpolation){
        for (int i = 0; i < param.num_shifts; i++){
	for (int j = i; j < param.num_shifts; j++){

	//this does all the product terms
        std::string filename_traces = "interp_traces_shifts_" + std::to_string(i) + std::to_string(j) + "_" + param.mass_label + "_" + param.probing_file;
        std::string filename_vars = "interp_variances_shifts_" + std::to_string(i) + std::to_string(j) + "_" + param.mass_label + "_" + param.probing_file;
        std::ifstream file_t;
        std::ifstream file_v;
        file_t.open(filename_traces, std::ios::in);
        file_v.open(filename_vars, std::ios::in);
        if (file_t && file_v){
        QDPIO::cout << "Trace and Variance checkpoint files exist. Reading in the values." << std::endl;
        file_t.close();
        file_v.close();
        restart = checkpoint_in(filename_traces, filename_vars, param.max_path_length);
        dbrs_mean[i][j] = restart.db_trace;
        dbrs_var[i][j] = restart.db_var;
        level_sum = restart.level_cost;
        } //file
	} //j
	
	//this does the single terms
	std::string filename_traces = "interp_traces_shifts_" + std::to_string(i) + "_" + param.mass_label + "_" + param.probing_file;
        std::string filename_vars = "interp_variances_shifts_" + std::to_string(i) + "_" + param.mass_label + "_" + param.probing_file;
        std::ifstream file_t;
        std::ifstream file_v;
        file_t.open(filename_traces, std::ios::in);
        file_v.open(filename_vars, std::ios::in);
        if (file_t && file_v){
        QDPIO::cout << "Trace and Variance checkpoint files exist. Reading in the values." << std::endl;
        file_t.close();
        file_v.close();
        restart = checkpoint_in(filename_traces, filename_vars, param.max_path_length);
        dbr_mean[i] = restart.db_trace;
        dbr_var[i] = restart.db_var;
	level_sum = restart.level_cost;
	}//file
        } //i
        rng_seed = restart_seed_in(seed_file);
        QDP::RNG::setrn(rng_seed);
        QDPIO::cout << "Seed " << rng_seed << " has been set" << std::endl;
      }else{
        restart.level = 0;
        restart.noise = 0;
      } */



      for (int noise = 0; noise < param.num_samples; noise++){
      //int noise;
      //for ( param.interp_restart ? noise = restart.noise+1 : noise = 0; noise < param.num_samples; noise++){ 
	  MatMap dbrs(param.num_shifts, VecMap(param.num_shifts));;
	  VecMap dbr(param.num_shifts);

	  QDPIO::cout << " Doing noise vector " << noise  << std::endl;
	
  
	//generate a random std::vector
        LatticeComplex vec ;
        LatticeReal rnd1, theta;
        random(rnd1);
        Real twopiN = Chroma::twopi / 4;
        theta = twopiN * floor(4*rnd1);
        vec = cmplx(cos(theta),sin(theta));


	const int N_rhs = (param.max_rhs_interp + Ns * Nc - 1) / Ns / Nc;
        for (int k1 = 0, dk = std::min(Nsrc, N_rhs); k1 < Nsrc ; k1 += dk, dk = std::min(Nsrc - k1, N_rhs)) {
	     X1dvector v_chi(Ns * Nc * dk);
	     for (int col=0; col<v_chi.size(); col++) v_chi[col].reset(new LatticeFermion);
	     X2dmatrix v_psi(param.num_shifts, X1dvector(Ns * Nc * dk));

	     for (int row = 0; row < param.num_shifts; row++){
		 for (int col = 0; col < v_psi[row].size(); col++){
		     v_psi[row][col].reset(new LatticeFermion);
		  }
	
	     }


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
		for (int row=0; row < param.num_shifts; row++) *v_psi[row][i_v * Ns * Nc + color_source * Ns + spin_source]  = zero;
	      }
	    }
	  }

	  for (int m1 = 0; m1 < param.num_shifts; m1++){
	      QDPIO::cout << "Solving for shift number " << m1 << " while generating the variances " <<  std::endl;
	      //modifyMass(param.prop.fermact.xml, tmpmass, param.shifts[m1]);
	      //modifyMass(param.prop.invParam.xml, tmpmass, param.shifts[m1]);
	      //std::istringstream xml_r(param.prop.fermact.xml);
	      //XMLReader fermacttop(xml_r);
              //Handle< FermionAction<T,P,Q> > S_r(TheFermionActionFactory::Instance().createObject(param.prop.fermact.id,
              //                                                 fermacttop,
              //                                                 param.prop.fermact.path));

	      //Handle< FermState<T,P,Q> > state_r(S_r->createState(u));

	      //Handle< SystemSolver<LatticeFermion> > PP = S_r->qprop(state_r, param.prop.invParam);
	      //std::vector<SystemSolverResults_t> res = (*PP)(v_psi[m1], std::vector<std::shared_ptr<const LatticeFermion>>(v_chi.begin(), v_chi.end()));
	      PP.shiftSolver(tmpmass+param.shifts[m1]);
	      StopWatch solve_time;
	      solve_time.start();
	      SB::doInversion(PP, v_psi[m1], std::vector<std::shared_ptr<const LatticeFermion>>(v_chi.begin(), v_chi.end()));
	      solve_time.stop();
	      level_sum[m1] += solve_time.getTimeInSeconds();
		//  for (int t = 0; t < res.size(); t++){
		 //     level_sum[m1] += 1.0 * res[t].n_count;
		 // }



	   }


	   //now do all the dot products, first does just the (D+sI)^{-1} terms
	   QDPIO::cout << "Entering dot products for product and single terms" << std::endl;
	   StopWatch swatch_dots;
       swatch_dots.start();
	   int clov_tracker = 0;
	   for (int m1 = 0; m1 < param.num_shifts; m1++){
	       assert(v_chi.size() == v_psi[m1].size());
	    if (param.use_hpe && do_hpe[m1] && m1 == which_clov[clov_tracker]){
		  std::vector<std::shared_ptr<LatticeFermion>> v_eta(Nc * Ns * dk);
		  for (int col = 0; col < v_eta.size(); col++){v_eta[col].reset(new LatticeFermion);}
		  LatticeFermion v_tmp;
		  //*v_eta[idk * Nc * Ns + col] = zero;
		  //for (int i = 0; i < v_psi[m1].size(); ++i){
		  for (int idk = 0; idk < dk; idk++){
		    for (int col = 0; col < Nc * Ns; col++){
		      *v_eta[idk * Nc * Ns + col] = *v_psi[m1][idk * Nc * Ns + col];
		      for (int p = 1; p < param.hpe_power+1; ++p){
			    D[clov_tracker].evenHoppingOp(v_tmp, *v_eta[idk * Nc * Ns + col], PLUS);
			    D[clov_tracker].oddHoppingOp(v_tmp, *v_eta[idk * Nc * Ns + col], PLUS);
			    *v_eta[idk * Nc * Ns + col] = v_tmp;
		      } //p
		  } //col
		} //idk
		do_disco(dbr[m1], v_chi, v_eta, ft, param.use_ferm_state_links ? state->getLinks() : u, param.max_path_length);
		clov_tracker++;
	    }else{ //do_hpe	 
		   do_disco(dbr[m1], v_chi, v_psi[m1], ft, param.use_ferm_state_links ? state->getLinks() : u, param.max_path_length);
	    } //do_hpe
		    
	    //now do all the dot products for the product combinations
	    //changing things up in the interpolation
	    //now need the terms with out the shift difference in front

	      for (int m2 = m1; m2 < param.num_shifts; m2++){
		  assert(v_psi[m1].size() == v_psi[m2].size());
		  std::vector<std::shared_ptr<LatticeFermion>> v_eta(Nc * Ns * dk);
		  for (int col = 0; col < v_eta.size(); col++){v_eta[col].reset(new LatticeFermion);}
		  
		  for (int idk = 0; idk < dk; idk++){
		    for (int col = 0; col < Nc * Ns; col++){
		      if (col == 2 || col == 3 || col == 6 || col == 7 || col == 10 || col == 11){
			 *v_eta[idk * Nc * Ns + col] = NOne * (Gamma(15) * *v_psi[m1][idk * Nc * Ns + col]);
		      }else{
			 *v_eta[idk * Nc * Ns + col] =  (Gamma(15) * *v_psi[m1][idk * Nc * Ns + col]); 
		     }
		} //col
	      } //idk
		  do_disco(dbrs[m1][m2], v_eta, v_psi[m2], ft, param.use_ferm_state_links ? state->getLinks() : u, param.max_path_length);
	    } //m2   

	    } // m1
	  swatch_dots.stop();
	  QDPIO::cout << "Computing inner products " << swatch_dots.getTimeInSeconds() << " secs" << std::endl;


	} //k1

	QDPIO::cout << "Updating the dbs " << std::endl;
        for (int m1 = 0; m1 < param.num_shifts; m1++){
	    for (int m2 = m1; m2 < param.num_shifts; m2++){
	    do_update(dbrs_mean[m1][m2], dbrs_var[m1][m2], dbrs[m1][m2], noise == 0);
	    //std::string filename_traces = "interp_traces_shifts_" + std::to_string(m1) + std::to_string(m2) + "_" + param.mass_label + "_" + param.probing_file;
	    //std::string filename_vars = "interp_variances_shifts_" + std::to_string(m1) + std::to_string(m2) + "_" + param.mass_label + "_" + param.probing_file;
	    //checkpoint_out(level, noise, rng_seed, param.num_shifts, level_sum, count, param.shifts, dbrs_mean[m1][m2], dbrs_var[m1][m2], filename_traces, filename_vars);
	    }
	    do_update(dbr_mean[m1], dbr_var[m1], dbr[m1], noise == 0);
	    //std::string filename_traces = "interp_traces_shifts_" + std::to_string(m1) + "_" + param.mass_label + "_" + param.probing_file;
	    //std::string filename_vars = "interp_variances_shifts_" + std::to_string(m1) + "_" + param.mass_label + "_" + param.probing_file;
	    //checkpoint_out(level, noise, rng_seed, param.num_shifts, level_sum, count, param.shifts, dbr_mean[m1], dbr_var[m1], filename_traces, filename_vars);
	}
        //QDP::RNG::savern(rng_seed);
        //restart_seed_out(seed_file, rng_seed);


	//need to divide the costs by the current number of noise vectors to get one sample variance
        for (int t = 0; t < costs.level_costs.size(); t++){
            costs.level_costs[t] = level_sum[t] / (noise + 1);
        }

	for (int k = 0; k < param.max_path_length+1; ++k){
        QDPIO::cout << "Retrieving product and single variances for displacement " << k << " and Gamma " << param.disp_gamma[1] <<  std::endl;
        //costs.rs_variances = retrieve_variances(dbrs_mean, dbrs_var, 1, noise+1, param.num_shifts, param.gamma_disp[0], param.gamma_disp[1]);
        costs.rs_variances = retrieve_variances(dbrs_mean, dbrs_var, 1, noise+1, param.num_shifts, k, param.disp_gamma[1]);
        //costs.r_variances = retrieve_variances(dbr_mean, dbr_var, 1, noise+1, param.num_shifts, param.gamma_disp[0], param.gamma_disp[1]);
        costs.r_variances = retrieve_variances(dbr_mean, dbr_var, 1, noise+1, param.num_shifts, k, param.disp_gamma[1]);


	//make sure all nodes get the same data
	   for (int i = 0; i < costs.r_variances.size(); i++){
	       QDPIO::cin >> costs.r_variances[i];
	       for (int j = 0; j < costs.rs_variances[i].size(); j++){
	       QDPIO::cin >> costs.rs_variances[i][j];
	       }
	    }

	//here is the call for interpolation
	if (noise > 0){
	mincosts = Interp::getMinShifts(costs, param.del_s, param.num_bcshifts, param.use_mg, k, param.disp_gamma[1]);
	//QDPIO::cout << "Cost for regular calculation of Disp = " << param.gamma_disp[0] << ", Gamma = " << param.gamma_disp[1] << "for noise vector " << noise << " is : " << mincosts[0].reg_costs << std::endl;
	QDPIO::cout << "Cost for regular calculation of Disp = " << k  << ", Gamma = " << param.disp_gamma[1] << "for noise vector " << noise << " is : " << mincosts[0].reg_costs << std::endl;
	}
      } //k
    
    //Need to report the cost of the original calculation (D^{-1}) to report speedups
    for (int k = 0; k < param.max_path_length+1; ++k){
	for (int g = 0; g < Ns * Ns; ++g){
	    std::vector<double> vs = retrieve_variances(dbr_mean, dbr_var, 1, noise+1, param.num_shifts, k, g);
	    QDPIO::cout << "Cost of D^{-1} with Disp = " << k << " and Gamma = " << g << "is : " << costs.level_costs[0]*vs[0] << std::endl;
	}
    }

    } //noise



    return mincosts;
    } // end getOptimalShifts
 
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
      
    const std::string name = "DISCO_FREQUENCY_SPLITTING_PROB_MM_SUPERB";

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


        //QDP::Seed rng_seed;
        //QDP::RNG::savern(rng_seed);
        //QDPIO::cout << "The seed before setup is  " << rng_seed << std::endl;

      // Calculate some gauge invariant observables just for info.
      MesPlq(xml_out, "Observables", u);

      

      //get the initial mass before anything!!
      Real m_qi = getMass(params.param.prop.fermact);
      //int nodelist = Layout::nodeNumber();


          modifyMass(params.param.prop.fermact.xml, m_qi, params.param.shifts[0]);
          modifyMass(params.param.prop.invParam.xml, m_qi, params.param.shifts[0]);
          std::istringstream  xml_s(params.param.prop.fermact.xml);
          XMLReader  fermactr(xml_s);
          Handle< FermionAction<T,P,Q> > S_s(TheFermionActionFactory::Instance().createObject(params.param.prop.fermact.id,
                                                               fermactr,
                                                               params.param.prop.fermact.path));
          Handle< FermState<T,P,Q> > state(S_s->createState(u));
	  //Handle< SystemSolver<LatticeFermion> > PP = S_s->qprop(state, params.param.prop.invParam);
	  SB::ShiftedChimeraSolver PP{params.param.prop.fermact, params.param.prop.fermact, u};


      std::vector<Interp::MinCosts_t> mincosts;
      if(params.param.use_interpolation){
      mincosts = getOptimalShifts(PP, params.param, u);
      
      std::vector<double> tmpcosts;
      tmpcosts.resize(mincosts.size());
      for (int i = 0; i < tmpcosts.size(); i++){
	  tmpcosts[i] = mincosts[i].min_costs;
      }
      std::vector<double>::iterator it = std::min_element(tmpcosts.begin(), tmpcosts.end());
      
      int pos = std::distance(tmpcosts.begin(), it);
      QDPIO::cout << "The minimum cost found with interpolation is " << *it << std::endl;

      int nshifts = mincosts[pos].shifts.size()+1;
      params.param.shifts.resize(nshifts);
      params.param.shifts[0] = 0.0;
      for (int i = 0; i < mincosts[pos].shifts.size(); i++){
	  params.param.shifts[i+1].elem() = mincosts[pos].shifts[i];
      }

      QDPIO::cout << "The minimum shifts are : " <<std::endl;
      for (int i = 1; i < params.param.shifts.size(); i++){
	  QDPIO::cout << params.param.shifts[i] << std::endl;
      } 
      
      }//use_interpolation

       

      if (params.param.do_trace_est){

      std::shared_ptr<Coloring> coloring;
      if (!params.param.probing_file.empty()) {
        QDPIO::cout << "Reading colors from file " << params.param.probing_file << std::endl;
        coloring.reset(new Coloring(params.param.probing_file));
      } else {
        QDPIO::cout << "Generating a " << params.param.probing_distance << "-distance coloring" << std::endl;
        coloring.reset(new Coloring(params.param.probing_distance, params.param.probing_power));
      }

      //Reset the mass
      modifyMass(params.param.prop.fermact.xml, m_qi, params.param.shifts[0]); 
	
	  
      // Initialize the slow Fourier transform phases
      int decay_dir           = Nd-1 ; // hadamard needs this for now
      //Initialize ft differently based on momentum list or max value.
      SftMom ft = params.param.use_p_list ? SftMom(params.param.p_list, decay_dir) : SftMom(params.param.p2_max, false, decay_dir);

      // number of colors
      int Nsrc = coloring->numColors();
      QDPIO::cout << "num colors " << Nsrc << std::endl;


      StopWatch swatch;
      swatch.start();

	  
      //reset the mass
      modifyMass(params.param.prop.fermact.xml, m_qi, params.param.shifts[0]);
      Complex NOne;
      Real a;
      a = -1.0;
      NOne = cmplx(a,0.0);
      int num_levels = params.param.shifts.size();
      QDPIO::cout << "Number of levels is " << num_levels << std::endl;
      int level_switch;

      StatHolder_t dbs;
      dbs.levels_var.resize(num_levels);
      dbs.levels_avg.resize(num_levels);
      std::vector<double> cl(num_levels);
      //std::vector<double> cl;

          //modifyMass(params.param.prop.fermact.xml, m_qi, params.param.shifts[0]);
          //modifyMass(params.param.prop.invParam.xml, m_qi, params.param.shifts[0]);
          //std::istringstream  xml_s(params.param.prop.fermact.xml);
          //XMLReader  fermactr(xml_s);
          //Handle< FermionAction<T,P,Q> > S_s(TheFermionActionFactory::Instance().createObject(params.param.prop.fermact.id,
          //                                                     fermactr,
          //                                                     params.param.prop.fermact.path));
          //Handle< FermState<T,P,Q> > state(S_s->createState(u));
	  //Handle< SystemSolver<LatticeFermion> > PP = S_s->qprop(state, params.param.prop.invParam);
	  //SB::ShiftedChimeraSolver PP{params.param.prop.fermact, params.param.prop.fermact, u};


      //will never be above 7 (for now) so I am ok with this for the time being
      //int max_num_levels = 7;

      //if using a checkpoint, this block will be used
      //this whole thing needs to be redone, but it works....
      /*Restart_t restart;
      QDP::Seed rng_seed;
      QDP::RNG::savern(rng_seed);
      std::string seed_file = "current_seed_" + params.param.mass_label + "_" + params.param.probing_file;
      if (params.param.trace_restart){
	for (int i = 0; i < max_num_levels; i++){
	std::string filename_traces = "traces_level_" + std::to_string(i) + "_" + params.param.mass_label + "_" + params.param.probing_file;
	std::string filename_vars = "variances_level_" + std::to_string(i) + "_" + params.param.mass_label + "_" + params.param.probing_file;
	std::ifstream file_t;
	std::ifstream file_v;
	file_t.open(filename_traces, std::ios::in);
	file_v.open(filename_vars, std::ios::in);
	if (file_t && file_v){
	QDPIO::cout << "Trace and Variance checkpoint files exist. Reading in the values." << std::endl;
	file_t.close();
	file_v.close();
	restart = checkpoint_in(filename_traces, filename_vars, params.param.max_path_length);
  
	//dbs.levels_avg[i] = restart.db_trace;
	//dbs.levels_var[i] = restart.db_var;
	dbs.levels_avg.push_back(restart.db_trace);
	dbs.levels_var.push_back(restart.db_var);

	//show_stats(dbs.levels_avg[i], dbs.levels_var[i], dbdet, 1, restart.noise+1, restart.level);
	cl = restart.level_cost;
	} //file
	} //i
	rng_seed = restart_seed_in(seed_file);
	num_levels = restart.num_levels;
	params.param.shifts.resize(num_levels);
	std::map<KeyOperator_t, ValOperator_t> db;
	for (int i = restart.level+1; i < restart.num_levels; i++){
	    dbs.levels_var.push_back(db);
	    dbs.levels_avg.push_back(db);
	}
	for (int i = 0; i < num_levels; i++){ params.param.shifts[i].elem() = restart.shifts[i]; }
	QDP::RNG::setrn(rng_seed);
	QDPIO::cout << "Seed " << rng_seed << " has been set" << std::endl;
      }else{
	restart.level = 0;
	restart.noise = 0;
	num_levels = params.param.shifts.size();
	cl.resize(num_levels);
	dbs.levels_var.resize(num_levels);
	dbs.levels_avg.resize(num_levels);
      } */

      //do the exact part of the HPE
      std::map< KeyOperator_t, ValOperator_t > dbdet;
      if(params.param.use_hpe){
      QDPIO::cout<<"Now computing the exact HPE contribution"<<std::endl;
      StopWatch swatch_det;
      swatch_det.start();
      exactHPE(dbdet, params.param, u);
      swatch_det.stop();
      QDPIO::cout << "Exact HPE contribution computed in time= " << swatch_det.getTimeInSeconds() << " secs" << std::endl;
      }


      //reset the mass
      modifyMass(params.param.prop.fermact.xml, m_qi, params.param.shifts[0]);
      //set the mass to the last shift
      modifyMass(params.param.prop.fermact.xml, m_qi, params.param.shifts[params.param.shifts.size()-1]);
      EvenOddPrecCloverLinOp D;
      if (params.param.use_hpe){
      std::istringstream  xml_s(params.param.prop.fermact.xml);
      XMLReader  fermacttop(xml_s);
      Handle< FermionAction<T,P,Q> > S_r(TheFermionActionFactory::Instance().createObject(params.param.prop.fermact.id,
                                                               fermacttop,
                                                               params.param.prop.fermact.path));
      Handle< FermState<T,P,Q> > state_r(S_r->createState(u));
      CloverFermActParams clov_params(fermacttop, params.param.prop.fermact.path);
      D.create(state_r, clov_params);
      }



      //for (int level = restart.level; level < num_levels; level++){
      for (int level = 0; level < num_levels; level++){
	
	if (level < num_levels-1){
	    level_switch = level;
	}else{
	    level_switch = level - 1;
	}
	QDPIO::cout << "On level : " << level << std::endl;
	//std::string filename_traces = "traces_level_" + std::to_string(level) + "_" + params.param.mass_label + "_" + params.param.probing_file;
	//std::string filename_vars = "variances_level_" + std::to_string(level) + "_" + params.param.mass_label + "_" + params.param.probing_file;
	QDPIO::cout << "Computing with shifts " << params.param.shifts[level_switch] << " and " << params.param.shifts[level_switch+1] << std::endl;

	double count_r = 0.0;
	double count_s = 0.0;
	std::vector<double> count(2);
	//if (params.param.trace_restart && level == restart.level){
	 //  count = restart.count;
	//}
	//count.resize(2);
	
	//int noise;
	//for (int noise = restart.noise ; noise < params.param.noise_vectors[level]; noise++) {
	//for ( (params.param.trace_restart && (level == restart.level)) ?  noise = restart.noise+1 :  noise = 0 ; noise < params.param.noise_vectors[level]; noise++) {
	for (int noise = 0; noise < params.param.noise_vectors[level]; noise++){
	std::map< KeyOperator_t, ValOperator_t > db;
        // doing a new noise vector
        QDPIO::cout << " Doing noise vector " << noise  << std::endl; 
        //QDPIO::cout << "The seed is " << rng_seed << " on level " << level << " and noise vector " << noise << std::endl;


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
	  X1dvector v_chi(Ns * Nc * dk);
          for (int col=0; col<v_chi.size(); col++) {v_chi[col].reset(new LatticeFermion);}

	  X2dmatrix v_psi(2, X1dvector(Ns * Nc * dk));
	  for (int row = 0; row < 2; row++){
	      for (int col = 0; col < v_psi[row].size(); col++){
		  v_psi[row][col].reset(new LatticeFermion);
	      }
	  }

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
                //*v_psi[i_v * Ns * Nc + color_source * Ns + spin_source]  = zero;
		for (int row=0; row < 2; row++) {*v_psi[row][i_v * Ns * Nc + color_source * Ns + spin_source]  = zero;}
              } 
            }
          }

	  if (level < num_levels-1){ //level control

	int sol = 0;
	for (int m = level_switch; m < level_switch+2; m++){
	  /*modifyMass(params.param.prop.fermact.xml, m_qi, params.param.shifts[m]);
	  modifyMass(params.param.prop.invParam.xml, m_qi, params.param.shifts[m]);
	  std::istringstream  xml_r(params.param.prop.fermact.xml);
	  XMLReader  fermactr(xml_r);
          Handle< FermionAction<T,P,Q> > S_r(TheFermionActionFactory::Instance().createObject(params.param.prop.fermact.id,
                                                               fermactr,
                                                               params.param.prop.fermact.path));
          Handle< FermState<T,P,Q> > state_r(S_r->createState(u));
	  Handle< SystemSolver<LatticeFermion> > PP = S_r->qprop(state_r, params.param.prop.invParam);
	  std::vector<SystemSolverResults_t> res_r = (*PP)(v_psi[sol], std::vector<std::shared_ptr<const LatticeFermion>>(v_chi.begin(), v_chi.end()));
          for (int t = 0; t < res_r.size(); t++){
          count[sol] += 1.0 * res_r[t].n_count;
          }*/

	  PP.shiftSolver(m_qi+params.param.shifts[m]);
	  StopWatch solve_time;
	  solve_time.start();
	  SB::doInversion(PP, v_psi[sol], std::vector<std::shared_ptr<const LatticeFermion>>(v_chi.begin(), v_chi.end()));
	  solve_time.stop();
	  count[sol] += solve_time.getTimeInSeconds();
          QDPIO::cout << "On level = " << level << " and noise vector " << noise << " and count = " << count[sol] << std::endl;
	  sol++;
	}


          // here the recursive call goes to compute 
          // the loops
          // result is ADDED to db
          StopWatch swatch_dots;
          swatch_dots.start();
	  assert(v_chi.size() == v_psi[0].size());
	  DComplex cmplxshifts;
	  cmplxshifts = cmplx(params.param.shifts[level_switch+1]-params.param.shifts[level_switch],0.0);
	  do_disco(db, v_psi[0], v_psi[1], ft, params.param.use_ferm_state_links ? state->getLinks() : u, params.param.max_path_length, cmplxshifts, dk);
          swatch_dots.stop();
          QDPIO::cout << "Computing inner products " << swatch_dots.getTimeInSeconds() << " secs"
                      << std::endl;

	  }else{ // level control
	  

          /*modifyMass(params.param.prop.fermact.xml, m_qi, params.param.shifts[level_switch+1]);
          modifyMass(params.param.prop.invParam.xml, m_qi, params.param.shifts[level_switch+1]);
          std::istringstream  xml_r(params.param.prop.fermact.xml);
          XMLReader  fermactr(xml_r);
          Handle< FermionAction<T,P,Q> > S_r(TheFermionActionFactory::Instance().createObject(params.param.prop.fermact.id,
                                                               fermactr,
                                                               params.param.prop.fermact.path));
          Handle< FermState<T,P,Q> > state_r(S_r->createState(u));
          Handle< SystemSolver<LatticeFermion> > PP = S_r->qprop(state_r, params.param.prop.invParam);
          std::vector<SystemSolverResults_t> res_r = (*PP)(v_psi[0], std::vector<std::shared_ptr<const LatticeFermion>>(v_chi.begin(), v_chi.end()));
          for (int t = 0; t < res_r.size(); t++){
          count[0] += 1.0 * res_r[t].n_count;
          } */

	  PP.shiftSolver(m_qi+params.param.shifts[level_switch+1]);
	  StopWatch solve_time;
	  solve_time.start();
	  SB::doInversion(PP, v_psi[0], std::vector<std::shared_ptr<const LatticeFermion>>(v_chi.begin(), v_chi.end()));
	  solve_time.stop();
	  count[0] += solve_time.getTimeInSeconds(); 
          QDPIO::cout << "On level = " << level << " and noise vector " << noise << " and count = " << count[0] << std::endl;

          StopWatch swatch_dots2;
          swatch_dots2.start();
          assert(v_chi.size() == v_psi[0].size());
	  if (params.param.use_hpe){
		  LatticeFermion v_eta;
		  LatticeFermion v_tmp;
          for (int idk = 0; idk < dk; idk++){
            for (int col = 0; col < Nc * Ns; col++){
		   v_eta = *v_psi[0][idk * Nc * Ns + col];
		   for (int p = 1; p < params.param.hpe_power+1; ++p){
			D.evenHoppingOp(v_tmp, v_eta, PLUS);
			D.oddHoppingOp(v_tmp, v_eta, PLUS);
			v_eta = v_tmp;
		   } //p
		   *v_psi[0][idk * Nc * Ns + col] = v_eta;
		  } //col
	    } // idk
		do_disco(db, v_chi, v_psi[0], ft, params.param.use_ferm_state_links ? state->getLinks() : u, params.param.max_path_length);
	  } else { //use_hpe
		do_disco(db, v_chi, v_psi[0], ft, params.param.use_ferm_state_links ? state->getLinks() : u, params.param.max_path_length);  
	  }
	  swatch_dots2.stop();
      QDPIO::cout << "Computing inner products " << swatch_dots2.getTimeInSeconds() << " secs" << std::endl;

	} //level control
		   
    } // for k1


// fix below!!! //
// comment out below for debugging...

       // Update the mean and the average
        do_update(dbs.levels_avg[level], dbs.levels_var[level], db, noise == 0);

	//checkpoint_out(level, noise, rng_seed, num_levels, cl, count, params.param.shifts, dbs.levels_avg[level], dbs.levels_var[level], filename_traces, filename_vars);
	//QDP::RNG::savern(rng_seed);
	//restart_seed_out(seed_file, rng_seed);

        // Show stats
        if(params.param.display_level_stats){
        show_stats(dbs.levels_avg[level], dbs.levels_var[level], dbdet, 1, noise+1, level);
	}
	
    } // noise

    if (level < num_levels - 1){
      cl[level] = (count[0] + count[1])/ (1.0 * params.param.noise_vectors[level]);
      QDPIO::cout << "Cost on level " << level << " is " << cl[level] << std::endl;
    }else{
      cl[level] = (count[0])/(1.0 * params.param.noise_vectors[level]);
      QDPIO::cout << "Cost on level " << level << " is " << cl[level] << std::endl;
    }



  } //level

 for (int k = 0; k < params.param.max_path_length+1; k++){
  for (int g = 0; g < Ns * Ns; g++){
  std::vector<double> target_vars(num_levels);
for (int level = 0; level < num_levels; level++){
    target_vars[level] = retrieve_variances(dbs.levels_avg[level], dbs.levels_var[level], 1, params.param.noise_vectors[level], k , g);
}
    QDPIO::cout << "Cost on each level is : " << std::endl;
    std::vector<double> total_cost(num_levels);
    for (int level = 0; level < num_levels; level++){
	total_cost[level] = std::sqrt(cl[level] * target_vars[level]);
	QDPIO::cout << "Level " << level << ": " << total_cost[level] << std::endl;
    }
    //total_cost = std::pow(total_cost, 2.0);
    QDPIO::cout << "Total cost of the calculation for Displacement : " << k << " and Gamma : " << g << " is " << std::pow(std::accumulate(total_cost.begin(), total_cost.end(), 0.0), 2.0) << std::endl;
    QDPIO::cout << "Total variance of the calculation for Displacement : " << k << " and Gamma : " << g << " is " << std::accumulate(target_vars.begin(), target_vars.end(), 0.0) << std::endl;
//comment out below for debugging

  } //g
 } //k


      // Normalize the traces
     for (int level = 0; level < num_levels; level++){
      for(std::map< KeyOperator_t, ValOperator_t >::iterator it=dbs.levels_avg[level].begin();it != dbs.levels_avg[level].end(); it++){
        for(int k=0;k<it->second.op.size();k++){
          it->second.op[k] = it->second.op[k]/toDouble(params.param.noise_vectors[level]);
        }
      }
     }

    
      //add the traces of the levels
      for (int level = 0; level < num_levels; level++){
        for(std::map< KeyOperator_t, ValOperator_t >::iterator it=dbs.total_avg.begin();it != dbs.total_avg.end(); it++){
            it->second.op += dbs.levels_avg[level][it->first].op;
        }
      }


      // Add the deterministic part to the traces
      for(std::map< KeyOperator_t, ValOperator_t >::iterator it=dbs.total_avg.begin();it != dbs.total_avg.end(); it++)
        it->second.op += dbdet[it->first].op;

      swatch.stop();
      QDPIO::cout << "Traces were  computed: time= " 
		  << swatch.getTimeInSeconds() 
		  << " secs" << std::endl;

      // write out the results
      
      // DB storage         
      BinaryStoreDB<SerialDBKey<KeyOperator_t>,SerialDBData<ValOperator_t> > qdp_db;
      
      // Open the file, and write the meta-data and the binary for this operator
      {
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

	//qdp_db.open(params.named_obj.sdb_file, O_RDWR | O_CREAT, 0664);
	//Slightly modify code to account for changes from multifile write.
	//Be consistent with old mode of filename write.
	std::string file_name = params.named_obj.sdb_file;
	qdp_db.open(file_name, O_RDWR | O_CREAT, 0664);

	qdp_db.insertUserdata(file_str);
      }
   
      SerialDBKey <KeyOperator_t> key ;
      SerialDBData<ValOperator_t> val ;
      std::map< KeyOperator_t, ValOperator_t >::iterator it;
      // Store all the data
      for(it=dbs.total_avg.begin();it!=dbs.total_avg.end();it++){
	key.key()  = it->first  ;
        key.key().mass_label = params.param.mass_label;
	val.data().op.resize(it->second.op.size()) ;
	for(int i(0);i<it->second.op.size();i++)
          val.data().op[i] = it->second.op[i];
	qdp_db.insert(key,val);
      }


      pop(xml_out);  // close last tag


      snoop.stop();
      QDPIO::cout << name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << std::endl;
      
      QDPIO::cout << name << ": ran successfully" << std::endl;
      
      END_CODE(); 

      } //debug
  
  
    } //func


  }// namespace

} // namespace Chroma
// vim: sw=2 sts=2
