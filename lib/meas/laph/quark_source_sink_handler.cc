#include "quark_source_sink_handler.h"
#include "meas/sources/zN_src.h"
#include "util/ft/sftmom.h"
#include "fermact.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "util/ferm/diractodr.h"


namespace Chroma {
  namespace LaphEnv {
   

QuarkSourceSinkHandler::Key::~Key() {}

QuarkSourceSinkHandler::Key::Key(XMLReader& xml_in) : noise(xml_in)
{
 if (xml_tag_count(xml_in,"QuarkSourceSinkKey")!=1){
    QDPIO::cerr << "Bad XML input to QuarkSourceSinkKey"<<endl;
    QDPIO::cerr << "Expected one <QuarkSourceSinkKey> tag"<<endl;
    QDP_abort(1);}
 XMLReader xmlr(xml_in, "./descendant-or-self::QuarkSourceSinkKey");
 xmlread(xmlr,"SourceTime",source_time,"QuarkSourceSinkKey");
 xmlread(xmlr,"DilutionIndex",dilution_index,"QuarkSourceSinkKey");
}

QuarkSourceSinkHandler::Key::Key(const LaphNoiseInfo& in_noise, int in_time, 
                                 int in_dil_ind)
    : noise(in_noise), source_time(in_time), dilution_index(in_dil_ind) {}


QuarkSourceSinkHandler::Key::Key(const QuarkSourceSinkHandler::Key& in)
    : noise(in.noise), source_time(in.source_time), 
      dilution_index(in.dilution_index) {}


QuarkSourceSinkHandler::Key& QuarkSourceSinkHandler::Key::operator=(
        const QuarkSourceSinkHandler::Key& in)
{
 noise=in.noise;
 source_time=in.source_time;
 dilution_index=in.dilution_index;
 return *this;
}

bool QuarkSourceSinkHandler::Key::operator<(const QuarkSourceSinkHandler::Key& rhs) const
{
 return ((source_time<rhs.source_time)||((source_time==rhs.source_time)
       &&(dilution_index<rhs.dilution_index)||((dilution_index==rhs.dilution_index)
       &&(noise<rhs.noise))));
}


void QuarkSourceSinkHandler::Key::output(XMLWriter& xmlout) const
{
 push(xmlout,"QuarkSourceSinkKey");
 noise.output(xmlout);
 write(xmlout,"SourceTime",source_time);
 write(xmlout,"DilutionIndex",dilution_index);
 pop(xmlout);
}


// *************************************************************************



QuarkSourceSinkHandler::QuarkSourceSinkHandler()
          : dilPtr(0), qactionPtr(0), invertPtr(0), 
            fileMode(0), maxFileNumber(0), m_serpar(QDPIO_SERIAL)
{
 create_handlers();
}

QuarkSourceSinkHandler::QuarkSourceSinkHandler(XMLReader& xml_in)
          : dilPtr(0), qactionPtr(0), invertPtr(0), 
            fileMode(0), maxFileNumber(0), m_serpar(QDPIO_SERIAL)
{
 create_handlers();
 set_info(xml_in);
}


void QuarkSourceSinkHandler::setInfo(XMLReader& xml_in)
{
 clear();
 set_info(xml_in);
}

void QuarkSourceSinkHandler::set_info(XMLReader& xml_in)
{
 if (xml_tag_count(xml_in,"QuarkSourceSinkInfo")!=1){
    QDPIO::cerr << "Bad XML input to QuarkSourceSinkHandler"<<endl;
    QDPIO::cerr << "Expected one <QuarkSourceSinkInfo> tag"<<endl;
    QDP_abort(1);}
 XMLReader xml_info(xml_in, "./descendant-or-self::QuarkSourceSinkInfo");

 xmlread(xml_info,"FileNameStub",fileStub,"QuarkSourceSinkHandler");
 fileStub=tidyString(fileStub);
 if (fileStub.empty()){
    QDPIO::cerr << "Blank file name in QuarkSourceSinkHandler"<<endl;
    QDP_abort(1);}
 xmlread(xml_info,"MaxFileNumber",maxFileNumber,"QuarkSourceSinkHandler");
 fileMode = 0;  // protect mode
 if (xml_tag_count(xml_info,"FileMode")==1){
    string fmode;
    xmlread(xml_info,"FileMode",fmode,"QuarkSourceSinkHandler");
    fmode=tidyString(fmode);
    if (fmode=="overwrite") fileMode=1;}
 m_serpar=QDPIO_SERIAL;
 if (xml_tag_count(xml_info,"IOMode")==1){
    string iomode;
    xmlread(xml_info,"IOMode",iomode,"QuarkSourceSinkHandler");
    iomode=tidyString(iomode);
    if (iomode=="parallel") m_serpar=QDPIO_PARALLEL;}

 invertPtr = new InverterInfo(xml_info);
 if  ((xml_tag_count(xml_info,"StoutLaphSmearing")==1)
    &&(xml_tag_count(xml_info,"GaugeConfigurationInfo")==1)
    &&(xml_tag_count(xml_info,"QuarkInfo")==1)
    &&(xml_tag_count(xml_info,"LaphDilutionScheme")==1))
    set_info_helper(xml_info);
 else
    set_info_from_file(xml_info);                        

 setup_file_map();

 QDPIO::cout << "Info set in QuarkSourceSinkHandler"<<endl;
 QDPIO::cout << outputInfo() <<endl;
}

void QuarkSourceSinkHandler::set_info_helper(XMLReader& xml_info)
{
   // one GaugeConfigurationInfo tag will initialize uPtr and part of smearPtr
 smearPtr->setInfo(xml_info); cout << "smearing info done"<<endl;
 uPtr->setInfo(xml_info);  cout << "gauge info done"<<endl;  
 try{
    dilPtr = new DilutionSchemeInfo(xml_info); cout << "dil done"<<endl;
    qactionPtr = new QuarkInfo(xml_info,uPtr->getInfo()); cout << "quark action"<<endl;}
 catch(...){
    QDPIO::cerr << "allocation problem in QuarkSourceSinkHandler"<<endl;
    QDP_abort(1);}
}


void QuarkSourceSinkHandler::set_info_from_file(XMLReader& xml_info)
{
 for (int suffix=0;suffix<=maxFileNumber;suffix++){
    string filename=make_file_name(suffix);
    if (fileExists(filename)){
       XMLReader header_xml;
       QDPFileReader fsource(header_xml,filename,m_serpar);
       set_info_helper(header_xml);
       fsource.close();
       break;}}
 if (!isInfoSet()){
    QDPIO::cerr << "could not get info from files in QuarkSourceSinkHandler"<<endl;
    QDP_abort(1);}
}


QuarkSourceSinkHandler::~QuarkSourceSinkHandler()
{
 destroy();
}

void QuarkSourceSinkHandler::clear()
{
 destroy();
 create_handlers();
 QDPIO::cout << "QuarkSourceSinkHandler cleared"<<endl;
}

void QuarkSourceSinkHandler::destroy()
{
 try{
    delete dilPtr;
    delete qactionPtr;
    delete invertPtr;}
 catch(...){ QDP_abort(1);}
 dilPtr=0;
 qactionPtr=0;
 invertPtr=0;
 fileStub.clear();
 maxFileNumber=0;
 fileMap.clear();
 fileMode=0;
 clearData();
 destroy_handlers();
}

bool QuarkSourceSinkHandler::isInfoSet() const
{
 return ((uPtr->isInfoSet())&&(smearPtr->isInfoSet())
        &&(dilPtr!=0) && (qactionPtr!=0) && (invertPtr!=0)
        && (!fileStub.empty()));
}

string QuarkSourceSinkHandler::outputInfo() const
{
 if (isInfoSet()){
    ostringstream oss;
    oss << "<QuarkSourceSinkInfo>"<<endl;
    oss << smearPtr->outputInfo() << endl;  // also outputs uPtr-> info
    oss << dilPtr->output() << endl;
    oss << qactionPtr->output() << endl;
    oss << "</QuarkSourceSinkInfo>";
    return oss.str();}
 return "";
}

string QuarkSourceSinkHandler::getHeader() const
{
 ostringstream oss;
 oss << "<QuarkSourceSinkHeader>"<<endl;
 oss << smearPtr->outputInfo() << endl; // also outputs uPtr-> info
 oss << dilPtr->output() << endl;
 oss << qactionPtr->output() << endl;
 oss << "</QuarkSourceSinkHeader>";
 return oss.str();
}

void QuarkSourceSinkHandler::outputInfo(XMLWriter& xmlout) const
{
 if (isInfoSet()){
    push(xmlout,"QuarkSourceSinkInfo");
    smearPtr->outputInfo(xmlout);  // also outputs uPtr-> info
    dilPtr->output(xmlout);
    qactionPtr->output(xmlout);
    pop(xmlout);}
}

void QuarkSourceSinkHandler::getHeader(XMLWriter& xmlout) const
{
 if (isInfoSet()){
    push(xmlout,"QuarkSourceSinkHeader");
    smearPtr->outputInfo(xmlout);  // also outputs uPtr-> info
    dilPtr->output(xmlout);
    qactionPtr->output(xmlout);
    pop(xmlout);}
}

void QuarkSourceSinkHandler::getFileMap(XMLWriter& xmlout) const
{
 if (isInfoSet()){
    push(xmlout,"FileMap");
    for (map<Key,int>::const_iterator it=fileMap.begin();
         it!=fileMap.end();it++){
       push(xmlout,"Entry");
       it->first.output(xmlout);
       write(xmlout,"Suffix",it->second);
       pop(xmlout);}
    pop(xmlout);}
}

const GaugeConfigurationInfo& QuarkSourceSinkHandler::getGaugeConfigurationInfo() const 
{
 if (!isInfoSet()){
    QDPIO::cerr << "error in QuarkSourceSinkHandler:"<<endl;
    QDPIO::cerr << "  must setInfo before calling getGaugeConfigurationInfo"<<endl;
    QDP_abort(1);}
 return uPtr->getInfo();
}

const FieldSmearingInfo& QuarkSourceSinkHandler::getFieldSmearingInfo() const
{
 if (!isInfoSet()){
    QDPIO::cerr << "error in QuarkSourceSinkHandler:"<<endl;
    QDPIO::cerr << "  must setInfo before calling getFieldSmearingInfo"<<endl;
    QDP_abort(1);}
 return smearPtr->getFieldSmearingInfo();
}

const DilutionSchemeInfo& QuarkSourceSinkHandler::getDilutionSchemeInfo() const 
{
 if (!isInfoSet()){
    QDPIO::cerr << "error in QuarkSourceSinkHandler:"<<endl;
    QDPIO::cerr << "  must setInfo before calling getDilutionSchemeInfo"<<endl;
    QDP_abort(1);}
 return *dilPtr;
}

const QuarkInfo& QuarkSourceSinkHandler::getQuarkInfo() const 
{
 if (!isInfoSet()){
    QDPIO::cerr << "error in QuarkSourceSinkHandler:"<<endl;
    QDPIO::cerr << "  must setInfo before calling getQuarkInfo"<<endl;
    QDP_abort(1);}
 return *qactionPtr;
}

const InverterInfo& QuarkSourceSinkHandler::getInverterInfo() const 
{
 if (!isInfoSet()){
    QDPIO::cerr << "error in QuarkSourceSinkHandler:"<<endl;
    QDPIO::cerr << "  must setInfo before calling getInverterInfo"<<endl;
    QDP_abort(1);}
 return *invertPtr;
}

int QuarkSourceSinkHandler::getNumberOfDilutionProjectors() const 
{
 if (!isInfoSet()){
    QDPIO::cerr << "error in QuarkSourceSinkHandler:"<<endl;
    QDPIO::cerr << "  must setInfo before calling getNumberOfDilutionProjectors"<<endl;
    QDP_abort(1);}
 return dilPtr->getNumberOfProjectors(smearPtr->getFieldSmearingInfo());
}



        // compute for all dilution indices and dump out to file
        // specified by "file_index";  compute sources for all source
        // times, sinks for just one source time but all sink times

void QuarkSourceSinkHandler::computeSource(const LaphNoiseInfo& noise)
{
 if (!isInfoSet()){
    QDPIO::cerr << "cannot compute in QuarkSourceSinkHandler until info set"<<endl;
    QDP_abort(1);}

     // set up the Laph eigenvectors
     
 const multi1d<LatticeColorVector>& Vs = set_up_laph_eigenvectors();
 
     // set up the dilution projectors

 vector<DilutionSchemeInfo::Projector> dilProjs;
 dilPtr->getProjectors(smearPtr->getFieldSmearingInfo(),dilProjs);
 
 QDPIO::cout << "Quark source computation for all dilutions, one noise, beginning"<<endl;
 START_CODE();
 StopWatch rolex;
 rolex.start();

     // set the noise vector (save current RNG seed, reset it, make noise, put
     // original seed back into RNG)
 
 Seed curr_seed;
 QDP::RNG::savern(curr_seed);
 QDP::RNG::setrn(noise.getSeed(uPtr->getInfo()));

 int Textent = uPtr->getInfo().getTimeExtent();
 int Tdir = uPtr->getInfo().getTimeDir();
 int Nspin = QDP::Ns;
 int nEigs = smearPtr->getFieldSmearingInfo().getNumberOfLaplacianEigenvectors();
 int Zn = noise.getZNGroup();
  
 multi3d<DComplex> laph_noise(Textent,Nspin,nEigs);
 for (int t=0;t<Textent;t++)
    for (int s=0;s<Nspin;s++)
       for (int v=0;v<nEigs;v++)
          laph_noise(t,s,v) = zN_rng(Zn);
 QDP::RNG::setrn(curr_seed);      //Return the seed to its previous value

 SftMom phases(0, false, Tdir);    // needed for time dilution (masks onto a time slice)
                                   // 0 = max momentum squared, true = avg over equiv mom 

     // loop over dilutions
 
 for (int dil=0;dil<dilProjs.size();dil++){
 
    QDPIO::cout << "doing dilution "<<dil<<endl;

    Key kval(noise,Textent,dil);    // Textent signals a source (not a sink)
    if ((fileMap.find(kval)!=fileMap.end())&&(fileMode!=1)){   // already computed!!
       QDPIO::cout << "warning: computeQuarkSource already computed..."
                   << "skip re-computing since fileMode is not overwrite"<<endl;}
    else{
    
        // get the lists of which spins and which eigenvectors are
        // "on" for this dilution projector
        
    int findex=first_available_suffix();
    const list<int>& on_spins=dilProjs[dil].onSpinIndices();
    const list<int>& on_eigs=dilProjs[dil].onEigvecIndices();
 
        //  initialize output field for sources
    LatticeFermion source = zero;
    for (list<int>::const_iterator vmask= on_eigs.begin(); vmask!=on_eigs.end(); vmask++){
    QDPIO::cout << "new vmask"<<endl;
       LatticeSpinVector sv = zero;
       for (list<int>::const_iterator smask= on_spins.begin(); smask!=on_spins.end(); smask++){
       QDPIO::cout << "new smask"<<endl;
          LatticeComplex temp = zero;
          for (int t0=0;t0<Textent;t0++){
             QDPIO::cout << "t0 = "<<t0<<endl;
             temp[phases.getSet()[t0]] = laph_noise(t0,*smask,*vmask);}
          pokeSpin(sv,temp,*smask);}
       source += sv * Vs[*vmask];
       }
    QDPIO::cout << "done computation: writing to file"<<endl;

        // output this dilution, all t0, to file

    string fileName=make_file_name(findex);
    XMLBufferWriter recordheader;
    kval.output(recordheader);
    XMLBufferWriter fileheader;
    getHeader(fileheader);

    if (filewrite(fileName,fileheader,recordheader,source)){
       QDPIO::cout << "write to file done"<<endl;
       fileMap.insert(make_pair(kval,findex));}
    }}

 rolex.stop();
 QDPIO::cout << "computeQuarkSource: one noise, all dilutions total time = "
             << rolex.getTimeInSeconds() << " secs" << endl;
 QDPIO::cout << "ran successfully" << endl;
}

 
   // ***************************************************************
   
void QuarkSourceSinkHandler::computeSink(const LaphNoiseInfo& noise, 
                                         int source_time)
{
 if (!isInfoSet()){
    QDPIO::cerr << "cannot compute in QuarkSourceSinkHandler until info set"<<endl;
    QDP_abort(1);}

 int Textent = uPtr->getInfo().getTimeExtent();
 if ((source_time<0)||(source_time>=Textent)){
    QDPIO::cerr << "invalid source time "<<source_time<<" for compute in QuarkSourceSinkHandler"<<endl;
    QDPIO::cerr << " skipping this computation"<<endl;
    return;}
 
     // set up the Laph eigenvectors
     
 const multi1d<LatticeColorVector>& Vs = set_up_laph_eigenvectors();

     // set up the dilution projectors

 vector<DilutionSchemeInfo::Projector> dilProjs;
 dilPtr->getProjectors(smearPtr->getFieldSmearingInfo(),dilProjs);
 
 QDPIO::cout << "Quark sink computation for all dilutions, one noise,"
             << " one source time beginning"<<endl;
 QDPIO::cout << " Source time = "<<source_time<<endl;
 START_CODE();
 StopWatch rolex;
 rolex.start();

     // set the noise vector (save current RNG seed, reset it, make noise, put
     // original seed back into RNG)
 
 Seed curr_seed;
 QDP::RNG::savern(curr_seed);
 QDP::RNG::setrn(noise.getSeed(uPtr->getInfo()));

 int Tdir = uPtr->getInfo().getTimeDir();
 int Nspin = QDP::Ns;
 int nEigs = smearPtr->getFieldSmearingInfo().getNumberOfLaplacianEigenvectors();
 int Zn = noise.getZNGroup();
  
 multi3d<DComplex> laph_noise(Textent,Nspin,nEigs);
 for (int t=0;t<Textent;t++)
    for (int s=0;s<Nspin;s++)
       for (int v=0;v<nEigs;v++)
          laph_noise(t,s,v) = zN_rng(Zn);
 QDP::RNG::setrn(curr_seed);      //Return the seed to its previous value

 SftMom phases(0, false, Tdir);  // needed for time dilution (masks onto a time slice)
                                 // 0 = max momentum squared, true = avg over equiv mom 
                                                     // rotate to DeGrand-Rossi, then 
 SpinMatrix SrcRotate = Gamma(8) * DiracToDRMat();   //  multiply by gamma_4
 SpinMatrix SnkRotate = adj(DiracToDRMat());    // rotate back to Dirac-Pauli

 string fermact_xml = qactionPtr->output();
 QDPIO::cout << "fermact_xml = "<<fermact_xml<<endl;
 string fermact_id = qactionPtr->getActionId();

   // Typedefs to save typing
 typedef LatticeFermion               T;
 typedef multi1d<LatticeColorMatrix>  P;
 typedef multi1d<LatticeColorMatrix>  Q;

 GroupXML_t solverInfo;
 solverInfo.xml =  invertPtr->output();
 solverInfo.id = invertPtr->getId();
 solverInfo.path = "//InvertParam";
 
 QDPIO::cout << "inverter xml = "<<solverInfo.xml<<endl;
 
   // Initialize fermion action
   
 istringstream xml_s(fermact_xml);
 XMLReader fermacttop0(xml_s);
 XMLReader fermacttop(fermacttop0,"/");   // due to XMLReader bug
 QDPIO::cout << "FermAct = " << fermact_id << endl;

 Handle< FermionAction<T,P,Q> > S_f; 
 Handle< FermState<T,P,Q> > state;
 Handle< SystemSolver<LatticeFermion> > PP;

 try{
    QDPIO::cout << "createObject"<<endl;
    S_f=TheFermionActionFactory::Instance().createObject(fermact_id,
                                               fermacttop,"//FermionAction");
    QDPIO::cout << "createState"<<endl;         
    state=S_f->createState(uPtr->getData());
    PP = S_f->qprop(state,solverInfo); QDPIO::cout << "createPP"<<endl;}
 catch(const string& err){
    QDPIO::cerr << " Fermion action and inverter could not be initialized"
                << " in QuarkSourceSinkHandler"<<endl;}
 QDPIO::cout << "Suitable factory found: compute all the quark props" << endl;
 
     // loop over dilutions
 
 for (int dil=0;dil<dilProjs.size();dil++){

    QDPIO::cout << "Starting dilution "<<dil<<endl;
    Key kval(noise,source_time,dil);
    if ((fileMap.find(kval)!=fileMap.end())&&(fileMode!=1)){   // already computed!!
       QDPIO::cout << "warning: computeQuarkSink already computed..."
                   << "skip re-computing since fileMode not overwrite"<<endl;}
    else{
    
        // get the lists of which spins and which eigenvectors are
        // "on" for this dilution projector
        
    int findex=first_available_suffix();
    const list<int>& on_spins=dilProjs[dil].onSpinIndices();
    const list<int>& on_eigs=dilProjs[dil].onEigvecIndices();

        //  initialize output field for sources
    LatticeFermion source = zero;
    for (list<int>::const_iterator vmask= on_eigs.begin(); vmask!=on_eigs.end(); vmask++){
       LatticeSpinVector sv = zero;
       for (list<int>::const_iterator smask= on_spins.begin(); smask!=on_spins.end(); smask++){
          LatticeComplex temp = zero;
          temp[phases.getSet()[source_time]] = laph_noise(source_time,*smask,*vmask);
          pokeSpin(sv,temp,*smask);}
       source[phases.getSet()[source_time]] += sv * Vs[*vmask];
       }
       
    source = SrcRotate * source;  // rotate to DeGrand-Rossi, mult by gamma_4
    LatticeFermion phi;
    QDPIO::cout << "source created...starting inversion"<<endl;
             // now do the inversion
    SystemSolverResults_t res = (*PP)(phi, source);

//   int res.n_count   Real res.resid

                   // sink = Vs * Vs^dagger * phi
    LatticeFermion sink;
    for (int s=0;s<Nspin;s++){
       LatticeColorVector phi_s = peekSpin(phi, s);
       LatticeColorVector sink_s = zero;
       for (int n=0;n<nEigs;n++){
          LatticeComplex tmp = localInnerProduct( Vs[n], phi_s );
          multi2d<DComplex> t_sum = phases.sft(tmp);
          for (int t=0;t<Textent;t++){
             sink_s[phases.getSet()[t]] += Vs[n] * t_sum[0][t];}}
       pokeSpin(sink, sink_s, s);}
    sink = SnkRotate * sink;         // rotate back to Dirac-Pauli

        // output this dilution to file

    string fileName=make_file_name(findex);
    XMLBufferWriter recordheader;
    kval.output(recordheader);
    XMLBufferWriter fileheader;
    getHeader(fileheader);

    if (filewrite(fileName,fileheader,recordheader,sink)){
       QDPIO::cout << "write to file done"<<endl;
       fileMap.insert(make_pair(kval,findex));}
    }}

 rolex.stop();
 QDPIO::cout << "computeQuarkSink: one noise, all dilutions, one source time, total time = "
             << rolex.getTimeInSeconds() << " secs" << endl;
 QDPIO::cout << "ran successfully" << endl;
}    

 
void QuarkSourceSinkHandler::computeSink(XMLReader& xml_in)
{
 LaphNoiseInfo noise(xml_in);
 int source_time;
 xmlread(xml_in,"source_time",source_time,"QuarkSourceSinkHandler");
 computeSink(noise,source_time);
}


// **********************************************************************


void QuarkSourceSinkHandler::setSources(const LaphNoiseInfo& noise, 
                                        int dilution_index)
{
 int Textent = uPtr->getInfo().getTimeExtent();
 setSink(noise,Textent,dilution_index);
}

void QuarkSourceSinkHandler::setSink(const LaphNoiseInfo& noise, 
                                     int source_time, int dilution_index)
{
 if (!isInfoSet()){
    QDPIO::cerr << "cannot set data in QuarkSourceSinkHandler until info set"<<endl;
    QDP_abort(1);}

 Key kval(noise,source_time,dilution_index);
 if (m_storage.find(kval)!=m_storage.end()){
    QDPIO::cout << "this element was already set in QuarkSourceSinkHandler"<<endl;
    return;}

 map<Key,int>::const_iterator it=fileMap.find(kval);
 if (it==fileMap.end()){
    QDPIO::cerr << "cannot setSink/Source in QuarkSourceSinkHandler since"
                << " not in any of the files"<<endl;
    QDP_abort(1);}
 int findex=it->second;
 string filename=make_file_name(findex);
 LatticeFermion *dataptr;
 try{
    dataptr=new LatticeFermion;}
 catch(...){
    QDPIO::cerr << "could not allocation memory for set in QuarkSourceSinkHandler"<<endl;
    QDP_abort(1);}
 XMLReader fileheader,recordheader;
 fileread(filename,fileheader,recordheader,*dataptr);
 m_storage.insert(make_pair(kval,dataptr));
}




const LatticeFermion& QuarkSourceSinkHandler::getSources(const LaphNoiseInfo& noise,
                                                 int dilution_index)
{
 int Textent = uPtr->getInfo().getTimeExtent();
 return getSink(noise,Textent,dilution_index);
}

const LatticeFermion& QuarkSourceSinkHandler::getSink(const LaphNoiseInfo& noise,
                              int source_time, int dilution_index)
{
 Key kval(noise,source_time,dilution_index);
 map<Key,LatticeFermion*>::const_iterator it=m_storage.find(kval);
 if (it==m_storage.end()){
    setSink(noise,source_time,dilution_index);
    it=m_storage.find(kval);
    if (it==m_storage.end()){
       QDPIO::cerr << "get in QuarkSourceSinkHandler failed...not in any of files"<<endl;
       QDP_abort(1);}}
 return *(it->second);    
}

        // remove from internal memory

void QuarkSourceSinkHandler::removeSourceData(const LaphNoiseInfo& noise, 
                                              int dilution_index)
{
 int Textent = uPtr->getInfo().getTimeExtent();
 removeSinkData(noise,Textent,dilution_index);
}

void QuarkSourceSinkHandler::removeSinkData(const LaphNoiseInfo& noise, 
                                    int source_time, int dilution_index)
{
 Key kval(noise,source_time,dilution_index);
 map<Key,LatticeFermion*>::iterator it=m_storage.find(kval);
 if (it!=m_storage.end()){
    delete it->second;
    m_storage.erase(kval);}
}


void QuarkSourceSinkHandler::clearData()
{
 for (map<Key,LatticeFermion*>::iterator it=m_storage.begin();
                   it!=m_storage.end();it++)
    delete it->second;
 m_storage.clear();
}


void QuarkSourceSinkHandler::create_handlers()
{
 try{
    uPtr = new GaugeConfigurationHandler;
    smearPtr = new FieldSmearingHandler;}
 catch(...){
    QDPIO::cerr << "allocation problem in QuarkSourceSinkHandler"<<endl;
    QDP_abort(1);}
}

void QuarkSourceSinkHandler::destroy_handlers()
{
 try {
    delete uPtr;
    delete smearPtr;}
 catch(...) {QDP_abort(1);} 
 uPtr=0;
 smearPtr=0;
}


    //  Existing files are checked for consistent header
    //  information.  The fileMap is then created.

void QuarkSourceSinkHandler::setup_file_map()
{
 fileMap.clear();
 XMLBufferWriter header_xml;
 getHeader(header_xml);
 string header_info(header_xml.str());

 for (int suffix=0;suffix<=maxFileNumber;suffix++){
    string filename=make_file_name(suffix);
    if (fileExists(filename)){
       XMLReader fheader_xml;
       QDPFileReader fsource(fheader_xml,filename,m_serpar);
       stringstream hstr;
       fheader_xml.print(hstr);
       string hdinf(hstr.str());
       if (!headerMatch(hdinf,header_info)){
          QDPIO::cerr << "header info in file "<<filename
                      << " does not match info in QuarkSourceSinkHandler"<<endl;
          QDP_abort(1);}
       XMLReader rec_xml;
       fsource.read(rec_xml);
       Key kval(rec_xml);
       map<Key,int>::iterator it=fileMap.find(kval);
       if (it!=fileMap.end()){
           QDPIO::cerr << "duplicate keys in fileMap in QuarkSourceSinkHandler"<<endl;
           QDPIO::cerr << " ... too confusing to continue"<<endl;
           QDPIO::cerr << "suffix "<<suffix<<" and suffix "<<it->second
                       <<" have same key"<<endl;
           XMLBufferWriter xmlout;
           kval.output(xmlout);
           QDPIO::cerr << " Key is "<<xmlout.str()<<endl;
           QDP_abort(1);}
       fileMap.insert(make_pair(kval,suffix));
       fsource.close();}}
}


const multi1d<LatticeColorVector>& QuarkSourceSinkHandler::set_up_laph_eigenvectors()
{
 if (!(smearPtr->isLaphEigenvectorsSet())){
    smearPtr->setLaphEigenvectors();  // try getting from named object space first
    if (!smearPtr->isLaphEigenvectorsSet()){
       smearPtr->computeLaphEigenvectors();   // if not in named object space, then compute!!
       if (!smearPtr->isLaphEigenvectorsSet()){
          QDPIO::cerr << "could not compute or set Laph eigenvectors so" << endl
                      << " so QuarkSourceSinkComputation is not possible"<<endl;
          QDP_abort(1);}}}
 return smearPtr->getLaphEigenvectors();
}

bool QuarkSourceSinkHandler::filewrite(const std::string& fileName, 
                                       XMLBufferWriter& fileHeader,
                                       XMLBufferWriter& recordHeader,
                                       const LatticeFermion& data)
{ 
 QDPFileWriter fout(fileHeader,fileName,QDPIO_SINGLEFILE,m_serpar);
 fout.write(recordHeader,data);
 bool flag=true;
 if (fout.bad()){
    QDPIO::cerr << "Error occurred while writing QuarkSourceSink to file"
                << fileName << endl;
    flag=false;}
 fout.close();
 return flag;
}

 
void QuarkSourceSinkHandler::fileread(const std::string& fileName, 
                                      XMLReader& fileHeader,
                                      XMLReader& recordHeader)
{
 QDPFileReader fin(fileHeader,fileName,m_serpar);
 fin.read(recordHeader);
 bool bad=fin.bad();
 fin.close();
 if (bad){
    QDPIO::cerr << "Error occurred while reading QuarkSourceSink from file"
                << fileName << endl;
    QDP_abort(1);}
}

void QuarkSourceSinkHandler::fileread(const std::string& fileName, 
                                      XMLReader& fileHeader,
                                      XMLReader& recordHeader, 
                                      LatticeFermion& data)
{
 QDPFileReader fin(fileHeader,fileName,m_serpar);
 fin.read(recordHeader,data);
 bool bad=fin.bad();
 fin.close();
 if (bad){
    QDPIO::cerr << "Error occurred while reading QuarkSourceSink from file"
                << fileName << endl;
    QDP_abort(1);}
}

void QuarkSourceSinkHandler::setParallelIO() 
{
 m_serpar=QDPIO_PARALLEL;
}

void QuarkSourceSinkHandler::setSerialIO()
{
 m_serpar=QDPIO_SERIAL;
}

string QuarkSourceSinkHandler::make_file_name(int suffix)
{
 stringstream fs;
 fs << fileStub << "." << suffix;
 return fs.str();
}

int QuarkSourceSinkHandler::first_available_suffix()
{
 for (int suffix=0;suffix<=maxFileNumber;suffix++){
    string filename=make_file_name(suffix);
    if (!fileExists(filename)) return suffix;}
 QDPIO::cerr << "no suffix numbers are available for writing"<<endl;
 QDPIO::cerr << " ... increase maxFilenumber"<<endl;
 QDP_abort(1);
}

// ***************************************************************
  }
}
 
