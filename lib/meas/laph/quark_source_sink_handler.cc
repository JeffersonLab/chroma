#include "quark_source_sink_handler.h"
#include "meas/sources/zN_src.h"
#include "util/ft/sftmom.h"
#include "fermact.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "util/ferm/diractodr.h"


namespace Chroma {
  namespace LaphEnv {

//  use of BinaryStoreDB in  qdp++/include/qdp_db_imp.h
//  with SerialDBKey and SerialDBData in chroma/lib/util/ferm/key_val_db.h

// ****************************************************************
// *                                                              *
// *  "QuarkSourceSinkHandler" handles access to the smeared gauge  *
// *  field and the eigenvectors of the Laplacian used in the     *
// *  quark field Laplacian Heaviside (Laph) smearing.            *
// *  See declarations below for available member functions.      *
// *                                                              *
// *  Synopsis of important tasks:                                *
// *     -- computes smeared gauge field (locally stored)         *
// *     -- computes Laplacian eigenvectors (NamedObjMap stored)  *
// *                                                              *
// *                                                              *
// *  Internally maintains a QuarkSourceSinkInfo.  Requires         *
// *  connection with an externally maintained                    *
// *  GaugeConfigurationHandler.                                  *
// *                                                              *
// ****************************************************************


   // use of a pointer to LaphNoiseInfo is dangerous!!  But this
   // class is only used internally so criminals have no access to it.
   
QuarkSourceSinkHandler::Key::Key()
    : noise(0), source_time(0), dilution_index(0) {}

QuarkSourceSinkHandler::Key::~Key()
{
 try {delete noise;} catch(...){QDP_abort(1);}
}

QuarkSourceSinkHandler::Key::Key(const LaphNoiseInfo& in_noise, int in_time, 
                                 int in_dil_ind)
    : source_time(in_time), dilution_index(in_dil_ind)
{
 try{
    noise=new LaphNoiseInfo(in_noise);}
 catch(...){
     QDPIO::cerr << "allocation error in QuarkSourceSinkHandler::Key"<<endl;
     QDP_abort(1);}
}

QuarkSourceSinkHandler::Key::Key(const QuarkSourceSinkHandler::Key& in)
    : source_time(in.source_time), 
      dilution_index(in.dilution_index)
{
 try{
    noise=new LaphNoiseInfo(*(in.noise));}
 catch(...){
     QDPIO::cerr << "allocation error in QuarkSourceSinkHandler::Key"<<endl;
     QDP_abort(1);}
}

QuarkSourceSinkHandler::Key& QuarkSourceSinkHandler::Key::operator=(
        const QuarkSourceSinkHandler::Key& in)
{
 if (noise!=in.noise){
    try {delete noise;} catch(...){QDP_abort(1);}
    try{
       noise=new LaphNoiseInfo(*(in.noise));}
    catch(...){
        QDPIO::cerr << "allocation error in QuarkSourceSinkHandler::Key"<<endl;
        QDP_abort(1);}}
 source_time=in.source_time;
 dilution_index=in.dilution_index;
 return *this;
}

bool QuarkSourceSinkHandler::Key::operator<(const QuarkSourceSinkHandler::Key& rhs) const
{
 return ((source_time<rhs.source_time)||((source_time==rhs.source_time)
       &&(dilution_index<rhs.dilution_index)||((dilution_index==rhs.dilution_index)
       &&(*noise<*rhs.noise))));
}


void QuarkSourceSinkHandler::Key::binaryWrite(BinaryWriter& out) const
{
 try{
    noise->binaryWrite(out);
    write(out,source_time);
    write(out,dilution_index);}
 catch(...){
    QDPIO::cerr << "binary write error in QuarkSourceSinkHandler::Key"<<endl;
    QDP_abort(1);}
}

void QuarkSourceSinkHandler::Key::binaryRead(BinaryReader& in)
{
 try{
    noise->binaryRead(in);
    read(in,source_time);
    read(in,dilution_index);}
 catch(...){
    QDPIO::cerr << "binary read error in QuarkSourceSinkHandler::Key"<<endl;
    QDP_abort(1);}
}


void write(BinaryWriter& out, const QuarkSourceSinkHandler::Key& key_out )
{
 key_out.binaryWrite(out);
}

void read(BinaryReader& in, QuarkSourceSinkHandler::Key& key_in )
{
 key_in.binaryRead(in);
}


// *************************************************************************



QuarkSourceSinkHandler::QuarkSourceSinkHandler()
          : uPtr(0), smearPtr(0), dilPtr(0), qactionPtr(0),
            invertPtr(0), fileMode(0) {}

QuarkSourceSinkHandler::QuarkSourceSinkHandler(XMLReader& xml_info)
{
 setInfo(xml_info);
}


void QuarkSourceSinkHandler::setInfo(XMLReader& xml_info)
{
 clear();
 set_info(xml_info);
 setFileNames(xml_info);
 setup_file_map();
 QDPIO::cout << "Info set in QuarkSourceSinkHandler"<<endl;
 QDPIO::cout << outputInfo() <<endl;
}

void QuarkSourceSinkHandler::set_info(XMLReader& xml_info)
{
 create_handlers();
 smearPtr->setInfo(xml_info);
 uPtr->setInfo(xml_info);
 try{
    dilPtr = new DilutionSchemeInfo(xml_info);
    qactionPtr = new QuarkActionInfo(xml_info);
    invertPtr = new InverterInfo(xml_info);}
 catch(...){
    QDPIO::cerr << "allocation problem in QuarkSourceSinkHandler"<<endl;
    QDP_abort(1);}
}

void QuarkSourceSinkHandler::setInfoFromFile(XMLReader& xml_info)
{
 clear();
 setFileNames(xml_info);
 for (int k=0;k<fileNames.size();k++)
    if (fileExists(fileNames[k])){
       BinaryStoreDB<SerialDBKey<Key>,SerialDBData<LatticeFermion> > fsource;
       fsource.open(fileNames[k], O_RDONLY , 0664);
       string headerInfo;
       fsource.getUserdata(headerInfo);
       stringstream oss; oss << headerInfo;
       XMLReader xmlr(oss);
       set_info(xmlr);
       fsource.close();
       break;}
 if (!isInfoSet()){
    QDPIO::cerr << "could not get info from files in QuarkSourceSinkHandler"<<endl;
    QDP_abort(1);}
 setup_file_map();
 QDPIO::cout << "Info set in QuarkSourceSinkHandler"<<endl;
 QDPIO::cout << outputInfo() <<endl;
}


QuarkSourceSinkHandler::~QuarkSourceSinkHandler()
{
 clear();
}

void QuarkSourceSinkHandler::clear()
{
 try{
    delete dilPtr;
    delete qactionPtr;
    delete invertPtr;}
 catch(...){ QDP_abort(1);}
 dilPtr=0;
 qactionPtr=0;
 invertPtr=0;
 fileNames.clear();
 fileMap.clear();
 fileMode=0;
 clearData();
 destroy_handlers();
 QDPIO::cout << "QuarkSourceSinkHandler cleared"<<endl;
}

bool QuarkSourceSinkHandler::isInfoSet() const
{
 return ((uPtr->isInfoSet())&&(smearPtr->isInfoSet())
        &&(dilPtr!=0) && (qactionPtr!=0) && (invertPtr!=0));
}

string QuarkSourceSinkHandler::outputInfo() const
{
 ostringstream oss;
 oss << uPtr->outputInfo() << endl;
 oss << smearPtr->outputInfo() << endl;
 oss << dilPtr->output() << endl;
 oss << qactionPtr->output() << endl;
 return oss.str();
}

        // compute for all dilution indices and dump out to file
        // specified by "file_index";  compute sources for all source
        // times, sinks for just one source time but all sink times

void QuarkSourceSinkHandler::computeSource(XMLReader& xml_in, int file_index)
{
 computeSource(LaphNoiseInfo(xml_in),file_index);
}



void QuarkSourceSinkHandler::computeSource(const LaphNoiseInfo& noise, 
                                           int file_index)
{
 if (!isInfoSet()){
    QDPIO::cerr << "cannot compute in QuarkSourceSinkHandler until info set"<<endl;
    QDP_abort(1);}
 if (fileMode==0){
    QDPIO::cerr << "cannot compute in QuarkSourceSinkHandler with read-only file mode"<<endl;
    QDP_abort(1);}

 int findex = file_index_helper(file_index);

     // set up the Laph eigenvectors
     
 const multi1d<LatticeColorVector>& Vs = set_up_laph_eigenvectors();
 
     // set up the dilution projectors

 vector<DilutionSchemeInfo::Projector> dilProjs;
 dilPtr->getProjectors(smearPtr->getInfo(),dilProjs);
 
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
 int nEigs = smearPtr->getInfo().getNumberOfLaplacianEigenvectors();
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

    Key kval(noise,Textent,dil);    // Textent signals a source (not a sink)
    if (fileMap.find(kval)!=fileMap.end()){   // already computed!!
       QDPIO::cout << "warning: computeQuarkSource already computed...skip re-computing"<<endl;}
    else{
    
        // get the lists of which spins and which eigenvectors are
        // "on" for this dilution projector
        
    const list<int>& on_spins=dilProjs[dil].onSpinIndices();
    const list<int>& on_eigs=dilProjs[dil].onEigvecIndices();
 
        //  initialize output field for sources
    LatticeFermion source = zero;
    for (int t0=0;t0<Textent;t0++){
       for (list<int>::const_iterator smask= on_spins.begin(); smask!=on_spins.end(); smask++)
       for (list<int>::const_iterator vmask= on_eigs.begin(); vmask!=on_eigs.end(); vmask++){
          Complex curr_n = laph_noise(t0,*smask,*vmask);
          LatticeFermion addto = zero;
          pokeSpin(addto[phases.getSet()[t0]], curr_n * Vs[*vmask], *smask);
          source[phases.getSet()[t0]] += addto;}
       }
       
        // output this dilution, all t0, to file

    BinaryStoreDB<SerialDBKey<Key>,SerialDBData<LatticeFermion> > fsource;
    fsource.open(fileNames[findex], O_RDWR , 0664);
    SerialDBKey<Key> dbkey;
    dbkey.key() = kval;
    SerialDBData<LatticeFermion> dbdata;
    dbdata.data() = source;
    fsource.insert(dbkey,dbdata);

    fileMap.insert(make_pair(kval,findex));
    fsource.close();}
    }

 rolex.stop();
 QDPIO::cout << "computeQuarkSource: one noise, all dilutions total time = "
             << rolex.getTimeInSeconds() << " secs" << endl;
 QDPIO::cout << "ran successfully" << endl;
}    
 
 
   // ***************************************************************
   
   
 
void QuarkSourceSinkHandler::computeSink(const LaphNoiseInfo& noise, int source_time, 
                                         int file_index)
{
 if (!isInfoSet()){
    QDPIO::cerr << "cannot compute in QuarkSourceSinkHandler until info set"<<endl;
    QDP_abort(1);}
 if (fileMode==0){
    QDPIO::cerr << "cannot compute in QuarkSourceSinkHandler with read-only file mode"<<endl;
    QDP_abort(1);}

 int findex = file_index_helper(file_index);
 int Textent = uPtr->getInfo().getTimeExtent();
 if ((source_time<0)||(source_time>=Textent)){
    QDPIO::cerr << "invalid source time "<<source_time<<" for compute in QuarkSourceSinkHandler"<<endl;
    QDPIO::cerr << " skipping this computation"<<endl;
    return;}

     // set up the Laph eigenvectors
     
 const multi1d<LatticeColorVector>& Vs = set_up_laph_eigenvectors();
 
     // set up the dilution projectors

 vector<DilutionSchemeInfo::Projector> dilProjs;
 dilPtr->getProjectors(smearPtr->getInfo(),dilProjs);
 
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
 int nEigs = smearPtr->getInfo().getNumberOfLaplacianEigenvectors();
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
 string fermact_id = qactionPtr->getActionId();

   // Typedefs to save typing
 typedef LatticeFermion               T;
 typedef multi1d<LatticeColorMatrix>  P;
 typedef multi1d<LatticeColorMatrix>  Q;

 GroupXML_t solverInfo;
 solverInfo.xml =  invertPtr->output();
 solverInfo.id = invertPtr->getId();
 solverInfo.path = "//InvertParam";
 
   // Initialize fermion action
   
 istringstream xml_s(fermact_xml);
 XMLReader fermacttop(xml_s);
 QDPIO::cout << "FermAct = " << fermact_id << endl;

 Handle< FermionAction<T,P,Q> > S_f; 
 Handle< FermState<T,P,Q> > state;
 Handle< SystemSolver<LatticeFermion> > PP;

 try{
    S_f=TheFermionActionFactory::Instance().createObject(fermact_id,
                                   fermacttop,"//FermionAction");
    state=S_f->createState(uPtr->getData());
    PP = S_f->qprop(state,solverInfo);}
 catch(const string& err){
    QDPIO::cerr << " Fermion action and inverter could not be initialized"
                << " in QuarkSourceSinkHandler"<<endl;}
 QDPIO::cout << "Suitable factory found: compute all the quark props" << endl;
 
     // loop over dilutions
 
 for (int dil=0;dil<dilProjs.size();dil++){

    Key kval(noise,source_time,dil);
    if (fileMap.find(kval)!=fileMap.end()){   // already computed!!
       QDPIO::cout << "warning: computeQuarkSink already computed...skip re-computing"<<endl;}
    else{
    
        // get the lists of which spins and which eigenvectors are
        // "on" for this dilution projector
        
    const list<int>& on_spins=dilProjs[dil].onSpinIndices();
    const list<int>& on_eigs=dilProjs[dil].onEigvecIndices();
 
        //  initialize output field for sources
    LatticeFermion source = zero;
    for (list<int>::const_iterator smask= on_spins.begin(); smask!=on_spins.end(); smask++)
    for (list<int>::const_iterator vmask= on_eigs.begin(); vmask!=on_eigs.end(); vmask++){
       Complex curr_n = laph_noise(source_time,*smask,*vmask);
       LatticeFermion addto = zero;
       pokeSpin(addto[phases.getSet()[source_time]], curr_n * Vs[*vmask], *smask);
       source[phases.getSet()[source_time]] += addto;}
       
    source = SrcRotate * source;  // rotate to DeGrand-Rossi, mult by gamma_4
    LatticeFermion phi;
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

        // output this dilution, all t0, to file

    BinaryStoreDB<SerialDBKey<Key>,SerialDBData<LatticeFermion> > fsink;
    fsink.open(fileNames[findex], O_RDWR , 0664);
    SerialDBKey<Key> dbkey;
    dbkey.key() = kval;
    SerialDBData<LatticeFermion> dbdata;
    dbdata.data() = sink;
    fsink.insert(dbkey,dbdata);

    fileMap.insert(make_pair(kval,findex));
    fsink.close();}
    }

 rolex.stop();
 QDPIO::cout << "computeQuarkSink: one noise, all dilutions, one source time, total time = "
             << rolex.getTimeInSeconds() << " secs" << endl;
 QDPIO::cout << "ran successfully" << endl;
}    
 
 
void QuarkSourceSinkHandler::computeSink(XMLReader& xml_in, int file_index)
{
 LaphNoiseInfo noise(xml_in);
 int source_time;
 read(xml_in,"//source_time",source_time);
 computeSink(noise,source_time,file_index);
}


// **********************************************************************

    // Database does not appear to do erase yet.

/*
void QuarkSourceSinkHandler::eraseSink(XMLReader& xml_in)
{
 LaphNoiseInfo noise(xml_in);
 int source_time;
 read(xml_in,"//source_time",source_time);
 eraseSink(noise,source_time);
}

void QuarkSourceSinkHandler::eraseSink(const LaphNoiseInfo& noise, int source_time)
{
 if (fileMode==0){
    QDPIO::cerr << "cannot erase in QuarkSourceSinkHandler...read-only"<<endl;
    return;}
 int Textent = uPtr->getInfo().getTimeExtent();
 if ((source_time<0)||(source_time>=Textent)){
    QDPIO::cerr << "invalid source time "<<source_time
                <<" to erase in QuarkSourceSinkHandler"<<endl;
    return;}

 int ndil=dilPtr->getNumberOfProjectors(smearPtr->getInfo());

 for (int dil=0;dil<ndil;dil++){

    Key kval(noise,source_time,dil);
    map<Key,int>::iterator it=fileMap.find(kval);
    if (it!=fileMap.end()){
       int findex=it->second;
       BinaryStoreDB<SerialDBKey<Key>,SerialDBData<LatticeFermion> > fsink;
       fsink.open(fileNames[findex], O_RDWR , 0664);
       SerialDBKey<Key> dbkey;
       dbkey.key() = kval;
       //fsink.erase(dbkey);
       fileMap.erase(kval);}
    }

 QDPIO::cout << "QuarkSourceSinkHandler: erased an entry" << endl;
}
 


void QuarkSourceSinkHandler::eraseSource(XMLReader& xml_in)
{
 eraseSource(LaphNoiseInfo(xml_in));
}
 
void QuarkSourceSinkHandler::eraseSource(const LaphNoiseInfo& noise)
{
 int Textent = uPtr->getInfo().getTimeExtent();
 eraseSink(noise,Textent);
}

*/
// ****************************************************************


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
    QDPIO::cerr << "cannot setSink in QuarkSourceSinkHandler since"
                << " not in any of the files"<<endl;
    QDP_abort(1);}
 int findex=it->second;
 BinaryStoreDB<SerialDBKey<Key>,SerialDBData<LatticeFermion> > fdata;
 fdata.open(fileNames[findex], O_RDONLY , 0664);
 SerialDBKey<Key> dbkey;
 dbkey.key() = kval;
 SerialDBData<LatticeFermion> *dataptr;
 try{
    dataptr=new SerialDBData<LatticeFermion>;}
 catch(...){
    QDPIO::cerr << "could not allocation memory for set in QuarkSourceSinkHandler"<<endl;
    QDP_abort(1);} 
 fdata.get(dbkey,*dataptr);
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
 map<Key,SerialDBData<LatticeFermion>*>::const_iterator it=m_storage.find(kval);
 if (it==m_storage.end()){
    setSink(noise,source_time,dilution_index);
    it=m_storage.find(kval);
    if (it==m_storage.end()){
       QDPIO::cerr << "get in QuarkSourceSinkHandler failed...not in any of files"<<endl;
       QDP_abort(1);}}
 return it->second->data();    
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
 map<Key,SerialDBData<LatticeFermion>*>::iterator it=m_storage.find(kval);
 if (it!=m_storage.end()){
    delete it->second;
    m_storage.erase(kval);}
}


void QuarkSourceSinkHandler::clearData()
{
 for (map<Key,SerialDBData<LatticeFermion>*>::iterator it=m_storage.begin();
                   it!=m_storage.end();it++)
    delete it->second;
 m_storage.clear();
}


       // storage and/or references to internal data



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


  //  This sets the fileName vector from the XML input.
  //  It checks that all file names are different (no repeats allowed)
  //  and aborts if no file names are given

void QuarkSourceSinkHandler::setFileNames(XMLReader& xml_in)
{
     // read list of filenames
 list<string> flist;
 fileMode = 0;  // read only
 try{
    XMLReader xmlr(xml_in,"//file_name_list");
    if (xmlr.count("//file_mode")==1){
       string fmode;
       read(xmlr,"//file_mode",fmode);
       fmode=tidyString(fmode);
       if (fmode=="must_exist") fileMode=1;
       else if (fmode=="create_if_not_there") fileMode=2;}
    int ntags=xmlr.count("//file_name");
    string fname;
    for (int i=1;i<=ntags;i++){
       ostringstream element_xpath;
       element_xpath << "//filename[" << i << "]";
       read(xmlr, element_xpath.str(),fname);
       flist.push_back(tidyString(fname));}
    }
 catch(const string& err){
    QDPIO::cerr << "could not read filenames in QuarkSourceSinkHandler"<<endl;
    QDP_abort(1);}

      // check for repeated file names

 fileNames.resize(flist.size());
 copy(flist.begin(),flist.end(),fileNames.begin());

 flist.sort();
 flist.unique();
 if (flist.size()!=fileNames.size()){
    QDPIO::cerr << "repeated filenames in QuarkSourceSinkHandler"<<endl;
    QDP_abort(1);}

 if (fileNames.empty()){
    QDPIO::cerr << "empty file list in QuarkSourceSinkHandler"<<endl;
    QDP_abort(1);}
}


    //  First does some checks on the files in fileNames based on fileMode.
    //  If a file does not exist and fileMode==2, it is created and header
    //  info inserted.  Existing files are checked for consistent header
    //  information.  The fileMap is then created.

void QuarkSourceSinkHandler::setup_file_map()
{
 string header_info;
 header_info=outputInfo();

      // make the fileMap from the files that exist; create files
      // if they don't exist and fileMode = 2 specified

 fileMap.clear();
 for (int k=0;k<fileNames.size();k++){
 
    if (!fileExists(fileNames[k])){
                           // file does not exist    
       if (fileMode==2){
                  // create new file since fileMode == 2
          BinaryStoreDB<SerialDBKey<Key>,SerialDBData<LatticeFermion> > fcreate;
          fcreate.open(fileNames[k], O_RDWR | O_CREAT, 0664);
          fcreate.insertUserdata(header_info);
          fcreate.close();}
       else{
              // non-existent file and fileMode != 2 -> die
          QDPIO::cerr << "file "<<fileNames[k]<<" could not be opened"
                      << " and must_exist was specified"<<endl;
          QDP_abort(1);}}

    else{   // file exists

       BinaryStoreDB<SerialDBKey<Key>,SerialDBData<LatticeFermion> > fsource;
       fsource.open(fileNames[k], O_RDONLY , 0664);
       string hdinf;
       fsource.getUserdata(hdinf);
       if (hdinf!=header_info){
          QDPIO::cerr << "header info in file "<<fileNames[k]
                      << " does not match info in QuarkSourceSinkHandler"<<endl;
          QDP_abort(1);}
       vector<SerialDBKey<Key> > keys;
       fsource.keys(keys);
       for (vector<SerialDBKey<Key> >::const_iterator kt=keys.begin();kt!=keys.end();kt++){
          if (fileMap.find(kt->key())!=fileMap.end()){
              QDPIO::cerr << "duplicate keys in fileMap in QuarkSourceSinkHandler"<<endl;
              QDPIO::cerr << " ... too confusing to continue"<<endl;
              QDP_abort(1);}
          fileMap.insert(make_pair(kt->key(),k));}
       fsource.close();}
    }
}

int QuarkSourceSinkHandler::file_index_helper(int file_index)
{
 if (file_index==-1) return fileNames.size()-1;  // last file
 if ((file_index<-1)||(file_index>=fileNames.size())){
    QDPIO::cerr << "bad file index in QuarkSourceSinkHandler"<<endl;
    QDP_abort(1);}
 return file_index;
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


// ***************************************************************
  }
}
 
