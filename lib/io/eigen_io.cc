#include "chromabase.h"
#include "io/eigen_io.h"
#include "io/readszinferm_w.h"
#include <iostream>
#include <iomanip>

namespace Chroma {

void read(XMLReader& xml, const string& path, RitzParams_t& header)
{


  try { 
    XMLReader paramtop(xml, path);
    read(paramtop, "Neig", header.Neig);
    read(paramtop, "RsdR", header.RsdR);
    read(paramtop, "RsdA", header.RsdA);
    read(paramtop, "RsdZero", header.RsdZero);
    read(paramtop, "ProjApsiP", header.ProjApsiP);
    read(paramtop, "Ndummy",    header.Ndummy);
    read(paramtop, "GammaFactor", header.GammaFactor);
    read(paramtop, "MaxKS", header.MaxKS);
    read(paramtop, "MaxCG", header.MaxCG);
    read(paramtop, "MinKSIter", header.MinKSIter);
    read(paramtop, "MaxKSIter", header.MaxKSIter);
    read(paramtop, "Nrenorm", header.Nrenorm);

  }
  catch( const string& error ) { 
    QDPIO::cerr << "Caught Exception " << error << endl;
    QDP_error_exit("Exiting\n");
  }
}

void write(XMLWriter& xml, const string& path, const RitzParams_t& header)
{
  push(xml, path);

  write(xml, "Neig", header.Neig);
  write(xml, "RsdR", header.RsdR);
  write(xml, "RsdA", header.RsdA);
  write(xml, "RsdZero", header.RsdZero);
  write(xml, "ProjApsiP", header.ProjApsiP);
  write(xml, "Ndummy",    header.Ndummy);
  write(xml, "GammaFactor", header.GammaFactor);
  write(xml, "MaxKS", header.MaxKS);
  write(xml, "MaxCG", header.MaxCG);
  write(xml, "MinKSIter", header.MinKSIter);
  write(xml, "MaxKSIter", header.MaxKSIter);
  write(xml, "Nrenorm", header.Nrenorm);
  pop(xml);
}

void read(XMLReader& xml, const string& path, EigenIO_t& io_header)
{

  try { 
    XMLReader paramtop(xml, path);
    read(paramtop, "eigen_file_stem", io_header.eigen_file);
    read(paramtop, "eigen_volfmt", io_header.eigen_volfmt);

    // If user specifies file format then read it
    if( paramtop.count("eigen_filefmt") == 1 ) {
      read(paramtop, "eigen_filefmt", io_header.eigen_filefmt);
    }
    else { 
      // Otherwise assume default file format -- SciDAC format
      io_header.eigen_filefmt = EVEC_TYPE_SCIDAC;
    }
  }
  catch( const string& error ) { 
    QDPIO::cerr << "Caught exception " << error << endl;
    QDP_error_exit("Exiting\n");
  }



}

void write(XMLWriter& xml, const string& path, const EigenIO_t& io_header)
{
  push(xml, path);
  write(xml,"eigen_filefmt", io_header.eigen_filefmt);
  write(xml,"eigen_file_stem", io_header.eigen_file);
  write(xml,"eigen_volfmt", io_header.eigen_volfmt);
  pop(xml);

}

void read(XMLReader &xml, const string& path, ChromaWilsonRitz_t& param)
{
  multi1d<int> seed_int(4);
  try { 
    XMLReader paramtop(xml, path);
    read(paramtop, "Param/version", param.version);

    // Read the fermion action
    XMLReader xml_tmp(paramtop, "Param/FermionAction");
    std::ostringstream os;
    xml_tmp.print(os);
    param.fermact = os.str();

    read(paramtop, "Param/nrow",     param.nrow);
    read(paramtop, "Param/rng",     param.seed);
    read(paramtop, "RitzParams", param.ritz_params);
    read(paramtop, "Cfg",       param.cfg);
    read(paramtop, "Eigen",      param.eigen_io_params);

    if( paramtop.count("StateInfo") == 1 ) { 
      XMLReader xml_state_info(paramtop, "StateInfo");
      std::ostringstream state_info_os;
      xml_state_info.print(state_info_os);
      param.state_info = state_info_os.str();
    }
    else {
      XMLBufferWriter s_i_xml;
      push(s_i_xml,"StateInfo");
      pop(s_i_xml);
      param.state_info = s_i_xml.str();
    }

  }
  catch( const string& error) { 
    QDPIO::cerr << "Caught exception " << error << endl;
    QDP_error_exit("Exiting \n");
  }
}

void write(XMLWriter &xml, const string& path, const ChromaWilsonRitz_t& param)
{
  push( xml, path );

  push( xml, "Param");
  write(xml, "version", param.version);
  write(xml, "FermionAction",    param.fermact);
  write(xml, "nrow",     param.nrow);
  write(xml, "rng",     param.seed);
  pop(xml);

  write(xml, "RitzParams", param.ritz_params);
  write(xml, "Cfg",       param.cfg);
  
  // Write out the state info struct, even if it is empty
  try { 
    // Make an istream from the XML
    std::istringstream s_i_i(param.state_info);

    // Make the XMLReader (toplevel)
    XMLReader r1(s_i_i);
    
    // Set context tag to "/StateInfo"
    XMLReader r2(r1, "/StateInfo");

    // Dump the StateInfo tag and descendets into an ostream
    ostringstream s_i_o;
    r2.print(s_i_o);

    // Dump the XML String without an extra header
    xml << s_i_o.str()  ;
  }
  catch( const string& e) { 
    QDPIO::cerr << "Caught exception : " << e << endl;
    QDP_abort(1);
  }
  write(xml, "Eigen",      param.eigen_io_params);
  
  pop(xml);
}

void writeEigen(const ChromaWilsonRitz_t& header, multi1d<Real>& lambda_lo,
		multi1d<LatticeFermion>& eigv_lo, Real& lambda_hi,
		QDP_serialparallel_t serpar)
{

  if(header.eigen_io_params.eigen_filefmt != EVEC_TYPE_SCIDAC ) { 
    QDPIO::cerr << "Writing Eigenvectors only supported in SciDAC format" << endl;
    QDP_abort(1);
  }


  XMLBufferWriter file_xml;
  int file_neig=1;

  push(file_xml, "WilsonRitzEigen");
  write(file_xml, "InputParams", header);
  write(file_xml, "lambda_hi", lambda_hi);
  write(file_xml, "Neig_file", file_neig);
  pop(file_xml);

  for(int lo=0; lo < header.ritz_params.Neig; lo++) { 
    ostringstream filename;
    filename << header.eigen_io_params.eigen_file << "_" << setw(3) << setfill('0') << lo;

    XMLBufferWriter record_eval_xml;
    push(record_eval_xml, "record_eval");
    write(record_eval_xml, "eval_num", lo);
    pop(record_eval_xml);

    XMLBufferWriter record_evec_xml;
    push(record_evec_xml, "record_evec");
    write(record_evec_xml, "evec_num", lo);
    pop(record_evec_xml);

    QDPIO::cout << "Opening eigenvector file: "<< filename.str() << endl;
    QDPFileWriter outfile(file_xml, filename.str(), 
			  header.eigen_io_params.eigen_volfmt, 
			  serpar, QDPIO_OPEN);

    multi1d<Real> lambda_lo_arr(1);
    lambda_lo_arr[0] = lambda_lo[lo];
    write(outfile, record_eval_xml, lambda_lo_arr);
    write(outfile, record_evec_xml, eigv_lo[lo]);
    close(outfile);
  }
}

void readEigenPair(Real& lambda_lo, int& eig_index,
		   LatticeFermion& eigv, 
		   const string& filename,
		   QDP_serialparallel_t serpar,
		   XMLReader& file_xml)
  /* XMLReader& record_xml) */
{
  
  XMLReader record_eval_xml;
  XMLReader record_evec_xml;

  QDPIO::cout << "Attempting to read from " << filename << endl;
  QDPFileReader from(file_xml, filename, serpar);
  

  // Now read the first eigenvalue/vector 
  multi1d<Real> lambda_lo_aux(1);
  read(from, record_eval_xml, lambda_lo_aux);
  lambda_lo = lambda_lo_aux[0];
  read(from, record_evec_xml, eigv);  
  /*
  // Get the lambda_lo
  try {
    read( record_xml, "/record_evalue/lambda", lambda_lo);
  }
  catch( const string& e) { 
    QDPIO::cerr << "Caught exception reading e-value " << e << endl;
    QDP_error_exit("Exiting\n");
  }
  */
  // Get the eval_num
  int eval_num; 
  int evec_num;
  try {
    read( record_eval_xml, "/record_eval/eval_num", eval_num);
  }
  catch( const string& e) { 
    QDPIO::cerr << "Caught exception reading eval_num " << e << endl;
    QDP_error_exit("Exiting\n");
  }

  try {
    read( record_evec_xml, "/record_evec/evec_num", evec_num);
  }
  catch( const string& e) { 
    QDPIO::cerr << "Caught exception reading evec_num " << e << endl;
    QDP_error_exit("Exiting\n");
  }

  if( eval_num != evec_num) { 
    QDP_error_exit("Eval num and evec_num are different\n");
  }

  eig_index = eval_num;
  // Done with file
  close(from);
}
    
	       
void readEigen(ChromaWilsonRitz_t& header, multi1d<Real>& lambda_lo,
	       multi1d<LatticeFermion>& eigv_lo, Real& lambda_hi,
	       const string& filename_stem, 
	       int Neig,
	       QDP_serialparallel_t serpar)	
{

  int neig_to_load = Neig;

  for(int lo=0; lo < Neig && lo < neig_to_load; lo++) { 
    XMLReader file_xml;
    
    int eig_index;
    
    Real lambda_lo_tmp;
    LatticeFermion eigv_tmp;
    
    ostringstream filename;
    
    filename << filename_stem << "_" << setw(3) << setfill('0') << lo;
    readEigenPair(lambda_lo_tmp,  eig_index, eigv_tmp,
		  filename.str(), serpar, 
		  file_xml);
    
    
    // Stuff we need to do once
    if( lo == 0 ) { 
      // Read header -- should be the same for all but I wont check for that
      try { 
	read(file_xml, "/WilsonRitzEigen/InputParams", header);
      }
      catch ( const string& e ) { 
	QDPIO::cerr << "Caught Exception reading header: " << e << endl;
	QDP_error_exit("Exiting\n");
      }
      
      
      // Get the highest e-value
      try { 
	read(file_xml, "/WilsonRitzEigen/lambda_hi", lambda_hi);
      }
      catch ( const string& e ) { 
	QDPIO::cerr << "Caught Exception reading lambda_hi: " << e << endl;
	QDP_error_exit("Exiting\n");
      }
      
      // Check how many e-values there are
      if( header.ritz_params.Neig < Neig ) { 
	QDPIO::cout << "Requested " << Neig << " eigenpairs but only " << header.ritz_params.Neig << " were computed. Will only read " << header.ritz_params.Neig << " pairs" << endl;
	neig_to_load = header.ritz_params.Neig;
      }
      else { 
	header.ritz_params.Neig = Neig;
      }
      
      // Resize arrays to accomodate
      lambda_lo.resize(neig_to_load);
      eigv_lo.resize(neig_to_load);
    }
    
    // Check index
    if( eig_index != lo) { 
      QDPIO::cerr << "Error: index and eig_index dont match lo = " << lo << " eig_index = " << eig_index << endl;
      QDP_error_exit("Exiting\n");
    }
    
    // Put things in their place
    lambda_lo[eig_index] = lambda_lo_tmp;
    eigv_lo[eig_index] = eigv_tmp;
  }
}


// Read SZIN eigenvalues. 
// Expects the following:
//  filename_stem.xml -- contains the Eigenvalues under tag Eigenvalues
//  filename_stem_XXX contains the SZIN eigenvectors
void readEigenSzin(multi1d<Real>& lambda_lo,
		   multi1d<LatticeFermion>& eigv_lo, Real& lambda_hi,
		   const int Neig, 
		   const string& filename_stem)
{

  ostringstream xml_filename;
  xml_filename << filename_stem << ".xml" ;

  QDPIO::cout << "Attempting to open " << xml_filename.str() << endl;
  // Try and get the e-values.
  XMLReader reader(xml_filename.str());

  try { 
    XMLReader paramtop(reader, "/Eigenvalues");
    read(paramtop, "lambda_lo", lambda_lo);
    read(paramtop, "lambda_hi", lambda_hi);
  }
  catch( const string& e) { 
    QDPIO::cerr << "Caught exception : " << e << endl;
    QDP_abort(1);
  }

  if ( lambda_lo.size() != Neig ) { 
    QDPIO::cerr << "Mismatch in no of low eigenvalues. Neig = " << Neig << 
      "but read " << lambda_lo.size() << endl;
    QDP_abort(1);
  }

  eigv_lo.resize(lambda_lo.size());

  for(int evec = 0; evec < lambda_lo.size(); evec++) { 

    // Create filename
    ostringstream filename;
    filename << filename_stem << "_" <<  setw(3) << setfill('0') << evec;
    QDPIO::cout << "Attempting to open " << filename.str() << endl;
    // read the evec
    try { 
      
      readSzinFerm(eigv_lo[evec], filename.str());
    }
    catch (const string& e) { 
      QDPIO::cerr << "Caught exception " << e << endl;
      QDP_abort(1);
    }


  }

  
}

}  // end namespace Chroma
  
