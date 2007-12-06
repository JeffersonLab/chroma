/***************************************************************************
 * make_test_props.cc
 * Driver routine to write some simple test sources/sinks to lime files 
 * *************************************************************************/



#include "chroma.h"
#include "meas/inline/make_xml_file.h"

#include "stdio.h"
#include "string.h"

using namespace QDP;
using namespace Chroma;


struct Param_t
{
	
	multi1d<int> layout; //Lattice dimensions 
	int decay_dir; //time direction

};


struct Files_t
{
  multi1d<std::string>   coeff_files; //Files of all the ops to make 
	multi1d<std::string>   elem_op_files;  //Files containing elemental ops. 
																			 //Each entry contains root of where
																			 //All elemental ops for a given config. 
																			 //may be found. 
																			 //elem_op_files.size() = Nbins
	
	multi1d<std::string>  op_files; 			 //Roots of output files
																			 //one entry for each bin


};

//! Mega-structure of all input
struct MakeOpsInput_t
{
  Param_t            param;
  Files_t             files;
};

// Reader for input parameters
void read(XMLReader& xml, const string& path, Param_t& param)
{
  XMLReader paramtop(xml, path);

  int version;
  read(paramtop, "version", version);

  switch (version) 
  {
  case 1:
    
		read(paramtop, "Layout", param.layout );
		read(paramtop, "Decay_dir", param.decay_dir);
	break;

  default :
    /**************************************************************************/

    cerr << "Input parameter version " << version << " unsupported." << endl;
    exit(1);
  }
}

//! Propagator parameters
void read(XMLReader& xml, const string& path, Files_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "Coeff_files", input.coeff_files);
  read(inputtop, "ElemOp_files", input.elem_op_files);
  read(inputtop, "OpOutput_files", input.op_files);
	
}

// Reader for input parameters
void read(XMLReader& xml, const string& path, MakeOpsInput_t& input)
{
  XMLReader inputtop(xml, path);

  // Read all the input groups
  try
  {
    // Read program parameters
    read(inputtop, "Param", input.param);

    // Read in the propagator(s) info
    read(inputtop, "Files", input.files);
  }
  catch (const string& e) 
  {
    cerr << "Error reading make_ops data: " << e << endl;
    exit(1);
  }
}


//! Baryon operator
    struct BaryonOperator_t
    {
      //! Quark orderings within a baryon operator
      struct Orderings_t
      {
	//! Baryon operator time slices
	struct TimeSlices_t
	{
	  //! Baryon operator dilutions
	  struct Dilutions_t
	  {
	    //! Momentum projected correlator
	    struct Mom_t
	    {
	      multi1d<int>       mom;       /*!< D-1 momentum of this operator */
	      multi1d<DComplex>  op;        /*!< Momentum projected operator */
	    };

	    multi1d<Mom_t> mom_projs;       /*!< Holds momentum projections of the operator */
	  };

	  int                  t0;          /*!< Source time location */
	  multi3d<Dilutions_t> dilutions;   /*!< Hybrid list indices */
	};
multi1d<int> perm;                  /*!< This particular permutation of quark orderings */
	multi1d<TimeSlices_t> time_slices;  /*!< Time slices of the lattice that are used */
      };

      multi1d< multi1d<int> > perms;   /*!< Permutations of quark enumeration */

			std::string   smearing;          /*!< String holding quark smearing xml */

      Seed          seed_l;            /*!< Id of left quark */
      Seed          seed_m;            /*!< Id of middle quark */
      Seed          seed_r;            /*!< Id of right quark */

      //int           operator_num;      /*!< Operator number within file */
      std::string   id;                /*!< Tag/ID used in analysis codes */

      int           mom2_max;          /*!< |\vec{p}|^2 */
      int           decay_dir;         /*!< Direction of decay */
      multi1d<Orderings_t> orderings;  /*!< Array is over quark orderings */
    };



struct GroupBaryonOperator_t
{
	struct Term_t
	{
		std::string el_op_name;
		DComplex coeff;
	};

	multi1d<Term_t> term;
	
	std::string name; 
	
};


//! BaryonOperator header reader
void read(XMLReader& xml, const string& path, BaryonOperator_t& param)
{
  XMLReader paramtop(xml, path);

  try
  {
    int version;
  	  read(paramtop, "version", version);
      read(paramtop, "id", param.id);
      //read(paramtop, "operator_num", param.operator_num);
      read(paramtop, "mom2_max", param.mom2_max);
      read(paramtop, "decay_dir", param.decay_dir);
      read(paramtop, "seed_l", param.seed_l);
      read(paramtop, "seed_m", param.seed_m);
      read(paramtop, "seed_r", param.seed_r);
      read(paramtop, "perms", param.perms);
  
			//read the smearing??
	}
  catch (const std::string& e) 
  {
    cerr << "BaryonOperator: Error reading: " << e << endl;
    exit(1);
  }
}

//! BaryonOperator header writer 
void write(XMLWriter& xml, const string& path, BaryonOperator_t& param)
{
  push(xml, path);

  try
  {
    
  	 // write(xml, "version", version);
      write(xml, "id", param.id);
      //read(paramtop, "operator_num", param.operator_num);
      write(xml, "mom2_max", param.mom2_max);
      write(xml, "decay_dir", param.decay_dir);
      write(xml, "seed_l", param.seed_l);
      write(xml, "seed_m", param.seed_m);
      write(xml, "seed_r", param.seed_r);
      write(xml, "perms", param.perms);
  
			//read the smearing??
	}
  catch (const std::string& e) 
  {
    cerr << "BaryonOperator: Error writing: " << e << endl;
    exit(1);
  }

	pop(xml);

}


//GroupBaryonOperator Writer
void write(XMLWriter &xml, const std::string &path, const GroupBaryonOperator_t::Term_t &param)
{
	push(xml, path);

	write(xml, "ElementalName" , param.el_op_name);
	write(xml, "Coefficient" , param.coeff);

	pop(xml);

}	

//GroupBaryonOperator Writer
void write(XMLWriter &xml, const std::string &path, const GroupBaryonOperator_t &param)
{
	push(xml, path);

	write(xml, "Name" , param.name);
	write(xml, "Terms" , param.term);

	pop(xml);

}	
//! BaryonOperator binary reader
void read(BinaryReader& bin, BaryonOperator_t::Orderings_t::TimeSlices_t::Dilutions_t::Mom_t& param)
    {
			read(bin, param.mom);
			read(bin, param.op);
		}

//! BaryonOperator binary writer
void read(BinaryReader& bin, BaryonOperator_t::Orderings_t::TimeSlices_t::Dilutions_t& param)
{
	read(bin, param.mom_projs);
}

//! BaryonOperator binary reader 
void read(BinaryReader& bin, BaryonOperator_t::Orderings_t::TimeSlices_t& param)
{
	read(bin, param.dilutions);
}


//! BaryonOperator binary reader
void read(BinaryReader& bin, BaryonOperator_t::Orderings_t& param)
{
	
	read(bin, param.perm);
	read(bin, param.time_slices);
}


//read BaryonOperator_t struct from binary bufffer reader.
void read(BinaryReader &bin, BaryonOperator_t &param)
{

  read(bin, param.seed_l);
	read(bin, param.seed_m);
  read(bin, param.seed_r);
	
	read(bin, param.mom2_max);
  read(bin, param.decay_dir);
	

	read(bin, param.perms);
 	
	read(bin, param.orderings);

}

//! BaryonOperator binary writer
void write(BinaryWriter& bin, const BaryonOperator_t::Orderings_t::TimeSlices_t::Dilutions_t::Mom_t& param)
    {
      write(bin, param.mom);
      write(bin, param.op);
    }

//! BaryonOperator binary writer
void write(BinaryWriter& bin, const BaryonOperator_t::Orderings_t::TimeSlices_t::Dilutions_t& param)
{
	write(bin, param.mom_projs);
}

//! BaryonOperator binary writer 
void write(BinaryWriter& bin, const BaryonOperator_t::Orderings_t::TimeSlices_t& param)
{
	write(bin, param.dilutions);
}


//! BaryonOperator binary writer
void write(BinaryWriter& bin, const BaryonOperator_t::Orderings_t& param)
{
	write(bin, param.perm);
  write(bin, param.time_slices);
}


//write BaryonOperator_t struct to binary bufffer writer.
void write(BinaryWriter &bin, const BaryonOperator_t &param)
{

  write(bin, param.seed_l);
	write(bin, param.seed_m);
  write(bin, param.seed_r);
  write(bin, param.mom2_max);
  write(bin, param.decay_dir);
 	write(bin, param.perms);
 	write(bin, param.orderings);
}


//This routine goes through (possibly several) coeff files and fills 
//the operator array
void readCoeffFiles(multi1d<GroupBaryonOperator_t> &arr, const multi1d<std::string> coeff_files) 
{

	int nops = 0;

	//First determine how many ops total 
	for (int i = 0 ; i < coeff_files.size() ; ++i)
	{

		TextFileReader reader(coeff_files[i]);

		int op;

		reader >> op;

		reader.close();

		nops += op; 
		
	//	QDPIO::cout<< "Nops = "<<nops<<endl;

	}

	arr.resize(nops);

	//Now read the coeffs

	for (int i = 0 ; i < coeff_files.size() ; ++i)
	{

		int op;

		TextFileReader reader(coeff_files[i]);
		reader >> op;
		
		for (int l = 0 ; l < op ; ++l)
		{
			int nelem;
			std::string name; 
			reader >> nelem >> name;

			arr[ i + l ].name = name;

			arr[ i + l ].term.resize( nelem );

			for (int m = 0 ; m < nelem ; ++m)
			{
				int a,b,c,i,j,k;
				Real re,im;
				char lparen,comma,rparen;

				reader >> a >> b >> c >> i >> j >> k >> lparen >> re >> comma >> im >> rparen;

				arr[ i + l].term[m].coeff = cmplx( re, im );

				//form the name of the elem. op file
				std:stringstream strm;
				
				strm << "Nuc_" << a << b << c << '_'; 		
				
				//displacements may be negative 
				if (i < 0)
				{
					strm<< 'm';
				}
				strm << abs(i);
	
				if (j < 0)
				{
					strm<< 'm';
				}
				strm << abs(j);

				if (k < 0)
				{
					strm<< 'm';
				}
				strm << abs(k);



				arr[i +l ].term[m].el_op_name = strm.str();

		} //m
			
	} //l 

	 		
	reader.close();

	} //i

}//void


//Fill the operator info from first elem. op so it isn't done for every 
//elem op
void initOp (BaryonOperator_t &oper , const BaryonOperator_t &elem_oper ) 
{
	oper.mom2_max= elem_oper.mom2_max;
	oper.decay_dir = elem_oper.decay_dir;
	oper.seed_l = elem_oper.seed_l;
	oper.seed_m = elem_oper.seed_m;
	oper.seed_r = elem_oper.seed_r;
	oper.perms = elem_oper.perms;	

	int Nord = elem_oper.perms.size();
	int Nt = elem_oper.orderings[0].time_slices.size();

	oper.orderings.resize(Nord);

	for (int p = 0 ; p < Nord ; ++p)
	{
		oper.orderings[p].perm = elem_oper.orderings[p].perm;

		oper.orderings[p].time_slices.resize(Nt);

		for (int t_0 = 0 ; t_0 < Nt ; ++t_0)
		{
			oper.orderings[p].time_slices[t_0].t0 = elem_oper.orderings[p].time_slices[t_0].t0;
				
			//Dilution sizes for each quark
			int Ni = elem_oper.orderings[p].time_slices[0].dilutions.size1();
			int Nj = elem_oper.orderings[p].time_slices[0].dilutions.size2();
			int Nk = elem_oper.orderings[p].time_slices[0].dilutions.size3();

			oper.orderings[p].time_slices[t_0].dilutions.resize(Ni, Nj, Nk);

			for(int i = 0 ; i < Ni ; ++i)   	
				for(int j = 0 ; j < Nj ; ++j)   	
					for(int k = 0 ; k < Nk ; ++k)   	
					{
						int Nmom =  elem_oper.orderings[p].time_slices[t_0].dilutions(i,j,k).mom_projs.size();

						//QDPIO::cout<<"Nmom = " <<Nmom<<endl; 
						oper.orderings[p].time_slices[t_0].dilutions(i,j,k).mom_projs.resize(Nmom);

						for (int m = 0 ; m < Nmom ; ++m)
						{
							oper.orderings[p].time_slices[t_0].dilutions(i,j,k).mom_projs[m].mom = 
							elem_oper.orderings[p].time_slices[t_0].dilutions(i,j,k).mom_projs[m].mom;

							//zero operator
							int Lop = elem_oper.orderings[p].time_slices[t_0].dilutions(i,j,k).mom_projs[m].op.size();

						//	QDPIO::cout<<"Lop = "<<Lop<<endl;
							oper.orderings[p].time_slices[t_0].dilutions(i,j,k).mom_projs[m].op.resize(Lop);

							for (int t = 0 ; t < Lop ; ++t)
							{
								//DComplex zro = cmplx(Real(0), Real(0));
								oper.orderings[p].time_slices[t_0].dilutions(i,j,k).mom_projs[m].op[t] = zero;

							}

						} //m
					}//ijk
		}//t_0
	}//p

} //void



//Add the elemental op to the final operator
void addTo(BaryonOperator_t &oper, const BaryonOperator_t &elem_oper, 
		const DComplex& coeff)
{

//	oper.perms = elem_oper.perms;

	int Nord = elem_oper.perms.size();
	int Nt = elem_oper.orderings[0].time_slices.size();

	for (int p = 0 ; p < Nord ; ++p)
	{
	//	oper.orderings[p].perm = elem_oper.orderings[p].perm;

		for (int t_0 = 0 ; t_0 < Nt ; ++t_0)
		{
		//	oper.orderings[p].time_slices[t_0].t0 = elem_oper.orderings[p].time_slices[t_0].t0;

			for(int i = 0 ; i < elem_oper.orderings[p].time_slices[0].dilutions.size1() ; ++i)   	
				for(int j = 0 ; j < elem_oper.orderings[p].time_slices[0].dilutions.size2() ; ++j)   	
					for(int k = 0 ; k < elem_oper.orderings[p].time_slices[0].dilutions.size3() ; ++k)   	
					{
						for (int m = 0 ; m < elem_oper.orderings[p].time_slices[t_0].dilutions(i,j,k).mom_projs.size() ; ++m)
						{
						//	oper.orderings[p].time_slices[t_0].dilutions(i,j,k).mom_projs[m].mom = 
							//	elem_oper.orderings[p].time_slices[t_0].dilutions(i,j,k).mom_projs[m].mom;

							oper.orderings[p].time_slices[t_0].dilutions(i,j,k).mom_projs[m].op += 
								elem_oper.orderings[p].time_slices[t_0].dilutions(i,j,k).mom_projs[m].op * coeff;
						} //m
					}//ijk
		}//t_0
	}//p

} //void


//fill a time slice of an op
void fillElem( BaryonOperator_t &fullOp , const BaryonOperator_t &timeslice_op, const int timeSlice)
{

	int Nord = timeslice_op.orderings.size();

	for (int p = 0 ; p < Nord ; ++p)
	{
	//	QDPIO::cout<<"Lt = "<<fullOp.orderings[p].time_slices.size()<<endl;

		fullOp.orderings[p].time_slices[timeSlice] = 
			timeslice_op.orderings[p].time_slices[0];
	}

}


int main(int argc, char **argv)
{

	

	// Put the machine into a known state
  Chroma::initialize(&argc, &argv);
  //  linkageHack();

  // Put this in to enable profiling etc.
  START_CODE();

//Read Input params from xml
	MakeOpsInput_t input;

	XMLReader xml_in;
	
	StopWatch swatch, snoop;




	try
  {
    xml_in.open(Chroma::getXMLInputFileName());
    read(xml_in, "/MakeOps", input);
  }
  catch(const std::string& e) 
  {
    QDPIO::cerr << "MAKEOPS: Caught Exception reading XML: " << e << endl;
    QDP_abort(1);
  }
  catch(std::exception& e) 
  {
    QDPIO::cerr << "MAKEOPS: Caught standard library exception: " << e.what() << endl;
    QDP_abort(1);
  }
  catch(...)
  {
    QDPIO::cerr << "MAKEOPS: caught generic exception reading XML" << endl;
    QDP_abort(1);
  }

  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "MakeOps");

  // Write out the input
  write(xml_out, "Input", xml_in);
		

  // Boiler plate magic to lay out our lattice
  // Happens in every QDP++ program more or less this
  // way.

  // Define the Lattice Size in the layout
  // Note: this expects a type of multi1d<int>
  Layout::setLattSize(input.param.layout);
  Layout::create();   // Setup the layout

	
	multi1d<GroupBaryonOperator_t> op;

	QDPIO::cout<< "Reading Coeff Files" << endl;

	//Read coeff files 
	readCoeffFiles(op, input.files.coeff_files);
	

	int time_dir = input.param.decay_dir;
	int Nt = input.param.layout[time_dir]; 
	int Nops = op.size();
	
	int Nbins = input.files.elem_op_files.size();

	//Write oplist to output xml
	push(xml_out, "Ops");
	for (int l = 0 ; l < Nops ; ++l)
	{
		write(xml_out, "Name", op[l].name );
	}

	pop(xml_out); //Ops

	QDPIO::cout << " MAKE_OPS: construct baryon operators" << endl;


	//loop over configurations
	for (int i = 0 ; i < Nbins ; ++i)
	{
		QDPIO::cout << "Forming Ops: Bin "<< i << endl; 

		for (int l = 0 ; l < Nops ; ++l)
		{
		 	
			snoop.reset();
			snoop.start();
				
			BaryonOperator_t source, sink;

			source.id = op[l].name; 
			sink.id = op[l].name;
			
			QDPIO::cout<< "Making Op " << source.id << endl;

			int Nterms = op[l].term.size();

			bool init = false;

			for (int m = 0 ; m < Nterms ; ++m)
			{
				
				BaryonOperator_t elem_source, elem_sink;

				bool initElem = false;
				
				//Loop over timeslice files for a single elemental
				for (int t_0 = 0 ; t_0 < Nt ; ++t_0)
				{
					BaryonOperator_t elem_source_t, elem_sink_t;
					
					std::stringstream stm;

					stm << input.files.elem_op_files[i] << op[l].term[m].el_op_name 
						<< "_" << t_0 << ".lime";

					std::string filename = stm.str();

					XMLReader file_xml, src_record_xml, sink_record_xml;
					BinaryBufferReader src_recordBin, sink_recordBin; 

					QDPFileReader elem_file(file_xml, filename, QDPIO_SERIAL);

					QDPIO::cout<<"Reading Source"<<endl;
					
					swatch.reset();
					swatch.start();

					//Read source 
					read(elem_file, src_record_xml, src_recordBin);
					read(src_recordBin, elem_source_t);

					swatch.stop();
					QDPIO::cout<<"Read Source: "<< swatch.getTimeInSeconds() << " secs"<< endl;


					QDPIO::cout<<"Reading Sink" <<endl;

					swatch.reset();
					swatch.start();


					//Read sink 
					read(elem_file, sink_record_xml, sink_recordBin);
					read(sink_recordBin, elem_sink_t);

					swatch.stop();
					QDPIO::cout<<"Read Sink: "<< swatch.getTimeInSeconds() << " secs"<< endl;

					elem_file.close();
			
					swatch.reset();
					swatch.start();
					
					if (!initElem)
						{
							int Nord = elem_sink_t.orderings.size();
							QDPIO::cout << "Nord = "<<Nord<<endl;


							//Initialize elemental ops
							elem_sink.orderings.resize(Nord);
							elem_source.orderings.resize(Nord);

							for (int p = 0 ; p < Nord ; ++p)
							{
								elem_source.orderings[p].perm = elem_source_t.orderings[p].perm;
								elem_source.orderings[p].time_slices.resize(Nt);
			
								elem_sink.orderings[p].perm = elem_sink_t.orderings[p].perm;
								elem_sink.orderings[p].time_slices.resize(Nt);

							}
							
							//QDPIO::cout<<"reading xml's"<<endl;
							//Snarf operator info that is same for source and sink 
							read(src_record_xml, "/BaryonCreationOperator", elem_source);   
							read(sink_record_xml, "/BaryonAnnihilationOperator", elem_sink);   
	
							initElem = true; 

						} //if !initElem
					
						//QDPIO::cout<<"Filling timeslices"<<endl;
						//fill a timeslice of the elem. op 
						fillElem(elem_sink, elem_sink_t, t_0);
						fillElem(elem_source, elem_source_t, t_0);

				} //t_0
				
				swatch.stop();
				QDPIO::cout<<	"Elem Op filled: " << swatch.getTimeInSeconds() << 
				" secs" << endl;
			/*
				QDPIO::cout<< "Elem Op "<<  op[l].term[m].el_op_name << endl ;
				QDPIO::cout << "Source: (t_0 = 3 , p = 2, dil= 0,0,0, t = 0  ) ="
					<<elem_source.orderings[2].time_slices[3].dilutions(0,0,0).mom_projs[0].op[0] << endl;


				QDPIO::cout << "Sink: (t_0 = 3 , p = 2, dil= 0,0,0, t = 0  ) ="
					<<elem_sink.orderings[2].time_slices[3].dilutions(0,0,0).mom_projs[0].op[0] << endl;
*/			
				
				swatch.reset();
				swatch.start();

				if (!init)
				{

					QDPIO::cout<<"init. group baryon op source"<<endl;				
					//Set all the headers, perms, mom_projs, t_0's only once 
					initOp(source , elem_source);
					QDPIO::cout<<"init. group baryon op sink"<<endl;				
					initOp(sink , elem_sink);
			
					QDPIO::cout<<"Finished init ops"<<endl;
					init = true;
				}
				QDPIO::cout<<"Adding elemental to group baryon op"<<endl;
				//exit(0);

				//add the elem op to the final op
				addTo(source, elem_source, op[l].term[m].coeff);
				addTo(sink, elem_sink, op[l].term[m].coeff);
				
			} //m
	
			swatch.stop();
			QDPIO::cout<<"Source and sink ops constructed: "<< swatch.getTimeInSeconds() << " secs " 
				<< endl;




			
			QDPIO::cout<< "Final Op "<<  op[l].name << endl ;
			
			int Nord = source.orderings.size();

			for (int t_0 = 0 ; t_0 < Nt ; ++t_0)
				for (int p = 0 ; p < Nord ; ++p)
				{
			QDPIO::cout << "Source: (t_0 = "<< t_0 << " , p = " << p <<
				 " , dil= 0,0,0, t = 0  ) =" << source.orderings[p].time_slices[t_0].dilutions(0,0,0).mom_projs[0].op[0] << endl;

	QDPIO::cout << "Sink: (t_0 = "<< t_0 << " , p = " << p <<
				 " , dil= 0,0,0, t = 0  ) =" << sink.orderings[p].time_slices[t_0].dilutions(0,0,0).mom_projs[0].op[0] << endl;
				}


			swatch.reset();
			swatch.start();
			
			//Write the source and sink op to file; 
			XMLBufferWriter file_xml; 
			//Put some stuff in the file xml; 
			write(file_xml, "GroupBaryonOperator", op[l]);

			BinaryBufferWriter source_bin, sink_bin; 
			
			std::string outfilename = input.files.op_files[i] + op[l].name;

			QDPIO::cout<<"Filename = "<<outfilename<<endl;

			QDPFileWriter outfile(file_xml, outfilename, QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);

			//write Source op
			XMLBufferWriter source_xml;
			//put some stuff in the source xml
			write(source_xml, "CreationOperator" , source);

			write(source_bin, source);
			write(outfile, source_xml, source_bin);


			//write Sink op
			XMLBufferWriter sink_xml;
			//put some stuff in sink xml
			write(sink_xml, "AnnihilationOperator", sink);

			write(sink_bin, sink);
			write(outfile,sink_xml, sink_bin);

			outfile.close();

			swatch.stop();
			QDPIO::cout<< "source and sink written: " << swatch.getTimeInSeconds() << " secs" << endl;
		
			snoop.stop();
			QDPIO::cout << "Op " << l << " constructed: " << snoop.getTimeInSeconds()
				<< " secs " << endl;

		} //l

		
		
	}//i		

	pop(xml_out); //MakeOps 

 // Clean up QDP++
  QDP_finalize();
  exit(0);  // Normal exit
}


