// $Id: t_io.cc,v 3.1 2007-06-12 13:59:14 edwards Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace Chroma;

const std::string foo_xml="<?xml version='1.0'?><foo><coeffs1>ZOLOTAREV</coeffs1><coeffs2>TANH</coeffs2></foo>";

int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {2,2,2,2};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml_out("t_io.xml");
  push(xml_out, "t_io");

  {
    LatticeReal a;
    Double d = 17;
    random(a);

    push(xml_out,"xml_test");
    write(xml_out,"a", a);
    write(xml_out,"d", d);
    pop(xml_out);

    BinaryFileWriter tobinary("t_io.bin");
    write(tobinary, a);
    write(tobinary, d);
    tobinary.close();
  }
 
  {
    TextFileWriter totext("t_io.txt");
    totext << int(17) << "\n";
  }

  {
    float x;
    QDPIO::cout << "Read some data from t_io.txt\n";
    TextFileReader from("t_io.txt");
    from >> x;
    from.close();
    
    QDPIO::cout << "you entered :" << x << ":\n";
  }
  
  {
    multi1d<LatticeColorMatrix> u(Nd);
    for(int i=0; i < u.size(); ++i)
      gaussian(u[i]);

    XMLBufferWriter file_xml, record_xml;
    push(file_xml, "The_file_xml");
    write(file_xml, "mary", int(3));
    pop(file_xml);

    push(record_xml, "The_record_xml");
    write(record_xml, "lamb", int(7));
    pop(record_xml);

    push(xml_out, "QIO_write_test");
    write(xml_out, "file_xml", file_xml);
    write(xml_out, "record_xml", record_xml);
    write(xml_out, "u", u);
    pop(xml_out);

    writeGauge(file_xml, record_xml, u, "t_io.cfg", 
	       QDPIO_SINGLEFILE, QDPIO_SERIAL);
  }

  {
    multi1d<LatticeColorMatrix> u(Nd);
    XMLReader file_xml, record_xml;

    readGauge(file_xml, record_xml, 
	      u, "t_io.cfg", 
	      QDPIO_SERIAL);

    push(xml_out, "QIO_read_test");
    write(xml_out, "file_xml", file_xml);
    write(xml_out, "record_xml", record_xml);
    write(xml_out, "u", u);
    pop(xml_out);
  }


  CoeffType t1;
  CoeffType t2;
  istringstream foo_xml_is(foo_xml);
  XMLReader fooreader(foo_xml_is);

  try { 
    read(fooreader, "/foo/coeffs1", t1);
    read(fooreader, "/foo/coeffs2", t2);
  }
  catch(const string& e) { 
    QDPIO::cerr << "Caught error : " << e << endl;
    QDP_abort(1);
  }

  try { 
    write(xml_out, "coeffs1", t1);
    write(xml_out, "coeffs2", t2);
  }
  catch(const string& e) { 
    QDPIO::cerr << "Caught error : " << e << endl;
    QDP_abort(1);
  }

  pop(xml_out);
  
  
  // Time to bolt
  Chroma::finalize();

  exit(0);
}
