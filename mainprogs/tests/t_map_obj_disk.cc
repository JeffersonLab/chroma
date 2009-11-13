#include "chroma.h"
#include <cstdlib>

#include "util/ferm/map_obj.h"
#include "util/ferm/map_obj/map_obj_disk.h"

using namespace QDP;
using namespace Chroma;


void fail()
{
  QDPIO::cout << "FAIL" << endl;
  Chroma::finalize();
  exit(1);
}


int main(int argc, char *argv[])
{

  Chroma::initialize(&argc, &argv);

  // Setup the QDP++
  const int foo[] = {2,2,2,2};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  // Params to create a map object disk
  MapObjectDiskParams p("./mapObjTest");

  // Make a disk map object -- keys are ints, data floats
  MapObjectDisk<char,float> the_map(p);

  // Open the map for 'filling'
  the_map.openWrite();
  QDPIO::cout << "Inserting (key,value): " << endl;
  // store a quadratic
  for(char i=0; i < 10; i++) { 
    float val = (float)i;
    val *= val;
    char key = 'a'+i; 
    QDPIO::cout << "  ( " << key <<", "<< val <<" )" << endl;
    the_map.insert(key, val);
  }
  
  the_map.closeWrite();

  /* Now reopen - random access */
  the_map.openRead();

  QDPIO::cout << "Forward traversal test: " << endl;
  // Traverse in sequence

  QDPIO::cout << "Looking up key: " ;
  for(char i=0; i<10; i++) {
    char key = 'a'+i;
    float val;
  
    QDPIO::cout << " " << key;

    the_map.lookup(key, val);
    float refval = i*i;
    float diff = fabs(refval-val);
    if ( diff > 1.0e-8 ) {
      fail();
    }
  }

  QDPIO::cout << endl << "OK" << endl;


  QDPIO::cout << "Reverse traversal check" << endl;
  QDPIO::cout << "Looking up key: " ;
  // Traverse in reverse

  for(char i=9; i >= 0; i--) {
    char key='a'+i; 
    float val;
    QDPIO::cout << " " << key;
    the_map.lookup(key, val);
    float refval = i*i;
    float diff = fabs(refval-val);
    if ( diff > 1.0e-8 ) {
      fail();
    }
  }
  QDPIO::cout << endl << "OK" << endl;


  // 20 'random' reads
  QDPIO::cout << "Random access... 20 reads" << endl;
  QDPIO::cout << "Looking up key: " ;
  for(int j=0; j < 20; j++){ 
    char key = 'a'+std::rand() % 10; // Pick randomly in 0..10
    float val;
    QDPIO::cout << " " << key ;
    the_map.lookup(key, val);

    float refval = (key-'a')*(key-'a');

    float diff = fabs(refval-val);
    if ( diff > 1.0e-8 ) {
      fail();
    }
  }
  QDPIO::cout << endl << "OK" << endl;

  the_map.closeRead();

  QDPIO::cout << "OK" << endl;

  Chroma::finalize();
  return 0;
}
