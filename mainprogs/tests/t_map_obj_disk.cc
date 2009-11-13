#include "chroma.h"

#include "util/ferm/map_obj.h"
#include "util/ferm/map_obj/map_obj_disk.h"

using namespace QDP;
using namespace Chroma;


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
  MapObjectDisk<int,float> the_map(p);

  // Open the map for 'filling'
  the_map.openWrite();
  
  // store a quadratic
  for(int i=0; i < 10; i++) { 
    float val = (float)i;
    val *= val;

    the_map.insert(i, val);
  }
  the_map.closeWrite();



  Chroma::finalize();
}
