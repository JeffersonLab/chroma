#include "chroma.h"

using namespace Chroma;

int main(int argc, char *argv[]) 
{ 
  Chroma::initialize(&argc, &argv);

  CircularBuffer<int> c(3);


  for(int i=0; i < 10; i++) { 
    c.push(i);

    QDPIO::cout << "Getting at elements: " << endl;
    QDPIO::cout << "Size is " << c.size() << endl;
    QDPIO::cout << "contents: " ;
    for(int j=0; j < c.size(); j++) { 
      int ith_elem;
      c.get(j, ith_elem);
      QDPIO::cout << " " << ith_elem;
    }
    QDPIO::cout << endl;
  }


  // Nuke buffer
  c.reset();

  try { 
    int foo;
    c.get(0,foo);
  }
  catch (const CircularBuffer<int>::OutOfBoundsException& e) {
    QDPIO::cout << "Caught deliberate exception: " << e.error_string << endl;
    QDPIO::cout << "Index requested : " << e.i << endl;
    QDPIO::cout << "Buffer size :" << e.size << endl;
  }


  for(int i=11; i < 20; i++) { 
    c.push(i);

    QDPIO::cout << "Getting at elements: " << endl;
    QDPIO::cout << "Size is " << c.size() << endl;
    QDPIO::cout << "contents: " ;
    for(int j=0; j < c.size(); j++) { 
      int ith_elem;
      c.get(j, ith_elem);
      QDPIO::cout << " " << ith_elem;
    }
    QDPIO::cout << endl;
  }

  try { 
    int foo;
     c.get(7, foo);
  }
  catch (const CircularBuffer<int>::OutOfBoundsException& e) {
    QDPIO::cout << "Caught deliberate exception: " << e.error_string << endl;
    QDPIO::cout << "Index requested : " << e.i << endl;
    QDPIO::cout << "Buffer size :" << e.size << endl;
  }


  Chroma::finalize();
  exit(0);
}
