#include "chroma.h"

using namespace QDP;
using namespace Chroma;

int main(int argc, char *argv[]) { 

  ChromaInitialize(&argc, &argv);

  CircularBuffer<int> c(3);


  for(int i=0; i < 10; i++) { 
    c.push(i);

    QDPIO::cout << "Getting at elements: " << endl;
    QDPIO::cout << "Size is " << c.size() << endl;
    QDPIO::cout << "contents: " ;
    for(int j=0; j < c.size(); j++) { 
      int ith_elem = c[j];
      QDPIO::cout << " " << ith_elem;
    }
    QDPIO::cout << endl;
  }


  // Nuke buffer
  c.reset();

  try { 
    int foo = c[0];
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
      int ith_elem = c[j];
      QDPIO::cout << " " << ith_elem;
    }
    QDPIO::cout << endl;
  }

  try { 
    int foo = c[7];
  }
  catch (const CircularBuffer<int>::OutOfBoundsException& e) {
    QDPIO::cout << "Caught deliberate exception: " << e.error_string << endl;
    QDPIO::cout << "Index requested : " << e.i << endl;
    QDPIO::cout << "Buffer size :" << e.size << endl;
  }


  ChromaFinalize();
}
