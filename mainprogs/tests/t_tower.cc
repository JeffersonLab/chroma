#include <chroma.h>
#include <tower.h>

using namespace QDP;
using namespace Chroma;

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

  Tower<float> a(5);
  a[0] = 1;
  a[1] = 1;
  a[2] = 1;
  a[3] = 1;
  a[4] = 1;

  Tower<float> b(5);
  b = a;

  Tower<float> c = a*b;

  for(int i=0; i < 5; i++) {
    QDPIO::cout << "c["<<i<<"]=" << c[i] << endl;
  }
  // Time to bolt
  Chroma::finalize();

  exit(0);
}
