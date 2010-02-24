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
  a[0] = 1; a[1]=1; a[2]=1; a[3]=1; a[4]=1;

  Tower<float> b(5);
  b = a;

  Tower<float> c = a*b;

  for(int i=0; i < 5; i++) {
    QDPIO::cout << "c["<<i<<"]=" << c[i] << endl;
  }

  Tower<float> d(3);
  d[0] = 1; d[1]=2; d[2]=3;

  Tower<float> e(3);
  e[0] = 4; e[1]=5; e[2]=6;

  Tower<float> f=d*e;

  QDPIO::cout << endl;
  for(int i=0; i < 3; i++) {
    QDPIO::cout << "f["<<i<<"]=" << f[i] << endl;
  }

  Tower<ComplexD> c1(3);
  Tower<ComplexD> c2(3);
  Tower<ComplexD> c3(3);

  random(c1);
  random(c2);
  random(c3);

  Tower<ComplexD> cf=(c1*c2)*c3;
  Tower<ComplexD> cf2=c1*(c2*c3);
  Tower<ComplexD> diff=cf-cf2;

  QDPIO::cout << "Diff of towers" << endl;
  print(diff);

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
