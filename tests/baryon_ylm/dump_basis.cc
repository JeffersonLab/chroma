#include "meas/inline/hadron/baryon_derivative_ylm.h"
#include "util/ferm/clebsch.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using namespace Chroma::BaryonDerivativeYlm;
void require(bool ok) { if (!ok) throw std::runtime_error("Baryon basis check failed"); }
int main() {
  // Half-integer tests catch accidental physical-spin arguments in this API.
  require(std::abs(Hadron::clebsch(1,1,1,-1,0,0)-1/std::sqrt(2.0))<1e-13);
  require(std::abs(Hadron::clebsch(1,-1,1,1,0,0)+1/std::sqrt(2.0))<1e-13);
  require(std::abs(Hadron::clebsch(1,1,1,1,2,2)-1)<1e-13);
  std::vector<Key> keys{{0}};
  for (int m=-1;m<=1;++m) keys.push_back({3,m});
  for (int f : {33,23}) for (int L=0;L<=2;++L) for (int m=-L;m<=L;++m) keys.push_back({f,L,m});
  require(expand(keys).size()==22);
  keys.push_back({0}); require(expand(keys).size()==22);
  for (const Key& k : std::vector<Key>{{},{3},{1,0},{33,3,0},{23,1,2},{0,0}}) {
    bool rejected=false; try {component(k);} catch (const std::invalid_argument&) {rejected=true;}
    require(rejected);
  }
  std::cout << std::setprecision(17);
  for (const auto& c : expand(keys)) {
    std::cout << "K"; for (int x:c.key) std::cout << " " << x; std::cout << "\n";
    for (const auto& term:c.paths) {
      for (int q=0;q<3;++q) {
        if (q) std::cout << "|";
        for (int x:term.first[q]) std::cout << x << " ";
      }
      std::cout << ": " << term.second.real() << " " << term.second.imag() << "\n";
    }
  }
}
