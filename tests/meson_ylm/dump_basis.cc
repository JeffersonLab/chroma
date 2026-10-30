#include "meas/inline/hadron/meson_derivative_ylm.h"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <stdexcept>
int main() {
  using namespace Chroma::MesonDerivativeYlm;
  auto all=expand({{0},{1},{2},{3}});
  assert(all.size()==40);
  assert(expand({{0},{1},{2},{3}},true).size()==26);
  assert(expand({{1},{1,0}}).size()==3);
  for(Key bad:std::vector<Key>{{},{4},{1,2},{2,3,0},{3,0,0,0},{3,3,3,0}}) {
    bool caught=false; try {expand({bad});} catch(const std::invalid_argument&) {caught=true;}
    assert(caught);
  }
  std::cout<<std::setprecision(17);
  for(const auto& c:all) {
    std::cout<<"K";for(int x:c.key)std::cout<<' '<<x;
    std::cout<<" : "<<adjointEta(c.key)<<'\n';
    for(const auto& t:c.paths) {
      std::cout<<"P";for(int x:t.first)std::cout<<' '<<x;
      std::cout<<" : "<<t.second.real()<<' '<<t.second.imag()<<'\n';
    }
  }
}
