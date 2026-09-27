#include "meas/inline/hadron/baryon_derivative_ylm.h"
#include "util/ferm/clebsch.h"
#include <cmath>
#include <set>
#include <stdexcept>
namespace Chroma { namespace BaryonDerivativeYlm {
namespace {
using Complex = std::complex<double>;
std::map<int, Complex> circular(int m) {
  // Redstar threequark_deriv_ops.cc, rightDerivT1; spatial paths are 1-based.
  if (m == 1) return {{1, Complex(0,0.5)}, {2, Complex(-0.5,0)}};
  if (m == 0) return {{3, Complex(0,-1/std::sqrt(2.0))}};
  if (m == -1) return {{1, Complex(0,-0.5)}, {2, Complex(-0.5,0)}};
  throw std::invalid_argument("Baryon YLM magnetic component out of range");
}
}
Component component(const Key& key) {
  Component c{key,{}};
  if (key == Key{0}) { c.paths[Path{}] = 1; return c; }
  if (key.size() == 2 && key[0] == 3 && std::abs(key[1]) <= 1) {
    for (const auto& a : circular(key[1])) {
      Path p{}; p[2] = {a.first}; c.paths[p] += a.second;
    }
    return c;
  }
  if (key.size() != 3 || (key[0] != 33 && key[0] != 23) ||
      key[1] < 0 || key[1] > 2 || std::abs(key[2]) > key[1])
    throw std::invalid_argument("Expected baryon YLM (0), (3,m), (33,L,M), or (23,L,M)");
  for (int a=-1; a<=1; ++a) for (int b=-1; b<=1; ++b) {
    // EVERY angular momentum and projection argument uses twice its physical value.
    const double cg = Hadron::clebsch(2,2*a,2,2*b,2*key[1],2*key[2]);
    if (std::abs(cg) < 1e-14) continue;
    for (const auto& x : circular(a)) for (const auto& y : circular(b)) {
      Path p{};
      if (key[0] == 33) p[2] = {y.first,x.first}; // d_a d_b: apply b first
      else {p[1] = {x.first}; p[2] = {y.first};}
      c.paths[p] += cg*x.second*y.second;
    }
  }
  for (auto it=c.paths.begin(); it!=c.paths.end(); ) {
    if (std::abs(it->second)<1e-14) it=c.paths.erase(it); else ++it;
  }
  return c;
}
std::vector<Component> expand(const std::vector<Key>& requests) {
  std::vector<Component> result;
  std::set<Key> seen;
  for (const auto& key : requests) {
    auto c=component(key);
    if (seen.insert(key).second) result.push_back(c);
  }
  if (result.empty()) throw std::invalid_argument("Empty baryon ylm_list");
  return result;
}
} }
