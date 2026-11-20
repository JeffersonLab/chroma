#ifndef CHROMA_BARYON_DERIVATIVE_YLM_H
#define CHROMA_BARYON_DERIVATIVE_YLM_H
#include <array>
#include <complex>
#include <map>
#include <vector>
namespace Chroma { namespace BaryonDerivativeYlm {
using Key = std::vector<int>;
using Path = std::array<std::vector<int>, 3>;
struct Component { Key key; std::map<Path, std::complex<double>> paths; };
// Explicit physical-integer descriptors: (0), (3,m), (33,L,M), (23,L,M).
Component component(const Key& key);
std::vector<Component> expand(const std::vector<Key>& requests);
} }
#endif
