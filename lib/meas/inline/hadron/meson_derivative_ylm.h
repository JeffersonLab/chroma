// -*- C++ -*-
#ifndef CHROMA_MESON_DERIVATIVE_YLM_H
#define CHROMA_MESON_DERIVATIVE_YLM_H

// Spin-independent Redstar derivative conventions. No QDP/Superb dependency.
#include <complex>
#include <map>
#include <vector>
namespace Chroma { namespace MesonDerivativeYlm {
using Complex = std::complex<double>;
using Key = std::vector<int>;
using Expansion = std::map<Key, Complex>;
struct Component { Key key; Expansion paths; };
Component component(const Key& key);
// Singleton derivative-order requests expand to every allowed coupling and m.
std::vector<Component> expand(const std::vector<Key>& requests, bool drop_negative_m=false);
// Phi_m(p)^dagger = eta (-1)^m Phi_-m(-p), for equal vector phasings.
int adjointEta(const Key& key);
} }
#endif
