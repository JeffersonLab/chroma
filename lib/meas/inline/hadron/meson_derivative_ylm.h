// -*- C++ -*-
#ifndef CHROMA_MESON_DERIVATIVE_YLM_H
#define CHROMA_MESON_DERIVATIVE_YLM_H

// Spin-independent Redstar derivative conventions. No QDP/Superb dependency.
#include <algorithm>
#include <cmath>
#include <complex>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

namespace Chroma
{
  namespace MesonDerivativeYlm
  {
    using Complex = std::complex<double>;
    using Key = std::vector<int>;
    using Expansion = std::map<Key, Complex>;
    struct Component { Key key; Expansion paths; };

    inline double factorial(int n) {
      double r = 1;
      for (int i = 2; i <= n; ++i) r *= i;
      return r;
    }
    inline int sign(int n) { return n % 2 == 0 ? 1 : -1; }
    // Condon--Shortley <a,m;b,n|J,M>, integer spins through three derivatives.
    inline double cg(int a, int b, int J, int m, int n, int M) {
      if (std::abs(m)>a || std::abs(n)>b || m+n!=M || std::abs(M)>J || J<std::abs(a-b) || J>a+b) return 0;
      double p = std::sqrt((2*J+1)*factorial(J+a-b)*factorial(J-a+b)*factorial(a+b-J)
			   /factorial(a+b+J+1)*factorial(J+M)*factorial(J-M)*factorial(a-m)*factorial(a+m)
			   *factorial(b-n)*factorial(b+n));
      double s=0;
      for (int k=0; k<=a+b; ++k) {
	int v[]={k,a+b-J-k,a-m-k,b+n-k,J-b+m+k,J-a-n+k};
	double d=1; bool valid=true;
	for (int q:v) { if(q<0) {valid=false; break;} d*=factorial(q); }
	if(valid) s+=sign(k)/d;
      }
      return p*s;
    }
    inline std::map<int,Complex> circular(int m) {
      const double h=1/std::sqrt(2.0);
      if(m==-1) return {{1,Complex(0,h)},{2,Complex(h,0)}};
      if(m==0) return {{3,Complex(0,1)}};
      if(m==1) return {{1,Complex(0,-h)},{2,Complex(h,0)}};
      throw std::invalid_argument("Invalid magnetic component");
    }
    inline void append(Expansion& out, const std::vector<int>& ms, double weight,
		       int pos, Key path, Complex coefficient) {
      if(pos<0) {out[path]+=weight*coefficient; return;}
      for(const auto& d:circular(ms[pos])) {
	Key next=path; next.push_back(d.first);
	append(out,ms,weight,pos-1,next,coefficient*d.second);
      }
    }
    inline Component component(const Key& key) {
      if(key.empty()) throw std::invalid_argument("Empty YLM descriptor");
      int n=key[0], J=0, j13=0, M=0;
      if(n==0 && key.size()==1) return {key,{{{},Complex(1,0)}}};
      if(n==1 && key.size()==2) {J=1; M=key[1];}
      else if(n==2 && key.size()==3) {J=key[1]; M=key[2]; if(J<0||J>2) throw std::invalid_argument("Invalid D2 J");}
      else if(n==3 && key.size()==4) {
	j13=key[1]; J=key[2]; M=key[3];
	if(j13<0||j13>2||J<std::abs(j13-1)||J>j13+1) throw std::invalid_argument("Invalid D3 coupling");
      } else throw std::invalid_argument("Expected (0), (1,m), (2,J,m), or (3,J13,J,m)");
      if(std::abs(M)>J) throw std::invalid_argument("YLM m outside [-J,J]");
      Component r{key,{}};
      for(int m1=-1;m1<=1;++m1) {
	if(n==1) {if(m1==M) append(r.paths,{m1},1,0,{},1); continue;}
	for(int m2=-1;m2<=1;++m2) {
	  if(n==2) {append(r.paths,{m1,m2},cg(1,1,J,m1,m2,M),1,{},1); continue;}
	  for(int m3=-1;m3<=1;++m3)
	    append(r.paths,{m1,m2,m3},cg(1,1,j13,m1,m3,m1+m3)*cg(j13,1,J,m1+m3,m2,M),2,{},1);
	}
      }
      for(auto it=r.paths.begin();it!=r.paths.end();) {
	if(std::abs(it->second)<1e-14) it=r.paths.erase(it); else ++it;
      }
      return r;
    }
    // A single integer requests all allowed couplings and m at that order.
    inline std::vector<Component> expand(const std::vector<Key>& requests, bool drop_negative_m=false) {
      std::vector<Component> result; std::set<Key> seen;
      auto add=[&](Key k) {
	Component c=component(k); // validate even components excluded by the storage policy
	if(drop_negative_m && k.size()>1 && k.back()<0) return;
	if(seen.insert(k).second) result.push_back(c);
      };
      for(const Key& k:requests) {
	if(k.size()!=1 || k[0]==0) {add(k); continue;}
	if(k[0]==1) for(int m=-1;m<=1;++m) add({1,m});
	else if(k[0]==2) for(int J=0;J<=2;++J) for(int m=-J;m<=J;++m) add({2,J,m});
	else if(k[0]==3) for(int j13=0;j13<=2;++j13) for(int J=std::abs(j13-1);J<=j13+1;++J)
						       for(int m=-J;m<=J;++m) add({3,j13,J,m});
	else throw std::invalid_argument("Supported derivative orders are 0 through 3");
      }
      if(result.empty()) throw std::invalid_argument("No YLM components requested");
      return result;
    }
    // Phi_m(p)^dagger = eta (-1)^m Phi_-m(-p), with equal vector phasings.
    inline int adjointEta(const Key& k) {
      component(k);
      return k[0]==3 ? sign(3-k[2]+k[1]) : 1;
    }
  }
}
#endif
