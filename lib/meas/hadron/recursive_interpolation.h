//file for interpolation of variances for Recursive Frequency Splitting
//

#ifndef __INCLUDE_RECURSIVE_INTERPOLATION__
#define __INCLUDE_RECURSIVE_INTERPOLATION__

#include "chromabase.h"
#include <complex>
#include <utility>
#include <numeric>
#include <boost/math/interpolators/pchip.hpp>
#include "pchip_matlab.hpp"

namespace Chroma {
//Interface for accessing interpolation data

namespace RecInterp {

//struct to hold the data
struct CostHolder_t {

    //the solver costs (iterations)
    std::vector<double> level_costs;

    //the variances of the \Gamma(D+sI)^{-1} terms
    std::vector<double> r_variances;

    //the variances of the (D+s_i I)^{-1}\Gamma(D+s_i I)^{-1} terms
    std::vector<std::vector<double>> rr_variances;

    //the variances of the (D+s_{i+1}I)^{-2}\Gamma(D+s_{i}I)^{-1} terms
    std::vector<std::vector<double>> rrs_variances;

    //the shifts
    std::vector<double> shifts;


};


//struct to hold the predicted data found from minimums
struct MinCosts_t {

    std::vector<double> shifts;
    double min_costs;
    double reg_costs;
    //the individual optimal level costs;
    std::vector<double> optimal_level_costs;
    
    //the individual optimal level variances;
    std::vector<double> optimal_level_variances;

};

//struct InterpPow_t {

//     std::vector<double> p_a;
//     std::vector<double> p_d;

//};

template<typename T>
std::vector<T> slice(std::vector<T> const &v, int m, int n)
{
	auto first = v.cbegin() + m;
	auto last = v.cbegin() + n + 1;
	std::vector<T> vec(first, last);
	return vec;
}

template<typename T>
std::vector<T> vecdiff(std::vector<T> const &x, std::vector<T> const &y)
{
	assert(x.size() == y.size());
	std::vector<T> vec;
	vec.resize(x.size());
	for (int i = 0; i < x.size(); i++){
		vec[i] = x[i]-y[i];
	}
	return vec;
}

template<typename T>
T vecnorm2(std::vector<T> const &x)
{
	T a;
	for (auto i : x){
		a+=i*i;
	}
	a = std::sqrt(a);
	return a;
}	

template<typename T>
std::vector<T> raiseToPow(std::vector<T> const &x, double p)
{
	std::vector<T> vec;
	vec.resize(x.size());
	for (int i = 0; i < x.size(); i++){
		vec[i] = std::pow(x[i], p);
	}
	return vec;
}

template<typename T>
std::vector<T> vecunion(std::vector<T> const &x, std::vector<T> const &y)
{
	std::vector<T> vec;
	vec.reserve(x.size() + y.size());
	vec.insert(vec.end(), x.begin(), x.end());
	vec.insert(vec.end(), y.begin(), y.end());
	sort(vec.begin(), vec.end());
	vec.erase(std::unique(vec.begin(), vec.end(), [](T l, T r){return std::abs(l-r) < 1e-12;}), vec.end());
	
	return vec;
}

template<typename T>
int find_shift(const std::vector<T>& x, const T& y)
{
    int i = 0;
    while( i!= x.size()){
        if (std::abs(x[i]-y) < 1e-12) {
            break;
        }
        i++;
    }
    return i;
}

template<typename T>
bool float_ismember(const std::vector<T>& x, const T& y)
{
	bool out = false;
	for (int i = 0; i < x.size(); i++){
		if (std::abs(x[i]-y) < 1e-12){
		out = true;
		}
	}
	return out;
}

template<typename T>
bool isequal(const T& x, const T& y)
{
	bool out = false;
	if (std::abs(x-y) < 1e-12){
	out = true;
	}
	return out;
}

/*std::vector<int> intersect(const std::vector<double>& x, const std::vector<double>& y)
{
	std::vector<int> z;
	for (int i = 0; i < x.size(); ++i){
	    for (int j = 0; j < y.size(); ++j){
		if ( std::abs(x[i]-y[j]) < 1e-12){ z.push_back(j);}
	    }
	}
	return z;
} */

//linear interpolation for when we need it with small shifts

std::vector<double> interpolationShifts(multi1d<double>& bcs, multi1d<int>& num);

double linearInterpolation(const double& y0, const double& y1, const double& x0, const double& x1, const double& x);

std::vector<int> intersect(const std::vector<double>& x, const std::vector<double>& y);

//MinCosts_t one_shift(MinCosts_t mincosts, const CostHolder_t& costs);
void one_shift(MinCosts_t& mincosts, const CostHolder_t& costs);

//MinCosts_t two_shifts(MinCosts_t mincosts, const CostHolder_t& costs);
void two_shifts(MinCosts_t& mincosts, const CostHolder_t& costs);

//MinCosts_t three_shifts(MinCosts_t mincosts, const CostHolder_t& costs);
void three_shifts(MinCosts_t& mincosts, const CostHolder_t& costs);

//MinCosts_t four_shifts(MinCosts_t mincosts, const CostHolder_t& costs);
void four_shifts(MinCosts_t& mincosts, const CostHolder_t& costs);

//MinCosts_t five_shifts(MinCosts_t mincosts, const CostHolder_t& costs);
void five_shifts(MinCosts_t& mincosts, const CostHolder_t& costs);

//MinCosts_t six_shifts(MinCosts_t mincosts, const CostHolder_t& costs);
void six_shifts(MinCosts_t& mincosts, const CostHolder_t& costs);

MinCosts_t findMinShifts(const int& num_shifts_to_calc, const CostHolder_t& costs);

std::vector<double> interpolate(const std::vector<double>& shifts, const std::vector<double>& vars, const std::vector<double>& intshifts);

//InterpPow_t findP(const CostHolder_t& costs);

CostHolder_t interpolate_variances(const CostHolder_t& costs, multi1d<double>& dels, multi1d<int>& num_bcshifts, const bool& use_mg);

std::vector<MinCosts_t> getMinShifts(const CostHolder_t& costs, multi1d<double>& dels, multi1d<int>& num_bcshifts, const bool& use_mg, int& disp, int& gamma);
					  




} //namespace RecInterp

}
#endif // __INCLUDE_RECURSIVE_INTERPOLATION__
