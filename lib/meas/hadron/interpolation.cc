//file for interpolating variances

#include "chromabase.h"
#include "meas/hadron/interpolation.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>


namespace Chroma {

namespace Interp {

std::vector<double> interpolationShifts(multi1d<double>& bcs, multi1d<int>& num)
{

	//bcs contain the power you want to raise 10 to in order
	//num is the number of points in the interval
	std::vector<double> shifts;
	int a = bcs.size();
	int b = num.size();
	assert( b == a - 1);
	int i = 0;
	double start;
	double stop;
	double tmp;
	double last;
	double step;
	while (i < a - 1)
	{
	start = 1.0 * bcs[i]; stop = 1.0 * bcs[i+1];
	step = std::abs(start - stop)/(1.0*num[i]-1);
	tmp = start;
	last = stop;
		while( tmp <= last){
		shifts.push_back(std::pow(10.0, tmp));
		tmp += step;
		}
	i += 1;
	}

	return shifts;

}


double linearInterpolation(const double& y0, const double& y1, const double& x0, const double& x1, const double& x)
{
	double y;
	return y = y0 + (x - x0) * ( ( y1 - y0) / (x1 - x0));
}

std::vector<int> intersect(const std::vector<double>& x, const std::vector<double>& y)
{
        std::vector<int> z;
        for (int i = 0; i < x.size(); ++i){
            for (int j = 0; j < y.size(); ++j){
                if ( std::abs(x[i]-y[j]) < 1e-12){ z.push_back(j);}
            }
        }
        return z;
}

void one_shift(MinCosts_t& mincosts, const CostHolder_t& costs){

        //QDPIO::cout << "Calculating the minimum 1 shift " << std::endl;
	//mincosts.shifts.resize(1);
	mincosts.reg_costs = costs.level_costs[0]*costs.r_variances[0];
	double tmpcost = std::numeric_limits<double>::max();
	double tmpmincost;
	std::vector<double> level_vals;
	std::vector<double> level_vars;
	level_vals.resize(2);
	level_vars.resize(2);
	for (int i = 1; i < costs.shifts.size(); i++){
		tmpmincost = std::sqrt((costs.level_costs[0] + costs.level_costs[i])*costs.rs_variances[0][i]) + std::sqrt(costs.level_costs[i]*costs.r_variances[i]);
		level_vars[0] =costs.rs_variances[0][i];
		level_vars[1] = costs.r_variances[i]; 
		level_vals[0] = std::sqrt((costs.level_costs[0] + costs.level_costs[i])*costs.rs_variances[0][i]);
		level_vals[1] = std::sqrt(costs.level_costs[i]*costs.r_variances[i]);
		if (std::pow(tmpmincost, 2.0) < tmpcost){
		 	tmpcost = std::pow(tmpmincost,2.0);
			mincosts.min_costs =  std::pow(tmpmincost,2.0);
			mincosts.optimal_level_costs = level_vals;
			mincosts.optimal_level_variances = level_vars;
			//mincosts.reg_costs = std::pow(tmpmincost,2);
			mincosts.shifts[0] = costs.shifts[i];
		}
	}
	//return mincosts;
}

void two_shifts(MinCosts_t& mincosts, const CostHolder_t& costs){

        //QDPIO::cout << "Calculating the minimum 2 shifts " << std::endl;
	//mincosts.shifts.resize(2);
	mincosts.reg_costs = costs.level_costs[0]*costs.r_variances[0];
	double tmpcost = std::numeric_limits<double>::max();
        double tmpmincost;
	std::vector<double> level_vals;
	std::vector<double> level_vars;
	level_vars.resize(3);
	level_vals.resize(3);
	for (int i = 1; i < costs.shifts.size(); i++){
		level_vals[0] = std::sqrt((costs.level_costs[0] + costs.level_costs[i])*costs.rs_variances[0][i]);
		level_vars[0] = costs.rs_variances[0][i];
		for (int j = i+1; j < costs.shifts.size(); j++){
		level_vals[1] = std::sqrt((costs.level_costs[i] + costs.level_costs[j])*costs.rs_variances[i][j]);
		level_vars[1] = costs.rs_variances[i][j];
		level_vals[2] = std::sqrt(costs.level_costs[j] * costs.r_variances[j]);
		level_vars[2] = costs.r_variances[j];
		tmpmincost = std::pow(std::accumulate(level_vals.begin(), level_vals.end(), 0.0),2.0);
		if (tmpmincost < tmpcost){
			tmpcost = tmpmincost;
			mincosts.min_costs = tmpmincost;
			mincosts.optimal_level_costs = level_vals;
			mincosts.optimal_level_variances = level_vars;
			mincosts.shifts[0] = costs.shifts[i];
			mincosts.shifts[1] = costs.shifts[j];
		}
	}
	}
	//return mincosts;

}

void three_shifts(MinCosts_t& mincosts, const CostHolder_t& costs){

        //QDPIO::cout << "Calculating the minimum 3 shifts " << std::endl;
        //mincosts.shifts.resize(3);
        mincosts.reg_costs = costs.level_costs[0]*costs.r_variances[0];
        double tmpcost = std::numeric_limits<double>::max();
        double tmpmincost;
        std::vector<double> level_vals;
	std::vector<double> level_vars;
        level_vals.resize(4);
	level_vars.resize(4);
        for (int i = 1; i < costs.shifts.size(); i++){
                level_vals[0] = std::sqrt((costs.level_costs[0] + costs.level_costs[i])*costs.rs_variances[0][i]);
		level_vars[0] = costs.rs_variances[0][i];
                for (int j = i+1; j < costs.shifts.size(); j++){
                level_vals[1] = std::sqrt((costs.level_costs[i] + costs.level_costs[j])*costs.rs_variances[i][j]);
		level_vars[1] = costs.rs_variances[i][j];
			for (int k = j+1; k < costs.shifts.size(); k++){
			level_vals[2] = std::sqrt( (costs.level_costs[j] + costs.level_costs[k]) * costs.rs_variances[j][k]);
			level_vars[2] =  costs.rs_variances[j][k];
                	level_vals[3] = std::sqrt(costs.level_costs[k] * costs.r_variances[k]);
			level_vars[3] = costs.r_variances[k];
                	tmpmincost = std::pow(std::accumulate(level_vals.begin(), level_vals.end(), 0.0),2.0);
                	if (tmpmincost < tmpcost){
                        	tmpcost = tmpmincost;
				mincosts.optimal_level_costs = level_vals;
				mincosts.optimal_level_variances = level_vars;
                        	mincosts.min_costs = tmpmincost;
                        	mincosts.shifts[0] = costs.shifts[i];
                        	mincosts.shifts[1] = costs.shifts[j];
				mincosts.shifts[2] = costs.shifts[k];
                	}
			}
              }
        }
        //return mincosts;

}

void four_shifts(MinCosts_t& mincosts, const CostHolder_t& costs){

        //QDPIO::cout << "Calculating the minimum 4 shifts " << std::endl;
        //mincosts.shifts.resize(4);
        mincosts.reg_costs = costs.level_costs[0]*costs.r_variances[0];
        double tmpcost = std::numeric_limits<double>::max();
        double tmpmincost;
        std::vector<double> level_vals;
	std::vector<double> level_vars;
        level_vals.resize(5);
	level_vars.resize(5);
        for (int i = 1; i < costs.shifts.size(); i++){
                level_vals[0] = std::sqrt((costs.level_costs[0] + costs.level_costs[i])*costs.rs_variances[0][i]);
		level_vars[0] = costs.rs_variances[0][i];
                for (int j = i+1; j < costs.shifts.size(); j++){
                level_vals[1] = std::sqrt((costs.level_costs[i] + costs.level_costs[j])*costs.rs_variances[i][j]);
		level_vars[1] = costs.rs_variances[i][j];
                        for (int k = j+1; k < costs.shifts.size(); k++){
                        level_vals[2] = std::sqrt( (costs.level_costs[j] + costs.level_costs[k]) * costs.rs_variances[j][k]);
			level_vars[2] = costs.rs_variances[j][k];
			for (int l = k+1; l < costs.shifts.size(); l++){
			level_vals[3] = std::sqrt( (costs.level_costs[k] + costs.level_costs[l])*costs.rs_variances[k][l]);
			level_vars[3] = costs.rs_variances[k][l];
                        level_vals[4] = std::sqrt(costs.level_costs[l] * costs.r_variances[l]);
			level_vars[4] = costs.r_variances[l];
                        tmpmincost = std::pow(std::accumulate(level_vals.begin(), level_vals.end(), 0.0),2.0);
                        if (tmpmincost < tmpcost){
                                tmpcost = tmpmincost;
				mincosts.optimal_level_costs = level_vals;
				mincosts.optimal_level_variances = level_vars;
                                mincosts.min_costs = tmpmincost;
                                mincosts.shifts[0] = costs.shifts[i];
                                mincosts.shifts[1] = costs.shifts[j];
                                mincosts.shifts[2] = costs.shifts[k];
				mincosts.shifts[3] = costs.shifts[l];
                        }
                        }
              }
        }
     }
     //return mincosts;

}

void five_shifts(MinCosts_t& mincosts, const CostHolder_t& costs){

        //QDPIO::cout << "Calculating the minimum 5 shifts " << std::endl;
        //mincosts.shifts.resize(5);
        mincosts.reg_costs = costs.level_costs[0]*costs.r_variances[0];
        double tmpcost = std::numeric_limits<double>::max();
        double tmpmincost;
        std::vector<double> level_vals;
	std::vector<double> level_vars;
        level_vals.resize(6);
	level_vars.resize(6);
        for (int i = 1; i < costs.shifts.size(); i++){
                level_vals[0] = std::sqrt((costs.level_costs[0] + costs.level_costs[i])*costs.rs_variances[0][i]);
		level_vars[0] = costs.rs_variances[0][i];
                for (int j = i+1; j < costs.shifts.size(); j++){
                level_vals[1] = std::sqrt((costs.level_costs[i] + costs.level_costs[j])*costs.rs_variances[i][j]);
		level_vars[1] = costs.rs_variances[i][j];
                        for (int k = j+1; k < costs.shifts.size(); k++){
                        level_vals[2] = std::sqrt( (costs.level_costs[j] + costs.level_costs[k]) * costs.rs_variances[j][k]);
			level_vars[2] =  costs.rs_variances[j][k];
                        for (int l = k+1; l < costs.shifts.size(); l++){
                        level_vals[3] = std::sqrt( (costs.level_costs[k] + costs.level_costs[l])*costs.rs_variances[k][l]);
			level_vars[3] = costs.rs_variances[k][l];
			for (int m = l+1; m < costs.shifts.size(); m++){
			level_vals[4] = std::sqrt( (costs.level_costs[l] + costs.level_costs[m]) * costs.rs_variances[l][m]);
			level_vars[4] = costs.rs_variances[l][m];
                        level_vals[5] = std::sqrt(costs.level_costs[m] * costs.r_variances[m]);
			level_vars[5] = costs.r_variances[m];
                        tmpmincost = std::pow(std::accumulate(level_vals.begin(), level_vals.end(), 0.0),2.0);
                        if (tmpmincost < tmpcost){
                                tmpcost = tmpmincost;
                                mincosts.min_costs = tmpmincost;
				mincosts.optimal_level_costs = level_vals;
				mincosts.optimal_level_variances = level_vars;
                                mincosts.shifts[0] = costs.shifts[i];
                                mincosts.shifts[1] = costs.shifts[j];
                                mincosts.shifts[2] = costs.shifts[k];
                                mincosts.shifts[3] = costs.shifts[l];
				mincosts.shifts[4] = costs.shifts[m];
                        }
                        }
			}
              }
        }
     }
     //return mincosts;

}

void six_shifts(MinCosts_t& mincosts, const CostHolder_t& costs){

	//QDPIO::cout << "Calculating the minimum 6 shifts " << std::endl;
        //mincosts.shifts.resize(6);
        mincosts.reg_costs = costs.level_costs[0]*costs.r_variances[0];
        double tmpcost = std::numeric_limits<double>::max();
        double tmpmincost;
        std::vector<double> level_vals;
	std::vector<double> level_vars;
        level_vals.resize(7);
	level_vars.resize(7);
        for (int i = 1; i < costs.shifts.size(); i++){
                level_vals[0] = std::sqrt((costs.level_costs[0] + costs.level_costs[i])*costs.rs_variances[0][i]);
		level_vars[0] = costs.rs_variances[0][i];
                for (int j = i+1; j < costs.shifts.size(); j++){
                level_vals[1] = std::sqrt((costs.level_costs[i] + costs.level_costs[j])*costs.rs_variances[i][j]);
		level_vars[1] = costs.rs_variances[i][j];
                        for (int k = j+1; k < costs.shifts.size(); k++){
                        level_vals[2] = std::sqrt( (costs.level_costs[j] + costs.level_costs[k]) * costs.rs_variances[j][k]);
			level_vars[2] = costs.rs_variances[j][k];
                        for (int l = k+1; l < costs.shifts.size(); l++){
                        level_vals[3] = std::sqrt( (costs.level_costs[k] + costs.level_costs[l])*costs.rs_variances[k][l]);
			level_vars[3] = costs.rs_variances[k][l];
                        for (int m = l+1; m < costs.shifts.size(); m++){
			level_vars[4] =  costs.rs_variances[l][m];
                        level_vals[4] = std::sqrt( (costs.level_costs[l] + costs.level_costs[m]) * costs.rs_variances[l][m]);
			for (int n = m+1; n < costs.shifts.size(); n++){
			level_vars[5] = costs.rs_variances[m][n];
			level_vals[5] = std::sqrt( (costs.level_costs[m] + costs.level_costs[n]) * costs.rs_variances[m][n]);
                        level_vals[6] = std::sqrt(costs.level_costs[n] * costs.r_variances[n]);
			level_vars[6] = costs.r_variances[n];
                        tmpmincost = std::pow(std::accumulate(level_vals.begin(), level_vals.end(), 0.0),2.0);
                        if (tmpmincost < tmpcost){
                                tmpcost = tmpmincost;
                                mincosts.min_costs = tmpmincost;
				mincosts.optimal_level_costs = level_vals;
				mincosts.optimal_level_variances = level_vars;
                                mincosts.shifts[0] = costs.shifts[i];
                                mincosts.shifts[1] = costs.shifts[j];
                                mincosts.shifts[2] = costs.shifts[k];
                                mincosts.shifts[3] = costs.shifts[l];
                                mincosts.shifts[4] = costs.shifts[m];
				mincosts.shifts[5] = costs.shifts[n];
			}
                        }
                        }
                        }
              }
        }
     }
     //return mincosts;

}

MinCosts_t findMinShifts(const int& num_shifts_to_calc, const CostHolder_t& costs)
{

	MinCosts_t mincosts;

	switch (num_shifts_to_calc) {
	case 1:
	if (costs.shifts.size() < 1){
	QDPIO::cout << "Not Enough Shifts" << std::endl;
	//return mincosts;
	break;
	}else{
	//mincosts.min_costs = 0.0;
	//mincosts.reg_costs = 0.0;
	mincosts.shifts.resize(1);
	mincosts.optimal_level_costs.resize(2);
	mincosts.optimal_level_variances.resize(2);
	//mincosts = one_shift(mincosts, costs);
	one_shift(mincosts, costs);
	QDPIO::cout << "The regular cost is : " << mincosts.reg_costs << std::endl;
	QDPIO::cout << "The minimum shifts(s) for " << num_shifts_to_calc << " shift(s) are : " << std::endl;
	for (int i = 0; i < mincosts.shifts.size(); i++){
	QDPIO::cout << mincosts.shifts[i] << std::endl;
	}
	QDPIO::cout << "with minimum cost : " << mincosts.min_costs << std::endl;
	QDPIO::cout << "and predicted speedup : " << mincosts.reg_costs/mincosts.min_costs << std::endl;
	QDPIO::cout << "and variances on each level : " << std::endl;
	for (int i = 0; i < mincosts.optimal_level_variances.size(); i++){
	QDPIO::cout << mincosts.optimal_level_variances[i] << std::endl;
	}
	QDPIO::cout << "and costs on each level : " << std::endl;
        for (int i = 0; i < mincosts.optimal_level_costs.size(); i++){
        QDPIO::cout << mincosts.optimal_level_costs[i] << std::endl;
        } 
	break;
	}
	case 2:
        if (costs.shifts.size() < 2){
        QDPIO::cout << "Not Enough Shifts" << std::endl;
        //return mincosts;
        break;
        }else{
        //mincosts.min_costs = 0.0;
        //mincosts.reg_costs = 0.0;
        mincosts.shifts.resize(2);
        mincosts.optimal_level_costs.resize(3);
        mincosts.optimal_level_variances.resize(3);
	//mincosts = two_shifts(mincosts, costs);
	two_shifts(mincosts, costs);
	QDPIO::cout << "The regular cost is : " << mincosts.reg_costs << std::endl;
        QDPIO::cout << "The minimum shifts(s) for " << num_shifts_to_calc << " shift(s) are : " << std::endl;
        for (int i = 0; i < mincosts.shifts.size(); i++){
        QDPIO::cout << mincosts.shifts[i] << std::endl;
        }
        QDPIO::cout << "with minimum cost : " << mincosts.min_costs << std::endl;
	QDPIO::cout << "and predicted speedup : " << mincosts.reg_costs/mincosts.min_costs << std::endl;
        QDPIO::cout << "and variances on each level : " << std::endl;
        for (int i = 0; i < mincosts.optimal_level_variances.size(); i++){
        QDPIO::cout << mincosts.optimal_level_variances[i] << std::endl;
        }
        QDPIO::cout << "and costs on each level : " << std::endl;
        for (int i = 0; i < mincosts.optimal_level_costs.size(); i++){
        QDPIO::cout << mincosts.optimal_level_costs[i] << std::endl;
        }
	break;
	}
	case 3:
        if (costs.shifts.size() < 3){
        QDPIO::cout << "Not Enough Shifts" << std::endl;
        //return mincosts;
        break;
        }else{
        //mincosts.min_costs = 0.0;
        //mincosts.reg_costs = 0.0;
        mincosts.shifts.resize(3);
        mincosts.optimal_level_costs.resize(4);
        mincosts.optimal_level_variances.resize(4);
	//mincosts = three_shifts(mincosts, costs);
	three_shifts(mincosts, costs);
	QDPIO::cout << "The regular cost is : " << mincosts.reg_costs << std::endl;
        QDPIO::cout << "The minimum shifts(s) for " << num_shifts_to_calc << " shift(s) are : " << std::endl;
        for (int i = 0; i < mincosts.shifts.size(); i++){
        QDPIO::cout << mincosts.shifts[i] << std::endl;
        }
        QDPIO::cout << "with minimum cost : " << mincosts.min_costs << std::endl;
	QDPIO::cout << "and predicted speedup : " << mincosts.reg_costs/mincosts.min_costs << std::endl;
        QDPIO::cout << "and variances on each level : " << std::endl;
        for (int i = 0; i < mincosts.optimal_level_variances.size(); i++){
        QDPIO::cout << mincosts.optimal_level_variances[i] << std::endl;
        }
        QDPIO::cout << "and costs on each level : " << std::endl;
        for (int i = 0; i < mincosts.optimal_level_costs.size(); i++){
        QDPIO::cout << mincosts.optimal_level_costs[i] << std::endl;
        }
	break;
	}
	case 4:
        if (costs.shifts.size() < 4){
        QDPIO::cout << "Not Enough Shifts" << std::endl;
        //return mincosts;
        break;
        }else{
        //mincosts.min_costs = 0.0;
        //mincosts.reg_costs = 0.0;
        mincosts.shifts.resize(4);
        mincosts.optimal_level_costs.resize(5);
        mincosts.optimal_level_variances.resize(5);
	//mincosts = four_shifts(mincosts, costs);
	four_shifts(mincosts, costs);
	QDPIO::cout << "The regular cost is : " << mincosts.reg_costs << std::endl;
        QDPIO::cout << "The minimum shifts(s) for " << num_shifts_to_calc << " shift(s) are : " << std::endl;
        for (int i = 0; i < mincosts.shifts.size(); i++){
        QDPIO::cout << mincosts.shifts[i] << std::endl;
        }
        QDPIO::cout << "with minimum cost : " << mincosts.min_costs << std::endl;
	QDPIO::cout << "and predicted speedup : " << mincosts.reg_costs/mincosts.min_costs << std::endl;
        QDPIO::cout << "and variances on each level : " << std::endl;
        for (int i = 0; i < mincosts.optimal_level_variances.size(); i++){
        QDPIO::cout << mincosts.optimal_level_variances[i] << std::endl;
        }
        QDPIO::cout << "and costs on each level : " << std::endl;
        for (int i = 0; i < mincosts.optimal_level_costs.size(); i++){
        QDPIO::cout << mincosts.optimal_level_costs[i] << std::endl;
        }
	break;
	}
	case 5:
        if (costs.shifts.size() < 5){
        QDPIO::cout << "Not Enough Shifts" << std::endl;
        //return mincosts;
        break;
        }else{
        //mincosts.min_costs = 0.0;
        //mincosts.reg_costs = 0.0;
        mincosts.shifts.resize(5);
        mincosts.optimal_level_costs.resize(6);
        mincosts.optimal_level_variances.resize(6);
	//mincosts = five_shifts(mincosts, costs);
	five_shifts(mincosts, costs);
	QDPIO::cout << "The regular cost is : " << mincosts.reg_costs << std::endl;
        QDPIO::cout << "The minimum shifts(s) for " << num_shifts_to_calc << " shift(s) are : " << std::endl;
        for (int i = 0; i < mincosts.shifts.size(); i++){
        QDPIO::cout << mincosts.shifts[i] << std::endl;
        }
        QDPIO::cout << "with minimum cost : " << mincosts.min_costs << std::endl;
	QDPIO::cout << "and predicted speedup : " << mincosts.reg_costs/mincosts.min_costs << std::endl;
        QDPIO::cout << "and variances on each level : " << std::endl;
        for (int i = 0; i < mincosts.optimal_level_variances.size(); i++){
        QDPIO::cout << mincosts.optimal_level_variances[i] << std::endl;
        }
        QDPIO::cout << "and costs on each level : " << std::endl;
        for (int i = 0; i < mincosts.optimal_level_costs.size(); i++){
        QDPIO::cout << mincosts.optimal_level_costs[i] << std::endl;
        }
	break;
	}
	case 6:
        if (costs.shifts.size() < 6){
        QDPIO::cout << "Not Enough Shifts" << std::endl;
        //return mincosts;
        break;
        }else{
        //mincosts.min_costs = 0.0;
        //mincosts.reg_costs = 0.0;
        mincosts.shifts.resize(6);
        mincosts.optimal_level_costs.resize(7);
        mincosts.optimal_level_variances.resize(7);
	//mincosts = six_shifts(mincosts, costs);
	six_shifts(mincosts, costs);
	QDPIO::cout << "The regular cost is : " << mincosts.reg_costs << std::endl;
        QDPIO::cout << "The minimum shifts(s) for " << num_shifts_to_calc << " shift(s) are : " << std::endl;
        for (int i = 0; i < mincosts.shifts.size(); i++){
        QDPIO::cout << mincosts.shifts[i] << std::endl;
        }
        QDPIO::cout << "with minimum cost : " << mincosts.min_costs << std::endl;
	QDPIO::cout << "and predicted speedup : " << mincosts.reg_costs/mincosts.min_costs << std::endl;
        QDPIO::cout << "and variances on each level : " << std::endl;
        for (int i = 0; i < mincosts.optimal_level_variances.size(); i++){
        QDPIO::cout << mincosts.optimal_level_variances[i] << std::endl;
        }
        QDPIO::cout << "and costs on each level : " << std::endl;
        for (int i = 0; i < mincosts.optimal_level_costs.size(); i++){
        QDPIO::cout << mincosts.optimal_level_costs[i] << std::endl;
        }
	break;
	}
	default :
	QDPIO::cout << "Only six shifts and less supported for interpolation. If you want to use more than six shifts" << std::endl;
	QDPIO::cout << "for trace estimation, set your desired shifts and set use_interpolation to false. WARNING: This " << std::endl;
	QDPIO::cout << "may increase the cost of your trace estimation." << std::endl;
	QDP_abort(1);
	}
	return mincosts;


}

std::vector<double> interpolate(const std::vector<double>& shifts, const std::vector<double>& vars, const std::vector<double>& intshifts, const std::string& dir)
{
	assert(shifts.size() == vars.size());
	std::vector<double> int_vals;
	std::vector<double> tmps = shifts;
	std::vector<double> tmpv = vars;
	int_vals.resize(intshifts.size());
	if (shifts.size() >= 3){
	   if (dir == "vertical"){
	   boost::math::interpolators::pchip<std::vector<double>> spline(std::move(tmps), std::move(tmpv));
           for (int i = 0; i < int_vals.size(); i++){
                int_vals[i] = spline(intshifts[i]);
           }
	   } else {
	   boost::math::interpolators::pchip_matlab<std::vector<double>> spline(std::move(tmps), std::move(tmpv));
           for (int i = 0; i < int_vals.size(); i++){
                int_vals[i] = spline(intshifts[i]);
           }
	   }
	} else if (shifts.size() == 3){
	   for (int i = 0; i < int_vals.size(); i++){
		if (intshifts[i] < shifts[1]){
		   int_vals[i] = linearInterpolation(vars[0], vars[1], shifts[0], shifts[1], intshifts[i]);
		}else{
		   int_vals[i] = linearInterpolation(vars[1], vars[2], shifts[1], shifts[2], intshifts[i]);
		}
	   }
	}else if (shifts.size() == 2) {
	   for (int i = 0; i < int_vals.size(); i++){
		int_vals[i] = linearInterpolation(vars[0], vars[1], shifts[0], shifts[1], intshifts[i]);
	   }
	}else{
	QDPIO::cout << "Cannot Interpolate with less than 2 shifts and 2 values. Aborting" << std::endl;
	QDP_abort(1);
	}
		
	return int_vals;

}


/*InterpPow_t findP(const CostHolder_t& costs)
{
	//QDPIO::cout << "In findP " << std::endl;
	//QDPIO::cout << "Using " << costs.shifts.size() << " total shifts " << std::endl;
	InterpPow_t int_p;
	//QDPIO::cout << "Initializing InterpPow struct " << std::endl;
	int_p.p_a.resize(test_shifts_a.size());
	int_p.p_d.resize(test_shifts_d.size());
	std::vector<double> pa;
	std::vector<double> pd;
	//using default values [1:.01:3]
	double tmp = 1.0;
	pa.push_back(tmp);
	//QDPIO::cout << "Initializing the power vector " << std::endl;
	while( tmp < 3.0){
	tmp += 0.01;
	pa.push_back(tmp);
	}
	pd.resize(pa.size());
	pd = pa;

	//first do the variances down.
	//QDPIO::cout << "Finding vertical p " << std::endl;
	std::vector<double> tvars;
	std::vector<double> tshifts;
	std::string dir = "vertical";
	for (int i = 1; i < costs.shifts.size(); i++){
	     tshifts = slice(costs.shifts, 0, i);
	     tvars.clear();
	     for (int j = 0; j < i+1; j++){
	     tvars.push_back(costs.rs_variances[j][i]);
	     }
	     std::vector<double> int_vals;
	     std::vector<double> norms;
	     norms.resize(pd.size());
	     std::vector<double> intshifts = tshifts;
	     //test point will be the last shift here
	     //old
	     //intshifts.insert(intshifts.end(), test_shifts_d[i-1]);
	     //new
	     intshifts.insert(intshifts.end()-1, test_shifts_d[i-1]);
	     std::vector<double> compare_vals = tvars;
	     //test point will be the last value here
	     //old
	     //compare_vals.insert(compare_vals.end(), test_vals_d[i-1]);
	     //new
	     compare_vals.insert(compare_vals.end()-1, test_vals_d[i-1]);
	     //QDPIO::cout << " On column " << i << " and shifts size is " << tshifts.size() << " with sample size " << tvars.size() << std::endl;
	     for (int j = 0; j < pd.size();  j++){
		tvars = raiseToPow(tvars, 1.0/pd[j]);
		int_vals = raiseToPow(interpolate(tshifts, tvars, intshifts, dir), pd[j]);
		tvars = raiseToPow(tvars, pd[j]);
		norms[j] = vecnorm2(vecdiff(int_vals,compare_vals));
	     }
	     std::vector<double>::iterator min_p = std::min_element(norms.begin(), norms.end());
	     int_p.p_d[i-1] = pd[std::distance(norms.begin(), min_p)];
	}		
	
	//QDPIO::cout << "Finding the horizontal p " << std::endl;	
	dir = "horizontal";
	for (int i = 0; i < costs.shifts.size()-1; i++){
	     tshifts = slice(costs.shifts,i,costs.shifts.size()-1);
	     tvars.clear();
	     for (int j = i; j < costs.shifts.size(); j++){
		tvars.push_back(costs.rs_variances[i][j]);
	     }
             std::vector<double> int_vals;
             std::vector<double> norms;
             norms.resize(pa.size());
             std::vector<double> intshifts = tshifts;
	     //test point will be there first shift here
	     //old
	     //intshifts.insert(intshifts.begin(), test_shifts_a[i]);
	     //new
	     intshifts.insert(intshifts.begin()+1, test_shifts_a[i]);
	     std::vector<double> compare_vals = tvars;
	     //test point will be the first value here
	     //old
	     //compare_vals.insert(compare_vals.begin(), test_vals_a[i]);
	     //new
	     compare_vals.insert(compare_vals.begin()+1, test_vals_a[i]);
	     //QDPIO::cout << " On row " << i << " and shifts size is " << tshifts.size() << " with sample size " << tvars.size() << std::endl;
             for (int j = 0; j < pa.size(); j++){
                tvars = raiseToPow(tvars, 1.0/pa[j]);
                int_vals = raiseToPow(interpolate(tshifts, tvars, intshifts, dir), pa[j]);
                tvars = raiseToPow(tvars, pa[j]);
                norms[j] = vecnorm2(vecdiff(int_vals,compare_vals));
             }
             std::vector<double>::iterator min_p = std::min_element(norms.begin(), norms.end());
             int_p.p_a[i] = pa[std::distance(norms.begin(), min_p)];
        }
	return int_p;

} */

CostHolder_t interpolate_variances(const CostHolder_t& costs, multi1d<double>& dels, multi1d<int>& num_bcshifts, const bool& use_mg)
{

	//QDPIO::cout << "In interpolate_variances " << std::endl;
	CostHolder_t int_vars;
	//for now just do a uniform shift-mesh discretization
	/*std::vector<double> intshifts;
	double tmp = 0.0;
	intshifts.push_back(tmp);
	while( tmp < costs.shifts[costs.shifts.size()-1]){
	tmp += dels;
	intshifts.push_back(tmp);
	}
	intshifts = vecunion(intshifts, costs.shifts);
	*/

	//for now, hardcoding in some small shifts...blehhhh
	multi1d<double> bcs;
	bcs.resize(2);
	bcs[0] = -5.0;
	bcs[1] = -2.0;
	multi1d<int> nums;
	nums.resize(1);
	nums[0] = 4;
	std::vector<double> lowshifts = interpolationShifts(bcs, nums);
	//QDPIO::cout << "The smallest shifts for interpolation are : " << std::endl;
	//for (auto i : lowshifts) {QDPIO::cout << i << std::endl; }

	std::vector<double> tmpintshifts = interpolationShifts(dels, num_bcshifts);
	std::vector<double> intshifts;
	for (int i = 0; i < lowshifts.size(); i++){intshifts.push_back(lowshifts[i]);}
	for (int i = 0; i < tmpintshifts.size(); i++){intshifts.push_back(tmpintshifts[i]);}
	tmpintshifts.clear();
	lowshifts.clear();
	
	intshifts = vecunion(intshifts, costs.shifts);
	int_vars.r_variances.resize(intshifts.size());	
	int_vars.rs_variances.resize(intshifts.size(), std::vector<double>(intshifts.size()));

	QDPIO::cout << "Printing the shifts used for interpolation " << std::endl;
	for (int i = 0; i < intshifts.size(); i++){
		QDPIO::cout << intshifts[i] << std::endl;
	} 

	//do the easy part first...interpolate the (D+sI)^{-1}
	std::vector<double> tshifts = costs.shifts;
	std::vector<double> tvars;
	for (int i = 0; i < costs.shifts.size(); ++i){
		tvars.push_back(std::log(costs.r_variances[i]));
	}
	boost::math::interpolators::pchip_matlab<std::vector<double>> spline(std::move(tshifts), std::move(tvars));	
	//assign the values
	QDPIO::cout << "Printing interpolated (D+sI)^{-1} values " << std::endl;
	for (int i = 0; i < intshifts.size(); i++){
		int_vars.r_variances[i] = std::exp(spline(intshifts[i]));
		QDPIO::cout << int_vars.r_variances[i] << std::endl;
	}

	QDPIO::cout << "Interpolating the boundary values " << std::endl;
	tvars.clear();
	//need to reassign the shifts because of move...
	tshifts = costs.shifts;
	for (int i = 0; i < costs.shifts.size(); ++i){
		QDPIO::cout << "On boundary " << i << std::endl;
		//tvars[i] = std::log(costs.rs_variances[i][i]);
		tvars.push_back(std::log(costs.rs_variances[i][i]));
	}
	QDPIO::cout << "Evaluating polynomial " << std::endl;
	boost::math::interpolators::pchip_matlab<std::vector<double>> diag_spline(std::move(tshifts), std::move(tvars));
	for (int i = 0; i < intshifts.size(); ++i){
		int_vars.rs_variances[i][i] = std::exp(diag_spline(intshifts[i]));
	}

	

	//QDPIO::cout << "Printing info from the interpolation: " <<std::endl;
	//QDPIO::cout << spline << std::endl;
	//QDPIO::cout << "Performed interpolation of (D+sI)^{-1} " << std::endl;

	//get the p values
	//InterpPow_t int_p = findP(costs, test_shifts_a, test_vals_a, test_shifts_d, test_vals_d);
 	//print the power values
 	/* QDPIO::cout << "Powers horizontally are : " << std::endl;
 	for (int i = 0; i < int_p.p_a.size(); i++){
		QDPIO::cout << int_p.p_a[i] << std::endl;
	}
	QDPIO::cout << "Powers vertically are : " << std::endl;	
        for (int i = 0; i < int_p.p_d.size(); i++){
                QDPIO::cout << int_p.p_d[i] << std::endl;
        } */

	//int_vars.rs_variances.resize(intshifts.size(), std::vector<double>(intshifts.size()));

	//here need to interpolate down first to get the values at some shift
	QDPIO::cout << "Interpolating vertically" << std::endl;
	std::vector<std::vector<double>> tmpvals_d;
	tmpvals_d.resize(costs.shifts.size()-1, std::vector<double>(costs.shifts.size()));
	//okay, I know that this is vertical interpolation, but we have switched things up so 
	//the matlab BC's are better than boosts, so the flag is horizontal
	std::string dir = "horizontal";
	for ( int i = 1; i < costs.shifts.size(); i++){
		int it = find_shift(intshifts, costs.shifts[i]);
		std::vector<double> s_intshifts = slice(intshifts, 0, it);
		std::vector<double> tshifts = slice(costs.shifts, 0, i);
		std::vector<double> tvars;
		//tvars.resize(tshifts.size());
		tmpvals_d[i-1].resize(s_intshifts.size());
		for (int j = 0; j < i+1; j++){
			//tvars.push_back(costs.rs_variances[j][i]);
			tvars.push_back(std::log(costs.rs_variances[j][i]));
		}
		//tvars = raiseToPow(tvars, 1.0/int_p.p_d[i-1]);
		//tmpvals_d[i-1] = raiseToPow(interpolate(tshifts, tvars, s_intshifts, dir), int_p.p_d[i-1]);
		tmpvals_d[i-1] = interpolate(tshifts, tvars, s_intshifts, dir);
		for (int j = 0; j < tmpvals_d[i-1].size(); j++){
			tmpvals_d[i-1][j] = std::exp(tmpvals_d[i-1][j]);
		}
		
	}

        /*QDPIO::cout << "Printing the intermediate products down " << std::endl;
        for (int i = 0; i < costs.shifts.size()-1; i++){
              QDPIO::cout << "Printing column " << i << std::endl;
              for (int j = 0; j < tmpvals_d[i].size(); j++){
              QDPIO::cout << tmpvals_d[i][j] << std::endl;
              }
        }*/

	QDPIO::cout << "Interpolating horizontally " << std::endl;
	int k = 0; int i = 0; int t = 0;
	std::vector<int> ix = intersect(slice(costs.shifts, 1, costs.shifts.size()-1), intshifts);
	//QDPIO::cout << "Shift intersection is : " << std::endl;
	for (auto ii : ix){ QDPIO::cout << ii << std::endl;}

	//std::vector<double> s_intshifts;
	//QDPIO::cout << "Starting loop " << std::endl;
	while (k < intshifts.size()-1){
	//QDPIO::cout << "Grabbing an iterator " << std::endl;
	std::vector<double>::iterator it = int_vars.rs_variances[k].begin();
	it += k;
	std::vector<double> tvars;
	//QDPIO::cout << "Setting the variances to be interpolated" << std::endl;
	if (!isequal(k, ix[i]))
	{
	tvars.push_back(std::log(int_vars.rs_variances[k][k]));
	for (int j = i; j < costs.shifts.size()-1; j++){
	    //QDPIO::cout << "On iteration " << k << " and on column " << j << ", not on a boundary " << std::endl;
	    tvars.push_back(std::log(tmpvals_d[j][k]));
	}
	} else {
        for (int j = i; j < costs.shifts.size()-1; j++){
	    //QDPIO::cout << "On iteration " << k << " and on column " << j << ", on a boundary " << std::endl;
            tvars.push_back(std::log(tmpvals_d[j][k]));
	}
	}
	
	std::vector<double> tshifts = slice(costs.shifts, t, costs.shifts.size()-1);
	//QDPIO::cout << "Setting the shifts to be interpolated" << std::endl;
	if (!float_ismember(tshifts, intshifts[k]))
	{
	tshifts.erase(tshifts.begin());
	tshifts.insert(tshifts.begin(), intshifts[k]);
	}
	std::vector<double> s_intshifts = slice(intshifts, k, intshifts.size()-1);
	//hopefully to control floating point errors...
	s_intshifts.erase(s_intshifts.begin());
	s_intshifts.insert(s_intshifts.begin(), tshifts[0]);
	//QDPIO::cout << "On iteration " << k << " with shifts : " << std::endl;
	//for (auto ii : tshifts){
	//	QDPIO::cout << ii << std::endl;
	//}
	//QDPIO::cout << "and variances : " << std::endl;
	//for (auto ii : tvars){
	//	QDPIO::cout << ii << std::endl;
	//}
	tvars = interpolate(tshifts, tvars, s_intshifts, dir);
	for (int j = 0; j < tvars.size(); j++){
		tvars[j] = std::exp(tvars[j]);
	}
	int_vars.rs_variances[k].insert(it, tvars.begin(), tvars.end());

	k += 1;
	if ( k > ix[i]){ i += 1;}
	if (isequal(k, ix[i])) {t += 1;}
	}
	int_vars.rs_variances[k].insert(int_vars.rs_variances[k].end(), costs.rs_variances[costs.shifts.size()-1][costs.shifts.size()-1]);


	//QDPIO::cout << "Printing the intermediate products down " << std::endl;
	//for (int i = 0; i < costs.shifts.size()-1; i++){
	//	QDPIO::cout << "Printing column " << i << std::endl;
	//	for (int j = 0; j < tmpvals_d[i].size(); j++){
	//	QDPIO::cout << tmpvals_d[i][j] << std::endl;
	//	}
	//}

	//now interopolate across
	//original
	/* QDPIO::cout << "Interpolating horizontally" << std::endl;
	int i = 0;
	int it = 1;
	dir = "horizontal";
	while ( it <= costs.shifts.size()-1){
		std::vector<double>::iterator itt = int_vars.rs_variances[i].begin();
		itt += i;
		std::vector<double> tmp;
		std::vector<double> tshifts = slice(costs.shifts,it,costs.shifts.size()-1);
		tshifts.insert(tshifts.begin(), intshifts[i]);
		std::vector<double> s_intshifts = slice(intshifts,i,intshifts.size()-1);
		tmp.insert(tmp.begin(), 0.0);
		for (int j = it; j < costs.shifts.size(); j++){
			tmp.push_back(tmpvals_d[j-1][i]);
		}
		tmp = raiseToPow(tmp, 1.0/int_p.p_a[it-1]);
		tmp = raiseToPow(interpolate(tshifts, tmp, s_intshifts, dir), int_p.p_a[it-1]);
		int_vars.rs_variances[i].insert(itt, tmp.begin(), tmp.end());
		i++;
		if ( intshifts[i] > costs.shifts[it-1]){
			it++;
		}
	} */

	//new
        /*QDPIO::cout << "Interpolating horizontally" << std::endl;
        int i = 0;
        int it = 1;
        dir = "horizontal";
	bool isin;
	while ( it < costs.shifts.size()-1){
                std::vector<double>::iterator itt = int_vars.rs_variances[i].begin();
                itt += i;
                std::vector<double> tmp;
                std::vector<double> tshifts = slice(costs.shifts,it,costs.shifts.size()-1);
		isin = float_ismember(tshifts, intshifts[i]);
		if (!isin){
			tmp.insert(tmp.begin(), 0.0);
		}		
                for (int j = it; j < costs.shifts.size(); j++){
                        tmp.push_back(tmpvals_d[j-1][i]);
                }
		
		if(!isin){
                tshifts.insert(tshifts.begin(), intshifts[i]);
		}
                std::vector<double> s_intshifts = slice(intshifts,i,intshifts.size()-1);
                //tmp.insert(tmp.begin(), 0.0);
                //for (int j = it; j < costs.shifts.size(); j++){
                //        tmp.push_back(tmpvals_d[j-1][i]);
                //}
                tmp = raiseToPow(tmp, 1.0/int_p.p_a[it-1]);
                tmp = raiseToPow(interpolate(tshifts, tmp, s_intshifts, dir), int_p.p_a[it-1]);
		//old
                //int_vars.rs_variances[i].insert(itt, tmp.begin(), tmp.end());
                //new
                int_vars.rs_variances[i].insert(itt, tmp.begin(), tmp.end());
                i++;
                if ( intshifts[i] > costs.shifts[it]){
                        it++;
                }
        }
	

	//the last section
	QDPIO::cout << "Interpolating last shift partition starting with row " << i <<  std::endl;
	
	//orginal
	//it--;
	
	//new, no decrementing it
	
	//old
	//int j = costs.shifts.size()-2;
	//new
	int j = it-1;
	QDPIO::cout << "It after horizontal is " << it << std::endl;
	while (i < intshifts.size()-1) {
	QDPIO::cout << "On iteration " << i << std::endl;
	std::vector<double>::iterator itt = int_vars.rs_variances[i].begin();
	itt += i;
	std::vector<double> tmp;
	QDPIO::cout << "Preparing variances " << std::endl;
	tmp.push_back(tmpvals_d[j][i]);
	tmp.insert(tmp.begin(), 0);
	QDPIO::cout << "Preparing shifts " << std::endl;
        std::vector<double> tshifts = slice(costs.shifts,it,costs.shifts.size()-1);
        tshifts.insert(tshifts.begin(), intshifts[i]);
	QDPIO::cout << "Preparing interpolation shifts " << std::endl;
        std::vector<double> s_intshifts = slice(intshifts,i,intshifts.size()-1);
        QDPIO::cout << "Size of shifts : " << tshifts.size() << std::endl;
	QDPIO::cout << "Size of interpolating shifts " << s_intshifts.size() << std::endl;
        QDPIO::cout << "Size of variancances : " << tmp.size() << std::endl;
	tmp = raiseToPow(tmp,1.0/int_p.p_a[it-1]);
	QDPIO::cout << "Interpolating the variances " << std::endl;
	tmp = raiseToPow(interpolate(tshifts, tmp, s_intshifts, dir), int_p.p_a[it-1]);
	QDPIO::cout << "Size of interpolated variances : " << tmp.size() << std::endl;
	QDPIO::cout << "Inserting the variances in the final array " << std::endl;
	//old
	//int_vars.rs_variances[i].insert(itt, tmp.begin(), tmp.end());
	//new
	int_vars.rs_variances[i].insert(itt, tmp.begin(), tmp.end());
        i++;
	//itt++;
	} */

	QDPIO::cout << "Multiplying by the shifts..." << std::endl;
	double sd;
	for (int i = 0; i < intshifts.size(); ++i){
		for (int j = i; j < intshifts.size(); ++j){
		sd = intshifts[j]-intshifts[i];
		int_vars.rs_variances[i][j] = std::pow(sd,2.0) * int_vars.rs_variances[i][j];
		}
	}

	

	//now do the iteration costs
	QDPIO::cout << "Setting the level costs" << std::endl;
	//This version below does "nearest neighbor interpolation"
	int_vars.level_costs.resize(intshifts.size());	
	if (use_mg){
	int j = 0;
	int it = 0;
	for (int k = 0; k < intshifts.size(); k++){
		int_vars.level_costs[k] = costs.level_costs[it];
		if (intshifts[k] >= costs.shifts[it]){
		it++;
		}
	}
	}else{
	//this version does pchip interpolation of the costs
	std::vector<double> tmpintcosts = costs.level_costs;
	std::vector<double> tmpintshifts = costs.shifts;
	int_vars.level_costs.resize(intshifts.size());
	boost::math::interpolators::pchip_matlab<std::vector<double>> costsspline(std::move(tmpintshifts), std::move(tmpintcosts));
	for (int ii = 0; ii < intshifts.size(); ++ii){
		int_vars.level_costs[ii] = costsspline(intshifts[ii]);
	}
	}


	QDPIO::cout << "Printing the level costs : " << std::endl;
	for (int ii = 0; ii < int_vars.level_costs.size(); ii++){
	QDPIO::cout << int_vars.level_costs[ii] << std::endl;
	}
	//assign the shifts;
	int_vars.shifts = intshifts;
	return int_vars;

}

std::vector<MinCosts_t> getMinShifts(const CostHolder_t& costs, multi1d<double>& dels, multi1d<int>& num_bcshifts, const bool& use_mg, int& disp, int& gamma)
{
	//QDPIO::cout << "In getMinShifts" << std::endl;
	std::vector<MinCosts_t> stats;
	//bleh, hardcoding 6 for now;
	stats.resize(6);
	
	CostHolder_t int_costs;
	int_costs = interpolate_variances(costs, dels, num_bcshifts, use_mg);

	for (int i = 0; i < int_costs.shifts.size(); i++){
	QDPIO::cout << "Printing row " << i << std::endl;
	for (int j = i; j < int_costs.shifts.size(); j++){
	QDPIO::cout << int_costs.rs_variances[i][j] << std::endl;
	}
	}
	QDPIO::cout << "Finding the minimum shifts for Displacement = " << disp << " and Gamma = " << gamma << std::endl; 
	for (int i = 0; i < 6; i++){
		stats[i] = findMinShifts(i+1,int_costs);
	}
	return stats; 		
}
/* std::vector<int> findGammaDisp(const CostHolder& costs){

	std::vector<int> gk;
	gk.resize(2);
	return gk;

} */


} //namespace

} //Chroma namespace
