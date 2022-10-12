//file for interpolating variances

#include "chromabase.h"
#include "meas/hadron/recursive_interpolation.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>


namespace Chroma {

namespace RecInterp{

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
	level_vals.resize(3);
	level_vars.resize(3);
	for (int i = 1; i < costs.shifts.size(); i++){
		tmpmincost = std::sqrt((costs.level_costs[0] + 2*costs.level_costs[i])*costs.rrs_variances[0][i]) + std::sqrt(costs.level_costs[i]*costs.rr_variances[0][i]) + std::sqrt(costs.level_costs[i]*costs.r_variances[i]);
		level_vars[0] = costs.rrs_variances[0][i];
		level_vars[1] = costs.rr_variances[0][i];
		level_vars[2] = costs.r_variances[i];
		level_vals[0] = std::sqrt((costs.level_costs[0] + 2*costs.level_costs[i])*costs.rrs_variances[0][i]);
		level_vals[1] = std::sqrt(costs.level_costs[i]*costs.rr_variances[0][i]);
		level_vals[2] = std::sqrt(costs.level_costs[i]*costs.r_variances[i]);
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
	level_vars.resize(5);
	level_vals.resize(5);
	for (int i = 1; i < costs.shifts.size(); i++){
		level_vals[0] = std::sqrt((costs.level_costs[0] + 2*costs.level_costs[i])*costs.rrs_variances[0][i]);
		level_vars[0] = costs.rrs_variances[0][i];
		level_vals[1] = std::sqrt(costs.level_costs[i] * costs.rr_variances[0][i]);
		level_vars[1] = costs.rr_variances[0][i];
		for (int j = i+1; j < costs.shifts.size(); j++){
		level_vals[2] = std::sqrt((costs.level_costs[i] + 2*costs.level_costs[j])*costs.rrs_variances[i][j]);
		level_vars[2] = costs.rrs_variances[i][j];
		level_vals[3] = std::sqrt( costs.level_costs[j] * costs.rr_variances[i][j]);
		level_vars[3] = costs.rr_variances[i][j];
		level_vals[4] = std::sqrt(costs.level_costs[j] * costs.r_variances[j]);
		level_vars[4] = costs.r_variances[j];
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
        level_vals.resize(7);
	level_vars.resize(7);
        for (int i = 1; i < costs.shifts.size(); i++){
                level_vals[0] = std::sqrt((costs.level_costs[0] + 2*costs.level_costs[i])*costs.rrs_variances[0][i]);
		level_vars[0] = costs.rrs_variances[0][i];
		level_vals[1] = std::sqrt( costs.level_costs[i] * costs.rr_variances[0][i] );
		level_vars[1] = costs.rr_variances[0][i];
                for (int j = i+1; j < costs.shifts.size(); j++){
                level_vals[2] = std::sqrt((costs.level_costs[i] + 2*costs.level_costs[j])*costs.rrs_variances[i][j]);
		level_vars[2] = costs.rrs_variances[i][j];
		level_vals[3] = std::sqrt( costs.level_costs[j] * costs.rr_variances[i][j]);
		level_vars[3] = costs.rr_variances[i][j];
			for (int k = j+1; k < costs.shifts.size(); k++){
			level_vals[4] = std::sqrt( (costs.level_costs[j] + 2*costs.level_costs[k]) * costs.rrs_variances[j][k]);
			level_vars[4] =  costs.rrs_variances[j][k];
			level_vals[5] = std::sqrt( costs.level_costs[k] * costs.rr_variances[j][k]);
			level_vars[5] = costs.rr_variances[j][k];
                	level_vals[6] = std::sqrt(costs.level_costs[k] * costs.r_variances[k]);
			level_vars[6] = costs.r_variances[k];
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
        level_vals.resize(9);
	level_vars.resize(9);
        for (int i = 1; i < costs.shifts.size(); i++){
                level_vals[0] = std::sqrt((costs.level_costs[0] + 2*costs.level_costs[i])*costs.rrs_variances[0][i]);
		level_vars[0] = costs.rrs_variances[0][i];
		level_vals[1] = std::sqrt( costs.level_costs[i] * costs.rr_variances[0][i]);
		level_vars[1] = costs.rr_variances[0][i];
                for (int j = i+1; j < costs.shifts.size(); j++){
                level_vals[2] = std::sqrt((costs.level_costs[i] + 2*costs.level_costs[j])*costs.rrs_variances[i][j]);
		level_vars[2] = costs.rrs_variances[i][j];
		level_vals[3] = std::sqrt( costs.level_costs[j] * costs.rr_variances[i][j]);
		level_vars[3] = costs.rr_variances[i][j];
                        for (int k = j+1; k < costs.shifts.size(); k++){
                        level_vals[4] = std::sqrt( (costs.level_costs[j] + 2*costs.level_costs[k]) * costs.rrs_variances[j][k]);
			level_vars[4] = costs.rrs_variances[j][k];
			level_vals[5] = std::sqrt( costs.level_costs[k] * costs.rr_variances[j][k]);
			level_vars[5] = costs.rr_variances[j][k];
			for (int l = k+1; l < costs.shifts.size(); l++){
			level_vals[6] = std::sqrt( (costs.level_costs[k] + 2*costs.level_costs[l])*costs.rrs_variances[k][l]);
			level_vars[6] = costs.rrs_variances[k][l];
			level_vals[7] = std::sqrt( costs.level_costs[l] * costs.rr_variances[k][l]);
			level_vars[7] = costs.rr_variances[k][l];
                        level_vals[8] = std::sqrt(costs.level_costs[l] * costs.r_variances[l]);
			level_vars[8] = costs.r_variances[l];
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
        level_vals.resize(11);
	level_vars.resize(11);
        for (int i = 1; i < costs.shifts.size(); i++){
                level_vals[0] = std::sqrt((costs.level_costs[0] + 2*costs.level_costs[i])*costs.rrs_variances[0][i]);
		level_vars[0] = costs.rrs_variances[0][i];
		level_vals[1] = std::sqrt( costs.level_costs[i] * costs.rr_variances[0][i]);
		level_vars[1] = costs.rr_variances[0][i];
                for (int j = i+1; j < costs.shifts.size(); j++){
                level_vals[2] = std::sqrt((costs.level_costs[i] + 2*costs.level_costs[j])*costs.rrs_variances[i][j]);
		level_vars[2] = costs.rrs_variances[i][j];
		level_vals[3] = std::sqrt( costs.level_costs[j] * costs.rr_variances[i][j]);
		level_vars[3] = costs.rr_variances[i][j];
                        for (int k = j+1; k < costs.shifts.size(); k++){
                        level_vals[4] = std::sqrt( (costs.level_costs[j] + 2*costs.level_costs[k]) * costs.rrs_variances[j][k]);
			level_vars[4] = costs.rrs_variances[j][k];
			level_vals[5] = std::sqrt( costs.level_costs[k] * costs.rr_variances[j][k]);
			level_vars[5] = costs.rr_variances[j][k];
                        for (int l = k+1; l < costs.shifts.size(); l++){
                        level_vals[6] = std::sqrt( (costs.level_costs[k] + costs.level_costs[l])*costs.rrs_variances[k][l]);
			level_vars[6] = costs.rrs_variances[k][l];
			level_vals[7] = std::sqrt( costs.level_costs[l] * costs.rr_variances[k][l]);
			level_vars[7] = costs.rr_variances[k][l];
			for (int m = l+1; m < costs.shifts.size(); m++){
			level_vals[8] = std::sqrt( (costs.level_costs[l] + 2*costs.level_costs[m]) * costs.rrs_variances[l][m]);
			level_vars[8] = costs.rrs_variances[l][m];
			level_vals[9] = std::sqrt( costs.level_costs[m] * costs.rr_variances[l][m]);
			level_vars[9] = costs.rr_variances[l][m];
                        level_vals[10] = std::sqrt(costs.level_costs[m] * costs.r_variances[m]);
			level_vars[10] = costs.r_variances[m];
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
        level_vals.resize(13);
	level_vars.resize(13);
        for (int i = 1; i < costs.shifts.size(); i++){
                level_vals[0] = std::sqrt((costs.level_costs[0] + 2*costs.level_costs[i])*costs.rrs_variances[0][i]);
		level_vars[0] = costs.rrs_variances[0][i];
		level_vals[1] = std::sqrt( costs.level_costs[i] * costs.rr_variances[0][i]);
		level_vars[1] = costs.rr_variances[0][i];
                for (int j = i+1; j < costs.shifts.size(); j++){
                level_vals[2] = std::sqrt((costs.level_costs[i] + 2*costs.level_costs[j])*costs.rrs_variances[i][j]);
		level_vars[2] = costs.rrs_variances[i][j];
		level_vals[3] = std::sqrt( costs.level_costs[j] * costs.rr_variances[i][j]);
		level_vars[3] = costs.rr_variances[i][j];
                        for (int k = j+1; k < costs.shifts.size(); k++){
                        level_vals[4] = std::sqrt( (costs.level_costs[j] + 2*costs.level_costs[k]) * costs.rrs_variances[j][k]);
			level_vars[4] = costs.rrs_variances[j][k];
			level_vals[5] = std::sqrt( costs.level_costs[k] * costs.rr_variances[j][k]);
			level_vars[5] = costs.rr_variances[j][k];
                        for (int l = k+1; l < costs.shifts.size(); l++){
                        level_vals[6] = std::sqrt( (costs.level_costs[k] + 2*costs.level_costs[l])*costs.rrs_variances[k][l]);
			level_vars[6] = costs.rrs_variances[k][l];
			level_vals[7] = std::sqrt( costs.level_costs[l] * costs.rr_variances[k][l]);
			level_vars[7] = costs.rr_variances[k][l];
                        for (int m = l+1; m < costs.shifts.size(); m++){
			level_vars[8] =  costs.rrs_variances[l][m];
                        level_vals[8] = std::sqrt( (costs.level_costs[l] + 2*costs.level_costs[m]) * costs.rrs_variances[l][m]);
			level_vals[9] = std::sqrt( costs.level_costs[m] * costs.rr_variances[l][m]);
			level_vars[9] = costs.rr_variances[l][m];
			for (int n = m+1; n < costs.shifts.size(); n++){
			level_vars[10] = costs.rrs_variances[m][n];
			level_vals[10] = std::sqrt( (costs.level_costs[m] + 2*costs.level_costs[n]) * costs.rrs_variances[m][n]);
			level_vals[11] = std::sqrt( costs.level_costs[n] * costs.rr_variances[m][n]);
			level_vars[11] = costs.rr_variances[m][n];
                        level_vals[12] = std::sqrt(costs.level_costs[n] * costs.r_variances[n]);
			level_vars[12] = costs.r_variances[n];
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
	break;
	}else{
	mincosts.shifts.resize(1);
	mincosts.optimal_level_costs.resize(3);
	mincosts.optimal_level_variances.resize(3);
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
        break;
        }else{
        mincosts.shifts.resize(2);
        mincosts.optimal_level_costs.resize(5);
        mincosts.optimal_level_variances.resize(5);
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
        break;
        }else{
        mincosts.shifts.resize(3);
        mincosts.optimal_level_costs.resize(7);
        mincosts.optimal_level_variances.resize(7);
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
        break;
        }else{
        mincosts.shifts.resize(4);
        mincosts.optimal_level_costs.resize(9);
        mincosts.optimal_level_variances.resize(9);
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
        break;
        }else{
        mincosts.shifts.resize(5);
        mincosts.optimal_level_costs.resize(11);
        mincosts.optimal_level_variances.resize(11);
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
        break;
        }else{
        mincosts.shifts.resize(6);
        mincosts.optimal_level_costs.resize(13);
        mincosts.optimal_level_variances.resize(13);
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
	/*} else if (shifts.size() == 3){
	   for (int i = 0; i < int_vals.size(); i++){
		if (intshifts[i] < shifts[1]){
		   int_vals[i] = linearInterpolation(vars[0], vars[1], shifts[0], shifts[1], intshifts[i]);
		}else{
		   int_vals[i] = linearInterpolation(vars[1], vars[2], shifts[1], shifts[2], intshifts[i]);
		}
	   }*/
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



CostHolder_t interpolate_variances(const CostHolder_t& costs, multi1d<double>& dels, multi1d<int>& num_bcshifts, const bool& use_mg)
{

	CostHolder_t int_vars;

	//for now, hardcoding in some small shifts...blehhhh
	//this needs to change!!! The small shifts should be user defined
	multi1d<double> bcs;
	bcs.resize(2);
	bcs[0] = -5.0;
	bcs[1] = -2.0;
	multi1d<int> nums;
	nums.resize(1);
	nums[0] = 4;
	std::vector<double> lowshifts = interpolationShifts(bcs, nums);
	//std::cout << "The smallest shifts for interpolation are : " << std::endl;
	//for (auto i : lowshifts) {std::cout << i << std::endl; }

	QDPIO::cout << "Setting the rest of the shifts" << std::endl;
	std::vector<double> tmpintshifts = interpolationShifts(dels, num_bcshifts);
	std::vector<double> intshifts;
	for (int i = 0; i < lowshifts.size(); i++){intshifts.push_back(lowshifts[i]);}
	for (int i = 0; i < tmpintshifts.size(); i++){intshifts.push_back(tmpintshifts[i]);}
	tmpintshifts.clear();
	lowshifts.clear();
	
	QDPIO::cout << "Finding the union of the shifts" << std::endl;
	//now have all the shifts used for evaluating the polynomials.
	intshifts = vecunion(intshifts, costs.shifts);

	//allocate the space for the predicted variances
	QDPIO::cout << "Allocating holders for the predicted variances" << std::endl;
	int_vars.r_variances.resize(intshifts.size());	
	int_vars.rrs_variances.resize(intshifts.size(), std::vector<double>(intshifts.size()));
	int_vars.rr_variances.resize(intshifts.size(), std::vector<double>(intshifts.size()));


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
	boost::math::interpolators::pchip_matlab<std::vector<double>> r_spline(std::move(tshifts), std::move(tvars));	
	//assign the values
	QDPIO::cout << "Printing interpolated (D+sI)^{-1} values " << std::endl;
	for (int i = 0; i < intshifts.size(); i++){
		int_vars.r_variances[i] = std::exp(r_spline(intshifts[i]));
		QDPIO::cout << int_vars.r_variances[i] << std::endl;
	}
	tvars.clear();

	//do the next easy part...interpolate the (D+sI)^{-1}\Gamma(D+sI)^{-1} values
	//reassign the shifts because of std::move
	tshifts = costs.shifts;
	for (int i = 0; i < costs.shifts.size(); ++i){
		tvars.push_back(std::log(costs.rr_variances[0][i]));
	}
	boost::math::interpolators::pchip_matlab<std::vector<double>> rr_spline(std::move(tshifts), std::move(tvars));
	//assign the values. Since the double inverse terms always include the same shifts, this is easy
	for (int i = 0; i < intshifts.size(); i++){
		for (int j = i; j < intshifts.size(); j++){
		int_vars.rr_variances[i][j] = std::exp(rr_spline(intshifts[j]));
		}
	}

	


	QDPIO::cout << "Interpolating the boundary values " << std::endl;
	tvars.clear();
	//need to reassign the shifts because of move...
	tshifts = costs.shifts;
	for (int i = 0; i < costs.shifts.size(); ++i){
		QDPIO::cout << "On boundary " << i << std::endl;
		//tvars[i] = std::log(costs.rs_variances[i][i]);
		tvars.push_back(std::log(costs.rrs_variances[i][i]));
	}
	QDPIO::cout << "Evaluating polynomial " << std::endl;
	boost::math::interpolators::pchip_matlab<std::vector<double>> diag_spline(std::move(tshifts), std::move(tvars));
	for (int i = 0; i < intshifts.size(); ++i){
		int_vars.rrs_variances[i][i] = std::exp(diag_spline(intshifts[i]));
	}

	


	//here need to interpolate down first to get the values at some shift
	QDPIO::cout << "Interpolating vertically" << std::endl;
	std::vector<std::vector<double>> tmpvals_d;
	tmpvals_d.resize(costs.shifts.size()-1, std::vector<double>(costs.shifts.size()));
	//okay, I know that this is vertical interpolation, but we have switched things up so 
	//the matlab BC's for the derivatives are better than boosts, so the flag is horizontal
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
			tvars.push_back(std::log(costs.rrs_variances[j][i]));
		}
		tmpvals_d[i-1] = interpolate(tshifts, tvars, s_intshifts, dir);
		for (int j = 0; j < tmpvals_d[i-1].size(); j++){
			tmpvals_d[i-1][j] = std::exp(tmpvals_d[i-1][j]);
		}
		
	}



	// this portion does the interpolation horizontally across the region defined by the sampled 
	// variances of the product terms
	QDPIO::cout << "Interpolating horizontally " << std::endl;
	int k = 0; int i = 0; int t = 0;
	std::vector<int> ix = intersect(slice(costs.shifts, 1, costs.shifts.size()-1), intshifts);
	//for (auto ii : ix){ QDPIO::cout << ii << std::endl;}

	while (k < intshifts.size()-1){
	std::vector<double>::iterator it = int_vars.rrs_variances[k].begin();
	it += k;
	std::vector<double> tvars;
	if (!isequal(k, ix[i]))
	{
	tvars.push_back(std::log(int_vars.rrs_variances[k][k]));
	for (int j = i; j < costs.shifts.size()-1; j++){
	    tvars.push_back(std::log(tmpvals_d[j][k]));
	}
	} else {
        for (int j = i; j < costs.shifts.size()-1; j++){
            tvars.push_back(std::log(tmpvals_d[j][k]));
	}
	}
	
	std::vector<double> tshifts = slice(costs.shifts, t, costs.shifts.size()-1);
	if (!float_ismember(tshifts, intshifts[k]))
	{
	tshifts.erase(tshifts.begin());
	tshifts.insert(tshifts.begin(), intshifts[k]);
	}
	std::vector<double> s_intshifts = slice(intshifts, k, intshifts.size()-1);

	//control floating point errors
	s_intshifts.erase(s_intshifts.begin());
	s_intshifts.insert(s_intshifts.begin(), tshifts[0]);

	tvars = interpolate(tshifts, tvars, s_intshifts, dir);
	for (int j = 0; j < tvars.size(); j++){
		tvars[j] = std::exp(tvars[j]);
	}
	int_vars.rrs_variances[k].insert(it, tvars.begin(), tvars.end());

	k += 1;
	if ( k > ix[i]){ i += 1;}
	if (isequal(k, ix[i])) {t += 1;}
	}
	int_vars.rrs_variances[k].insert(int_vars.rrs_variances[k].end(), costs.rrs_variances[costs.shifts.size()-1][costs.shifts.size()-1]);



	QDPIO::cout << "Multiplying by the shifts..." << std::endl;
	double sd;
	for (int i = 0; i < intshifts.size(); ++i){
		for (int j = i; j < intshifts.size(); ++j){
		sd = intshifts[j]-intshifts[i];
		int_vars.rrs_variances[i][j] = std::pow(sd,4.0) * int_vars.rrs_variances[i][j];
		int_vars.rr_variances[i][j] = std::pow(sd,2.0) * int_vars.rr_variances[i][j];
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
	std::vector<MinCosts_t> stats;
	QDPIO::cout << "In recursive getMinShifts" << std::endl;
	//bleh, hardcoding 6 for now as nothing more is gained after 6 shifts;
	//this needs to change because this is likely dependedent on lattice size
	//but to find minimum cost all shifts have to be iterated through
	//and you have to know the number of shifts....
	stats.resize(6);
	
	CostHolder_t int_costs;
	int_costs = interpolate_variances(costs, dels, num_bcshifts, use_mg);

	for (int i = 0; i < int_costs.shifts.size(); i++){
	QDPIO::cout << "Printing row " << i << " of triple inverse term" << std::endl;
	for (int j = i; j < int_costs.shifts.size(); j++){
	QDPIO::cout << int_costs.rrs_variances[i][j] << std::endl;
	}
	}

	for (int i = 0; i < int_costs.shifts.size(); i++){
	QDPIO::cout << " Printing row " << i << " of double inverse term" << std::endl;
	for (int j = i; j < int_costs.shifts.size(); j++){
	QDPIO::cout << int_costs.rr_variances[i][j] << std::endl;
	}
	}
	
	QDPIO::cout << "Finding the minimum shifts for Displacement = " << disp << " and Gamma = " << gamma << std::endl; 
	for (int i = 0; i < 6; i++){
		stats[i] = findMinShifts(i+1,int_costs);
	}
	return stats; 		
}

} //RecInterp namespace

} //Chroma namespace
