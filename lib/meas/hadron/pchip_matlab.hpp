// Copyright Nick Thompson, 2020
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_INTERPOLATORS_PCHIP_MATLAB_HPP
#define BOOST_MATH_INTERPOLATORS_PCHIP_MATLAB_HPP
#include <memory>
#include <boost/math/interpolators/detail/cubic_hermite_detail.hpp>

namespace boost {
namespace math {
namespace interpolators {

template<class RandomAccessContainer>
class pchip_matlab {
public:
    using Real = typename RandomAccessContainer::value_type;

    pchip_matlab(RandomAccessContainer && x, RandomAccessContainer && y,
          Real left_endpoint_derivative = std::numeric_limits<Real>::quiet_NaN(),
          Real right_endpoint_derivative = std::numeric_limits<Real>::quiet_NaN())
    {
        using std::isnan;
        if (x.size() < 3)
        {
            throw std::domain_error("Must be at least four data points.");
        }
        RandomAccessContainer s(x.size(), std::numeric_limits<Real>::quiet_NaN());
        if (isnan(left_endpoint_derivative))
        {
            // O(h) finite difference derivative:
            // This, I believe, is the only derivative guaranteed to be monotonic:
            // This is great, but lets change it to non centered three point

            //s[0] = (y[1]-y[0])/(x[1]-x[0]);

	    //below is my addition
            Real h1 = (x[0]-x[1]);
	    Real h2 = (x[1]-x[2]);
	    Real del1 = (y[0]-y[1])/h1;
	    Real del2 = (y[1]-y[2])/h2;
	
	    s[0] = ( (2*h1 + h2)*del1 - h1*del2)/(h1 + h2);
	    if ( (s[0] > 0 && del1 < 0) || (s[0] < 0 && del1 > 0) )
	    {
		s[0] = 0;
	    } else if ( ((del1 > 0 && del2 < 0) || (del1 < 0 && del2 > 0)) && (std::abs(s[0]) > std::abs(3*del1)) ) {
		s[0] = 3*del1;
	    }
	    //end of my addition
        }
        else
        {
            s[0] = left_endpoint_derivative;
        }

        for (decltype(s.size()) k = 1; k < s.size()-1; ++k) {
            Real hkm1 = x[k] - x[k-1];
            Real dkm1 = (y[k] - y[k-1])/hkm1;

            Real hk = x[k+1] - x[k];
            Real dk = (y[k+1] - y[k])/hk;
            Real w1 = 2*hk + hkm1;
            Real w2 = hk + 2*hkm1;
            if ( (dk > 0 && dkm1 < 0) || (dk < 0 && dkm1 > 0) || dk == 0 || dkm1 == 0)
            {
                s[k] = 0;
            }
            else
            {
                s[k] = (w1+w2)/(w1/dkm1 + w2/dk);
            }

        }
        // Quadratic extrapolation at the other end:
        auto n = s.size();
        if (isnan(right_endpoint_derivative))
        {
            //s[n-1] = (y[n-1]-y[n-2])/(x[n-1] - x[n-2]);
            Real hn1 = x[n-1] - x[n-2];
	    Real hn2 = x[n-2] - x[n-3];
	    Real deln1 = (y[n-1] - y[n-2])/hn1;
	    Real deln2 = (y[n-2] - y[n-3])/hn2;
	    s[n-1] = ( (2*hn1 + hn2)*deln1 - hn1*deln2)/(hn1+hn2);
            if ( (s[n-1] > 0 && deln1 < 0) || (s[n-1] < 0 && deln1 > 0) )
            {
                s[n-1] = 0;
            } else if ( ((deln1 > 0 && deln2 < 0) || (deln1 < 0 && deln2 > 0)) && (std::abs(s[n-1]) > std::abs(3*deln1))) {
                s[n-1] = 3*deln1;
            }
            //end of my addition 
        }
        else
        {
            s[n-1] = right_endpoint_derivative;
        }
        impl_ = std::make_shared<detail::cubic_hermite_detail<RandomAccessContainer>>(std::move(x), std::move(y), std::move(s));
    }

    Real operator()(Real x) const {
        return impl_->operator()(x);
    }

    Real prime(Real x) const {
        return impl_->prime(x);
    }

    friend std::ostream& operator<<(std::ostream & os, const pchip_matlab & m)
    {
        os << *m.impl_;
        return os;
    }

    void push_back(Real x, Real y) {
        using std::abs;
        using std::isnan;
        if (x <= impl_->x_.back()) {
             throw std::domain_error("Calling push_back must preserve the monotonicity of the x's");
        }
        impl_->x_.push_back(x);
        impl_->y_.push_back(y);
        impl_->dydx_.push_back(std::numeric_limits<Real>::quiet_NaN());
        auto n = impl_->size();
        impl_->dydx_[n-1] = (impl_->y_[n-1]-impl_->y_[n-2])/(impl_->x_[n-1] - impl_->x_[n-2]);
        // Now fix s_[n-2]:
        auto k = n-2;
        Real hkm1 = impl_->x_[k] - impl_->x_[k-1];
        Real dkm1 = (impl_->y_[k] - impl_->y_[k-1])/hkm1;

        Real hk = impl_->x_[k+1] - impl_->x_[k];
        Real dk = (impl_->y_[k+1] - impl_->y_[k])/hk;
        Real w1 = 2*hk + hkm1;
        Real w2 = hk + 2*hkm1;
        if ( (dk > 0 && dkm1 < 0) || (dk < 0 && dkm1 > 0) || dk == 0 || dkm1 == 0)
        {
            impl_->dydx_[k] = 0;
        }
        else
        {
            impl_->dydx_[k] = (w1+w2)/(w1/dkm1 + w2/dk);
        }
    }

private:
    std::shared_ptr<detail::cubic_hermite_detail<RandomAccessContainer>> impl_;
};

}
}
}
#endif
