// $Id: remez_gmp.cc,v 3.1 2007-04-17 03:13:04 edwards Exp $
/*! \file
 *  \brief Remez algorithm for finding nth roots
 */

#include "update/molecdyn/monomial/remez_gmp.h"

#define JMAX 1000 //Maximum number of iterations of Newton's approximation

namespace Chroma
{

  // Constructor
  RemezGMP::RemezGMP()
  {
    START_CODE();

    alloc = false;

    END_CODE();
  }


  // Constructor
  RemezGMP::RemezGMP(const Real& lower, const Real& upper, long precision) 
  {
    START_CODE();

    prec = precision;
    bigfloat::setDefaultPrecision(prec);
    
    alloc = false;
    apstrt = lower;
    apend = upper;
    apwidt = apend - apstrt;
    QDPIO::cout << __func__ << ": approximation bounds are [" 
		<< (double)apstrt << " , " << (double)apend << "]"<< endl;

    n = 0;
    d = 0;

    tolerance = 1e-15;

    // Only require the approximation spread to be less than tolerance 
    if (sizeof(Real)==sizeof(double)) {
      tolerance = 1e-15;
    } else if (sizeof(Real)==sizeof(float)) {
      tolerance = 1e-8;
    }

    END_CODE();
  }

  // Free memory and reallocate as necessary
  void RemezGMP::allocate(int num_degree, int den_degree)
  {
    START_CODE();

    // Note use of new and delete in memory allocation - cannot run on qcdsp
    param.resize(num_degree+den_degree+1);
    roots.resize(num_degree);
    poles.resize(den_degree);
    xx.resize(num_degree+den_degree+3);
    mm.resize(num_degree+den_degree+2);

    alloc = true;

    END_CODE();
  }

  // Reset the bounds of the approximation
  void RemezGMP::setBounds(const Real& lower, const Real& upper)
  {
    START_CODE();

    apstrt = lower;
    apend = upper;
    apwidt = apend - apstrt;

    END_CODE();
  }

  // Generate the rational approximation x^(pnum/pden)
  Real RemezGMP::generateApprox(int degree, unsigned long pnum, unsigned long pden)
  {
    return generateApprox(degree, degree, pnum, pden);
  }

  // Generate the rational approximation x^(pnum/pden)
  Real RemezGMP::generateApprox(int num_degree, int den_degree, 
				unsigned long pnum, unsigned long pden)
  {
    START_CODE();

    // Reallocate arrays, since degree has changed
    if (num_degree != n || den_degree != d) allocate(num_degree,den_degree);

    step.resize(num_degree+den_degree+2);

    power_num = pnum;
    power_den = pden;
    spread = 1.0e37;
    iter = 0;

    n = num_degree;
    d = den_degree;
    neq = n + d + 1;

    initialGuess();
    stpini(step);

    while (spread > tolerance) 
    {
      //iterate until convergance

      if (iter++%100==0) 
	QDPIO::cout << __func__ << ": Iteration " << iter-1 
		    << " spread " << spread << " delta " << delta << endl;

      equations();

      if (delta < tolerance)
      {
	QDPIO::cerr << __func__ << ":Delta too small, try increasing precision" << endl;;
	QDP_abort(1);
      }

      search(step);

    }

    int sign;
    Real error = (Real)getErr(mm[0],sign);
    QDPIO::cout << __func__ << " Converged at " << iter << " iterations, error = " << error << endl;

    // Once the approximation has been generated, calculate the roots
    if (!root()) 
    {
      QDPIO::cerr << __func__ << ": Root finding failed" << endl;
      QDP_abort(1);
    }
  
    END_CODE();

    // Return the maximum error in the approximation
    return error;
  }

  // Return the partial fraction expansion of the approximation x^(pnum/pden)
  RemezCoeff_t RemezGMP::getPFE()
  {
    START_CODE();

    if (n!=d) 
    {
      QDPIO::cerr << __func__ << ": Cannot handle case: Numerator degree neq Denominator degree" << endl;
      QDP_abort(1);
    }

    if (!alloc) 
    {
      QDPIO::cerr << __func__ << ": Approximation not yet generated" << endl;
      QDP_abort(1);
    }

    RemezCoeff_t coeff;
    coeff.res.resize(n);
    coeff.pole.resize(d);

    multi1d<bigfloat> r(n);
    multi1d<bigfloat> p(d);
  
    for (int i=0; i<n; i++) r[i] = roots[i];
    for (int i=0; i<d; i++) p[i] = poles[i];
  
    // Perform a partial fraction expansion
    pfe(r, p, norm);

    // Convert to Real and return
    coeff.norm = (double)norm;
    for (int i=0; i<n; i++) coeff.res[i] = (double)r[i];
    for (int i=0; i<d; i++) coeff.pole[i] = (double)p[i];

    END_CODE();

    return coeff;
  }

  // Return the partial fraction expansion of the approximation x^(-pnum/pden)
  RemezCoeff_t RemezGMP::getIPFE()
  {
    START_CODE();

    if (!alloc) 
    {
      QDPIO::cerr << __func__ << ": Approximation not yet generated" << endl;
      QDP_abort(1);
    }

    RemezCoeff_t coeff;
    coeff.res.resize(n);
    coeff.pole.resize(d);

    multi1d<bigfloat> r(d);
    multi1d<bigfloat> p(n);
  
    // Want the inverse function
    for (int i=0; i<n; i++) {
      r[i] = poles[i];
      p[i] = roots[i];
    }

    // Perform a partial fraction expansion
    pfe(r, p, (bigfloat)1l/norm);

    // Convert to Real and return
    coeff.norm = (double)((bigfloat)1l/(norm));
    for (int i=0; i<n; i++) {
      coeff.res[i] = (double)r[i];
      coeff.pole[i] = (double)p[i];
    }

    END_CODE();

    return coeff;
  }

  // Initial values of maximal and minimal errors
  void RemezGMP::initialGuess() 
  {
    START_CODE();

    // Supply initial guesses for solution points
    long ncheb = neq;			// Degree of Chebyshev error estimate
    bigfloat a, r;

    // Find ncheb+1 extrema of Chebyshev polynomial

    a = ncheb;
    mm[0] = apstrt;
    for (long i = 1; i < ncheb; i++) {
      r = 0.5 * (1 - cos((M_PI * i)/(double) a));
      //r *= sqrt_bf(r);
      r = (exp((double)r)-1.0)/(exp(1.0)-1.0);
      mm[i] = apstrt + r * apwidt;
    }
    mm[ncheb] = apend;

    a = 2.0 * ncheb;
    for (long i = 0; i <= ncheb; i++) {
      r = 0.5 * (1 - cos(M_PI * (2*i+1)/(double) a));
      //r *= sqrt_bf(r); // Squeeze to low end of interval
      r = (exp((double)r)-1.0)/(exp(1.0)-1.0);
      xx[i] = apstrt + r * apwidt;
    }

    END_CODE();
  }

  // Initialise step sizes
  void RemezGMP::stpini(multi1d<bigfloat>& step) 
  {
    START_CODE();

    if (step.size() == 0)
      QDP_error_exit("%s: step not allocated", __func__);

    xx[neq+1] = apend;
    delta = 0.25;
    step[0] = xx[0] - apstrt;
    for (int i = 1; i < neq; i++) step[i] = xx[i] - xx[i-1];
    step[neq] = step[neq-1];

    END_CODE();
  }

  // Search for error maxima and minima
  void RemezGMP::search(multi1d<bigfloat>& step) 
  {
    START_CODE();

    if (step.size() == 0)
      QDP_error_exit("%s: step not allocated", __func__);

    bigfloat a, q, xm, ym, xn, yn, xx0, xx1;
    int i, j, meq, emsign, ensign, steps;

    meq = neq + 1;
    multi1d<bigfloat> yy(meq);

    bigfloat eclose = 1.0e30;
    bigfloat farther = 0l;

    j = 1;
    xx0 = apstrt;

    for (i = 0; i < meq; i++) {
      steps = 0;
      xx1 = xx[i]; // Next zero
      if (i == meq-1) xx1 = apend;
      xm = mm[i];
      ym = getErr(xm,emsign);
      q = step[i];
      xn = xm + q;
      if (xn < xx0 || xn >= xx1) {	// Cannot skip over adjacent boundaries
	q = -q;
	xn = xm;
	yn = ym;
	ensign = emsign;
      } else {
	yn = getErr(xn,ensign);
	if (yn < ym) {
	  q = -q;
	  xn = xm;
	  yn = ym;
	  ensign = emsign;
	}
      }
  
      while(yn >= ym) {		// March until error becomes smaller.
	if (++steps > 10) break;
	ym = yn;
	xm = xn;
	emsign = ensign;
	a = xm + q;
	if (a == xm || a <= xx0 || a >= xx1) break;// Must not skip over the zeros either side.
	xn = a;
	yn = getErr(xn,ensign);
      }

      mm[i] = xm;			// Position of maximum
      yy[i] = ym;			// Value of maximum

      if (eclose > ym) eclose = ym;
      if (farther < ym) farther = ym;

      xx0 = xx1; // Walk to next zero.
    } // end of search loop

    q = (farther - eclose);	// Decrease step size if error spread increased
    if (eclose != 0.0) q /= eclose; // Relative error spread
    if (q >= spread) delta *= 0.5; // Spread is increasing; decrease step size
    spread = q;

    for (i = 0; i < neq; i++) {
      q = yy[i+1];
      if (q != 0.0) q = yy[i] / q  - (bigfloat)1l;
      else q = 0.0625;
      if (q > (bigfloat)0.25) q = 0.25;
      q *= mm[i+1] - mm[i];
      step[i] = q * delta;
    }
    step[neq] = step[neq-1];
  
    for (i = 0; i < neq; i++) {	// Insert new locations for the zeros.
      xm = xx[i] - step[i];
      if (xm <= apstrt) continue;
      if (xm >= apend) continue;
      if (xm <= mm[i]) xm = (bigfloat)0.5 * (mm[i] + xx[i]);
      if (xm >= mm[i+1]) xm = (bigfloat)0.5 * (mm[i+1] + xx[i]);
      xx[i] = xm;
    }

    END_CODE();
  }

  // Solve the equations
  void RemezGMP::equations(void) 
  {
    START_CODE();

    bigfloat x, y, z;
    int i, j, ip;
    bigfloat *aa;

    multi1d<bigfloat> AA(neq*neq);
    multi1d<bigfloat> BB(neq);
  
    for (i = 0; i < neq; i++) {	// set up the equations for solution by simq()
      ip = neq * i;		// offset to 1st element of this row of matrix
      x = xx[i];			// the guess for this row
      y = func(x);		// right-hand-side vector

      z = (bigfloat)1l;
      aa = &AA[ip];
      for (j = 0; j <= n; j++) {
	*aa++ = z;
	z *= x;
      }

      z = (bigfloat)1l;
      for (j = 0; j < d; j++) {
	*aa++ = -y * z;
	z *= x;
      }
      BB[i] = y * z;		// Right hand side vector
    }

    // Solve the simultaneous linear equations.
    if (simq(AA, BB, param, neq))
    {
      QDPIO::cerr << __func__ << ": simq failed" << endl;
      QDP_abort(1);
    }

    END_CODE();
  }

  // Evaluate the rational form P(x)/Q(x) using coefficients
  // from the solution vector param
  bigfloat RemezGMP::approx(const bigfloat& x) 
  {
    START_CODE();

    bigfloat yn, yd;
    int i;

    // Work backwards toward the constant term.
    yn = param[n];		// Highest order numerator coefficient
    for (i = n-1; i >= 0; i--) yn = x * yn  +  param[i]; 
    yd = x + param[n+d];	// Highest degree coefficient = 1.0
    for (i = n+d-1; i > n; i--) yd = x * yd  +  param[i];

    END_CODE();

    return(yn/yd);
  }

  // Compute size and sign of the approximation error at x
  bigfloat RemezGMP::getErr(const bigfloat& x, int& sign) 
  {
    bigfloat e, f;

    f = func(x);
    e = approx(x) - f;
    if (f != 0) e /= f;
    if (e < (bigfloat)0.0) {
      sign = -1;
      e = -e;
    }
    else 
      sign = 1;
  
    return(e);
  }

  // Calculate function required for the approximation
  bigfloat RemezGMP::func(const bigfloat& x) 
  {
    START_CODE();

    bigfloat y,dy,f=1l,df;

    // initial guess to accelerate convergance
    y = (bigfloat)pow((double)x,(double)((bigfloat)1l/(bigfloat)power_den));
    while (abs_bf(f)>(bigfloat)1l/pow_bf((bigfloat)10,prec)) { // approx good to 10^(-prec)
      f = pow_bf(y,power_den) - x;
      df = (bigfloat)power_den*pow_bf(y,power_den-1);// need power_den-1 because of diff
      dy = f/df;
      y -= dy;
    }

    END_CODE();

    return pow_bf(y,power_num);
  }

  // Solve the system AX=B
  int RemezGMP::simq(multi1d<bigfloat>& A, multi1d<bigfloat>& B, multi1d<bigfloat>& X, int n) 
  {
    START_CODE();

    if (A.size() == 0 || B.size() == 0)
      QDP_error_exit("%s: A or B not allocated", __func__);

    int i, j, ij, ip, ipj, ipk, ipn;
    int idxpiv, iback;
    int k, kp, kp1, kpk, kpn;
    int nip, nkp, nm1;
    bigfloat em, q, rownrm, big, size, pivot, sum;
    bigfloat *aa;

    multi1d<int> IPS(neq);		// simq() work vector

    nm1 = n - 1;
    // Initialize IPS and X
  
    ij = 0;
    for (i = 0; i < n; i++) {
      IPS[i] = i;
      rownrm = 0.0;
      for(j = 0; j < n; j++) {
	q = abs_bf(A[ij]);
	if(rownrm < q) rownrm = q;
	++ij;
      }
      if (rownrm == (bigfloat)0l) 
      {
	QDPIO::cout << "simq rownrm=0\n" << endl;
	return(1);
      }
      X[i] = (bigfloat)1.0 / rownrm;
    }
  
    for (k = 0; k < nm1; k++) {
      big = 0.0;
      idxpiv = 0;
    
      for (i = k; i < n; i++) {
	ip = IPS[i];
	ipk = n*ip + k;
	size = abs_bf(A[ipk]) * X[ip];
	if (size > big) {
	  big = size;
	  idxpiv = i;
	}
      }
    
      if (big == (bigfloat)0l) 
      {
	QDPIO::cout << "simq big=0" << endl;
	return(2);
      }
      if (idxpiv != k) {
	j = IPS[k];
	IPS[k] = IPS[idxpiv];
	IPS[idxpiv] = j;
      }
      kp = IPS[k];
      kpk = n*kp + k;
      pivot = A[kpk];
      kp1 = k+1;
      for (i = kp1; i < n; i++) {
	ip = IPS[i];
	ipk = n*ip + k;
	em = -A[ipk] / pivot;
	A[ipk] = -em;
	nip = n*ip;
	nkp = n*kp;
	aa = &A[nkp+kp1];
	for (j = kp1; j < n; j++) {
	  ipj = nip + j;
	  A[ipj] = A[ipj] + em * *aa++;
	}
      }
    }
    kpn = n * IPS[n-1] + n - 1;	// last element of IPS[n] th row
    if (A[kpn] == (bigfloat)0l) 
    {
      QDPIO::cout << "simq A[kpn]=0" << endl;
      return(3);
    }

  
    ip = IPS[0];
    X[0] = B[ip];
    for (i = 1; i < n; i++) {
      ip = IPS[i];
      ipj = n * ip;
      sum = 0.0;
      for (j = 0; j < i; j++) {
	sum += A[ipj] * X[j];
	++ipj;
      }
      X[i] = B[ip] - sum;
    }
  
    ipn = n * IPS[n-1] + n - 1;
    X[n-1] = X[n-1] / A[ipn];
  
    for (iback = 1; iback < n; iback++) {
      //i goes (n-1),...,1
      i = nm1 - iback;
      ip = IPS[i];
      nip = n*ip;
      sum = 0.0;
      aa = &A[nip+i+1];
      for (j= i + 1; j < n; j++) 
	sum += *aa++ * X[j];
      X[i] = (X[i] - sum) / A[nip+i];
    }
  
    END_CODE();

    return(0);
  }

  // Calculate the roots of the approximation
  int RemezGMP::root() 
  {
    START_CODE();

    long i,j;
    bigfloat x,dx=0.05;
    bigfloat upper=1, lower=-100000;
    bigfloat tol = 1e-20;

    multi1d<bigfloat> poly(neq+1);
    
    // First find the numerator roots
    for (i=0; i<=n; i++) poly[i] = param[i];
    for (i=n-1; i>=0; i--) {
      roots[i] = rtnewt(poly,i+1,lower,upper,tol);
      if (roots[i] == 0.0) 
      {
	QDPIO::cerr << __func__ << ": Failure to converge on root" << endl;
	return 0;
      }
      poly[0] = -poly[0]/roots[i];
      for (j=1; j<=i; j++) poly[j] = (poly[j-1] - poly[j])/roots[i];
    }
  
    // Now find the denominator roots
    poly[d] = 1l;
    for (i=0; i<d; i++) poly[i] = param[n+1+i];
    for (i=d-1; i>=0; i--) {
      poles[i]=rtnewt(poly,i+1,lower,upper,tol);
      if (poles[i] == 0.0) {
	QDPIO::cerr << __func__ << ": Failure to converge on root" << endl;
	return 0;
      }
      poly[0] = -poly[0]/poles[i];
      for (j=1; j<=i; j++) poly[j] = (poly[j-1] - poly[j])/poles[i];
    }

    norm = param[n];
    QDPIO::cout << __func__ << ": Normalisation constant is " << (double)norm << endl;
    for (i=0; i<n; i++) 
      QDPIO::cout << "root[" << i << "] = " << (double)roots[i] << endl;
    for (i=0; i<d; i++) 
      QDPIO::cout << "pole[" << i << "] = " << (double)poles[i] << endl;

    END_CODE();

    return 1;
  }

  // Evaluate the polynomial
  bigfloat RemezGMP::polyEval(const bigfloat& x, const multi1d<bigfloat>& poly, long size) 
  {
    bigfloat f = poly[size];
    for (int i=size-1; i>=0; i--) f = f*x + poly[i];
    return f;
  }

  // Evaluate the differential of the polynomial
  bigfloat RemezGMP::polyDiff(const bigfloat& x, const multi1d<bigfloat>& poly, long size) 
  {
    bigfloat df = (bigfloat)size*poly[size];
    for (int i=size-1; i>0; i--) 
      df = df*x + (bigfloat)i*poly[i];
    return df;
  }

  // Newton's method to calculate roots
  bigfloat RemezGMP::rtnewt(const multi1d<bigfloat>& poly, long i, 
			    const bigfloat& x1, const bigfloat& x2, const bigfloat& xacc) 
  {
    START_CODE();

    int j;
    bigfloat df, dx, f, rtn;
  
    rtn=(bigfloat)0.5*(x1+x2);
    for (j=1; j<=JMAX;j++) 
    {
      f = polyEval(rtn, poly, i);
      df = polyDiff(rtn, poly, i);
      dx = f/df;
      rtn -= dx;
      if ((x1-rtn)*(rtn-x2) < (bigfloat)0.0)
	QDPIO::cerr << __func__ << ": Jumped out of brackets in rtnewt" << endl;
      if (abs_bf(dx) < xacc) return rtn;
    }
    QDPIO::cerr << __func__ << ": Maximum number of iterations exceeded in rtnewt" << endl;

    END_CODE();
    return 0.0;
  }

  // Evaluate the partial fraction expansion of the rational function
  // with res roots and poles poles.  Result is overwritten on input
  // arrays.
  void RemezGMP::pfe(multi1d<bigfloat>& res, multi1d<bigfloat>& poles, const bigfloat& norm) 
  {
    START_CODE();

//  QDPIO::cout << __func__ << " : enter" << endl;

    if (res.size() == 0 || poles.size() == 0)
      QDP_error_exit("%s: res or poles not allocated", __func__);

    int i,j,small;
    bigfloat temp;
    multi1d<bigfloat> numerator(n);
    multi1d<bigfloat> denominator(d);

    // Construct the polynomials explicitly 
    for (i=1; i<n; i++) {
      numerator[i] = 0l;
      denominator[i] = 0l;
    }
    numerator[0]=1l;
    denominator[0]=1l;

    for (j=0; j<n; j++) {
      for (i=n-1; i>=0; i--) {
	numerator[i] *= -res[j];
	denominator[i] *= -poles[j];
	if (i>0) {
	  numerator[i] += numerator[i-1];
	  denominator[i] += denominator[i-1];
	}
      }
    }

    // Convert to proper fraction form.
    // Fraction is now in the form 1 + n/d, where O(n)+1=O(d)
    for (i=0; i<n; i++) numerator[i] -= denominator[i];

    // Find the residues of the partial fraction expansion and absorb the
    // coefficients.
    for (i=0; i<n; i++) {
      res[i] = 0l;
      for (j=n-1; j>=0; j--) {
	res[i] = poles[i]*res[i]+numerator[j];
      }
      for (j=n-1; j>=0; j--) {
	if (i!=j) res[i] /= poles[i]-poles[j];
      }
      res[i] *= norm;
    }  

    // res now holds the residues
    j = 0;
    for (i=0; i<n; i++) poles[i] = -poles[i];

    // Move the ordering of the poles from smallest to largest
    for (j=0; j<n; j++) 
    {
      small = j;
      for (i=j+1; i<n; i++) {
	if (poles[i] < poles[small]) small = i;
      }
      if (small != j) {
	temp = poles[small];
	poles[small] = poles[j];
	poles[j] = temp;
	temp = res[small];
	res[small] = res[j];
	res[j] = temp;
      }
      QDPIO::cout << __func__ << ": Residue = " << (double)res[j] << " Pole = " << (double)poles[j] << endl;
    }

//  QDPIO::cout << __func__ << " : exit" << endl;

    END_CODE();
  }


  // Given a partial fraction expansion, evaluate it at x
  Real RemezGMP::evalPFE(const Real& x, const RemezCoeff_t& coeff)
  {
    if ((coeff.res.size() != coeff.pole.size()) || coeff.res.size() == 0)
    {
      QDPIO::cerr << __func__ << ": invalid res and pole" << endl;
      QDP_abort(1);
    }

    Real f = coeff.norm;
    for (int i=0; i < coeff.res.size(); ++i)
      f += coeff.res[i] / (x + coeff.pole[i]);

    return f;
  }

}  // end namespace Chroma

