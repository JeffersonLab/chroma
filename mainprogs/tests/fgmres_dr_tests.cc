#include "gtest/gtest.h"
#include "chromabase.h"
#include "qdp-lapack.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/invert/syssolver_fgmres_dr_params.h"
#include "actions/ferm/invert/syssolver_linop_fgmres_dr.h"
//using namespace std;
using namespace Chroma;

  const std::string xml_for_param = 
    "<?xml version='1.0'?> \
   <Params>					      \
     <FermionAction>				      \
        <FermAct>CLOVER</FermAct>		      \
        <Kappa>0.115</Kappa>			      \
        <clovCoeff>1.17</clovCoeff>		      \
        <clovCoeffR>0.91</clovCoeffR>		      \
        <clovCoeffT>1.07</clovCoeffT>		      \
        <AnisoParam>				      \
          <anisoP>true</anisoP>			      \
          <t_dir>3</t_dir>			      \
          <xi_0>2.464</xi_0>			      \
          <nu>0.95</nu>				      \
        </AnisoParam>				      \
        <FermState>				      \
          <Name>STOUT_FERM_STATE</Name>		      \
          <rho>0.22</rho>			      \
          <n_smear>2</n_smear>			      \
          <orthog_dir>3</orthog_dir>		      \
          <FermionBC>				      \
            <FermBC>SIMPLE_FERMBC</FermBC>	      \
            <boundary>1 1 1 -1</boundary>	      \
          </FermionBC>				      \
        </FermState>				      \
       </FermionAction>				      \
     <InvertParam>				      \
     <invType>FGMRESDR_INVERTER</invType>	      \
     <RsdTarget>1.0e-7</RsdTarget>		      \
     <NKrylov>5</NKrylov>			      \
     <NDefl>3</NDefl>				      \
     <MaxIter>130</MaxIter>			      \
     <PrecondParams>				      \
       <invType>MR_INVERTER</invType>		      \
       <RsdMR>0.1</RsdMR>			      \
       <MaxMR>10</MaxMR>			      \
     </PrecondParams>				      \
   </InvertParam>				      \
  </Params>";


class FGMRESDRTests : public ::testing::Test {
public:

  // Type aliases should be visible to all tests
  using T = LatticeFermion; 
  using Q = multi1d<LatticeColorMatrix>;
  using P = multi1d<LatticeColorMatrix>;


  FGMRESDRTests()
  {
    u.resize(Nd);
    for(int mu=0; mu < Nd; ++mu) {
      gaussian(u[mu]);
      reunit(u[mu]);
    }

    std::istringstream input(xml_for_param);
    XMLReader xml_in(input);
    
    S_f = dynamic_cast<FermAct4D<T,P,Q>*>(TheFermionActionFactory::Instance().createObject("CLOVER", 
											   xml_in, 
											   "FermionAction")
					  );
    state = S_f->createState(u);
    linop = S_f->linOp(state);

  }

  void TearDown()
  {
   
  }

  // Virtual destructor
  virtual
  ~FGMRESDRTests() {}


  multi1d<LatticeColorMatrix> u;
  Handle< FermAct4D<T,P,Q> > S_f;
  Handle< FermState<T,P,Q> > state;
  Handle< LinearOperator<T> > linop;
};

class FGMRESDRTestsFloatParams : public FGMRESDRTests,
				 public ::testing::WithParamInterface<float>
{

};

TEST_F(FGMRESDRTests, canConstructDefault)
{
  SysSolverFGMRESDRParams p;
  ASSERT_EQ(p.NKrylov, 0);
  ASSERT_EQ(p.NDefl, 0);
  ASSERT_EQ(p.MaxIter, 0);
  ASSERT_EQ(p.PrecondParams.path, "/PrecondParams");
}


TEST_F(FGMRESDRTests, canReadXML)
{
  std::istringstream input(xml_for_param);
  XMLReader xml_in(input);
  SysSolverFGMRESDRParams p( xml_in, "/Params/InvertParam" );
  ASSERT_EQ(p.NKrylov, 5);
  ASSERT_EQ(p.NDefl, 3);
  ASSERT_EQ(p.MaxIter, 130);
  ASSERT_EQ(p.PrecondParams.path, "/PrecondParams");
  ASSERT_EQ(p.PrecondParams.id, "MR_INVERTER");
}

TEST_F(FGMRESDRTests, canCreateFGMRESDRClassFromParam)
{
  std::istringstream input(xml_for_param);
  XMLReader xml_in(input);
  SysSolverFGMRESDRParams p( xml_in, "/Params/InvertParam" );
  LinOpSysSolverFGMRESDR fgmres_dr_solver(linop,state,p);
}

TEST_F(FGMRESDRTests, canCreateFGMRESDRClassFromFactory)
{
  std::istringstream input(xml_for_param);
  XMLReader xml_in(input);
  Handle< LinOpSystemSolver<T> > solver_handle =  TheLinOpFermSystemSolverFactory::Instance().createObject( "FGMRESDR_INVERTER", xml_in, std::string("/Params/InvertParam"), state, linop);
}

TEST_F(FGMRESDRTests, multi2d)
{
  int n_krylov = 5;
  int rows=n_krylov+1;
  int cols=n_krylov;
  multi2d<DComplex> H(cols,rows);

  for(int col=0; col < n_krylov; ++col) {
    for(int row=0; row < n_krylov+1; ++row)  {
      H(col,row) = DComplex(0);
    }
  }
  
  int dim = n_krylov;
  for(int row=0; row < dim+1; ++row)  {
    QDPIO::cout << "[ ";
    for(int col=0; col < dim; ++col) {
      QDPIO::cout << H(col, row) << " ";
    }
    QDPIO::cout << " ] " << std::endl;
  }
}

TEST_F(FGMRESDRTests, arnoldi5)
{
  std::istringstream input(xml_for_param);
  XMLReader xml_in(input);
  SysSolverFGMRESDRParams p( xml_in, "/Params/InvertParam" );
  LinOpSysSolverFGMRESDR sol(linop,state,p);


  const Subset& s = sol.subset();


    // Create a gaussian source 
  LatticeFermion rhs;
  gaussian(rhs, s);

  
  Real rsd_target(1.0e-6);
  rsd_target*=sqrt(norm2(rhs,s));
  int n_krylov=5;
  int n_deflate=0;
  multi2d<DComplex> H(n_krylov,n_krylov+1); // The H matrix
  multi2d<DComplex> R(n_krylov,n_krylov+1); // R = H diagonalized with Givens rotations
  multi1d<T> V(n_krylov+1);  // K(A)
  multi1d<T> Z(n_krylov+1);  // K(MA)
  multi1d< Handle<Givens> > givens_rots(n_krylov+1);
  multi1d<DComplex> g(n_krylov+1);
  multi2d<DComplex> Qk_Hk;
  multi1d<DComplex> Qk_Hk_taus;

  for(int col=0; col < n_krylov; ++col) {
    for(int row=0; row < n_krylov+1; ++row)  {
      H(col,row) = DComplex(0);
      R(col,row) = DComplex(0);
    }
  }
  

  Double beta=sqrt(norm2(rhs, s));
  for(int j=0; j < g.size(); ++j) { 
    g[j] = DComplex(0); 
  }
  g[0] = beta;
  // Set up initial V[0] = rhs / || r^2 || -- 
  Double beta_inv = Double(1)/beta;
  V[0][s] = beta_inv * rhs;
  int dim;

  sol.FlexibleArnoldi(n_krylov, n_deflate,
		      rsd_target,
		      V,
		      Z, 
		      H, 
		      R,
		      givens_rots,
		      g,
		      Qk_Hk,
		      Qk_Hk_taus,
		      dim);

  ASSERT_EQ(dim,n_krylov);

  // Instead of printing, Assert that Upper Hessenberg
  for(int row=0; row < dim+1; ++row) {
    for(int col = 0; col < row-1; ++col) { 

      // These are below the row below the diagonal -- should be 0

      double re = toDouble( real( H(col,row) ) );
      double im = toDouble( imag( H(col,row) ) );
      EXPECT_DOUBLE_EQ( re, 0);
      EXPECT_DOUBLE_EQ( im, 0);
    }

  }

    
  
  // Now the various V's should be orthonormal. Check this.
  for(int j=0; j < dim; ++j) { 
    EXPECT_DOUBLE_EQ( toDouble(sqrt(norm2(V[j],s))), 1.0);
  }

  for(int j=0; j < dim; ++j) { 
    for(int i=j+1; i < dim; ++i) {
      ASSERT_NEAR( 0, toDouble(norm2(innerProduct(V[i],V[j],s))),1.0e-16);
    }
  }
}

TEST_F(FGMRESDRTests, GivensF0G0)
{
  multi2d<DComplex> H(1,2);
  // Fill both 'f' and 'g' with 0
  for(int row=0; row < 2; ++row) {
    H(0,row) = DComplex(0);
  }
  Givens g(0,H);
  
  // Apply Givens to col 0. 
  g(0,H);


  // H(1,0) should be 0
  EXPECT_DOUBLE_EQ( 0, toDouble( real(H(0,1)) ) );
  EXPECT_DOUBLE_EQ( 0, toDouble( imag(H(0,1)) ) );


  EXPECT_DOUBLE_EQ( 0, toDouble( real(H(0,0)) ) );
  EXPECT_DOUBLE_EQ( 0, toDouble( imag(H(0,0)) ) );
}

TEST_F(FGMRESDRTests, GivensF0GNot0)
{
  multi2d<DComplex> H(1,2);
  // Fill both 'f' and 'g' with 0
  for(int row=0; row < 2; ++row) {
    H(0,row) = DComplex(0);
  }

  DComplex g;
  random(g);


  H(0,1) = g;  // Nonzero 'g'

  // Compute Givens
  Givens G(0,H);
  
  // Apply Givens to col 0. 
  G(0,H);
  // H(1,0) should be 0
  EXPECT_DOUBLE_EQ( 0, toDouble( real(H(0,1)) ) );
  EXPECT_DOUBLE_EQ( 0, toDouble( imag(H(0,1)) ) );

  EXPECT_DOUBLE_EQ( toDouble(real(g)), toDouble( real(H(0,0)) ) );
  EXPECT_DOUBLE_EQ( toDouble(imag(g)), toDouble( imag(H(0,0)) ) );
}

TEST_F(FGMRESDRTests, GivensFNot0G0)
{
  multi2d<DComplex> H(1,2);
  // Fill both 'f' and 'g' with 0
  for(int row=0; row < 2; ++row) {
    H(0,row) = DComplex(0);
  }

  DComplex f;
  random(f);


  H(0,0) = f;  // Nonzero 'f'

  // Compute Givens
  Givens G(0,H);
  
  // Apply Givens to col 0. 
  G(0,H);
  // H(1,0) should be 0
  EXPECT_DOUBLE_EQ( 0, toDouble( real(H(0,1)) ) );
  EXPECT_DOUBLE_EQ( 0, toDouble( imag(H(0,1)) ) );

  EXPECT_DOUBLE_EQ( toDouble(real(sqrt(norm2(f)))), toDouble( real(H(0,0)) ) );
  EXPECT_DOUBLE_EQ( toDouble(imag(sqrt(norm2(f)))), toDouble( imag(H(0,0)) ) );
}

TEST_F(FGMRESDRTests, GivensFNot0GNot0)
{
  multi2d<DComplex> H(2,2);
  // Fill both 'f' and 'g' with 0
  for(int row=0; row < 2; ++row) {
    H(0,row) = DComplex(0);
  }

  DComplex f,g;
  random(f);
  random(g);


  H(0,0) = f;  // Nonzero f& g
  H(0,1) = g;

  H(1,0) = f;
  H(1,1) = g;

  // Compute Givens
  Givens G(0,H); // Col 0 -- same as used to compute (get r automatically)

  G(0,H);  // Col 0 -- same as used to compute (get r automatically)
  G(1,H);  // Col 1 -- same as used to compute (get r automatically)

  // H(0,1) should be 0
  EXPECT_DOUBLE_EQ( 0, toDouble( real(H(0,1)) ) );
  EXPECT_DOUBLE_EQ( 0, toDouble( imag(H(0,1)) ) );

  // H(1,1) should be 0
  EXPECT_DOUBLE_EQ( 0, toDouble( real(H(1,1)) ) );
  EXPECT_DOUBLE_EQ( 0, toDouble( imag(H(1,1)) ) );

  // H(0,0) and H(1,0) should be the same
  EXPECT_DOUBLE_EQ( toDouble( real(H(0,0)) ) , toDouble( real(H(1,0)) ) );
  //EXPECT_DOUBLE_EQ( toDouble( imag(H(0,0)) ) , toDouble( imag(H(1,0)) ) );


 
}
		


// TODO:!!!!
// Get rid of original arnoldi test.
// Fold its elements into this test
// Parameterize the tests on nkrylov, and rsd_target	   
// So we can engineer situations where we converge in the cycle, 
// and ones where we do not.

TEST_P(FGMRESDRTestsFloatParams, arnoldi5GivensRot)
{
  std::istringstream input(xml_for_param);
  XMLReader xml_in(input);
  SysSolverFGMRESDRParams p( xml_in, "/Params/InvertParam" );
  LinOpSysSolverFGMRESDR sol(linop,state,p);


  const Subset& s = sol.subset();


    // Create a gaussian source 
  LatticeFermion rhs;
  gaussian(rhs, s);

  double rsd_target_in=GetParam();
  Real rsd_target(rsd_target_in);
  rsd_target *= sqrt(norm2(rhs,s));

  int n_krylov=5;
  int n_deflate=0;

  // Work spaces
  multi2d<DComplex> H(n_krylov, n_krylov+1); // The H matrix
  multi2d<DComplex> R(n_krylov, n_krylov+1); // R = H diagonalized with Givens rotations
  for(int col = 0; col < n_krylov; ++col) { 
    for(int row = 0; row < n_krylov+1; ++row) {
      
      H(col,row) = DComplex(0);
      R(col,row) = DComplex(0);
    }
  }
  multi1d<T> V(n_krylov+1);  // K(A)
  multi1d<T> Z(n_krylov+1);  // K(MA)
  multi1d< Handle<Givens> > givens_rots(n_krylov+1);
  multi1d<DComplex> g(n_krylov+1);
  multi2d<DComplex> Qk;
  multi1d<DComplex> Qk_tau;

  // Assume zero initial guess
  Double beta=sqrt(norm2(rhs, s));

  for(int j=0; j < g.size(); ++j) { 
    g[j] = DComplex(0); 
  }
  g[0] = beta;
  // Set up initial V[0] = rhs / || r^2 || -- 
  Double beta_inv = Double(1)/beta;
  V[0][s] = beta_inv * rhs;

  int dim;
  sol.FlexibleArnoldi(n_krylov, 
		      n_deflate,
		      rsd_target,
		      V,
		      Z, 
		      H, 
		      R,
		      givens_rots,
		      g,
		      Qk,
		      Qk_tau,
		      dim);
  
  ASSERT_LE(dim,n_krylov);



  // Instead of printing, Assert that Upper Hessenberg
  for(int row=0; row < dim+1; ++row) {
    // These are below the row below the diagonal -- should be 0
    for(int col = 0; col < row-1; ++col) { 
      double re = toDouble( real( H(col,row) ) );
      double im = toDouble( imag( H(col,row) ) );
      EXPECT_DOUBLE_EQ( re, 0);
      EXPECT_DOUBLE_EQ( im, 0);
    }

  }

    
  
  // Now the various V's should be orthonormal. Check this.
  for(int j=0; j < dim+1; ++j) { 
    EXPECT_DOUBLE_EQ( toDouble(sqrt(norm2(V[j],s))), 1.0);
  }

  for(int j=0; j < dim+1; ++j) { 
    for(int i=j+1; i < dim+1; ++i) {
      ASSERT_NEAR( 0, toDouble(norm2(innerProduct(V[i],V[j],s))),1.0e-16);
    }
  }

  multi1d< Handle<Givens> > GivensRots(dim+1);

  multi1d< DComplex > resvec(dim+1);
  for(int row=0; row < dim+1; ++row) {
    resvec[row] = DComplex(0);
  }
  resvec[0] = beta;

  // For each column
  for(int j=0; j < dim; ++j) {

    // For column 0 this would apply 0 rotations
    // For column 1 this would apply that of column 0 etc.
    // For column 2 this would apply that of col 0 & 1
    for(int i=0; i < j; ++i) {
      (*GivensRots[i])(j, H);
    }
    // Make a givens rotation for my column
    GivensRots[j] = new Givens(j, H);
    (*GivensRots[j])(j,H);
    (*GivensRots[j])(resvec); // Apply also to g
  }

  // At this point H should be Upper Triangular
  for(int row=1; row < dim+1; ++row) {
    for(int col=0; col < row; ++col)   {
      EXPECT_DOUBLE_EQ( 0, toDouble(real(H(col,row))));
      EXPECT_DOUBLE_EQ( 0, toDouble(imag(H(col,row))));
    }
  }

#if 0 
  // R is not computed anymore 
  // The transformed H should agree with R
  for(int row=0; row < dim+1; ++row) {
    for(int col=row; col < dim; ++col)   {
      EXPECT_DOUBLE_EQ( toDouble(real(R(col,row))), toDouble(real(H(col,row))));
      EXPECT_DOUBLE_EQ( toDouble(imag(R(col,row))), toDouble(imag(H(col,row))));
    }
  }
#endif

  // The resvec we generated should agree with what
  for(int col=0; col < dim+1; ++col) { 
    QDPIO::cout << "col = " << col << " g=" <<g[col] << " resvec" << resvec[col] << std::endl;
    EXPECT_DOUBLE_EQ(toDouble( real(g[col]) ), toDouble( real(resvec[col])) );
    EXPECT_DOUBLE_EQ(toDouble( imag(g[col]) ), toDouble( imag(resvec[col])) );

  }


}





TEST_F(FGMRESDRTests, testLeastSquares)
{
  std::istringstream input(xml_for_param);
  XMLReader xml_in(input);
  SysSolverFGMRESDRParams p( xml_in, "/Params/InvertParam" );
  LinOpSysSolverFGMRESDR sol(linop,state,p);

  int n_rows=4;
  int n_cols=3;

  multi2d<DComplex> R(n_cols,n_rows);
  multi1d<DComplex> b(n_rows);
  multi1d<DComplex> x(n_rows);
  for(int row=0; row < n_rows; ++row) { 
    for(int col=0; col < n_cols; ++col) { 
      R(col,row)=DComplex(0);
    }
    random(b[row]);
    x[row] = zero;
  }

  // Fill upper triangular with Random junk
  for(int row=0; row < n_rows-1; ++row) { 
    for(int col=row; col < n_cols;  ++col) { 
      random(R(col,row));
    }
  }

  sol.LeastSquaresSolve(R,b,x,n_cols);


  for(int row=0; row < n_cols; ++row) {
    DComplex prod=DComplex(0);
    for(int col=row; col < n_cols; ++col) { 
      prod += R(col,row)*x[col];
    }
    EXPECT_DOUBLE_EQ( toDouble(real(prod)), toDouble(real(b[row])));
    EXPECT_DOUBLE_EQ( toDouble(imag(prod)), toDouble(imag(b[row])));
  }
}


TEST_P(FGMRESDRTestsFloatParams, testFullSolverNoDeflate)
{
  std::istringstream input(xml_for_param);
  XMLReader xml_in(input);
  SysSolverFGMRESDRParams p( xml_in, "/Params/InvertParam" );
  p.NDefl = 0;
  // Reset target residuum as per input
  float rsd_target_in = GetParam();
  p.RsdTarget = Real(rsd_target_in);

  // Construct the solver 
  LinOpSysSolverFGMRESDR sol(linop,state,p);

  // Grab the subset
  const Subset& s =  linop->subset();

  // Gaussian RHS
  LatticeFermion rhs; 
  gaussian(rhs,s);

  // zero initial guess
  LatticeFermion x = zero;
  
  // Solve system
  (sol)(x,rhs);

  // Compute residuum

  LatticeFermion r= zero;
  (*linop)(r,x,PLUS);   // r = Ax
  r[s]-=rhs;        // r = Ax - b = -(b-Ax);

  Double resid_rel = sqrt( norm2(r,s)/norm2(rhs,s) ); // sqrt( || r ||^2/||b||^2 )
  ASSERT_LE( toDouble(resid_rel), toDouble(p.RsdTarget) );
  
}
  

TEST_P(FGMRESDRTestsFloatParams, testFullSolverDeflate)
{
  std::istringstream input(xml_for_param);
  XMLReader xml_in(input);
  SysSolverFGMRESDRParams p( xml_in, "/Params/InvertParam" );
  p.NDefl = 3;
  p.NKrylov = 6;
  // Reset target residuum as per input
  float rsd_target_in = GetParam();
  p.RsdTarget = Real(rsd_target_in);

  // Construct the solver 
  LinOpSysSolverFGMRESDR sol(linop,state,p);

  // Grab the subset
  const Subset& s =  linop->subset();

  // Gaussian RHS
  LatticeFermion rhs; 
  gaussian(rhs,s);

  // zero initial guess
  LatticeFermion x = zero;
  
  // Solve system
  (sol)(x,rhs);

  // Compute residuum

  LatticeFermion r= zero;
  (*linop)(r,x,PLUS);   // r = Ax
  r[s]-=rhs;        // r = Ax - b = -(b-Ax);

  Double resid_rel = sqrt( norm2(r,s)/norm2(rhs,s) ); // sqrt( || r ||^2/||b||^2 )
  ASSERT_LE( toDouble(resid_rel), toDouble(p.RsdTarget) );
  
}
  



INSTANTIATE_TEST_CASE_P(FGMRESDRTests,
			FGMRESDRTestsFloatParams,
			testing::Values(1.0e-3,1.0e-9));


TEST_F(FGMRESDRTests, testQDPLapackZGETRFZGETRS)
{
  std::istringstream input(xml_for_param);
  XMLReader xml_in(input);
  SysSolverFGMRESDRParams p( xml_in, "/Params/InvertParam" );
  LinOpSysSolverFGMRESDR sol(linop,state,p);


  const Subset& s = sol.subset();


    // Create a gaussian source 
  LatticeFermion rhs;
  gaussian(rhs, s);

  Real rsd_target(1.0e-9);
  rsd_target *= sqrt(norm2(rhs,s));

  int n_krylov=5;
  int n_deflate=0;

  // Work spaces
  multi2d<DComplex> H(n_krylov,n_krylov+1); // The H matrix
  multi2d<DComplex> R(n_krylov,n_krylov+1); // R = H diagonalized with Givens rotations
  for(int row = 0; row < n_krylov+1; ++row) {
    for(int col = 0; col < n_krylov; ++col) { 
      H(col,row) = DComplex(0);
      R(col,row) = DComplex(0);
    }
  }
  multi1d<T> V(n_krylov+1);  // K(A)
  multi1d<T> Z(n_krylov+1);  // K(MA)
  multi1d< Handle<Givens> > givens_rots(n_krylov+1);
  multi1d<DComplex> g(n_krylov+1);
  multi2d<DComplex> Qk;
  multi1d<DComplex> Qk_tau;
  // Assume zero initial guess
  Double beta=sqrt(norm2(rhs, s));

  for(int j=0; j < g.size(); ++j) { 
    g[j] = DComplex(0); 
  }
  g[0] = beta;
  // Set up initial V[0] = rhs / || r^2 || -- 
  Double beta_inv = Double(1)/beta;
  V[0][s] = beta_inv * rhs;

  int dim;
  sol.FlexibleArnoldi(n_krylov, 
		      n_deflate,
		      rsd_target,
		      V,
		      Z, 
		      H, 
		      R,
		      givens_rots,
		      g,
		      Qk,
		      Qk_tau,
		      dim);
  
  // NOw want to solve System H^\dagger f_m = h_m
  multi1d<DComplex> f_m(n_krylov);
  multi1d<DComplex> h_m(n_krylov);
  multi2d<DComplex> h1(n_krylov, n_krylov);

  // h1 is the copy I will pass to LAPACK 
  // I actually want to solve with the dagger of H
  // so I should transpose as well as conjugate
  // however LAPACK expects Fortran order: rows run fastest, columns slower
  // multi2d has Fortran like indexing (row first, column second)
  // but underneath the hood, columns still run fastest.
  // So in principle we need to transpose in and out of Fortran
  // however as I want to solve with the Transpose conjugate
  // all I have to do is leave the indices untransposed, and just conjugate
  for(int col=0; col < n_krylov; ++col) {
    for(int row =0; row < n_krylov; ++row) {
      h1(col,row) = H(col,row);
    }
  }

  for(int col=0; col < n_krylov; ++col) {
    h_m[col] = conj(H(col, n_krylov)); // Keeping this for checking
    f_m[col] = h_m[col]; // Solution will overwrite f_m
  }

  multi1d<int> ipiv(n_krylov);
  int info;
  
  // Factor H
  QDPIO::cout << "CALLING ZGETRF to factor H" << std::endl;
  QDPLapack::zgetrf(n_krylov,n_krylov,h1,n_krylov,ipiv,info);
  ASSERT_EQ(info,0);

  // Now solve using ZGETRS
  QDPIO::cout << "CALLING ZGETRS to solve for f" << std::endl;
  char trans='C'; // Decompose the Hermitian Conjugate
 
  QDPLapack::zgetrs(trans, n_krylov,1,h1,n_krylov,ipiv,f_m,n_krylov,info);
  ASSERT_EQ(info,0);
  
  // Check Solution
  multi1d<DComplex> diff(n_krylov);
  
  for(int row=0; row < n_krylov; ++row) { 
    diff[row] = zero;
  }

  // Check by multiplying back with H^H
  trans='C'; // Hermitian conjugate H, Leading Dimension of H is n_krylov+1
             // since now we are storing column major order. (Rows run fastest) 
  QDPLapack::zgemv(trans,n_krylov, n_krylov, H, n_krylov+1, f_m, diff);
  for(int row=0; row < n_krylov; ++row) { 
    diff[row] -= h_m[row];
    EXPECT_LT( toDouble(real(diff[row])), 1.0e-14 );
    EXPECT_LT( toDouble(imag(diff[row])), 1.0e-14 );
  }

  for(int row=0; row < n_krylov;++row) { 
    QDPIO::cout << "diff["<<row<<"] = " << diff[row] << std::endl;
  }

  // Now form Hm + f_m h_m^H
  // We will send this to Fortran, and we don't want to transpose
  // or conjugate. So just transpose indices on H_tilde now, to ensure
  // correct column major order for Fortran
  multi2d<DComplex> H_tilde(n_krylov,n_krylov);
  for(int row=0; row < n_krylov;++row) { 
    for(int col=0; col< n_krylov; ++col) { 
      H_tilde(col,row) = H(col,row) + f_m[row]*conj(h_m[col]);
    }
  }

  // Space for eigenvectors and eivenvalues
  multi1d<DComplex> evals(n_krylov);  // Evalues come back in this
  multi2d<DComplex> evecs(n_krylov, n_krylov); // Evecs come back in this in Fortran order

  QDPIO::cout << "CALLING ZGEEV" << std::endl;
  QDPLapack::zgeev(n_krylov, H_tilde, evals, evecs);
  QDPIO::cout << "Checking evecs" << std::endl;

  // Recreate H_tilde as the evec routine overwrites it
  for(int row=0; row < n_krylov;++row) { 
    for(int col=0; col< n_krylov; ++col) { 
      H_tilde(col,row) = H(col,row) + f_m[row]*conj(h_m[col]);
    }
  }

  // For each evec check   H_tilde v - lambda v is small
  //
  for(int evec=0; evec < n_krylov; ++evec) { 
  
    // Extract evec from 'evecs'
    multi1d<DComplex> my_evec(n_krylov);
    for(int row=0; row<n_krylov; ++row) {
      my_evec[row] = evecs(evec,row);
    }

    // Multiply: H_tilde my_evec -> diff
    trans='N'; 
    QDPLapack::zgemv(trans, n_krylov, n_krylov, H_tilde, n_krylov, my_evec, diff);
    
    // Subtract lambda my_evec
    for(int row=0; row < n_krylov; ++row) { 
      diff[row] -= evals[evec]*evecs(evec,row);

      // check result is small
      EXPECT_LT( toDouble(real(diff[row])), 1.0e-14 );
      EXPECT_LT( toDouble(imag(diff[row])), 1.0e-14 );
 
      QDPIO::cout << "eigenvec: " << evec << " diff[" << row <<"]=" << diff[row] << "  lambda=" << evals[evec] << std::endl;
    }
  }
}

TEST_F(FGMRESDRTests, testGetEigenvector)
{
  std::istringstream input(xml_for_param);
  XMLReader xml_in(input);
  SysSolverFGMRESDRParams p( xml_in, "/Params/InvertParam" );
  LinOpSysSolverFGMRESDR sol(linop,state,p);


  const Subset& s = sol.subset();


    // Create a gaussian source 
  LatticeFermion rhs;
  gaussian(rhs, s);

  Real rsd_target(1.0e-9);
  rsd_target *= sqrt(norm2(rhs,s));

  int n_krylov=5;
  int n_deflate=0;

  // Work spaces
  multi2d<DComplex> H(n_krylov,n_krylov+1); // The H matrix
  multi2d<DComplex> R(n_krylov,n_krylov+1); // R = H diagonalized with Givens rotations
  for(int row = 0; row < n_krylov+1; ++row) {
    for(int col = 0; col < n_krylov; ++col) { 
      H(col,row) = DComplex(0);
      R(col,row) = DComplex(0);
    }
  }
  multi1d<T> V(n_krylov+1);  // K(A)
  multi1d<T> Z(n_krylov+1);  // K(MA)
  multi1d< Handle<Givens> > givens_rots(n_krylov+1);
  multi1d<DComplex> g(n_krylov+1);
  multi2d<DComplex> Qk;
  multi1d<DComplex> Qk_tau;
  // Assume zero initial guess
  Double beta=sqrt(norm2(rhs, s));

  for(int j=0; j < g.size(); ++j) { 
    g[j] = DComplex(0); 
  }
  g[0] = beta;
  // Set up initial V[0] = rhs / || r^2 || -- 
  Double beta_inv = Double(1)/beta;
  V[0][s] = beta_inv * rhs;

  int dim;
  sol.FlexibleArnoldi(n_krylov, 
		      n_deflate,
		      rsd_target,
		      V,
		      Z, 
		      H, 
		      R,
		      givens_rots,
		      g,
		      Qk,
		      Qk_tau,
		      dim);
  
  // NOw want to solve System H^\dagger f_m = h_m
  multi1d<DComplex> f_m(n_krylov);
  multi2d<DComplex> evecs(n_krylov,n_krylov);
  multi1d<DComplex> evals(n_krylov);
  multi1d<int> order_array(n_krylov);

  sol.GetEigenvectors(n_krylov,
		      H,
		      f_m,
		      evecs,
		      evals,
		      order_array);

  // Check that the order array is correct:
  DComplex prev_eval=evals[order_array[0]];
  QDPIO::cout << "evals[0] = " << prev_eval << " || evals[0] || = " << norm2(prev_eval) << std::endl;
  for(int i=1; i < n_krylov; i++) { 
    DComplex curr_eval=evals[order_array[i]];
    QDPIO::cout << "evals["<<i<<"] = " << curr_eval << " ||evals["<<i<<"] || = " << norm2(curr_eval) <<std::endl;

    ASSERT_LE( toDouble(norm2(prev_eval)), toDouble(norm2(curr_eval)));
    prev_eval = curr_eval;
  }

  // Check the eigenvalue stuff
  multi1d<DComplex> h_m(n_krylov);

  // h1 is the copy I will pass to LAPACK 
  // I actually want to solve with the dagger of H
  // so I should transpose as well as conjugate
  // however LAPACK expects Fortran order: rows run fastest, columns slower
  // multi2d has Fortran like indexing (row first, column second)
  // but underneath the hood, columns still run fastest.
  // So in principle we need to transpose in and out of Fortran
  // however as I want to solve with the Transpose conjugate
  // all I have to do is leave the indices untransposed, and just conjugate
  for(int col=0; col < n_krylov; ++col) {
    h_m[col] = conj(H(col, n_krylov)); // Keeping this for checking
  }
  // Now form Hm + f_m h_m^H
  // We will send this to Fortran, and we don't want to transpose
  // or conjugate. So just transpose indices on H_tilde now, to ensure
  // correct column major order for Fortran
  multi2d<DComplex> H_tilde(n_krylov,n_krylov);
  for(int row=0; row < n_krylov;++row) { 
    for(int col=0; col< n_krylov; ++col) { 
      H_tilde(col,row) = H(col,row) + f_m[row]*conj(h_m[col]);
    }
  }

  // For each evec check   H_tilde v - lambda v is small
  multi1d<DComplex> diff(n_krylov);
  for(int evec=0; evec < n_krylov; ++evec) { 
  
    // Extract evec from 'evecs'
    multi1d<DComplex> my_evec(n_krylov);
    for(int row=0; row<n_krylov; ++row) {
      my_evec[row] = evecs(evec,row);
    }

    // Multiply: H_tilde my_evec -> diff
    char trans='N'; 
    QDPLapack::zgemv(trans, n_krylov, n_krylov, H_tilde, n_krylov, my_evec, diff);
    
    // Subtract lambda my_evec
    for(int row=0; row < n_krylov; ++row) { 
      diff[row] -= evals[evec]*evecs(evec,row);

      // check result is small
      EXPECT_LT( toDouble(real(diff[row])), 1.0e-14 );
      EXPECT_LT( toDouble(imag(diff[row])), 1.0e-14 );
 
      QDPIO::cout << "eigenvec: " << evec << " diff[" << row <<"]=" << diff[row] << "  lambda=" << evals[evec] << std::endl;
    }
  }
}


TEST_F(FGMRESDRTests, testQR)
{
  std::istringstream input(xml_for_param);
  XMLReader xml_in(input);
  SysSolverFGMRESDRParams p( xml_in, "/Params/InvertParam" );
  LinOpSysSolverFGMRESDR sol(linop,state,p);


  const Subset& s = sol.subset();


    // Create a gaussian source 
  LatticeFermion rhs;
  gaussian(rhs, s);

  Real rsd_target(1.0e-9);
  rsd_target *= sqrt(norm2(rhs,s));

  int n_krylov=p.NKrylov;
  int n_deflate=p.NDefl;
  int total_dim = n_krylov+n_deflate;

  // Work spaces
  multi2d<DComplex> H(total_dim, total_dim + 1); // The H matrix
  multi2d<DComplex> R(total_dim, total_dim + 1); // R = H diagonalized with Givens rotations

  for(int col = 0; col < total_dim; ++col) { 
    for(int row = 0; row < total_dim+1; ++row) {
      H(col,row) = DComplex(0);
      R(col,row) = DComplex(0);
    }
  }

  multi1d<T> V(total_dim+1);  // K(A)
  multi1d<T> Z(total_dim+1);  // K(MA)
  multi1d< Handle<Givens> > givens_rots(total_dim+1);
  multi1d<DComplex> g(total_dim + 1);
  multi1d<DComplex> c(total_dim + 1);
  multi2d<DComplex> Qk;
  multi1d<DComplex> Qk_tau;
  // Assume zero initial guess
  Double beta=sqrt(norm2(rhs, s));

  for(int j=0; j < g.size(); ++j) { 
    g[j] = DComplex(0); 
    c[j] = DComplex(0);
  }
  g[0] = beta;
  c[0] = beta;

  // Set up initial V[0] = rhs / || r^2 || -- 
  Double beta_inv = Double(1)/beta;
  V[0][s] = beta_inv * rhs;

  int dim;
  sol.FlexibleArnoldi(n_krylov, 
		      0, // Guaranteed no deflation
		      rsd_target,
		      V,
		      Z, 
		      H, 
		      R,
		      givens_rots,
		      g,
		      Qk,
		      Qk_tau,
		      dim);


  multi1d<DComplex> eta(n_krylov);
  sol.LeastSquaresSolve(R,g,eta,dim); // Now we can form c - H \eta

  
  // NOw want to solve System H^\dagger f_m = h_m
  multi1d<DComplex> f_m(n_krylov);
  multi2d<DComplex> evecs(n_krylov,n_krylov);
  multi1d<DComplex> evals(n_krylov);
  multi1d<int> order_array(n_krylov);

  // Get the EVs
  sol.GetEigenvectors(n_krylov,
		      H,
		      f_m,
		      evecs,
		      evals,
		      order_array);

 

  // This is where we will store G_{k} = [ g_1 | .. | g_k ]
  //
  // G_{k+1} = [ G_k c - H \eta ]
  //           [  0             ]
  QDPIO::cout << "Setting up G_k and G_{k+1}" << std::endl << std::flush; 
  multi2d<DComplex> Qkplus1(n_deflate+1, n_krylov+1);

  // First copy in the eigenvectors
  for(int col=0; col < n_deflate; ++col) { 
    for(int row=0; row < n_krylov; ++row) { 
      // Copy in the ones corresponding to the lowest k
      // eigenvectors
      Qkplus1(col,row) = evecs( order_array[col], row); 
    }
    Qkplus1(col,n_krylov) = zero;
  }

  for(int row=0; row < n_krylov+1; ++row) { 
    Qkplus1(n_deflate,row) = c[row];
    for(int col=0; col < n_krylov; ++col) { 
      Qkplus1(n_deflate,row) -= H(col,row)*eta(col);
    }
  }

  // NB: Lapack will overwrite both G_{k} and G_{k+1} with the factorizations.

  // QR reflection coefficients.
  multi1d<DComplex> tau_kplus1;
  QDPIO::cout << "QR Decomposing Gk+1" << std::endl;
  QDPLapack::zgeqrf(n_krylov+1, n_deflate+1, Qkplus1, tau_kplus1);

  // Now I want to compute H^{pr}_m Q_k
  // I can do this using ZUNMQR but it will overwrite the 'H' I use, 
  // which at this point is OK in the algorithm, however, here, I will make a copy.
  // to check my tests.

  multi2d<DComplex> H_copy(n_krylov, n_krylov+1);
  multi2d<DComplex> H_copy2(n_krylov, n_krylov+1);
  for(int col=0; col < n_krylov; ++col) {
    for(int row=0; row < n_krylov+1; ++row) { 
      H_copy(col,row) = H(col,row);
      H_copy2(col,row) = H(col,row);
    }
  }

  // Now the total dim here is actually n_krylov the other data in H is invalid
  // So first lets Hit H_copy1 from the right with Q_k from G_k
  char side='R';
  char trans='N';
  QDPIO::cout << "Post multiplying with Q_k using ZUNMQR 1" << std::endl;

  // Here, M and N are the dimensions of the final storage. and 'K' is the order of the decomposition, 
  // E.g. to right Multiply by Q I can just set it to n_deflate rather than n_deflate+1
  
  // !!!! NB: The dimensions here are still the original dimensions of H_copy2,
  // !!!! even tho at the end of this result only the (rows=n_krylov+1)x(cols=n_deflate) portion is 
  // !!!! valid
  QDPLapack::zunmqr2(side,trans, n_krylov+1, n_krylov, n_deflate, Qkplus1, tau_kplus1, H_copy);


  QDPIO::cout << "Pre multiply with Q^{H}_{k+1} usign ZUNMQR2" << std::endl;
  trans='C'; // Multiply with Herm Conjugate
  side='L';  // Multiply from Left
  
  // !!!! NB: The dimensions here are still the original dimensions of H_copy2 
  // !!!! Even tho when we are done here, really only the (rows=k+1)x(cols=k) portion of it 
  // !!!! is what is relevant
  QDPLapack::zunmqr2(side,trans, n_krylov+1, n_krylov, n_deflate+1, Qkplus1, tau_kplus1, H_copy);

  // Now we should check this. One way is by explicitly forming Q_plus_1 -- I need this
  // anyhow for transforming the V's and the Z's.

  // Recover Q_plus1 -- this reconstructs it explicitly. We will need Qkplus1
  // and Qk to reconstruct V_{k+1} and Z_{k}
  QDPLapack::zungqr(n_krylov+1,n_deflate+1,n_deflate+1, Qkplus1, tau_kplus1);


  // NOw try and form  H_k Q_k using ZGEMM
  // 
  // H_k has n_krylov+1 rows, and n_krylov columns
  // Q_k has n_krylov rows, and n_deflate columns
  // 
  // so H_k Q_k has n_krylov+1 rows and n_deflate columns
  multi2d<DComplex> HKQK(n_deflate, n_krylov+1);
  
  char transa = 'N';
  char transb = 'N'; 
  DComplex calpha(1);
  DComplex cbeta = zero;
  int m = n_krylov+1;
  int k = n_krylov;
  int n = n_deflate;

  // This will evaluate:  H_m Q_k
  // by beating it into the form alpha * H_m * Q_k + 0 * H_copy1 => H_copy1
  // NB: I don't form Q_k, but just use Q_{k+1}, and take care of the fact 
  // that it is bigger than I want by using the leading dimension of Q 
  // and by setting the m,n,k dimensions appropriately.
  QDPLapack::zgemm(transa, transb,
		   m,n,k,
		   calpha,
		   H_copy2,
		   m,
		   Qkplus1,
		   m,      // The leading dimension of Qkplus1 is n_krylov+1 i.e. 'm' even if we multiply only with the kxn portion
		   cbeta,
		   HKQK,
		   m);

  multi2d<DComplex> QK1HKQK(n_deflate, n_deflate+1);  

  transa = 'C';
  transb = 'N'; 
  
  m = n_deflate+1; // rows of Q_{k+1}^H
  k = n_krylov+1;  // cols of Q_{k+1}^H
  n = n_deflate; // cols of Q^{H}_{k+1} H_{k} Q_{k} 

  // This will evaluate:  H_m Q_k
  // by beating it into the form alpha * H_m * Q_k + 0 * H_copy1 => H_copy1
  // NB: I don't form Q_k, but just use Q_{k+1}, and take care of the fact 
  // that it is bigger than I want by using the leading dimension of Q 
  // and by setting the m,n,k dimensions appropriately.
  QDPLapack::zgemm(transa, transb,
		   m,n,k,
		   calpha,
		   Qkplus1,
		   k,
		   HKQK,
		   k,      // The leading dimension of Qkplus1 is n_krylov+1 i.e. 'm' even if we multiply only with the kxn portion
		   cbeta,
		   QK1HKQK,
		   m);

  QDPIO::cout << "QK1HKQK - H_copy" << std::endl;

  for(int row = 0; row < n_deflate+1; ++row) {
    for(int col=0; col < n_deflate; ++col) { 
      QK1HKQK(col,row) -= H_copy(col,row);
      QDPIO::cout <<  QK1HKQK(col,row) << "\t" ;
      EXPECT_LT( toDouble(sqrt(norm2(QK1HKQK(col,row)))), 1.0e-15 );     
    }
    QDPIO::cout << std::endl;
  }

  // OK Now I need to form: V_{k+1} = V^{pr}_{m+1} Q_{k+1}
  multi1d<LatticeFermion> new_V(n_deflate+1);
  for(int i=0; i < n_deflate+1; ++i) {
    new_V[i] = zero;
    for(int j=0; j < n_krylov+1; ++j) { 
      new_V[i] += V[j]*Qkplus1(i,j);
    }
  }

  // Check Columns are Orthonormal
  for(int i=0; i < n_deflate+1; ++i) {
    for(int j=0; j < n_deflate+1; ++j ) { 
      DComplex iprod = innerProduct(new_V[i], new_V[j], s);
      QDPIO::cout << "< V[" << i <<"],V["<<j<<"] > = " << iprod ;
      
      if ( i==j ) { 
	// iprod should be 1
	Double diff_r = Double(1) - sqrt(real(iprod));
	Double diff_i = fabs(imag(iprod));
	QDPIO::cout << " ABS DIFF=" << fabs(diff_r) << " ABS DIFF IM=" << diff_i << std::endl;
	EXPECT_LT( toDouble(diff_r), 1.0e-11);
	EXPECT_LT( toDouble(diff_i), 1.0e-12);

      }
      else { 
	Double diff_r = fabs(real(iprod));
	Double diff_i = fabs(imag(iprod));
	QDPIO::cout << " ABS DIFF RE=" << diff_r << " ABS DIFF IM=" << diff_i << std::endl;
	EXPECT_LT( toDouble(diff_r), 1.0e-9);
	EXPECT_LT( toDouble(diff_i), 1.0e-9);
      }
    }
  }

  // 
  multi1d<LatticeFermion> new_Z(n_deflate);
  for(int i=0; i < n_deflate; ++i) {
    new_Z[i] = zero;
    for(int j=0; j < n_krylov; ++j) { 
      new_Z[i] += Z[j]*Qkplus1(i,j);
    }
  }



}
