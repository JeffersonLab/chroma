// $Id: mesq.cc,v 3.0 2006-04-03 04:58:58 edwards Exp $

#error "NOT FULLY CONVERTED"

/* Return the value of the trivial geometric definition of */
/* the topological Charge */

/* Only for U(1) in 2 dim (Nc=1,Nd=2)*/

/* u -- gauge field (Read) */
/* q -- topological charge (Write) */

include(types.mh)

SUBROUTINE(MesQ, u, q)

multi1d<LatticeColorMatrix> u(Nd);
Double q;
{
  include(COMMON_DECLARATIONS)
  
  LatticeColorMatrix tmp_0;
  LatticeColorMatrix tmp_1;
  LatticeComplex wplaq_tmp;
  LatticeReal tmp_r;
  LatticeReal tmp_x;
  LatticeReal tmp_y;
  LatticeReal phi;

  Double tmp;
  Double tmp2;
  Double tmp3;

  int cb;
  int mu;
  int nu;

  START_CODE();
   
               
  phi = 0;

  if ( Nd != 2 ) 
    QDP_error_exit("only support 2 dimensions");

  if ( Nc != 1 ) 
    QDP_error_exit("only support for U[1]);
  
  mu=0;
  nu=1;

  for(cb=0; cb < 2; ++cb)
  {
    /* tmp_0 = u(x+mu,nu)*u_dag(x+nu,mu) */
    tmp_0 = shift(u[nu][1-cb], cb, FORWARD, mu) * shift(adj[u[mu][1-cb]], cb, FORWARD, nu);
    
    /* tmp_1 = tmp_0*u_dag(x,nu)=u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu) */
    tmp_1 = tmp_0 * adj(u[nu][cb]);
    
    /* wplaq_tmp = tr(u(x,mu)*tmp_1=u(x,mu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu)) */
    wplaq_tmp = trace(u[mu][cb] * tmp_1);
    
    tmp_x = real(wplaq_tmp);
    tmp_y = imag(wplaq_tmp);
    tmp_r = atan2(tmp_y, tmp_x);
    phi += tmp_r;
  }
 
  tmp = -sum(phi);

              
  tmp2 = 8*atan(1);
  q = tmp / tmp2;

  END_CODE();
}
