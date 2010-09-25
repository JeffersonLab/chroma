// $Id: wilslp.cc,v 3.3 2006-07-08 04:34:55 edwards Exp $
/*! \file
 *  \brief Calculate Wilson loops
 */

#include "chromabase.h"
#include "meas/glue/wilslp.h"
#include "meas/gfix/axgauge.h"

namespace Chroma 
{

  //! Calculate Wilson loops
  /*!
   * \ingroup glue
   *
   * Calculates, depending on option, (1) "space-like" planar   
   * Wilson loops in the directions perpendicular to j_decay    
   * that have equal length, (2) "time-like" planar Wilson loops
   * with time direction j_decay and space directions the       
   * perpendicular ones that have equal length and (3) off-axis 
   * "time-like" Wilson loops along 3 paricular paths in the    
   * space directions that have equal length.                   
   *
   * \param u          gauge field (Read)                                              
   * \param j_decay    decay direction (Read)                                     
   * \param t_dir      time direction (Read)                                     
   * \param kind       binary-combined YES/NO [1/0] of the three options (Read)      
   *                   e.g. kind = 2 gives planar t-like, kind=6 is 
   *                   planar + off-axis: sqrt(2), sqrt(5), sqrt(3)
   */

  void wilslp(const multi1d<LatticeColorMatrix>& u, 
	      int j_decay, int t_dir, int kind,
	      XMLWriter& xml, const string& xml_group)
  {
    START_CODE();

    multi1d<int> nrow(Nd);
    nrow = Layout::lattSize();
    
    multi1d<LatticeColorMatrix> ug(Nd);
    multi1d<int> space_dir(Nd);
    LatticeColorMatrix   u_t;
    LatticeColorMatrix   u_tmp;
    LatticeColorMatrix   u_diag;
    LatticeColorMatrix   u_space;
    LatticeColorMatrix tmp_2;
    LatticeColorMatrix tmp_3;
    LatticeReal wl_trace;
    LatticeInteger t_coord;
    LatticeBoolean btmp;
    Double dummy;
    Real ftmp;

    int lengthr;
    int lengtht;
    int lsizet;
    int lsizer;
    int nspace;
    int mu;
    int nu;
    int rho;
    int r;
    int t;
    int tt;
    int nr;
    int i;
    int j;
    int r_off;
  
    lsizet  = nrow[j_decay];
    lengtht = lsizet / 2;
  
    /* Determine 'lattice shape' */
    space_dir = 1; 
    nr = 1;
    for(mu = 0;mu  <= ( Nd-2); ++mu )
      if ( mu != j_decay && mu != t_dir )
      {
	lsizer = nrow[mu];
	for(nu = mu+1;nu  < ( Nd); ++nu )
	  if ( nu != j_decay && nu != t_dir )
	  {
	    r = nrow[nu];
	    if( r == lsizer )
	      space_dir[mu] += nr;
	  }
      }

    nspace = space_dir[0];
    nu = 0;
    for(mu = 1;mu  <= ( Nd-2); ++mu )
      if ( mu != j_decay && mu != t_dir )
      {
	nr = space_dir[mu];
	if( nr > nspace )
	{
	  nspace = nr;
	  nu = mu;
	}
      }

    lsizer = nrow[nu];
    nr = 0;
    space_dir[nr] = nu;
    for(mu = nu+1;mu  < ( Nd); ++mu )
      if ( mu != j_decay && mu != t_dir )
      {
	r = nrow[mu];
	if( r == lsizer )
	{
	  nr = nr + 1;
	  space_dir[nr] = mu;
	}
      }

    QDPIO::cout << "nspace = " << nspace
		<< "  space_dir = " 
		<< space_dir[0] << " " 
		<< space_dir[1] << " " 
		<< space_dir[2] << " " 
		<< space_dir[3] << endl;

    if( nr != (nspace-1) )
    {
      QDPIO::cerr << __func__ 
		  << ": Trouble with space dimensions for Wilson loops: nr = " << nr
		  << "  nspace = " << nspace
		  << "  nu = " << nu << endl;
      QDP_abort(1);
    }

//    lengthr = lsizer/2;    // Old style with measurement on L/2
    lengthr = lsizer;      // Now use all L

    /* Compute "space-like" planar Wilson loops, if desired */
    if ( (kind & 1) != 0 )
    {
      if( nspace < 2 )
	for(mu = 0;mu  < ( Nd); ++mu )
	  QDP_error_exit("Wrong lattice size for space-like Wilson loops: ", 
			 mu, nrow[mu]);
                    
      multi2d<Double> wils_loop1(lengthr, lengthr);
      wils_loop1 = 0;			/* initialize the Wilson loops */

      for(j = 1;j  < ( nspace); ++j )
      {
	nu = space_dir[j];

	/*+ */
	/* Fix to axial gauge in nu-direction */
	/*- */
	ug = u;
	axGauge (ug, nu);
	t_coord = Layout::latticeCoordinate(nu);

	for(i = 0;i  < ( j); ++i )
	{
	  mu = space_dir[i];

	  for(r = 0;r  < ( lengthr); ++r )
	  {
	    /* Gather t-links (i.e. nu-links) from r-direction (i.e. mu) */
	    /* and make "r-links" */

	    if ( r == 0 )
	    {
	      u_t = shift(ug[nu], FORWARD, mu);
	      u_space = ug[mu];
	    }
	    else
	    {
	      tmp_2 = shift(u_t, FORWARD, mu);
	      u_t = tmp_2;

	      tmp_2 = shift(u_space, FORWARD, mu);
	      u_space = ug[mu] * tmp_2;
	    }

	    for(t = 0;t  < ( lengthr); ++t )
	    {
	      /* Gather r-link from t-direction (i.e. nu) */

	      if ( t == 0 )
	      {
		u_tmp = shift(u_space, FORWARD, nu);
	      }
	      else
	      {
		tmp_2 = shift(u_tmp, FORWARD, nu);
		u_tmp = tmp_2;
	      }

	      tt = lsizer - t - 1;
	      btmp = t_coord < tt;

	      tmp_2 = ug[nu] * u_tmp;
	      tmp_3 = tmp_2 * adj(u_t);
	      copymask(tmp_3, btmp, u_tmp);
	      tmp_2 = tmp_3 * adj(u_space);
	      wl_trace = real(trace(tmp_2));
	      wils_loop1[r][t] += sum(wl_trace);
	    }   /* end t loop */
	  }     /* end r loop */
	}       /* end i loop (for mu) */
      }         /* end j loop (for nu) */

      dummy = 2.0 / double (Layout::vol()*Nc*nspace*(nspace-1)) ;

      push(xml, "wils_loop1"); // XML tag for wils_wloop1
      write(xml, "lengthr", lengthr);
      push(xml, "wloop1");
                  
      multi1d<Double> wloop1(lengthr);

      for(r = 0; r < lengthr; ++r)
      {
	for(t = 0; t < lengthr; ++t)
	{
	  wils_loop1[t][r] = wils_loop1[t][r] * dummy;
	  wloop1[t]      = wils_loop1[t][r];
	}
	push(xml, "elem");
	write(xml, "r", r);
	write(xml, "loop", wloop1);   // write out wils_wloop1
	pop(xml); // elem
      } // end for r

      pop(xml);               // XML end tag for wloop1
      pop(xml);               // XML end tag for wils_wloop1
      QDPIO::cout << "wils_loop1 data written to .xml file " << endl;        

    }           /* end of option "space-like planar Wilson loops" */

    /* Fix to axial gauge */
    ug = u;
    axGauge (ug, j_decay);
    t_coord = Layout::latticeCoordinate(j_decay);

    /* Compute "time-like" planar Wilson loops, if desired */
    if ( (kind & 2) != 0 )
    {
      multi2d<Double> wils_loop2(lengtht, lengthr);
      wils_loop2 = 0;                       /* initialize the Wilson loops */

      QDPIO::cout << "computing time-like Wilson loops" << endl;

      for(i = 0;i  < ( nspace); ++i )
      {
	mu = space_dir[i];
	for(r = 0;r  < ( lengthr); ++r )
	{
	  /* Gather t-links from space-direction (i.e. mu) */
	  /* and make "space-links" */

	  if ( r == 0 )
	  {
	    u_t = shift(ug[j_decay], FORWARD, mu);
	    u_space = ug[mu];
	  }
	  else
	  {
	    tmp_2 = shift(u_t, FORWARD, mu);
	    u_t   = tmp_2;

	    tmp_2 = shift(u_space, FORWARD, mu);
	    tmp_3 = ug[mu] * tmp_2;
	    u_space = tmp_3;
	  }

	  for(t = 0;t  < ( lengtht); ++t )
	  {
	    /* Gather space-link from t-direction (i.e. j_decay) */

	    if ( t == 0 )
	    {
	      u_tmp = shift(u_space, FORWARD, j_decay);
	    }
	    else
	    {
	      tmp_2 = shift(u_tmp, FORWARD, j_decay);
	      u_tmp = tmp_2;
	    }
	    tt = lsizet - t - 1;
	    btmp = t_coord < tt;

	    tmp_2 = ug[j_decay] * u_tmp;
	    tmp_3 = tmp_2 * adj(u_t);
	    copymask(tmp_3, btmp, u_tmp);
	    tmp_2 = tmp_3 * adj(u_space);
	    wl_trace = real(trace(tmp_2));
	    wils_loop2[t][r] += sum(wl_trace);

	  }    /* end t loop */
	}      /* end r loop */
      }        /* end i loop (for mu) */

      dummy = 1.0 / double (Layout::vol()*Nc*nspace) ;

      push(xml,"wils_loop2"); // XML tag for wils_wloop2
      write(xml, "lengthr", lengthr);
      write(xml, "lengtht", lengtht);
      push(xml, "wloop2");

      multi1d<Double> wloop2(lengtht);

      for(r = 0; r < lengthr; ++r)
      {
	for(t = 0; t < lengtht; ++t)
	{
	  wils_loop2[t][r] *= dummy;
	  wloop2[t]         = wils_loop2[t][r];
	}
	push(xml, "elem");
	write(xml, "r", r);
	write(xml, "loop", wloop2);   // write out wils_loop2
	pop(xml); // elem
      } // end for r

      pop(xml);               // XML end tag for loop2
      pop(xml);               // XML end tag for wils_loop2
      QDPIO::cout << "wils_loop2 data written to .xml file " << endl;  
    }          /* end on option "time-like planar Wilson loops" */

    /* Compute "time-like" off-axis Wilson loops, if desired */
    if ( (kind & 4) != 0 )
    {
      int length;
      if( nspace < 2 )
      {
	for(mu = 0;mu  < ( Nd); ++mu )
	  QDP_error_exit("Wrong lattice size for off-axis Wilson loops: ", 
			 mu, nrow[mu]);
      }
      else if ( nspace == 2 )
      {
	length = lengthr + lengthr/2;
      }
      else
	length = 2 * lengthr + lengthr/2;

      QDPIO::cout << "wils_loop3: lengtht=" << lengtht << "  length=" << length 
		  << "  nspace=" << nspace << "   j_decay=" << j_decay
		  << "  t_dir=" << t_dir << endl;
      multi2d<Double> wils_loop3(lengtht,  length);  
      wils_loop3 = 0;     /* initialize the non-planar Wilson loops */

      for(j = 1;j  < ( nspace); ++j )
      {
	nu = space_dir[j];

	for(i = 0;i  < ( j); ++i )
	{
	  mu = space_dir[i];

	  /*+ */
	  /* Do off-axis "sqrt(2)" loops in (mu,nu)-plane */
	  /*- */

	  /* Make the "corner link" in the (mu,nu) direction */
	  ftmp = 1.0 / 2.0 ;
	  tmp_2 = shift(ug[nu], FORWARD, mu);
	  u_diag = ug[mu] * tmp_2;
	  tmp_2 = shift(ug[mu], FORWARD, nu);
	  u_diag += ug[nu] * tmp_2;
	  u_diag = u_diag * ftmp;

	  for(r = 0;r  < ( lengthr); ++r )
	  {
	    /* Gather t-links from space-directions (i.e. mu and nu) */
	    /* and make "space-links" */

	    if ( r == 0 )
	    {
	      tmp_2 = shift(ug[j_decay], FORWARD, mu);
	      u_t   = shift(tmp_2, FORWARD, nu);
	      u_space = u_diag;
	    }
	    else
	    {
	      tmp_2 = shift(u_t, FORWARD, mu);
	      u_t   = shift(tmp_2, FORWARD, nu);
	      tmp_2 = shift(u_space, FORWARD, mu);
	      tmp_3 = shift(tmp_2, FORWARD, nu);
	      u_space = u_diag * tmp_3;
	    }

	    for(t = 0;t  < ( lengtht); ++t )
	    {
	      /* Gather space-link from t-direction (i.e. j_decay) */
	      if ( t == 0 )
	      {
		u_tmp = shift(u_space, FORWARD, j_decay);
	      }
	      else
	      {
		tmp_2 = shift(u_tmp, FORWARD, j_decay);
		u_tmp = tmp_2;
	      }

	      tt = lsizet - t - 1;
	      btmp = t_coord < tt;
	      tmp_2 = ug[j_decay] * u_tmp;
	      tmp_3 = tmp_2 * adj(u_t);
	      copymask(tmp_3, btmp, u_tmp);
	      tmp_2 = tmp_3 * adj(u_space);
	      wl_trace = real(trace(tmp_2));
	      wils_loop3[t][r] += sum(wl_trace);
	    }        /* end t loop */
	  }          /* end r loop */

	  /*+ */
	  /* Do off-axis "sqrt(2)" loops in (mu,-nu)-plane */
	  /*- */

	  /* Make the "corner link" in the (mu,-nu) direction */
	  ftmp = 1.0 / 2.0;
	  tmp_2 = adj(ug[nu]) * ug[mu];
	  u_diag = shift(tmp_2, BACKWARD, nu);
	  tmp_2 = shift(ug[nu], FORWARD, mu);
	  tmp_3 = shift(tmp_2, BACKWARD, nu);
	  u_diag += ug[mu] * adj(tmp_3);
	  u_diag = u_diag * ftmp;

	  for(r = 0;r  < ( lengthr); ++r )
	  {
	    /* Gather t-links from space-directions (i.e. mu and -nu) */
	    /* and make "space-links" */

	    if ( r == 0 )
	    {
	      tmp_2 = shift(ug[j_decay], FORWARD, mu);
	      u_t   = shift(tmp_2, BACKWARD, nu);
	      u_space = u_diag;
	    }
	    else
	    {
	      tmp_2 = shift(u_t, FORWARD, mu);
	      u_t   = shift(tmp_2, BACKWARD, nu);
	      tmp_2 = shift(u_space, FORWARD, mu);
	      tmp_3 = shift(tmp_2, BACKWARD, nu);
	      u_space = u_diag * tmp_3;
	    }

	    for(t = 0;t  < ( lengtht); ++t )
	    {
	      /* Gather space-link from t-direction (i.e. j_decay) */
	      if ( t == 0 )
	      {
		u_tmp = shift(u_space, FORWARD, j_decay);
	      }
	      else
	      {
		tmp_2 = shift(u_tmp, FORWARD, j_decay);
		u_tmp = tmp_2;
	      }

	      tt = lsizet - t - 1;
	      btmp = t_coord < tt;
	      tmp_2 = ug[j_decay] * u_tmp;
	      tmp_3 = tmp_2 * adj(u_t);
	      copymask(tmp_3, btmp, u_tmp);
	      tmp_2 = tmp_3 * adj(u_space);
	      wl_trace = real(trace(tmp_2));
	      wils_loop3[t][r] += sum(wl_trace);

	    }        /* end t loop */
	  }          /* end r loop */

	  /*+ */
	  /* Do off-axis "sqrt(5)" loops in (2*mu,nu)-plane */
	  /*- */
	  r_off = lengthr;

	  /* Make the "corner link" in the (mu,nu) direction */
	  ftmp = 1.0 / 2.0;

	  tmp_2 = shift(ug[nu], FORWARD, mu);
	  tmp_3 = ug[mu] * tmp_2;
	  tmp_2 = shift(tmp_3, FORWARD, mu);
	  u_diag = ug[mu] * tmp_2;

	  tmp_2 = shift(ug[mu], FORWARD, mu);
	  tmp_3 = ug[mu] * tmp_2;
	  tmp_2 = shift(tmp_3, FORWARD, nu);
	  u_diag += ug[nu] * tmp_2;

	  u_diag = u_diag * ftmp;

	  for(r = 0;r  < ( lengthr/2); ++r )
	  {
	    /* Gather t-links from space-directions (i.e. 2*mu and nu) */
	    /* and make "space-links" */

	    if ( r == 0 )
	    {
	      tmp_2 = shift(ug[j_decay], FORWARD, mu);
	      tmp_3 = shift(tmp_2, FORWARD, mu);
	      u_t   = shift(tmp_3, FORWARD, nu);
	      u_space = u_diag;
	    }
	    else
	    {
	      tmp_2 = shift(u_t, FORWARD, mu);
	      tmp_3 = shift(tmp_2, FORWARD, mu);
	      u_tmp = shift(tmp_3, FORWARD, nu);
	      u_t = u_tmp;

	      tmp_2 = shift(u_space, FORWARD, mu);
	      tmp_3 = shift(tmp_2, FORWARD, mu);
	      tmp_2 = shift(tmp_3, FORWARD, nu);
	      u_tmp = u_diag * tmp_2;
	      u_space = u_tmp;
	    }

	    for(t = 0;t  < ( lengtht); ++t )
	    {
	      /* Gather space-link from t-direction (i.e. j_decay) */
	      if ( t == 0 )
	      {
		u_tmp = shift(u_space, FORWARD, j_decay);
	      }
	      else
	      {
		tmp_2 = shift(u_tmp, FORWARD, j_decay);
		u_tmp = tmp_2;
	      }

	      tt = lsizet - t - 1;
	      btmp = t_coord < tt;
	      tmp_2 = ug[j_decay] * u_tmp;
	      tmp_3 = tmp_2 * adj(u_t);
	      copymask(tmp_3, btmp, u_tmp);
	      tmp_2 = tmp_3 * adj(u_space);
	      wl_trace = real(trace(tmp_2));
	      wils_loop3[t][r_off+r] += sum(wl_trace);

	    }        /* end t loop */
	  }          /* end r loop */

	  /*+ */
	  /* Do off-axis "sqrt(5)" loops in (2*mu,-nu)-plane */
	  /*- */

	  /* Make the "corner link" in the (2*mu,-nu) direction */
	  ftmp = 1.0 / 2.0;

	  tmp_2 = shift(ug[mu], FORWARD, mu);
	  tmp_3 = ug[mu] * tmp_2;
	  tmp_2 = adj(ug[nu]) * tmp_3;
	  u_diag = shift(tmp_2, BACKWARD, nu);

	  tmp_2 = shift(ug[nu], FORWARD, mu);
	  tmp_3 = shift(tmp_2, BACKWARD, nu);
	  tmp_2 = ug[mu] * adj(tmp_3);
	  tmp_3 = shift(tmp_2, FORWARD, mu);
	  u_diag += ug[mu] * tmp_3;

	  u_diag = u_diag * ftmp;

	  for(r = 0;r  < ( lengthr/2); ++r )
	  {
	    /* Gather t-links from space-directions (i.e. 2*mu and -nu) */
	    /* and make "space-links" */

	    if ( r == 0 )
	    {
	      tmp_2 = shift(ug[j_decay], FORWARD, mu);
	      tmp_3 = shift(tmp_2, FORWARD, mu);
	      u_t   = shift(tmp_3, BACKWARD, nu);
	      u_space = u_diag;
	    }
	    else
	    {
	      tmp_2 = shift(u_t, FORWARD, mu);
	      tmp_3 = shift(tmp_2, FORWARD, mu);
	      u_tmp = shift(tmp_3, BACKWARD, nu);
	      u_t = u_tmp;

	      tmp_2 = shift(u_space, FORWARD, mu);
	      tmp_3 = shift(tmp_2, FORWARD, mu);
	      tmp_2 = shift(tmp_3, BACKWARD, nu);
	      u_tmp = u_diag * tmp_2;
	      u_space = u_tmp;
	    }

	    for(t = 0;t  < ( lengtht); ++t )
	    {
	      /* Gather space-link from t-direction (i.e. j_decay) */
	      if ( t == 0 )
	      {
		u_tmp = shift(u_space, FORWARD, j_decay);
	      }
	      else
	      {
		tmp_2 = shift(u_tmp, FORWARD, j_decay);
		u_tmp = tmp_2;
	      }


	      tt = lsizet - t - 1;
	      btmp = t_coord < tt;

	      tmp_2 = ug[j_decay] * u_tmp;
	      tmp_3 = tmp_2 * adj(u_t);
	      copymask(tmp_3, btmp, u_tmp);
	      tmp_2 = tmp_3 * adj(u_space);
	      wl_trace = real(trace(tmp_2));
	      wils_loop3[t][r_off+r] += sum(wl_trace);
	    }        /* end t loop */
	  }          /* end r loop */


	  /*+ */
	  /* Do off-axis "sqrt(5)" loops in (2*nu,mu)-plane */
	  /*- */

	  /* Make the "corner link" in the (2*nu,mu) direction */
	  ftmp = 1.0 / 2.0;

	  tmp_2 = shift(ug[mu], FORWARD, nu);
	  tmp_3 = ug[nu] * tmp_2;
	  tmp_2 = shift(tmp_3, FORWARD, nu);
	  u_diag = ug[nu] * tmp_2;

	  tmp_2 = shift(ug[nu], FORWARD, nu);
	  tmp_3 = ug[nu] * tmp_2;
	  tmp_2 = shift(tmp_3, FORWARD, mu);
	  u_diag += ug[mu] * tmp_2;

	  u_diag = u_diag * ftmp;

	  for(r = 0;r  < ( lengthr/2); ++r )
	  {
	    /* Gather t-links from space-directions (i.e. 2*nu and mu) */
	    /* and make "space-links" */

	    if ( r == 0 )
	    {
	      tmp_2 = shift(ug[j_decay], FORWARD, nu);
	      tmp_3 = shift(tmp_2, FORWARD, nu);
	      u_t   = shift(tmp_3, FORWARD, mu);
	      u_space = u_diag;
	    }
	    else
	    {
	      tmp_2 = shift(u_t, FORWARD, nu);
	      tmp_3 = shift(tmp_2, FORWARD, nu);
	      u_tmp = shift(tmp_3, FORWARD, mu);
	      u_t = u_tmp;

	      tmp_2 = shift(u_space, FORWARD, nu);
	      tmp_3 = shift(tmp_2, FORWARD, nu);
	      tmp_2 = shift(tmp_3, FORWARD, mu);
	      u_tmp = u_diag * tmp_2;
	      u_space = u_tmp;
	    }

	    for(t = 0;t  < ( lengtht); ++t )
	    {
	      /* Gather space-link from t-direction (i.e. j_decay) */
	      if ( t == 0 )
	      {
		u_tmp = shift(u_space, FORWARD, j_decay);
	      }
	      else
	      {
		tmp_2 = shift(u_tmp, FORWARD, j_decay);
		u_tmp = tmp_2;
	      }

	      tt = lsizet - t - 1;
	      btmp = t_coord < tt;

	      tmp_2 = ug[j_decay] * u_tmp;
	      tmp_3 = tmp_2 * adj(u_t);
	      copymask(tmp_3, btmp, u_tmp);
	      tmp_2 = tmp_3 * adj(u_space);
	      wl_trace = real(trace(tmp_2));
	      wils_loop3[t][r_off+r] += sum(wl_trace);
	    }        /* end t loop */
	  }          /* end r loop */


	  /*+ */
	  /* Do off-axis "sqrt(5)" loops in (2*nu,-mu)-plane */
	  /*- */

	  /* Make the "corner link" in the (2*nu,-mu) direction */
	  ftmp = 1.0 / 2.0;
	  tmp_2 = shift(ug[nu], FORWARD, nu);
	  tmp_3 = ug[nu] * tmp_2;
	  tmp_2 = adj(ug[mu]) * tmp_3;
	  u_diag = shift(tmp_2, BACKWARD, mu);

	  tmp_2 = shift(ug[mu], FORWARD, nu);
	  tmp_3 = shift(tmp_2, BACKWARD, mu);
	  tmp_2 = ug[nu] * adj(tmp_3);
	  tmp_3 = shift(tmp_2, FORWARD, nu);
	  u_diag += ug[nu] * tmp_3;

	  u_diag = u_diag * ftmp;

	  for(r = 0;r  < ( lengthr/2); ++r )
	  {
	    /* Gather t-links from space-directions (i.e. 2*nu and -mu) */
	    /* and make "space-links" */

	    if ( r == 0 )
	    {
	      tmp_2 = shift(ug[j_decay], FORWARD, nu);
	      tmp_3 = shift(tmp_2, FORWARD, nu);
	      u_t   = shift(tmp_3, BACKWARD, mu);
	      u_space = u_diag;
	    }
	    else
	    {
	      tmp_2 = shift(u_t, FORWARD, nu);
	      tmp_3 = shift(tmp_2, FORWARD, nu);
	      u_tmp = shift(tmp_3, BACKWARD, mu);
	      u_t   = u_tmp;

	      tmp_2 = shift(u_space, FORWARD, nu);
	      tmp_3 = shift(tmp_2, FORWARD, nu);
	      tmp_2 = shift(tmp_3, BACKWARD, mu);
	      u_tmp = u_diag * tmp_2;
	      u_space = u_tmp;
	    }

	    for(t = 0;t  < ( lengtht); ++t )
	    {
	      /* Gather space-link from t-direction (i.e. j_decay) */
	      if ( t == 0 )
	      {
		u_tmp = shift(u_space, FORWARD, j_decay);
	      }
	      else
	      {
		tmp_2 = shift(u_tmp, FORWARD, j_decay);
		u_tmp = tmp_2;
	      }

	      tt = lsizet - t - 1;
	      btmp = t_coord < tt;

	      tmp_2 = ug[j_decay] * u_tmp;
	      tmp_3 = tmp_2 * adj(u_t);
	      copymask(tmp_3, btmp, u_tmp);
	      tmp_2 = tmp_3 * adj(u_space);
	      wl_trace = real(trace(tmp_2));
	      wils_loop3[t][r_off+r] += sum(wl_trace);

	    }        /* end t loop */
	  }          /* end r loop */
	}            /* end i loop (for mu) */
      }              /* end j loop (for nu) */

      dummy = 1.0 / double (Layout::vol()*Nc*nspace*(nspace-1)) ;
      for(t = 0;t  < ( lengtht); ++t )
	for(r = 0;r  < ( lengthr); ++r )
	  wils_loop3[t][r] = wils_loop3[t][r] * dummy;

      dummy = 1.0 / double (Layout::vol()*Nc*2*nspace*(nspace-1)) ;
      for(t = 0;t  < ( lengtht); ++t )
	for(r = r_off;r  < ( r_off+lengthr/2); ++r )
	  wils_loop3[t][r] = wils_loop3[t][r] * dummy;

      if ( nspace > 2 )
      {
	int k;
	/*+ */
	/* Do off-axis "sqrt(3)" loops */
	/*- */
	r_off = lengthr + lengthr/2;

	for(k = 2;k  < ( nspace); ++k )
	{
	  rho = space_dir[k];

	  for(j = 1;j  < ( k); ++j )
	  {
	    nu = space_dir[j];

	    for(i = 0;i  < ( j); ++i )
	    {
	      mu = space_dir[i];

	      /*+ */
	      /* Do off-axis "sqrt(3)" loops in (mu,nu,rho)-plane */
	      /*- */

	      /* Make the "corner link" in the (mu,nu,rho) direction */
	      ftmp = 1.0 / 6.0;
	      tmp_2 = shift(ug[rho], FORWARD, nu);
	      u_tmp = ug[nu] * tmp_2;

	      tmp_2 = shift(ug[nu], FORWARD, rho);
	      u_tmp += ug[rho] * tmp_2;

	      tmp_2 = shift(u_tmp, FORWARD, mu);
	      u_diag = ug[mu] * tmp_2;

	      tmp_2 = shift(ug[rho], FORWARD, mu);
	      u_tmp = ug[mu] * tmp_2;

	      tmp_2 = shift(ug[mu], FORWARD, rho);
	      u_tmp += ug[rho] * tmp_2;

	      tmp_2 = shift(u_tmp, FORWARD, nu);
	      u_diag += ug[nu] * tmp_2;

	      tmp_2 = shift(ug[nu], FORWARD, mu);
	      u_tmp = ug[mu] * tmp_2;

	      tmp_2 = shift(ug[mu], FORWARD, nu);
	      u_tmp += ug[nu] * tmp_2;

	      tmp_2 = shift(u_tmp, FORWARD, rho);
	      u_diag += ug[rho] * tmp_2;

	      u_diag = u_diag * ftmp;

	      for(r = 0;r  < ( lengthr); ++r )
	      {
		/* Gather t-links from space-directions (i.e. mu, nu and rho) */
		/* and make "space-links" */

		if ( r == 0 )
		{
		  tmp_2 = shift(ug[j_decay], FORWARD, mu);
		  tmp_3 = shift(tmp_2, FORWARD, nu);
		  u_t   = shift(tmp_3, FORWARD, rho);

		  u_space = u_diag;
		}
		else
		{
		  tmp_2 = shift(u_t, FORWARD, mu);
		  tmp_3 = shift(tmp_2, FORWARD, nu);
		  u_tmp = shift(tmp_3, FORWARD, rho);
		  u_t   = u_tmp;

		  tmp_2 = shift(u_space, FORWARD, mu);
		  tmp_3 = shift(tmp_2, FORWARD, nu);
		  tmp_2 = shift(tmp_3, FORWARD, rho);
		  u_tmp = u_diag * tmp_2;
		  u_space = u_tmp;
		}

		for(t = 0;t  < ( lengtht); ++t )
		{
		  /* Gather space-link from t-direction (i.e. j_decay) */
		  if ( t == 0 )
		  {
		    u_tmp = shift(u_space, FORWARD, j_decay);
		  }
		  else
		  {
		    tmp_2 = shift(u_tmp, FORWARD, j_decay);
		    u_tmp = tmp_2;  
		  }

		  tt = lsizet - t - 1;
		  btmp = t_coord < tt;

		  tmp_2 = ug[j_decay] * u_tmp;
		  tmp_3 = tmp_2 * adj(u_t);
		  copymask(tmp_3, btmp, u_tmp);
		  tmp_2 = tmp_3 * adj(u_space);
		  wl_trace = real(trace(tmp_2));
		  wils_loop3[t][r_off+r] += sum(wl_trace);
		}    /* end t loop */
	      }      /* end r loop */

	      /*+ */
	      /* Do off-axis "sqrt(3)" loops in (mu,nu,-rho)-plane */
	      /*- */

	      /* Make the "corner link" in the (mu,nu,-rho) direction */
	      ftmp = 1.0 / 6.0;
	      tmp_2 = shift(ug[nu], FORWARD, mu);
	      u_tmp = ug[mu] * tmp_2;

	      tmp_2 = shift(ug[mu], FORWARD, nu);
	      u_tmp += ug[nu] * tmp_2;

	      tmp_2 = adj(ug[rho]) * u_tmp;
	      u_diag = shift(tmp_2, BACKWARD, rho);

	      tmp_2 = adj(ug[rho]) * ug[nu];
	      u_tmp = shift(tmp_2, BACKWARD, rho);

	      tmp_2 = shift(ug[rho], FORWARD, nu);
	      tmp_3 = shift(tmp_2, BACKWARD, rho);
	      u_tmp += ug[nu] * adj(tmp_3);

	      tmp_2 = shift(u_tmp, FORWARD, mu);
	      u_diag += ug[mu] * tmp_2;

	      tmp_2 = adj(ug[rho]) * ug[mu];
	      u_tmp = shift(tmp_2, BACKWARD, rho);

	      tmp_2 = shift(ug[rho], FORWARD, mu);
	      tmp_3 = shift(tmp_2, BACKWARD, rho);
	      u_tmp += ug[mu] * adj(tmp_3);

	      tmp_2 = shift(u_tmp, FORWARD, nu);
	      u_diag += ug[nu] * tmp_2;

	      u_diag = u_diag * ftmp;

	      for(r = 0;r  < ( lengthr); ++r )
	      {
		/* Gather t-links from space-directions (i.e. mu, nu and -rho) */
		/* and make "space-links" */

		if ( r == 0 )
		{
		  tmp_2 = shift(ug[j_decay], FORWARD, mu);
		  tmp_3 = shift(tmp_2, FORWARD, nu);
		  u_t   = shift(tmp_3, BACKWARD, rho);
		  u_space = u_diag;
		}
		else
		{
		  tmp_2 = shift(u_t, FORWARD, mu);
		  tmp_3 = shift(tmp_2, FORWARD, nu);
		  u_tmp = shift(tmp_3, BACKWARD, rho);
		  u_t = u_tmp;

		  tmp_2 = shift(u_space, FORWARD, mu);
		  tmp_3 = shift(tmp_2, FORWARD, nu);
		  tmp_2 = shift(tmp_3, BACKWARD, rho);
		  u_tmp = u_diag * tmp_2;
		  u_space = u_tmp;
		}

		for(t = 0;t  < ( lengtht); ++t )
		{
		  /* Gather space-link from t-direction (i.e. j_decay) */
		  if ( t == 0 )
		  {
		    u_tmp = shift(u_space, FORWARD, j_decay);
		  }
		  else
		  {
		    tmp_2 = shift(u_tmp, FORWARD, j_decay);
		    u_tmp = tmp_2;
		  }

		  tt = lsizet - t - 1;
		  btmp = t_coord < tt;

		  tmp_2 = ug[j_decay] * u_tmp;
		  tmp_3 = tmp_2 * adj(u_t);
		  copymask(tmp_3, btmp, u_tmp);
		  tmp_2 = tmp_3 * adj(u_space);
		  wl_trace = real(trace(tmp_2));
		  wils_loop3[t][r_off+r] += sum(wl_trace);

		}    /* end t loop */
	      }      /* end r loop */

	      /*+ */
	      /* Do off-axis "sqrt(3)" loops in (mu,-nu,rho)-plane */
	      /*- */

	      /* Make the "corner link" in the (mu,-nu,rho) direction */
	      ftmp = 1.0 / 6.0;
	      tmp_2 = shift(ug[rho], FORWARD, mu);
	      u_tmp = ug[mu] * tmp_2;

	      tmp_2 = shift(ug[mu], FORWARD, rho);
	      u_tmp += ug[rho] * tmp_2;

	      tmp_2 = adj(ug[nu]) * u_tmp;
	      u_diag = shift(tmp_2, BACKWARD, nu);

	      tmp_2 = adj(ug[nu]) * ug[rho];
	      u_tmp = shift(tmp_2, BACKWARD, nu);

	      tmp_2 = shift(ug[nu], FORWARD, rho);
	      tmp_3 = shift(tmp_2, BACKWARD, nu);
	      u_tmp += ug[rho] * adj(tmp_3);

	      tmp_2 = shift(u_tmp, FORWARD, mu);
	      u_diag += ug[mu] * tmp_2;

	      tmp_2 = adj(ug[nu]) * ug[mu];
	      u_tmp = shift(tmp_2, BACKWARD, nu);

	      tmp_2 = shift(ug[nu], FORWARD, mu);
	      tmp_3 = shift(tmp_2, BACKWARD, nu);
	      u_tmp += ug[mu] * adj(tmp_3);

	      tmp_2 = shift(u_tmp, FORWARD, rho);
	      u_diag += ug[rho] * tmp_2;

	      u_diag = u_diag * ftmp;

	      for(r = 0;r  < ( lengthr); ++r )
	      {
		/* Gather t-links from space-directions (i.e. mu, -nu and rho) */
		/* and make "space-links" */

		if ( r == 0 )
		{
		  tmp_2 = shift(ug[j_decay], FORWARD, mu);
		  tmp_3 = shift(tmp_2, BACKWARD, nu);
		  u_t   = shift(tmp_3, FORWARD, rho);
		  u_space = u_diag;
		}
		else
		{
		  tmp_2 = shift(u_t, FORWARD, mu);
		  tmp_3 = shift(tmp_2, BACKWARD, nu);
		  u_tmp = shift(tmp_3, FORWARD, rho);
		  u_t   = u_tmp;

		  tmp_2 = shift(u_space, FORWARD, mu);
		  tmp_3 = shift(tmp_2, BACKWARD, nu);
		  tmp_2 = shift(tmp_3, FORWARD, rho);
		  u_tmp = u_diag * tmp_2;
		  u_space = u_tmp;
		}

		for(t = 0;t  < ( lengtht); ++t )
		{
		  /* Gather space-link from t-direction (i.e. j_decay) */
		  if ( t == 0 )
		  {
		    u_tmp = shift(u_space, FORWARD, j_decay);
		  }
		  else
		  {
		    tmp_2 = shift(u_tmp, FORWARD, j_decay);
		    u_tmp = tmp_2;
		  }

		  tt = lsizet - t - 1;
		  btmp = t_coord < tt;

		  tmp_2 = ug[j_decay] * u_tmp;
		  tmp_3 = tmp_2 * adj(u_t);
		  copymask(tmp_3, btmp, u_tmp);
		  tmp_2 = tmp_3 * adj(u_space);
		  wl_trace = real(trace(tmp_2));
		  wils_loop3[t][r_off+r] += sum(wl_trace);

		}    /* end t loop */
	      }      /* end r loop */

	      /*+ */
	      /* Do off-axis "sqrt(3)" loops in (mu,-nu,-rho)-plane */
	      /*- */

	      /* Make the "corner link" in the (mu,-nu,-rho) direction */
	      ftmp = 1.0 / 6.0;

	      tmp_2 = adj(ug[nu]) * ug[mu];
	      u_tmp = shift(tmp_2, BACKWARD, nu);

	      tmp_2 = shift(ug[nu], FORWARD, mu);
	      tmp_3 = shift(tmp_2, BACKWARD, nu);
	      u_tmp += ug[mu] * adj(tmp_3);

	      tmp_2 = adj(ug[rho]) * u_tmp;
	      u_diag = shift(tmp_2, BACKWARD, rho);

	      tmp_2 = adj(ug[rho]) * ug[mu];
	      u_tmp = shift(tmp_2, BACKWARD, rho);

	      tmp_2 = shift(ug[rho], FORWARD, mu);
	      tmp_3 = shift(tmp_2, BACKWARD, rho);
	      u_tmp += ug[mu] * adj(tmp_3);

	      tmp_2 = adj(ug[nu]) * u_tmp;
	      tmp_3 = shift(tmp_2, BACKWARD, nu);
	      u_diag += tmp_3;

	      tmp_2 = shift(ug[rho], BACKWARD, rho);
	      tmp_3 = adj(ug[nu]) * adj(tmp_2);
	      u_tmp = shift(tmp_3, BACKWARD, nu);

	      tmp_2 = shift(ug[nu], BACKWARD, nu);
	      tmp_3 = adj(ug[rho]) * adj(tmp_2);
	      tmp_2 = shift(tmp_3, BACKWARD, rho);
	      u_tmp += tmp_2;

	      tmp_2 = shift(u_tmp, FORWARD, mu);
	      u_diag += ug[mu] * tmp_2;

	      u_diag = u_diag * ftmp;

	      for(r = 0;r  < ( lengthr); ++r )
	      {
		/* Gather t-links from space-directions (i.e. mu, -nu and -rho) */
		/* and make "space-links" */

		if ( r == 0 )
		{
		  tmp_2 = shift(ug[j_decay], FORWARD, mu);
		  tmp_3 = shift(tmp_2, BACKWARD, nu);
		  u_t   = shift(tmp_3, BACKWARD, rho);
		  u_space = u_diag;
		}
		else
		{
		  tmp_2 = shift(u_t, FORWARD, mu);
		  tmp_3 = shift(tmp_2, BACKWARD, nu);
		  u_tmp = shift(tmp_3, BACKWARD, rho);
		  u_t   = u_tmp;

		  tmp_2 = shift(u_space, FORWARD, mu);
		  tmp_3 = shift(tmp_2, BACKWARD, nu);
		  tmp_2 = shift(tmp_3, BACKWARD, rho);
		  u_tmp = u_diag * tmp_2;
		  u_space = u_tmp;
		}

		for(t = 0;t  < ( lengtht); ++t )
		{
		  /* Gather space-link from t-direction (i.e. j_decay) */
		  if ( t == 0 )
		  {
		    u_tmp = shift(u_space, FORWARD, j_decay);
		  }
		  else
		  {
		    tmp_2 = shift(u_tmp, FORWARD, j_decay);
		    u_tmp = tmp_2;
		  }

		  tt = lsizet - t - 1;
		  btmp = t_coord < tt;

		  tmp_2 = ug[j_decay] * u_tmp;
		  tmp_3 = tmp_2 * adj(u_t);
		  copymask(tmp_3, btmp, u_tmp);
		  tmp_2 = tmp_3 * adj(u_space);
		  wl_trace = real(trace(tmp_2));
		  wils_loop3[t][r_off+r] += sum(wl_trace);

		}    /* end t loop */
	      }      /* end r loop */
	    }        /* end i loop (for mu) */
	  }          /* end j loop (for nu) */
	}            /* end k loop (for rho) */




	dummy = 3.0 / double (Layout::vol()*Nc*2*nspace*(nspace-1)*(nspace-2)) ;
	for(t = 0;t  < ( lengtht); ++t )
	  for(r = r_off;r  < ( r_off+lengthr); ++r )
	    wils_loop3[t][r] = wils_loop3[t][r] * dummy;

      }              /* end if ( nspace > 2 ) */

      push(xml,"wils_loop3"); // XML tag for wils_wloop3
      write(xml, "lengthr", lengthr);
      write(xml, "lengtht", lengtht);
      write(xml, "length", length);
      push(xml, "wloop3");

      multi1d<Double> wloop3(wils_loop3.size2());

      for(r = 0; r < wils_loop3.size1(); ++r)
      {
        for(t = 0; t < wils_loop3.size2(); ++t)
        {
          wloop3[t] = wils_loop3[t][r];
        }
	push(xml, "elem");
        write(xml, "r", r);
        write(xml, "loop", wloop3);   // write out wils_wloop3
	pop(xml); // elem
      } // end for r

      pop(xml);        // XML end tag for wloop3
      pop(xml);        // XML end tag for wils_wloop3
      QDPIO::cout << "wils_loop3 data written to .xml file " << endl;  
    }      /* end of option "off-axis Wilson loops" */

    QDPIO::cout << "All wils_loop data written to .xml file " << endl;  

    END_CODE();
  } // end of wilslp

}  // end namespace Chroma

