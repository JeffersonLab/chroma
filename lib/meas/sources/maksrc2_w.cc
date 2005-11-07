/*   $Id: maksrc2_w.cc,v 1.2 2005-11-07 06:28:15 edwards Exp $ ($Date: 2005-11-07 06:28:15 $) */

#error "DO NOT USE. ONLY FOR REFERENCE."

namespace Chroma {


/* This routine is specific to Wilson fermions! */

/* Create the source for the quark propagator at the site with cartesian  */
/* coordinates t_source and colour and spin component colour_source and */
/* spin_source. The result will appear in "chi". */

/* chi -- source for the quark propagator ( Write ) */
/* colour_source -- color component of the source ( Read ) */
/* spin_source -- spin component of the source ( Read ) */
/* t_source -- cartesian coordinates of the source point ( Read ) */
/* j_decay -- the direction for the perpendicular hyperplane to be filled ( Read ) */
/* wave_fn_0 -- wave function, for "shell" source (cb=0) ( Read ) */
/* wave_fn_1 -- wave function, for "shell" source (cb=1) ( Read ) */
/* wvf_type -- wave function type, for "shell" source ( Read ) */
/* wvf_param -- wvf_param of wave function ( Read ) */
/* WvfIntPar -- number of iterations to approximate Gaussian */
/*              or terminate CG inversion for Wuppertal smearing ( Read ) */
/* Z3_flag  -- 2^Z3_flag random sources/dir ( Read ) */
/* Z3_wght -- Z3 weights ( Read ) */
/* src_type -- the type of source to fill ( Read ) */
/* src_color_vec -- color vector for (gauge invariant) source ( Modify ) */
/* new_source -- flag to indicate new source ( Read ) */

include(types.mh)


  SUBROUTINE(maksrc2, u, chi, colour_source, spin_source, t_source, j_decay,
	     wave_fn_0, wave_fn_1, wvf_type, wvf_param, WvfIntPar,
	     Z3_flag, Z3_wght, src_type, src_color_vec, new_source)

  multi1d<LatticeColorMatrix> u(Nd);
LatticeFermion chi;
LatticeColorVector src_color_vec;
LatticeComplex wave_fn_0;
LatticeComplex wave_fn_1;
int colour_source;
int spin_source;
multi1d<int> t_source(Nd);
int j_decay;
int wvf_type;
Real wvf_param;
int WvfIntPar;
int src_type;
multi1d<Complex> Z3_wght(Z3_num);
int Z3_flag;
int new_source;

{ /* Local variables */
  include(COMMON_DECLARATIONS)

    LatticeFermion lrtmp;
  Complex Z3_tmp;
  Real Z3_f;
  Real one;
  Real RsdCG;
  multi1d<int> tmp_source(Nd);
  multi1d<int> inc_source(Nd);
  multi1d<int> hyp_size(Nd);
  int Z3_num;
  int idiv;
  int sum;
  int cb_source;
  multi1d<int> t_cb_source(Nd);
  int point_source;
  int shell_source;
  int gauge_inv_shell_source;
  
  /* Loop Counters */
  int cb;
  int n;
  int mu;
  int i;
  
  START_CODE("subroutine");;
  
  if ( j_decay == 0 )
    QDP_error_exit("j_decay == 0 not supported");
  
  /* For a point source, create the source at site t_cb_source with spin and colour */
  /* equal to colour_source and spin_source. */
  /* For a wall or shell source, create it in "time" slice t_source(j_decay). */
  chi = 0;
  
  /* First pass to sort things */
  switch (src_type)
  {
  case OPTION(POINT_SOURCE):
    point_source = YES;
    shell_source = NO;
    gauge_inv_shell_source = NO;

    break;
  case OPTION(WALL_SOURCE):
    point_source = NO;
    shell_source = NO;
    gauge_inv_shell_source = NO;
    /* Fill on BOTH checkerboards! */
    for(cb = 0; cb < 2; ++cb)
      walfil (chi[cb], t_source[j_decay], j_decay, colour_source, spin_source);

    break;
  case OPTION(SHELL_SOURCE):
    point_source = NO;
    shell_source = YES;

    break;
  case OPTION(BNDST_SOURCE):
    point_source = NO;
    shell_source = YES;
    
    break;
  case OPTION(POINT_AND_BNDST_SOURCE):
    point_source = NO;
    shell_source = YES;
    
    break;
  case OPTION(SHELL_AND_BNDST_SOURCE):
    point_source = NO;
    shell_source = YES;
    
    break;
  case OPTION(POINT_AND_SHELL_AND_BNDST_SOURCE):
    point_source = NO;
    shell_source = YES;
    
    break;
  default:
    QDP_error_exit("invalid type of source", src_type);
  }

  /* Now deal with point or various shell sources */
  if( shell_source == YES )
  {
    if( wvf_type == OPTION(GAUGE_INV_GAUSSIAN_WVF) ||
        wvf_type == OPTION(WUPPERTAL_WVF) )
    {
      gauge_inv_shell_source = YES;
    }
    else
    {
      gauge_inv_shell_source = NO;
      wvffil (chi[0], t_source[j_decay], j_decay, wave_fn_0, colour_source, spin_source);
      wvffil (chi[1], t_source[j_decay], j_decay, wave_fn_1, colour_source, spin_source);
    }
  }

  if( point_source == YES || 
      ( new_source == YES && gauge_inv_shell_source == YES ))
  {
    if ( Z3_flag > 0 )
    {
      Z3_num = INTEGER_LSHIFT_FUNCTION(1,Z3_flag*(Nd-1)); /* Number of Z3 sources */

                              
      idiv = INTEGER_LSHIFT_FUNCTION(1,Z3_flag);

      for(mu = 0; mu < Nd; ++mu)
      {
	if ( mu != j_decay )
	{
	  FILL(hyp_size(mu),idiv);
	  inc_source[mu] = nrow[mu] / idiv;
	}
	else
	{
	  hyp_size[mu] = 1;
	  inc_source[mu] = 0;
	}
      }

      Z3_f = FLOAT(Z3_num);
      one = 1;
      Z3_f = one / Z3_f;

      for(n = 0; n < Z3_num; ++n)
      {
	/* Move to the next site of the crystal. Mod back into the zero crystal. */
	crtesn (n, tmp_source, hyp_size, Nd);
	tmp_source = tmp_source * inc_source;
	tmp_source += t_source;
	tmp_source = mod(tmp_source, nrow);

	/* Compute checkerboard and checkerboard coordinates of the source point */
	t_cb_source[0] = tmp_source[0]/2;
	sum = tmp_source[0];
	for(i = 1; i < Nd; ++i)
	{
	  t_cb_source[i] = tmp_source[i];
	  sum = sum + tmp_source[i];
	}
	cb_source = INTEGER_MOD_FUNCTION( sum, 2 );

	/* Construct the wave function centered at this point */
	srcfil (lrtmp, t_cb_source, colour_source, spin_source);

	/* Add the wavefunction weighted by  Z3_wght(n)/Z3_num */
	Z3_tmp = Z3_wght[n] * Z3_f;
	chi[cb_source] += lrtmp * Z3_tmp;
      }

    }
    else
    {
      /* Compute checkerboard and checkerboard coordinates of the source point */
      
      t_cb_source[0] = t_source[0]/2;
      sum = t_source[0];
      for(i = 1;i  < ( Nd); ++i )
      {
	t_cb_source[i] = t_source[i];
	sum = sum + t_source[i];
      }
      cb_source = INTEGER_MOD_FUNCTION( sum, 2 );

      srcfil (chi[cb_source], t_cb_source, colour_source, spin_source);
    }
  }

  if( gauge_inv_shell_source == YES )
  {
    if( new_source == YES )
    {
      for(cb = 0; cb < 2; cb++ )
        CvToFerm (src_color_vec[cb], chi[cb], spin_source, BACKWARD);

      switch (wvf_type)
      {
      case OPTION(GAUGE_INV_GAUSSIAN_WVF):
        gaus_smear (u, src_color_vec, wvf_param, WvfIntPar, j_decay);
        break;
      case OPTION(WUPPERTAL_WVF):
        FILL(RsdCG,FUZZ);
        wupp_smear (u, src_color_vec, wvf_param, WvfIntPar, j_decay, RsdCG);
        break;
      default:
        QDP_error_exit("Unknown gauge invariant wave function", wvf_type);
      }
    }

    for(cb = 0; cb < 2; cb++ )
      CvToFerm (src_color_vec[cb], chi[cb], spin_source, FORWARD);
  }
  
  END_CODE("subroutine");;
}

}  // end namespace Chroma
