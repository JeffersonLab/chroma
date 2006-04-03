// $Id: topol.cc,v 3.0 2006-04-03 04:58:58 edwards Exp $

#error "NOT FULLY CONVERTED"

/*  TOPOL -  Run the cooling and naive topological charge routines on a */
/*           gauge-field configuration. */

include(types.mh)

SUBROUTINE(topol, u, TopAccu, ActAccu, NumCool, NumTop)

multi1d<LatticeColorMatrix> u(Nd);
Real TopAccu;
Real ActAccu;
int NumCool;		/* Total number of cooling sweeps */
int NumTop;		/* Frequency of charge measurement */
{

include(COMMON_DECLARATIONS)

/* local variables */
multi1d<Double> qtop(NumTop);		/* Topological charge */
multi1d<Double> S_ratio(NumTop);		/* Action/continuum instanton action */

Double qtop_conv;
Double S_ratio_conv;
int conviter;
int i_cool;
int n;
Double ftmp1;
Double ftmp2;

START_CODE();


/* Do topology measurement on the initial "hot" configuration */
qnaive (u, ftmp1, ftmp2);
qtop[0] = ftmp1;
S_ratio[0] = ftmp2;

conviter = 0;
qtop_conv = 0;
S_ratio_conv = 0;

/* Do the cooling sweeps and topology measurements */
for(n=1; n<NumTop; ++n)
{
  for(i_cool=0; i_cool<NumCool; ++i_cool)
    cool (u);

  qnaive (u, ftmp1, ftmp2);
  qtop[n] = ftmp1;
  S_ratio[n] = ftmp2;

  FPRINTF(trm_out,"k= %d  qtop= %g  S_ratio= %g\n",n,qtop[n],S_ratio[n]);

  ftmp1 = qtop[n] - qtop[n-1];
  ftmp1 = fabs(ftmp1);
  ftmp2 = S_ratio[n] - S_ratio[n-1];
  ftmp2 = fabs(ftmp2);

  if ( (conviter == 0) && (ftmp1 < TopAccu) && (ftmp2 < ActAccu) )
  {
    conviter = n;
    qtop_conv = qtop[n];
    S_ratio_conv = S_ratio[n];   
    FPRINTF(trm_out,"Convergence at k= %d\n",n);
  }
}

push(xml_out,"Cooled_Topology");
write(xml_out, "NumTop", NumTop);
write(xml_out, "NumCool", NumCool);
write(xml_out, "conviter", conviter);
write(xml_out, "qtop_conv", qtop_conv);
write(xml_out, "S_ratio_conv", S_ratio_conv);
write(xml_out, "qtop", qtop);
write(xml_out, "S_ratio", S_ratio);
pop(xml_out);


END_CODE();
}
