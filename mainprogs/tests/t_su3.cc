// t_su3.cc, 2004/11/16 velytsky
// velytski@csit.fsu.edu
// gluodynamics (su(3)) heatbath updating with measurements
#include <iostream>
#include <cstdio>

#include "chroma.h"
#include "update/heatbath/hb_params.h"
#include "update/heatbath/su3_hb_sweep.h"

using namespace std;


HBParams::HBParams(int x1, double x2, double x3, bool x4)
{
	NmaxHB=x1;
	BetaMC=x2;
	XiBare=x3;
	xi02=XiBare*XiBare;
	AnisoP=x4;
}

int main(int argc, char *argv[]) {
	// Put the machine into a known state
	QDP_initialize(&argc, &argv);

	//set parameters: (NmaxHB, BetaMC,XiBare,AnisoP)
	HBParams hbp(-5,5.6,1.0,false);

	// Setup the layout
	const int foo[]={4,4,4,4};
	multi1d<int> nrow(Nd);
	nrow=foo;
	Layout::setLattSize(nrow);
	Layout::create();

	multi1d<LatticeColorMatrix> u(Nd);
	u=1.0; // ordered start u=diag(1,...1)

	for(int i=0;i<2;i++) {
		su3_hb_sweep(u,hbp); //one hb sweep
		//measure plaquettes and link
		Double w_plaq, s_plaq, t_plaq, link;
		MesPlq(u,w_plaq,s_plaq,t_plaq,link);
		//measure polyakov loops
		DComplex poly_loop;
		polylp(u,poly_loop,Nd-1);
		QDPIO::cout << i <<" "
                    <<w_plaq<<" "<<s_plaq<<" "<<t_plaq << " " 
		    <<link<<" "
                    << real(poly_loop) << " " <<imag(poly_loop)
		    <<endl;
	}
	
	Layout::destroy();
	QDP_finalize();
	return 1;
}
