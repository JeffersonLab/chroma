#ifndef HB_PARAMS_H
#define HB_PARAMS_H

namespace Chroma {

class HBParams 
{
	public:
		HBParams(int, double, double, bool);
		int nmax() { return NmaxHB; }
		Double beta() { return BetaMC; }
		Double xi() { return XiBare; }
		Double xi2() { return xi02; }
		bool aniso() {return AnisoP; }
		
	private:
		/**************************************************
		 * number of maximum HB tries for Creutz or KP a_0, 
		 * negative or zero value - update every single link
		 * (try infinitely long)
		 **************************************************/
		int NmaxHB;
		// MC SU(N) Beta
		Double BetaMC;
		//the bare anisotropy
		Double XiBare;
		Double xi02;
		bool AnisoP;
};

}  // end namespace Chroma

#endif
