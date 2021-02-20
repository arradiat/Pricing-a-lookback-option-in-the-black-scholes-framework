#ifndef PARAM_H_
#define PARAM_H_
#include <vector>
namespace lookback{



class Param {
	double T_;// maturité
	double sigma_; //volatilité"
	double T0_;  //* date de calul*/
	double opt_;  /*0 pour un call 1 pour un put*/
	double S0_;  /*prix a l'emission de l'actif sous jacent*/
	double r_; /* taux d'interet constant*/
	double N_; /* discretisation en Temps*/
	double M_; /*Nombre de répétition pour monte carlo */
	 /* parametre mu*/

public:
	/* constructeur surchargÈ crÈant
	 * une Èquation avec les paramËtres*/
	Param(double TO,double T,double opt_type,double S0,double r, double sigma, double N, double M);
	/**methode permettant de recupere le membre de donnees r*/
	double get_r();
	/**methode permettant de recupere le membre de donnees T*/
	double get_T();
	/**methode permettant de recupere le membre de donnees N*/
	double get_N();
	/**methode permettant de recupere le membre de donnees sigma*/
	double get_sigma();
	/**methode permettant de recupere le membre de donnees N*/

	/**methode permettant de recupere le membre de donnees M*/
	double get_M();
	/**methode permettant de recupere le membre de donnees TO*/
	double get_T0();
	/**methode permettant de recupere le membre de donnees S0*/
	double get_S0();
	/**methode permettant de recupere le membre de donnees opt_type*/
	double get_opt_type();

	//double get_mu();
	/**methode permettant de recupere le membre de donnees pay_off_*/

};
double mini(double * arr, int n);
double maxi(double * arr, int n);
double* simulate_S(Param p);
double monte_carlo_call(Param p);
double monte_carlo_put(Param p);
double pricer(Param p);
Param lire_csv();
double delta(Param p);
double gamma(Param p);
double theta(Param p);
double rho (Param p);
double vega(Param p);
double * evo_prix(Param p);
double * evo_delta(Param p);
void write_csv1();
void write_csv2();
double delta1(Param p);
double delta2(Param p);
double ATM_call(Param p) ;
}

#endif
