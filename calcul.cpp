#include "calcul.h"

#include<stdlib.h>
#include<ostream>
#include <utility> // std::pair
#include <stdexcept>
#include <ctime>
#include<iostream>  // For console output
#include <vector>   // For storing random numbers
#include <cmath>    // For pow() function and other math
#include <random>   //
#include <ctime>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <string.h>
#include <cstdlib>

namespace lookback{


Param::Param(double T0,double T,double opt_type,double S0,double r, double sigma, double N, double M){

	M_=M;
	opt_=opt_type;
	r_=r;
	T_=T;
	T0_=T0;
	sigma_=sigma;
	N_=N;
	r_=r;
	S0_=S0;


}

/**methode permettant de recupere le membre de donnees r_*/
//double Param::get_mu(){
//	return (mu_);
//}
double Param::get_r(){
	return(r_);
}
/**methode permettant de recupere le membre de donnees T_*/
double Param::get_T(){
	return(T_);
}
/**methode permettant de recupere le membre de donnees N_*/
double Param::get_N(){
	return(N_);
}
/**methode permettant de recupere le membre de donnees sigma_*/
double Param::get_sigma(){
	return(sigma_);
}

/**methode permettant de recupere le membre de donnees M_*/
double Param::get_M(){
	return(M_);
}/**methode permettant de recupere le membre de donnees S0*/
double Param::get_S0(){
	return(S0_);
}/**methode permettant de recupere le membre de donnees T0*/
double Param::get_T0(){
	return(T0_);
}
double Param::get_opt_type(){
	return(opt_);

}


using namespace std;
double mini(double * arr, int n){
   /** We are assigning the first array element to
    * the temp variable and then we are comparing
    * all the array elements with the temp inside
    * loop and if the element is smaller than temp
    * then the temp value is replaced by that. This
    * way we always have the smallest value in temp.
    * Finally we are returning temp.
    */
   double temp = arr[0];
   for(int i=0; i<n; i++) {
      if(temp>arr[i]){

         temp=arr[i];
      }
   }
   return temp;
}

double maxi(double * arr, int n){
   /* We are assigning the first array element to
    * the temp variable and then we are comparing
    * all the array elements with the temp inside
    * loop and if the element is bigger than temp
    * then the temp value is replaced by that. This
    * way we always have the smallest value in temp.
    * Finally we are returning temp.
    */
   double temp = arr[0];
   for(int i=0; i<n; i++) {
      if(temp<arr[i]) {
         temp=arr[i];
      }
   }
   return temp;
}
double* simulate_S(Param p){
	double S0=p.get_S0();
	double sigma=p.get_sigma();
	//double mu=p.get_mu();
	double r=p.get_r();
	int N= (int) p.get_N();
	double T=p.get_T();
	double T0=p.get_T0();
    double * S=new double[N];
    S[0]=S0;

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0,1.0);
    double* eps=new double[N];
    for (int i=0; i<N; i++) {
        eps[i] = distribution(generator);
        }
    double h=(T-T0)/N;
    for (int i=1;i<N;i++){
        S[i]=S[i-1]*exp(sigma*sqrt(h)*eps[i]+(r-0.5*sigma*sigma)*h);
    }

    return(S);
    }

/*
double monte_carlo_call(Param p){

		int N= p.get_N();
		double T=p.get_T();
		double T0=p.get_T();
		int M=p.get_M();
		double r= p.get_r();
	double * prix_vec=new double[M];
	for (int i=0;i<M;i++){
		double *S=simulate_S(p);
		double k=mini(S,N);
		prix_vec[i]=S[N-1]-k;


	}
	double a=0;
	for (int i=0;i<M;i++){
		a=a+prix_vec[i];
	}

	double disc=exp(-r*(T-T0));
	a=a*disc;
	return(a/M);
}*/

double monte_carlo_call(Param p){
	int N= (int)p.get_N();
			double T=p.get_T();
			double T0=p.get_T0();
			int M=(int) p.get_M();
			double r= p.get_r();
    double * prix_vec=new double[M];
    for (int i=0;i<M;i++){
        double *S=simulate_S(p);
        double k=mini(S,N);
        prix_vec[i]=S[N-1]-k;

        delete[] S;

    }
    double a=0;
    for (int i=0;i<M;i++){
        a=a+prix_vec[i];
    }
    delete[] prix_vec;
    double disc=exp(-r*(T-T0));
    a=a*disc;
    return(a/M);
}
double monte_carlo_put(Param p ){

			int N= p.get_N();
			double T=p.get_T();
			double T0=p.get_T0();
			int M=p.get_M();
			double r= p.get_r();
	double * prix_vec=new double[M];
	for (int i=0;i<M;i++){
		double *S=simulate_S(p);
		double k=maxi(S,N);
		//double k=2;
		prix_vec[i]=k-S[N-1];

	}
	double a=0;
	for (int i=0;i<M;i++){
		a=a+prix_vec[i];
	}
	//delete prix_vec;
	double disc=exp(-r*(T-T0));
	a=a*disc;

	return(a/M);
}
double pricer (Param p ){
	double prix;

			int opt_type=p.get_opt_type();
	if( opt_type==0){
		// l'option est un call
		prix=monte_carlo_call( p );

	}

	else{
		//l'option est un put
		prix=monte_carlo_put(p );
	}

	return (prix);
}

Param lire_csv(){
	using namespace std;
	double * param=new double[8];
	   ifstream file("/Users/princessemame/PAP/lookback_opt/Book1.csv");
	   string line;
	   vector<double> vec;
	   while (getline(file, line)){
	       //cout << line << endl;
	       istringstream ss(line);
	       string word;
	       while (getline(ss, word,';')){
	    	   double value = strtod(word.c_str(), NULL); // cnversion string to double
	    	   vec.push_back(value); //remplissage du vecteur
	       }
	       //std::cout<<vec.size()<<endl;


	   }
	   for(int i=8;i<16;i++){
	   	    	   param[i-8]=vec.at(i);// recuperation des valeurs (9 derniers elements de vec)
	   	       }
	   Param p(param[3],param[4],param[7],param[2],param[1],param[0],param[5],param[6]);
	   vec.clear();
	   file.close();
	   return(p);
}
double delta(Param p){
	double h= 0.001;
	Param p1(p.get_T0(),p.get_T(),p.get_opt_type(),p.get_S0()+h,p.get_r(),p.get_sigma(),p.get_N(),p.get_M());

	Param p2(p.get_T0(),p.get_T(),p.get_opt_type(),p.get_S0()-h,p.get_r(),p.get_sigma(),p.get_N(),p.get_M());
	double num=pricer(p1)-pricer(p2);
	double delta= num/(2*h);

	return(delta);
}
double gamma(Param p){
	double h= 0.001;
	Param p1(p.get_T0(),p.get_T(),p.get_opt_type(),p.get_S0()+h,p.get_r(),p.get_sigma(),p.get_N(),p.get_M());

	Param p2(p.get_T0(),p.get_T(),p.get_opt_type(),p.get_S0()-h,p.get_r(),p.get_sigma(),p.get_N(),p.get_M());
	double num=delta(p1)-delta(p2);
	double gamma= num/(2*h);
	return(gamma);
}


double theta(Param p){
	double h= 0.001;
	Param p1(p.get_T0(),p.get_T()+h,p.get_opt_type(),p.get_S0(),p.get_r(),p.get_sigma(),p.get_N(),p.get_M());

	Param p2(p.get_T0(),p.get_T()-h,p.get_opt_type(),p.get_S0(),p.get_r(),p.get_sigma(),p.get_N(),p.get_M());
	double num=pricer(p1)-pricer(p2);

	double theta= num/(2*h);
	return(theta);
}
double rho(Param p){
	double h= 0.001;
	Param p1(p.get_T0(),p.get_T(),p.get_opt_type(),p.get_S0(),p.get_r()+h,p.get_sigma(),p.get_N(),p.get_M());

	Param p2(p.get_T0(),p.get_T(),p.get_opt_type(),p.get_S0(),p.get_r()-h,p.get_sigma(),p.get_N(),p.get_M());
	double num=pricer(p1)-pricer(p2);

	double rho= num/(2*h);
	return(rho);
}
double vega(Param p){
	double h= 0.001;
	Param p1(p.get_T0(),p.get_T(),p.get_opt_type(),p.get_S0(),p.get_r(),p.get_sigma()+h,p.get_N(),p.get_M());

	Param p2(p.get_T0(),p.get_T(),p.get_opt_type(),p.get_S0(),p.get_r(),p.get_sigma()-h,p.get_N(),p.get_M());
	double num=pricer(p1)-pricer(p2);

	double vega= num/(2*h);
	return(vega);
}
double * evo_prix(Param p){
	double * prix=new double[10];
	for( int i=0; i<10;i++){
		Param p1(p.get_T0(),p.get_T(),p.get_opt_type(),p.get_S0()+i*10,p.get_r(),p.get_sigma(),p.get_N(),p.get_M());
		prix[i]=pricer(p1);
	}
	return(prix);

}
double * evo_delta(Param p){
	double * delt=new double[10];
	for( int i=0; i<10;i++){
		//Param p1(p.get_T0(),p.get_T(),p.get_opt_type(),p.get_S0()+i*10,p.get_r(),p.get_sigma(),p.get_mu(),p.get_N(),p.get_M());
		Param p1(p.get_T0(),p.get_T(),p.get_opt_type(),p.get_S0()+i*10,p.get_r(),p.get_sigma(),p.get_N(),p.get_M());
		delt[i]=delta(p1);

	}
	return(delt);

}



void write_csv1(){
	 std::ofstream myfile;
	 Param p= lire_csv();
	 double prix=pricer(p);
	 double delta1=delta(p);
	 double gamma1=gamma(p);
	 double vega1=vega(p);
	 double rho1=rho(p);
	 double theta1=theta(p);
	   myfile.open ("resulta1.csv");
	   myfile << "T0;T;opt_type;S0;r;sigma;N;M;prix;delta; gamma; rho; theta;vega; \n";
	   myfile << p.get_T0()<<";";
	   myfile<<p.get_T()<<";";
	   myfile<<p.get_opt_type()<<";";
	   myfile<<p.get_S0()<<";";
	   myfile<<p.get_r()<<";";
	   myfile<<p.get_sigma()<<";";
	  // myfile<<p.get_mu()<<";";
	   myfile<<p.get_N()<<";";
	   myfile<<p.get_M()<<";";
	   myfile<<prix<<";";
	   myfile << delta1<<";";
	   myfile << gamma1<<";";
	   myfile << rho1<<";";
	   myfile << theta1<<";";
	   myfile << vega1<<";";

	   myfile.close();


}
void write_csv2(){
	 std::ofstream myfile;
	 Param p= lire_csv();
	 double* prix=evo_prix(p);
	 double* delta1=evo_delta(p);

	   myfile.open ("resulta3.csv");
	   myfile << "S0;prix;delta \n";
	   for (int i=0;i<10;i++){
		   myfile << p.get_S0()+i*10<<";";
		   myfile << prix[i]<<";";
		   myfile << delta1[i]<<"\n";
	   	   }
	   myfile.close();

}


}

