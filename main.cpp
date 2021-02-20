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
#include "calcul.h"

using namespace lookback;
using namespace std;
int main(){
	
	/*Param p(0,1,0,100,0.1,0.01,1000,1000);
	double prix=pricer(p);

	std::cout<<prix<<endl;*/

	std::cout<<"le prix du lookback put est "<<endl;

	double prix1=pricer(lire_csv());
	std::cout<<prix1<<endl;
	//double delt=delta(lire_csv());
	std::cout<<"le delta du lookback Call est "<<endl;
	std::cout<<delta(lire_csv())<<endl;

	//std::cout<<delt<<endl;
	double gam=gamma(lire_csv());
	std::cout<<"le gamma du lookback Call est "<<endl;
	std::cout<<gam<<endl;
	/*double rh=rho(lire_csv());
	std::cout<<"le rho du lookback Call est "<<endl;
	std::cout<<rh<<endl;
	double thet=theta(lire_csv());
	std::cout<<"le theta du lookback Call est "<<endl;
	std::cout<<thet<<endl;
	double veg=vega(lire_csv());
	std::cout<<"le vega du lookback Call est "<<endl;
	std::cout<<veg<<endl;
	double * prixS =evo_prix(lire_csv());
	std::cout<<"evo prix"<<endl;
	std::cout<<prixS[5]<<endl;
	double * deltaS=evo_delta(lire_csv());
	std::cout<<"evo delta"<<endl;
	std::cout<<deltaS[5]<<endl;
	*/

	    // Write the vector to CSV
	   write_csv1();
	   write_csv2();
/*
	//double a= delta1(lire_csv());

	std::cout<<a<<endl;
	std::cout<<"le delta du lookback Call est avec formule exacte derive S_T par rapport a S_t "<<endl;
	//double b= delta2(lire_csv());

	std::cout<<b<<endl;
	Param p(0,1,0,100,0.1,0.01,1000,1000);
	//std::cout<<delta1(p)<<endl;
	Param p1(0,1,0,450,0.1,0.01,1000,1000);
	//
*/
	   //std::cout<<"le delta du lookback Call quand S_T  est proche du minimum de la trajectoire"<<endl;
	   //std::cout<<delta2(lire_csv())<<endl;






	return(0);

}
