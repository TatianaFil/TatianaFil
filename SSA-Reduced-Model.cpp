//                                     CODE for
//the Reduced Model in the paper:
//Modulation of nuclear and cytoplasmic mRNA fluctuations by
//time-dependent stimuli: analytical distributions
//***************************************************************************************
#include <iostream>
#include <fstream>
#include <stdio.h> 
#include <math.h>       
#include <random>
#include <cmath>
#include <boost/math/tools/roots.hpp>
using namespace std;
using boost::math::tools::bisect;
//***************************************************************************************
// Definition of parameters:
//***************************************************************************************
// If the number of mRNA life-cycle stages is larger than L=50, change the dimensions in all the arrays
 double f [50]; 
 double b [50];  
 int hist[50][5000];  
 int nmax[50], NMAXX; 
 double sum[50], sp;
 int S, Smax ,e, i, j, NR;
 double t, tmax, sb, su, r,  x;
 double  f0, r1,r2, tau;
 double  g1, omega, fi;
 int L;
 double k[50]; 
 double p_old[50]; 
 double p[50]; 
 double sumf0;
 double sumfs;
 double bb, PP;
 int mn;
 double Int1;
 int t1;
 int qq;
 int data [50];
 int NN;
 double An[50], Fn[50];
 double SumAn, SumAn2; 

// S is the number of iterations 
// T is the size of array for discretized time T=tmax/sp

//***************************************************************************************
// Termination condition for bisection method
//***************************************************************************************
struct TerminationCondition  {
  bool operator() (double min, double max)  {
    return abs(min - max) <= 1e-12;
  }
};
//***************************************************************************************
// bisection method to find the delay time, tau
//***************************************************************************************
double Fun(double x, double rr1, double tt, double sumfss)  
{
   int nn;
   double FSsum;

   FSsum=0;
   for (nn=1; nn<=NN; nn++) {  
     FSsum=FSsum+(An[nn]/(omega*nn))*(sin((x+tt)*omega*nn+Fn[nn])-sin(tt*omega*nn+Fn[nn]));
   }

  return sumfss*x+sb*FSsum-log(1/rr1);
}; 
//***************************************************************************************
//***************************************************************************************
//***************************************************************************************
//***************************************************************************************
//***************************************************************************************






//***************************************************************************************
//                             Gillespie algorithm 
//***************************************************************************************
int main () {
 

 
// Needed for random-number generator
/////////////////////////////////////
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> unidist(0, 1);


// Open "data.csv" file for output results
//////////////////////////////////////////
  ofstream myfile; 
  myfile.open ("data.csv");  
  



// Giving values to the model parameters:
/////////////////////////////////////////
  Smax = 100;   // Total number of trajectories for the algorithm (best S>=100000)
  S = 1;        // Initialising the number of trajectories



   
// Giving values to the model parameters and time
/////////////////////////////////////////////////

  tmax=8; // final time: this is the time at which we obtain the distributions


  L = 10; // number of mRNA life-cycle stages

  
  sb=5;        // gene activation rate constant ($\sigma$ in the paper)
  bb=20;       // mean burst size ($b$ in the paper)
  g1=0.5;      // this is the $\epsilon$ parameter in the paper 
  omega=0.5;   // this is the $\omega$ parameter in the paper
  fi= M_PI/2;  // this denotes the signal phase


  k[1]=4.2; // these are the hopping rates
  k[2]=2.6;
  k[3]=2.4;
  k[4]=3.7;
  k[5]=4.8;
  k[6]=2.7;
  k[7]=2.3;
  k[8]=3.7;
  k[9]=3.0;
  k[L]=4.6;
  
  
  NN=10; // number of terms of Fourier series

  for (i=1; i<=NN; i++) {
     An[i]=0; // Initial values: these are the Fourier coefficients $A_n$ in the paper
  }

  for (i=1; i<=NN; i++) {
     Fn[i]=0; // Initial values: there are the phases in the Fourier series, $\phi_n$ in the paper
  } 

  for (i=1; i<=NN; i=i+2) {                          
     An[i]=-g1*4/(M_PI*i) ; // Evaluation
  }
  for (i=1; i<=NN; i=i+2) {
     Fn[i]=M_PI/2; // Evaluation
  }
  


// Random generator for bursty production
  PP=1./(1.+bb); //for the geometric distribution with mean bb
  std::default_random_engine generator;
  std::geometric_distribution<int> distribution(PP);  
  


// Main loop for trajectories
/////////////////////////////
  while (S <= Smax) { 

  t = 0;      // initialise time
     

  for (i=1; i<=L; i++) {
     p_old[i] = 0; // initialise the (old) number of mRNA numbers at stage i=1,...,L
     p[i] = 0; // initialise the (new) number of mRNA numbers at stage i=1,...,L
  }



 
// Main loop for time
/////////////////////
/////////////////////
  while (t < tmax) {  

     for (i=1; i<=L; i++) {
         p_old[i] = p[i];  }   // save mRNA counts (old counts==new counts)



     NR=L+1;  //number of reactions
      
     for (i=1; i<=L; i++) {
     f[i] = k[i]*p[i];} // calculate propensities f_i for i=1,...,L
     f[L+1] = sb;       // calculate propensities f_{L+1}: this is not full expression, it is defined later
 


     sumfs=0;  // sum of all the time independent terms of propensities
     for (i=1; i<=NR; i++) {
        sumfs=sumfs+f[i];  }


 
// Generate uniform random number r1 
// for the random delay time for the next reaction
     do {r1 = unidist(gen);} 
     while ((r1==0)||(r1==1));
// Generate uniform random number r2 
// for the next reaction that takes place
     do {r2 = unidist(gen);} 
     while ((r2==0)||(r2==1));
   



//We need to use numerical bisection method in order to estimate the random time, tau
    double from = 0; 
    double to = 10;    

    auto x = boost::math::tools::bisect( // call the bisection function
            [](double x){ return Fun(x, r1, t, sumfs); },
            from, to,TerminationCondition()); 
    double root = (x.first + x.second) / 2; 
    tau=root;





// Sum of all the propensity functions
     sumf0=0;  
     for (i=1; i<=NR-1; i++) {
        sumf0=sumf0+f[i]; // time-independent terms
      }

     SumAn=0;
     for (i=1; i<=NN; i++) {
     SumAn=SumAn+An[i]*cos(omega*i*(t+tau)+Fn[i]); // time-dependent terms
     }

     f0 =sumf0+ f[NR]*(1.+SumAn); // sum of all the terms
  



   
// Selection of the next reaction 
     b[1] = f[1];
     for (i=2; i<=NR-1; i++) {
        b[i] = b[i-1] + f[i]; 
     }

     SumAn2=0;
     for (i=1; i<=NN; i++) {
     SumAn2=SumAn2+An[i]*cos(omega*i*(t+tau)+Fn[i]);
     }
     b[NR]=b[NR-1]+f[NR]*(1.+SumAn2);

     for (i=1; i<=NR; i++) {
        if (b[i] >= r2*f0) {j=i; break;} // the j-th reaction happens
     }

     


     t = t+tau; // Update the time





// Update the number of molecules after j-th reaction fires   
     if (t <= tmax){
            
         
            if (j <= L-1) {p [j] = p[j] - 1;  p [j+1] = p [j+1] + 1;}// mRNA stage change
            if (j == L) {p [L] = p[L] - 1; }  // mRNA degradation
            if (j==L+1){                    // bursty transcription from geometric distribution
                 mn = distribution(generator); 
                 p[1]=p[1]+mn; 
             }    
             
            for (qq=1; qq<=L; qq++) {     
              data[qq]=p [qq];   // save the new counts in an array
            }
       
     }// end if t<=tmax loop
     
     
     if (t > tmax){
        for (qq=1; qq<=L; qq++) {

                data[qq]=p_old[qq]; // in this case the new counts remain the same as the old ones

         }
     }/// end if t>tmax loop

  }// End of main loop for time
///////////////////////////////
///////////////////////////////


// saving data for the histograms of mRNA counts at each stage, qq=1,...,L
   for (qq=1; qq<=L; qq++) {   
        hist[qq][data[qq]+1]= hist[qq][data[qq]+1]+1;     
        if (data[qq] > nmax[qq]) {nmax[qq] = data[qq];} 
   }



// Update the number of trajectories  
    if (S % 20000==0 ){
    cout<< S<<" "<<endl;} //print the number of observed trajectrory
    S = S+1;

 }// End of main loop for trajectories
//////////////////////////////////////






// Now we analyse the data obtained after the simulations
/////////////////////////////////////////////////////////

   


      
   
// Normalisation condition for the distributions of mRNA numbers
   for (qq=1; qq<=L; qq++) {
        sum[qq]=0;
        for (i=1; i<=nmax[qq]+1; i++) {
             sum[qq]=sum[qq]+hist[qq][i]/(Smax*1.0);
        } 
        cout<<qq<<","<<sum[qq]<<" "<<endl; // prints the normalization condition for all: the output needs to be 1 
   }





// Choose the maximum number of bins for printing the histograms in the output file   
   NMAXX = nmax[qq];
      for (qq=1; qq<=L; qq++) {
         if(nmax[qq]>NMAXX){NMAXX =nmax[qq];}
      }




      
// Saving the histograms for all the mRNA stages 1,2,3,...,10 
// This needs to be changed accordingly to the chosen number L
   for (i=1; i<=NMAXX+1; i++) {
       myfile<<i-1<<", "<<hist[1][i]/((Smax*1.0))<<", "<<hist[2][i]/((Smax*1.0))<<", "<<hist[3][i]/((Smax*1.0))<<", "<<hist[4][i]/((Smax*1.0))<<", "<<hist[5][i]/((Smax*1.0))<<","<<hist[6][i]/((Smax*1.0))<<", "<<hist[7][i]/((Smax*1.0))<<", "<<hist[8][i]/((Smax*1.0))<<", "<<hist[9][i]/((Smax*1.0))<<", "<<hist[10][i]/((Smax*1.0))<<"\n";
    }   
    
// The distributions can be easily processed and plotted by importing the output data.csv file in Mathematica

  myfile.close();
  return 0;
}
// END the code
//************************************************************
