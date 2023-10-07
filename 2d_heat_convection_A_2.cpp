#include<iostream>
#include<stdio.h>
#include<math.h>
#include <cmath>
#include<string.h>
#include<iomanip>
#include<fstream>
#include<string>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <vector>


#define NI 27 //nodes in x direction
#define NJ 27 //nodes in y direction

using namespace std;


//float CELL_ZI[NI+1][NJ+1],CELL_ZETA[NI+1][NJ+1],

/////////////////////
/* geometric data */
////////////////////

float L,H; 
int i,j;
double hA, hC;

vector<double> X1;
vector<double> Y1;
double X[NJ][NI], Y[NJ][NI];

double  Y_FACE_CENTER_X[NJ][NI], X_FACE_CENTER_Y[NJ][NI];
double  X_NODE[NJ+1][NI+1], Y_NODE[NJ+1][NI+1], DXF[4][NJ+1][NI+1];//DXF node distance

double  S[4][NJ][NI]; //surface area
double  VOL[NJ][NI];  //volume 
double DELX[NJ][NI], DELY[NJ][NI];
double fx[NJ][NI],  fy[NJ][NI]; //non uniform mesh factors
//double delx,dely, ST_FACTOR_X, ST_FACTOR_Y; 
//double Rx, Ry, rx, ry, Alphax, Alphay, delxs, delys; //strech mesh parameters
//delxs smallest delx
//delys smallest dely


//////////////////////////
/* heat convection data */
//////////////////////////

vector<double> u; //u given
vector<double> v; //v given
double UA, UB, UC, VD;
double FA, TA, TL ; 

double U[2][NJ+1][NI+1];//U-field given here U[0][j][i] means u and U[1][j][i] means v
double U_FACE[2][NJ+1][NI+1];//U-field given face velocity

double Pe[4][NJ+1][NI+1]; //peclet number (0 means west, 1 means south, 2 means east and 3 means north)
double Mw, Me, Ms, Mn ;  //
double Fw, Fe, Fs, Fn ; //convective coefficient 
double AW[NJ+1][NI+1], AS[NJ+1][NI+1], AE[NJ+1][NI+1], AN[NJ+1][NI+1], AP[NJ+1][NI+1];
double SP[NJ+1][NI+1]; 
double T_OLD[NJ+1][NJ+1],T[NJ+1][NI+1];

//ax and ay are convective heat flux in x and y direction
//qx and qy are conductive heat flux in x and y direction

double ax[NJ+1][NI+1], ay[NJ+1][NI+1], qx[NJ+1][NI+1],qy[NJ+1][NI+1];
double T_YFACE[NJ+1][NI+1], T_XFACE[NJ+1][NI+1];
double k, Cp, alpha, rho, Q_VOL_GEN;  //properties
double Su; //source term zero here


///////////////////////////
/* convergence and error */
///////////////////////////

double RESIDUE[NJ+1][NI+1];
double MAX_ERROR, ABS_ERROR, RMS_ERROR, RRESIDUE;
int ITER = 0;

void SET_GEOMETRY_UNIFORM();
void SET_GEOMETRY_NONUNIFORM(); //here non uniform mesh is given
void APPLY_IC();
void APPLY_BC();
void UPDATE();
void CALC_CONV(); //equation solved using gauss seidal
//double CALC_ABS_ERROR();
double CALC_CONV_ERROR();

//TDMA algorithm
double a[NJ+1][NI+1], b[NJ+1][NI+1], c[NJ+1][NI+1], d[NJ+1][NI+1];
double P[NJ+1][NI+1], Q[NJ+1][NI+1];
void CALC_CONV_VTDMA(); //equation solved using vertical TDMA
void CALC_CONV_HTDMA(); //equation solved using horizontal TDMA


int FILE_WRITE1(); //writing post processing file in .dat format 
//int FILE_WRITE2();
void FILE_READ(); //reading geometric and u and v data
void INTERPOLATE(); // interpolating u and v data from node to face center (convective part)

double A;

int main ()
{
	L = 2; H = 2;
	hA  = 0.068*H;
	hC  = 0.068*H;
	rho  = 1;
	alpha = 0.002 ;//(k/Cp);
    A = 2;

	UA =  1;
	UB = 0;
	UC = 1;
	VD = 0;
	TA = 20;
	TL = 10;
	Q_VOL_GEN = 0;
	
	FA = rho*UA*hA*(TA - TL); //inlet mass flux * temp. different
	
//	ST_FACTOR_Y = 1.1;
//	ST_FACTOR_X = 1.08;
	MAX_ERROR = 1e-04;	
	
    FILE_READ();	
	
	SET_GEOMETRY_NONUNIFORM();
	INTERPOLATE();
    APPLY_BC();
    APPLY_IC();
 
    
    UPDATE();
    ofstream file3("ITER_VS_RESIDUE_04_HTDMA.dat");
ITER:

   ++ITER;
       
//  CALC_CONV();
//  CALC_CONV_VTDMA();
    CALC_CONV_HTDMA();
    
    APPLY_BC();
    
    RRESIDUE = CALC_CONV_ERROR();
//  RRESIDUE = CALC_CONV_HTDMA_ERROR();

    
    file3<<ITER<<"\t"<<RRESIDUE<<"\n";
//    cout<<ITER<<" "<<RRESIDUE<<" "<<endl;

	if(RRESIDUE > MAX_ERROR){
		UPDATE();
		goto ITER;
	}
    
    file3.close();
   
	FILE_WRITE1();
     

}

void FILE_READ(){
	
	double xdata, ydata, udata, vdata;
	int k;
	k = 0;
	
	ifstream file1("xc.dat");
	ifstream file2("yc.dat");
	ifstream file3("u.dat");
	ifstream file4("v.dat");
	
	while ( file1 >> xdata)
	{
		X1.push_back(xdata);
	}
	
	while( file2 >> ydata)
	{
		Y1.push_back(ydata);
	}
	
	while ( file3 >> udata)
	{
		u.push_back(udata);
	}
	
	while ( file4 >> vdata)
	{
		v.push_back(vdata);
	}
	
/*	for(int i = 0; i <= X1.size()-1;i++)
	{
	//cout << X1[i] <<" "<<Y1[i]<< endl;
	}
*/		
    for(i=1;i<=NI;i++)
	{
	  for(j=1;j<=NJ;j++)
	  {
         U[0][j][i] = u[k];
         U[1][j][i] = v[k];
	//	 cout<<U[0][i][j]<<""<<" "<<U[1][i][j]<<" "<<k<<endl;
          k  = k + 1 ; 
	  }
    }


	file1.close();
	file2.close();
	file3.close();
	file4.close();
		
}


void SET_GEOMETRY_NONUNIFORM()
{
	
/*	                     // in compuational domain generates 1x1 squre domain.
	delx = (L/(NI-2));    // X is horizontal direction in physical domain
	dely = (H/(NJ-2));  // Y is vertical direction in physical domain
	delv = delx * dely ;
//	cout<<delv<<endl ;
//	cout<<delx<<" "<<dely<<endl;
*/
  for(i=1;i<=NI-1;i++)
	{
	for(j=1;j<=NJ-1;j++)
	{
	    X[j][i]= X1[i-1];
		Y[j][i]= Y1[j-1];
//	   cout<<X[i][j]<<" "<<Y[i][j]<<endl;
	    } 
	}	

  for(i=1;i<=NI-2;i++)
	{
	for(j=1;j<=NJ-2;j++)
	 {	
		DELX[j][i] = X[j][i+1]-X[j][i];
		DELY[j][i] = Y[j+1][i] - Y[j][i];
     }
   }
		
   for(i=2;i<=NI-1;i++)
	{
	  for(j=2;j<=NJ-1;j++)
	  {
		X_NODE[j][i] = (X[j][i]+X[j][i-1])/2;
		Y_NODE[j][i] = (Y[j][i]+Y[j-1][i])/2;       
	  } 
	}	


for(j=2;j<=NJ-1;j++) // left and right boundary
	  {
	  	X_NODE[j][1] = 0;
	  	Y_NODE[j][1] = Y_NODE[j][2];
	  	X_NODE[j][NI] = L;
	  	Y_NODE[j][NI] =Y_NODE[j][NI-1];
	  }
for(i=2;i<=NI-1;i++) //top and bottom boundary
	  {
	  	X_NODE[1][i] = X_NODE[2][i];
	  	Y_NODE[1][i] = 0;
	  	X_NODE[NJ][i] = X_NODE[NJ-1][i];
	  	Y_NODE[NJ][i] = H;
	  }
	  
	  X_NODE[1][1] = 0; Y_NODE[1][1] = 0;  //four corner of mesh
	  X_NODE[1][NI] = L; Y_NODE[1][NI]=0;
	  X_NODE[NJ][1] = 0; Y_NODE[NJ][1] = H; 
	  X_NODE[NJ][NI] = L; Y_NODE[NJ][NI] = H;

//face center
for(i=1;i<=NI-1;i++)
	{
	for(j=2;j<=NJ-1;j++)
	{
		Y_FACE_CENTER_X[j][i]= X[j][i];
	//	cout<<Y_FACE_CENTER_X[j][i]<<endl;
	    } 
	}	

for(i=2;i<=NI-1;i++)
  {	
   for(j=1;j<=NJ-1;j++)
	{
        X_FACE_CENTER_Y[j][i]= Y[j][i];	
    //  cout<<X_FACE_CENTER_Y[j][i]<<endl;
	    } 
	}	
	  
//NODE DISTANCE

for(i=2;i<=NI-1;i++)
 {
	for(j=2;j<=NJ-1;j++)
	{
		DXF[0][j][i] = (X_NODE[j][i]-X_NODE[j][i-1]); //WEST
		DXF[1][j][i] = (Y_NODE[j][i]-Y_NODE[j-1][i]); //SOUTH
		DXF[2][j][i] = (X_NODE[j][i+1]-X_NODE[j][i]); //EAST
		DXF[3][j][i] = (Y_NODE[j+1][i]-Y_NODE[j][i]); //NORTH
	  //  cout<<DXF[0][j][i]<<" "<<DXF[1][j][i]<<" "<<DXF[2][j][i]<<" "<<DXF[3][j][i]<<endl;
	  //  cout<<"\n"<<endl;
		}    
    }
		
//SURFACE AREA

for(i=2;i<=NI-1;i++)
 {
	for(j=2;j<=NJ-1;j++)
	{
		S[0][j][i] = (Y[j][i-1] - Y[j-1][i-1]);//WEST FACE
		S[1][j][i] = (X[j-1][i] - X[j-1][i-1]);	//SOUTH FACE 
	    S[2][j][i] = (Y[j][i] - Y[j-1][i]); //EAST FACE
	    S[3][j][i] = (X[j][i] - X[j][i-1]);    //NORTH FACE
	    
	    VOL[j][i] = S[2][j][i] * S[3][j][i] ;
	//    cout<<S[2][j][i]<<" "<<S[3][j][i]<<endl;     
        }    
    }
				

}  

void INTERPOLATE(){
	
		double m, mass[NJ+1][NI+1];
	double m1, m2, mass1[NI+1], mass2[NJ+1];

	for(i=2;i<=NI-1;i++)
 {
	for(j=2;j<=NJ-1;j++)
	{
       fx[j][i] = 0.5*DELX[j-1][i-1]/(X_NODE[j][i+1]-X_NODE[j][i]);
	   fy[j][i] = 0.5*DELY[j-1][i-1]/(Y_NODE[j+1][i]-Y_NODE[j][i]);
    }
 }
 
 //y-face velocity (in x direction)
    for(i=1;i<=NI-1;i++)
	{
	  for(j=2;j<=NJ-1;j++)
	  {	
	     if (i==1)
	     {
	     	U_FACE[0][j][1] = U[0][j][1] ;
		 }
		 else if (i==NI-1){
		 	U_FACE[0][j][NI-1] = U[0][j][NI] ;
		 }
		 else{
         U_FACE[0][j][i] = fx[j][i]*U[0][j][i] + (1-fx[j][i])*U[0][j][i+1];	  		
	     }
	  //   cout<<U_FACE[0][j][i]<<endl;
	  }
    }

//x-face velocity (in y direction)
    for(i=2;i<=NI-1;i++)
	{
	  for(j=1;j<=NJ-1;j++)
	  {	
	     if (j==1){
	     	U_FACE[1][1][i] = U[1][1][i];
		 }
		 else if (j==NJ-1){
		 	U_FACE[1][NJ-1][i] = U[1][NJ][i];
		 }
		 else{
         U_FACE[1][j][i] = fy[j][i]*U[1][j][i] + (1-fy[j][i])*U[1][j+1][i];	  		
	     }
	//     cout<<U_FACE[1][j][i]<<endl;
	  }
    }
    

}
void APPLY_BC(){
	
   	 
	 for(j=1;j<=NJ;j++) // Right boundary condition
	{
	   if((Y_NODE[j][NI] > hA) && (Y_NODE[j][NI] <= H))
		 {
		T[j][NI] = 10;//right 10
	     }
	   else
	   {
	   	T[j][NI] = T[j][NI-1] ;
		 }  
	}

	 for(j=1;j<=NJ;j++) // left boundary condition
	{
		if((Y_NODE[j][1] > (H-hA)) && (Y_NODE[j][1] <= H))
		 {
		 T[j][1] = 20; //left
//		 U[0][1][j] = UA; //left	
//		 U[1][1][j] = 0; //left	
	     }
	    else
		 {
		 T[j][1] = T[j][2]; //left
//		 U[0][1][j] = 0; //left	
//		 U[1][1][j] = 0; //left			 	
		  } 
	}

/*
	 for(j=1;hC<j<=(NJ);j++) // Right boundary condition
	{
		U[0][NI][j] = UC; //Right
		U[1][NI][j] = 0; //Right		
	}
	
	for(i=1;i<=NI;i++) // top and bottom boundary condition
	{
		U[1][i][NJ] = VD ; //top
	}
*/	
	for(i=1;i<=NI;i++) // top and bottom boundary condition
	{
	//	T[NJ][i] = T[NJ-1][i] ; //top
	   T[NJ][i] = 15 ; //top
		T[1][i] = T[2][i] ; //bottom
	}


}

void APPLY_IC(){
	
	for(i=2;i<=NI-1;i++)
	{
		for(j=2;j<=NJ-1;j++) // interior nodes
	   {	
			T[j][i] = 5;
		}
	}
}

void CALC_CONV(){

  for(i=2;i<=NI-1;i++)
   {
      for(j=2;j<=NJ-1;j++)
	{

/*
	  kw[i][j] = 5 * (1 + (100/L) *Y_FACE_CENTER_X[i-1][j]);
      ks[i][j] = 5 * (1 + (100/L) *0.5*(Y_FACE_CENTER_X[i][j] + Y_FACE_CENTER_X[i-1][j]));
      ke[i][j] = 5 * (1 + (100/L) *Y_FACE_CENTER_X[i][j]);
      kn[i][j] = 5 * (1 + (100/L) *0.5*(Y_FACE_CENTER_X[i][j] + Y_FACE_CENTER_X[i-1][j]));
*/


//Pe[0] is for west face
// Pe[1] is for south face
// Pe[2] is for east face
// Pe[3] is for north face
// here all for west, east, south and north face, peclet number is stored on particular nodes.
      Pe[0][j][i] = (rho*U_FACE[0][j][i-1]*DXF[0][j][i])/alpha; 
      Pe[1][j][i] = (rho*U_FACE[1][j-1][i]*DXF[1][j][i])/alpha; 
      Pe[2][j][i] = (rho*U_FACE[0][j][i]*DXF[2][j][i])/alpha; 
      Pe[3][j][i] = (rho*U_FACE[1][j][i]*DXF[3][j][i])/alpha; 

//here if particular peclet number is less than 2 then, CD scheme is used,
// Otherwise the upwind scheme is used.

      //west face
      if (abs(Pe[0][j][i]) <= 2){ 
      	Fw = 0.5*rho*U_FACE[0][j][i-1]*S[0][j][i] ;
	  }
	  
	  else{
	  	Fw = fmax((rho*U_FACE[0][j][i-1]*S[0][j][i]),0) ;
	  }
      
      //south face
      if (abs(Pe[1][j][i]) <= 2){
      	Fs = 0.5*rho*U_FACE[1][j-1][i]*S[1][j][i] ;
	  }
	  
	  else{
	  	Fs = fmax((rho*U_FACE[1][j-1][i]*S[1][j][i]),0) ;
	  }
      		
      //east face
      if (abs(Pe[2][j][i]) <= 2){
      	Fe = -0.5*rho*U_FACE[0][j][i]*S[2][j][i];
	  }
	  
	  else{
	  	Fe = fmax(0,(-rho*U_FACE[0][j][i]*S[2][j][i])) ;
	  }

      //north face
      if (abs(Pe[3][j][i]) <= 2){
      	Fn = -0.5*rho*U_FACE[1][j][i]*S[3][j][i] ;
	  }
	  
	  else{
	  	Fn = fmax(0,(-rho*U_FACE[1][j][i]*S[3][j][i])) ;
	  }
         
      Su = Q_VOL_GEN*VOL[j][i]; //here Q_VOL_GEN is zero 
	  SP[i][j] = Q_VOL_GEN*VOL[j][i]/(T_OLD[j][i]);
	    
	  AW[j][i] = (alpha * S[0][j][i] / (DXF[0][j][i])) + Fw;
	  AS[j][i] = (alpha * S[1][j][i] / (DXF[1][j][i])) + Fs;
	  AE[j][i] = (alpha *S[2][j][i] / (DXF[2][j][i])) + Fe;
	  AN[j][i] = (alpha * S[3][j][i] / (DXF[3][j][i])) + Fn;
	  AP[j][i] = AE[j][i] + AW[j][i] + AN[j][i] + AS[j][i] - SP[j][i];

     // cout<<kw[j][i]<<" "<<ks[j][i]<<" "<<ke[j][i]<<" "<<kn[j][i]<<endl;
//     cout<<AW[j][i]<<" "<<AS[j][i]<<" "<<AE[j][i]<<" "<<AN[j][i]<<" "<<AP[j][i]<<endl;
		}	
	}
		
  for(i=2;i<=NI-1;i++)
	{
	 for(j=2;j<=NJ-1;j++)
	 {
	T[j][i] = (AE[j][i] * T[j][i+1] + AW[j][i] * T[j][i-1] + AN[j][i] * T[j+1][i] + AS[j][i] * T[j-1][i] - Su) / AP[j][i];			
	 }
	}
		
}


void CALC_CONV_VTDMA(){

  for(i=2;i<=NI-1;i++)
   {
      for(j=2;j<=NJ-1;j++)
	{

/*
	  kw[i][j] = 5 * (1 + (100/L) *Y_FACE_CENTER_X[i-1][j]);
      ks[i][j] = 5 * (1 + (100/L) *0.5*(Y_FACE_CENTER_X[i][j] + Y_FACE_CENTER_X[i-1][j]));
      ke[i][j] = 5 * (1 + (100/L) *Y_FACE_CENTER_X[i][j]);
      kn[i][j] = 5 * (1 + (100/L) *0.5*(Y_FACE_CENTER_X[i][j] + Y_FACE_CENTER_X[i-1][j]));
*/

      Pe[0][j][i] = (rho*U_FACE[0][j][i-1]*DXF[0][j][i])/alpha; 
      Pe[1][j][i] = (rho*U_FACE[1][j-1][i]*DXF[1][j][i])/alpha; 
      Pe[2][j][i] = (rho*U_FACE[0][j][i]*DXF[2][j][i])/alpha; 
      Pe[3][j][i] = (rho*U_FACE[1][j][i]*DXF[3][j][i])/alpha; 

      //west face
      if (abs(Pe[0][j][i]) <= 2){
      	Fw = 0.5*rho*U_FACE[0][j][i-1]*S[0][j][i] ;
	  }
	  
	  else{
	  	Fw = fmax((rho*U_FACE[0][j][i-1]*S[0][j][i]),0) ;
	  }
      
      //south face
      if (abs(Pe[1][j][i]) <= 2){
      	Fs = 0.5*rho*U_FACE[1][j-1][i]*S[1][j][i] ;
	  }
	  
	  else{
	  	Fs = fmax((rho*U_FACE[1][j-1][i]*S[1][j][i]),0) ;
	  }
      		
      //east face
      if (abs(Pe[2][j][i]) <= 2){
      	Fe = -0.5*rho*U_FACE[0][j][i]*S[2][j][i];
	  }
	  
	  else{
	  	Fe = fmax(0,(-rho*U_FACE[0][j][i]*S[2][j][i])) ;
	  }

      //north face
      if (abs(Pe[3][j][i]) <= 2){
      	Fn = -0.5*rho*U_FACE[1][j][i]*S[3][j][i] ;
	  }
	  
	  else{
	  	Fn = fmax(0,(-rho*U_FACE[1][j][i]*S[3][j][i])) ;
	  }
         
      Su = Q_VOL_GEN*VOL[j][i];
	  SP[i][j] = Q_VOL_GEN*VOL[j][i]/(T_OLD[j][i]);
	    
	  AW[j][i] = (alpha * S[0][j][i] / (DXF[0][j][i])) + Fw;
	  AS[j][i] = (alpha * S[1][j][i] / (DXF[1][j][i])) + Fs;
	  AE[j][i] = (alpha *S[2][j][i] / (DXF[2][j][i])) + Fe;
	  AN[j][i] = (alpha * S[3][j][i] / (DXF[3][j][i])) + Fn;
	  AP[j][i] = AE[j][i] + AW[j][i] + AN[j][i] + AS[j][i] - SP[j][i];

     // cout<<kw[j][i]<<" "<<ks[j][i]<<" "<<ke[j][i]<<" "<<kn[j][i]<<endl;
//     cout<<AW[j][i]<<" "<<AS[j][i]<<" "<<AE[j][i]<<" "<<AN[j][i]<<" "<<AP[j][i]<<endl;
		}	
	}

//a,b,c and d are coefficient for TDMA algorithm

  for(i=2;i<=NI-1;i++)
	{
	 for(j=2;j<=NJ-1;j++)
	 {
	   a[j][i] = AP[j][i];
	   b[j][i] = AN[j][i];
	   c[j][i] = AS[j][i];
	   d[j][i] = (AE[j][i]*T[j][i+1] + AW[j][i]*T[j][i-1] + Su);
	
		 if(j==2) {
		 P[2][i] = b[2][i]/a[2][i] ;
		 Q[2][i]= (d[2][i] + c[2][i]*T[1][i])/a[2][i]; 
	     }
	     else
	     {
		 P[j][i] = b[j][i]/(a[j][i] - c[j][i]*P[j-1][i]);
		 Q[j][i] = (c[j][i] * Q[j-1][i] + d[j][i])/(a[j][i] - c[j][i]*P[j-1][i]); 
		 }
	    T[j][i] = P[j][i] * T[j+1][i] + Q[j][i];
     }
 }
	
}

void CALC_CONV_HTDMA(){

  for(j=2;j<=NJ-1;j++)
   {
      for(i=2;i<=NI-1;i++)
	{

/*
	  kw[i][j] = 5 * (1 + (100/L) *Y_FACE_CENTER_X[i-1][j]);
      ks[i][j] = 5 * (1 + (100/L) *0.5*(Y_FACE_CENTER_X[i][j] + Y_FACE_CENTER_X[i-1][j]));
      ke[i][j] = 5 * (1 + (100/L) *Y_FACE_CENTER_X[i][j]);
      kn[i][j] = 5 * (1 + (100/L) *0.5*(Y_FACE_CENTER_X[i][j] + Y_FACE_CENTER_X[i-1][j]));
*/

      Pe[0][j][i] = (rho*U_FACE[0][j][i-1]*DXF[0][j][i])/alpha; 
      Pe[1][j][i] = (rho*U_FACE[1][j-1][i]*DXF[1][j][i])/alpha; 
      Pe[2][j][i] = (rho*U_FACE[0][j][i]*DXF[2][j][i])/alpha; 
      Pe[3][j][i] = (rho*U_FACE[1][j][i]*DXF[3][j][i])/alpha; 

      //west face
      if (abs(Pe[0][j][i]) <= 2){
      	Fw = 0.5*rho*U_FACE[0][j][i-1]*S[0][j][i] ;
	  }
	  
	  else{
	  	Fw = fmax((rho*U_FACE[0][j][i-1]*S[0][j][i]),0) ;
	  }
      
      //south face
      if (abs(Pe[1][j][i]) <= 2){
      	Fs = 0.5*rho*U_FACE[1][j-1][i]*S[1][j][i] ;
	  }
	  
	  else{
	  	Fs = fmax((rho*U_FACE[1][j-1][i]*S[1][j][i]),0) ;
	  }
      		
      //east face
      if (abs(Pe[2][j][i]) <= 2){
      	Fe = -0.5*rho*U_FACE[0][j][i]*S[2][j][i];
	  }
	  
	  else{
	  	Fe = fmax(0,(-rho*U_FACE[0][j][i]*S[2][j][i])) ;
	  }

      //north face
      if (abs(Pe[3][j][i]) <= 2){
      	Fn = -0.5*rho*U_FACE[1][j][i]*S[3][j][i] ;
	  }
	  
	  else{
	  	Fn = fmax(0,(-rho*U_FACE[1][j][i]*S[3][j][i])) ;
	  }
         
      Su = Q_VOL_GEN*VOL[j][i];
	  SP[i][j] = Q_VOL_GEN*VOL[j][i]/(T_OLD[j][i]);
	    
	  AW[j][i] = (alpha * S[0][j][i] / (DXF[0][j][i])) + Fw;
	  AS[j][i] = (alpha * S[1][j][i] / (DXF[1][j][i])) + Fs;
	  AE[j][i] = (alpha *S[2][j][i] / (DXF[2][j][i])) + Fe;
	  AN[j][i] = (alpha * S[3][j][i] / (DXF[3][j][i])) + Fn;
	  AP[j][i] = AE[j][i] + AW[j][i] + AN[j][i] + AS[j][i] - SP[j][i];

     // cout<<kw[j][i]<<" "<<ks[j][i]<<" "<<ke[j][i]<<" "<<kn[j][i]<<endl;
//     cout<<AW[j][i]<<" "<<AS[j][i]<<" "<<AE[j][i]<<" "<<AN[j][i]<<" "<<AP[j][i]<<endl;
		}	
	}

//a,b,c and d are coefficient for TDMA algorithm


  for(j=2;j<=NJ-1;j++)
	{
	 for(i=2;i<=NI-1;i++)
	 {
	   a[j][i] = AP[j][i];
	   b[j][i] = AE[j][i];
	   c[j][i] = AW[j][i];
	   d[j][i] = (AN[j][i]*T[j+1][i] + AS[j][i]*T[j-1][i] + Su);

		 if(i==2) {
		 P[j][i] = b[j][2]/a[j][2] ;
		 Q[j][i] = (d[j][2] + c[j][2]*T[j][1])/a[j][2]; 
	     }
	     else
	     {
		 P[j][i] = b[j][i]/(a[j][i] - c[j][i]*P[j][i-1]);
		 Q[j][i] = (c[j][i] * Q[j][i-1] + d[j][i])/(a[j][i] - c[j][i]*P[j][i-1]); 
		 }
	    T[j][i] = P[j][i] * T[j][i+1] + Q[j][i];

 //       cout<<T_OLD[j][i]<<endl;
//      cout<<a[j][i]<<" "<<b[j][i]<<" "<<c[j][i]<<" "<<d[j][i]<<endl;	    
//	    cout<<T[j][i]<<" "<<P[j][i]<<" "<<Q[j][i]<<endl;
     }
 }
 

	
}

void UPDATE(){
  for (i=2;i<=NI-1;i++)
	{	
	for (j=2;j<=NJ-1;j++)
	{	
	   T_OLD[j][i] = T[j][i];
      }
    }
}

double CALC_ABS_ERROR()
{	
/*	for (j=1;j<=NI;j++)
	{
		for (i=1;i<=NI;i++)
		{	
		  RESIDUE[j][i] =(T[j][i] - T_OLD[j][i]);			
		}
	}
	
//	RTOT = sqrt(RTOT/(NI*NI));
//	return RTOT;
*/
}


double CALC_CONV_ERROR()
{	
    double RTOT = 0;
    double TA, TC;
    
   for (i=2;i<=NI-1;i++)
  {
	  for (j=2;j<=NJ-1;j++)
	{
	  RESIDUE[j][i] = abs(AE[j][i] * T[j][i+1] + AW[j][i] * T[j][i-1] + AN[j][i] * T[j+1][i] + AS[j][i] * T[j-1][i] + Su - AP[j][i]*T[j][i]); 			
      RTOT=RTOT+RESIDUE[j][i];			
	}
  }
	RTOT = RTOT/FA;
	return RTOT;
}



int FILE_WRITE1(){  // how to write the data in the file ...

	
	//flux in x-direction
	//qx is conductive heat flux
	//ax is convective heat flux
	
	for (i=1;i<=NI-1;i++)
	{
	   for (j=2;j<=NJ-1;j++)	
	{
	 	qx[j][i] = -alpha*(T[j][i+1] - T[j][i])/(X_NODE[j][i+1] - X_NODE[j][i]);
		ax[j][i] = rho*U_FACE[0][j][i]*(fx[j][i]*T[j][i] + (1-fx[j][i])*T[j][i+1]);  
	    }
     }
	
	//flux in y-direction
	//qy is conductive heat flux
	//ay is convective heat flux
	
	for (i=2;j<=NI-1;i++)	
	{
	   for (j=1;j<=NJ-1;j++)
		{
		 qy[j][i] = -alpha*(T[j+1][i] - T[j][i])/(Y_NODE[j+1][i] - Y_NODE[j][i]);	
	     ay[j][i] = rho*U_FACE[1][j][i]*(fy[j][i]*T[j][i] + (1-fy[j][i])*T[j+1][i]);  
	    }
    }

/*
	for (i=1;i<=NI-1;i++)
	{
	   for (j=2;j<=NJ-1;j++)	
	{
          T_YFACE[j][i] = 0.5*(T[j][i] + T[j][i-1]); 
	    }
     }
     
    for (i=2;j<=NI-1;i++)	
	{
	   for (j=1;j<=NJ-1;j++)
		{
           T_XFACE[j][i] = 0.5*(T[j][i] + T[j-1][i]); 
	    }
    }
 
*/

double sum, sum1, sum2, sum3, sum4, sumbb3, sumbb4;
double suma1[NI+1], suma2[NI+1], sumb1[NJ+1], sumb2[NJ+1], sumb3[NJ+1], sumb4[NJ+1];

cout<<"\n"<<endl;

	for (i=2;i<=NI-1;i++)
	{
       suma1[i] = rho*(-U_FACE[1][1][i]*T[1][i]*S[1][2][i] + U_FACE[1][NJ-1][i]*T[NJ][i]*S[3][NJ-1][i]) ;
	   suma2[i]= -alpha*((T[2][i] - T[1][i])/(Y_NODE[1][i] - Y_NODE[2][i]))*S[1][2][i] + alpha*((T[NJ][i] - T[NJ-1][i])/(Y_NODE[NJ][i] - Y_NODE[NJ-1][i]))*S[3][NJ-1][i] ;
//	   cout<<(alpha*((T[2][i] - T[1][i])/(Y_NODE[1][i] - Y_NODE[2][i]))) + (alpha*((T[NJ][i] - T[NJ-1][i])/(Y_NODE[NJ][i] - Y_NODE[NJ-1][i])))<<endl;
//     cout<<rho*(U[1][1][i]*T[1][i] + U[1][NJ][i]*T[NJ][i])<<endl;
//     suma3[i] = 
	   sum1 = sum1 + suma1[i];
	   sum2 = sum2 + suma2[i];
    }

	for (j=2;j<=NJ-1;j++)
	{   
       sumb1[j] = rho*(U_FACE[0][j][1]*T[j][1]*S[0][j][2] - U_FACE[0][j][NI-1]*T[j][NI]*S[2][j][NI-1]) ;
	   sumb2[j]= -alpha*((T[j][2] - T[j][1])/(X_NODE[j][1] - X_NODE[j][2])) *S[0][j][2] + alpha*((T[j][NI] - T[j][NI-1])/(X_NODE[j][NI] - X_NODE[j][NI-1]))*S[2][j][NI-1];
       sum3 = sum3 + sumb1[j] ;
       sum4 = sum4 + sumb2[j] ;
       
       if((Y_NODE[j][1] > (H-hA)) && (Y_NODE[j][1] <= H))
		 {
		   sumb3[j] = rho*U_FACE[0][j][1]*T[j][1] ;
		   sumbb3 = sumbb3 + sumb3[j];
	     }
	     
	  if((Y_NODE[j][1] >= 0) && (Y_NODE[j][1] <= hA))
		 {
		   sumb4[j] = rho*U_FACE[0][j][NI-1]*T[j][NI] ;
		   sumbb4 = sumbb4 + sumb4[j];
	     }   
    }
     sum = (sum1 + sum3) + (sum2 + sum4);
     cout<<"\n"<<endl;
	 cout<<"different heat flux through the boundaries :"<<endl;
	 cout<<"x-direction convective flux = "<<sum3<<" "<<"x-direction conductive flux = "<<sum4<<endl;
     cout<<"\n"<<endl;
	 cout<<"y-direction convective flux = "<<sum1<<" "<<"y-direction conductive flux = "<<sum2<<endl;
	 cout<<"\n"<<endl;
	 cout<<"total heat flux = "<<sum<<endl;
     cout<<"\n"<<endl;
     cout<<"inlet convective flux = "<<sumbb3<<endl;
	 cout<<"outlet convective flux = "<<sumbb4<<endl;
     
	//grid line
	ofstream file5("grid_data_04_HTDMA.dat");
    file5<<"VARIABLES = "<<'"'<<"X"<<'"'<<", "<<'"'<<"Y"<<'"'<<endl;
    file5<<"ZONE I="<<NI-1<<", J="<<NJ-1<<", F=POINT"<<endl;

	for (i=2;i<=NI-1;i++)
	{
		for (j=2;j<=NI-1;j++)	
		{
//       cout<<AW[j][i]<<" "<<AS[j][i]<<" "<<AE[j][i]<<" "<<AN[j][i]<<" "<<AP[j][i]<<endl;
		}
	}



	for (i=1;i<=NI-1;i++)
	{
		for (j=1;j<=NI-1;j++)	
		{
     		file5<<X[j][i]<<" "<<Y[j][i]<<endl;
		}
	}
	file5.close();

    //temperature values at CV and heat flux (face values)
	ofstream writer("structural_grid_04_HTDMA.dat");
	if (!writer){
		cout<<"error"<<endl;
	}
	else{
		writer<<"VARIABLES = "<<'"'<<"X"<<'"'<<", "<<'"'<<"Y"<<'"'<<", "<<'"'<<"T"<<'"'<<", "<<'"'<<"u"<<'"'<<", "<<'"'<<"v"<<'"'<<", "<<'"'<<"qx"<<'"'<<", "<<'"'<<"qy"<<'"'<<", "<<'"'<<"ax"<<'"'<<", "<<'"'<<"ay"<<'"'<<endl;
		writer<<"ZONE I="<<NI<<", J="<<NJ<<", F=POINT"<<endl;
     for (i=1;i<=NI;i++)	
		{
			for (j=1;j<=NJ;j++)
			{
	 			writer<<X_NODE[j][i]<<" "<<Y_NODE[j][i]<<" "<<T[j][i]<<" "<<U[0][j][i]<<" "<<U[1][j][i]<<" "<<qx[j][i]<<" "<<qy[j][i]<<" "<<ax[j][i]<<" "<<ay[j][i]<<" "<<endl;
			}
		}
			writer.close();
	}
}

/*
int FILE_WRITE2(int n){  
	int aa;
	char ch[80];
 
	char str_interface[80];

    aa = sprintf (ch, "%d",n);
			
		strcpy (str_interface,"2d_heat_conduction_with_heat_gen");
		strcat (str_interface,ch);
		strcat (str_interface,".dat");
		
  
         ofstream file(str_interface);
		file<<"VARIABLES = "<<'"'<<"X"<<'"'<<", "<<'"'<<"Y"<<'"'<<endl;
		file<<"ZONE I="<<NI<<", J="<<NI<<", F=POINT"<<endl;
     for(j=1;j<=NI;j++) 	
		{
			for(i=1;i<=NI;i++)
			{
				file<<X_NODE[j][i]<<" "<<Y_NODE[j][i]<<" "<<endl;
			}
		}
		//	file.close();
	}
*/
