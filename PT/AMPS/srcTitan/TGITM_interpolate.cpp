//$Id$
#include <iostream>
#include <fstream>
#include <sstream>
//#include <array>
#include "Titan.h"


using namespace std;

	const int nd = Titan::tgitm_exobase::ndes, nc = Titan::tgitm_exobase::nclmns,
	np = Titan::tgitm_exobase::nintp,on=1,off=0,nspc=PIC::nTotalSpecies;
	
	double Titan::tgitm_exobase::tgitm_grid[nd][nc],realT,
	Titan::tgitm_exobase::interp_val[np],Titan::tgitm_exobase::maxflx[PIC::nTotalSpecies],
	Titan::tgitm_exobase::totalflux[PIC::nTotalSpecies];
	
	const double numT=Titan::tgitm_exobase::Tnum;



/////####################################################################################################3
/*
 * read_tgitm(): This function reads in a flux s^-1 and temperature map (N2,CH4,H2) from the 
 * Bell GITM result at the exobase (lat and long) current the exobase is set at 4010430.7796m 
 * e.g. main_lib.cpp line 38. The current setup adjusts the N2 flux on each surface element by
 * sqrt(numT/realT) where numT is the gas temperature set for an imposed numerical distribution 
 * and realT is the real temperature of the element.
 * 
 * Additionally the total integrated surface flux is calculated in this function e.g. 
 * totalflux[PIC::nTotalSpecies] this  values are used in the gettotalproductionrate function.
*/
	void Titan::tgitm_exobase::read_tgitm() {
	const int nLatitudes = 72;

	
	ifstream fin("data/input/Titan/TGITM_EXOBASE.dat");
	
	if(fin.is_open())
	{
		//initialize toatal flux values to 0
		for(int i=0;i<nspc;i++)Titan::tgitm_exobase::totalflux[i]=0.0;
		
		for(int i=0;i<nd;i++)
		{
			for(int j=0;j<nc;j++)
			{
				//Read in flux and temperature map
				fin >> Titan::tgitm_exobase::tgitm_grid[i][j];
			}
			
//for loop below over species 2 - N2, 3 - CH4, 4 - H2
				for(int j=2;j<5;j++){
						
						//If Snum is on implement numerical distribution (e.g. weighted velocity distribution function)
						if(Titan::tgitm_exobase::Snum == on && j == 2){
							realT=Titan::tgitm_exobase::tgitm_grid[i][5]; //index 5 is column with temperature
							//This line scaled the flux from the local surface element
							Titan::tgitm_exobase::tgitm_grid[i][j] = Titan::tgitm_exobase::tgitm_grid[i][j]*((j==2) ? pow(numT/realT,2) : 1.0);
						}
						//Sum up total fluxes
						Titan::tgitm_exobase::totalflux[j-2] = Titan::tgitm_exobase::totalflux[j-2]+
						Titan::tgitm_exobase::tgitm_grid[i][j];
				}
			
		}
	}
	else
	{
		cout<<"THE FILE DID NOT OPEN"<<endl;
	}

	}
	
/////###########################################################################################################3
/*
 This function interpolates a generated particles coordinated to calculated the nummerical flux
 * from a surface element in the simulation the interpolated values are put in array interp_val[i]
 * The function does a 2D interpolation and is used in Titan.h. In Titan.h the a particle is accepted using
 * acceptance rejection based on the local maximum flux Titan.h (line 69). In this function we find the local
 * maximum using the four interpolation points.
 */
  void Titan::tgitm_exobase::tgitm_interpolate(double polar, double azimuth)
//Angles should be in Titan centered frame
{
	double lat, lng, fp[4];
	int i_lat, i_lng, i_cell, itp_point[4];
	const double tgitm_lb_lat=1.25, tgitm_hb_lat=178.75, 
	tgitm_lb_lng=2.5, tgitm_hb_lng=357.5, dlat=2.5, dlng=5.0;
	const int nLatitudes = 72;
	const double rad2dg=57.2957795;

	//Spherical coordinates conveted to lattitude and east longitude for interpolation
	lat=180.0-polar*rad2dg;
	lng=360.0-azimuth*rad2dg;
	
	if(lng < tgitm_lb_lng) lng=tgitm_lb_lng+dlng*0.25;
	if(lng > tgitm_hb_lng) lng=tgitm_hb_lng-dlng*0.25;
	if(lat > tgitm_hb_lat) lat=tgitm_hb_lat-dlat*0.25;
	if(lat < tgitm_lb_lat) lat=tgitm_lb_lat+dlat*0.25;
	
	i_lat=int((lat-tgitm_lb_lat)/dlat);
	i_lng=int((lng-tgitm_lb_lng)/dlng);
	i_cell=i_lng*(nLatitudes-1)+i_lat;

	itp_point[0]=i_cell+i_lng;
	itp_point[1]=itp_point[0]+1;
	itp_point[2]=itp_point[0]+nLatitudes;
	itp_point[3]=itp_point[2]+1;
	
	//The above translates particle's position to map consistent with the Bell result
	///////////
	
	for(int i=0; i<np; i++)
	{
		fp[0] = Titan::tgitm_exobase::tgitm_grid[itp_point[0]][i+2];
		fp[1] = Titan::tgitm_exobase::tgitm_grid[itp_point[1]][i+2];
		fp[2] = Titan::tgitm_exobase::tgitm_grid[itp_point[2]][i+2];
		fp[3] = Titan::tgitm_exobase::tgitm_grid[itp_point[3]][i+2];
		
		
	if(i < np-1){
		Titan::tgitm_exobase::maxflx[i]=fp[0];

		for (int j=0; j<4; j++)
		{ 
			//calculate local maxium in flux
			if (Titan::tgitm_exobase::maxflx[i] <fp[j] ) Titan::tgitm_exobase::maxflx[i]=fp[j];
		}
	}
	
	////Interpolation routine
		Titan::tgitm_exobase::interp_val[i]=fp[0]*(Titan::tgitm_exobase::tgitm_grid[itp_point[2]][0]-lng)*(Titan::tgitm_exobase::tgitm_grid[itp_point[3]][1]-lat)+
		
		fp[1]*(Titan::tgitm_exobase::tgitm_grid[itp_point[2]][0]-lng)*(lat -Titan::tgitm_exobase::tgitm_grid[itp_point[2]][1])+
		
		fp[2]*(lng-Titan::tgitm_exobase::tgitm_grid[itp_point[0]][0])*(Titan::tgitm_exobase::tgitm_grid[itp_point[1]][1]-lat)+
		
		fp[3]*(lng-Titan::tgitm_exobase::tgitm_grid[itp_point[0]][0])*(lat-Titan::tgitm_exobase::tgitm_grid[itp_point[0]][1]);
			
		Titan::tgitm_exobase::interp_val[i]=Titan::tgitm_exobase::interp_val[i]/(Titan::tgitm_exobase::tgitm_grid[itp_point[2]][0]-Titan::tgitm_exobase::tgitm_grid[itp_point[0]][0])
		/(Titan::tgitm_exobase::tgitm_grid[itp_point[1]][1]-Titan::tgitm_exobase::tgitm_grid[itp_point[0]][1]);
	}
}



