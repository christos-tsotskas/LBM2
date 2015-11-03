#include "global_defines.cuh"

void redistribute2(const int &lx , const int &ly, const int &lz, const FLOATING &Reynolds,
		const FLOATING &r_small, const FLOATING &nu, FLOATING *node_grid){

	/************************************************************************
	 *                                                                      *
	 *     density redistribution in first lattice column                   *
	 *                                                                      *
	 *                                                                      *
	 *     Last change: 04/05/2003                                          *
	 *                                                                      *
	 ************************************************************************/
	/*
				c
				c.......directed flow is induced by density redistribution in the first
				c       lattice column. This is not too clever, since the resulting
				c       reynolds number can not be controlled and reaching steady state
				c       takes quite some time, but it is simple and it works ...
				c       use this to start with no initial field
	 */

	/*
	creates u_n, u_squ and assigns the values in node[]

	 */
	const FLOATING t_0 =  1.0  / 3.0  ;
		const FLOATING t_1 =  1.0 /  18.0 ;
		const FLOATING t_2 =  1.0 / 36.0 ;
		const FLOATING c_squ = 1.0  / 3.0  ;
	int R_big;
	int baffle_position=59;
	FLOATING    mass_flow;

	//.....local variables

	int  y, z, i, x,  yc, zc, yr, zr;


	FLOATING   rho/*local density*/, u_avg,A_out,A_inn,A_anu,pi,u_xa,u_xs;
	FLOATING   u_x,u_y, u_z, u_n[19] , u_squ;


	yc= (ly+1)/2 -1;//CHANGED! ORIGINALLY IT WAS yc= (ly +1)/2; AND zc= (ly +1)/2;
	zc= (ly+1)/2 -1;
	R_big=35;
	mass_flow=0.05;
	pi = acos(0.0);
	u_avg =Reynolds*nu/(2*R_big);
	u_xa= (R_big*R_big) / (R_big*R_big-r_small*r_small)*u_avg/(1+mass_flow);
	A_out=pi*R_big*R_big;
	A_inn=pi*r_small*r_small;
	A_anu=A_out-A_inn;
	u_xs=A_anu*u_xa*mass_flow/A_inn;



	//.....compute weighting factors (depending on lattice geometry) for
	//     increasing/decreasing inlet densities


	//8etei se olo to domain thn idia taxuthta, thn opoia 8a allaksei meta
	//gia ton eswteriko swlhna

	//todo: vale sto katw for, to x na paizei metaksu timwn pou prosdiorizontai apo to rank!
	//px. gia x=0...1/rank... 1/rank+margin... 2/rank.... etc!


	for( z = 0; z< lz ; ++z){
		for( y = 0; y< ly ; ++y){
			for( x = 0; x< baffle_position+1 ; ++x){

				zr=z-zc;
				yr=y-yc;

				if(  yr*yr+zr*zr < r_small*r_small  ){

					rho=0.0;
					for( i=0;i<DENSITIES; ++i)
						rho+=node_grid[index4D(z,y,x,i)];

					u_x = u_xs;
					u_y =  0.0;
					u_z =  0.0;

					u_n[0]= 0.0; //SHOULD NEVER USED!
					u_n[1] =   u_x; //u_xa
					u_n[2] =         u_y;
					u_n[3] = - u_x;
					u_n[4] =       - u_y;
					u_n[5] =   u_z;
					u_n[6] =       - u_z;
					u_n[7] =   u_x + u_y;
					u_n[8] = - u_x + u_y;
					u_n[9] = - u_x - u_y;
					u_n[10] =   u_x - u_y;
					u_n[11] =   u_x - u_z;
					u_n[12] = - u_x - u_z;
					u_n[13] = - u_x + u_z;
					u_n[14] =   u_x + u_z;
					u_n[15] =   u_z + u_y;
					u_n[16] = - u_z + u_y;
					u_n[17] = - u_z - u_y;
					u_n[18] =   u_z - u_y;

					u_squ = u_x*u_x + u_y*u_y + u_z*u_z;

					node_grid[index4D(z,y,x,0)] = (FLOATING) (t_0  * rho *(1.0  - u_squ / (2.0  * c_squ)));

					//...........axis speeds (factor: t_1)
					for( i = 1 ; i < 7 ; ++i)
						node_grid[index4D(z,y,x,i)] = (FLOATING) (t_1 * rho *	 (1.0+ ( u_n[i]/c_squ ) +  0.5* ( (u_n[i]*u_n[i])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));

					//...........diagonal speeds (factor: t_2)
					for( i = 7 ; i < 19 ; ++i)
						node_grid[index4D(z,y,x,i)] =  (FLOATING) (t_2  * rho * (1.0+ ( u_n[i]/c_squ ) +  0.5* ( (u_n[i]*u_n[i])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));

				}
			}
		}
	}
}

void LBM::redistribute(){

	/************************************************************************
	 *                                                                      *
	 *     density redistribution in first lattice column                   *
	 *                                                                      *
	 *                                                                      *
	 *     Last change: 04/05/2003                                          *
	 *                                                                      *
	 ************************************************************************/
	/*
				c
				c.......directed flow is induced by density redistribution in the first
				c       lattice column. This is not too clever, since the resulting
				c       reynolds number can not be controlled and reaching steady state
				c       takes quite some time, but it is simple and it works ...
				c       use this to start with no initial field
	 */

	/*
	creates u_n, u_squ and assigns the values in node[]

	 */

	if(data_location==GPU)
		copy_data_from_device_to_host();

	int R_big;
	int baffle_position=59;
	FLOATING    mass_flow;

	//.....local variables

	int  y, z, x,  yc, zc, yr, zr;


	FLOATING   rho/*local density*/, u_avg,A_out,A_inn,A_anu,pi,u_xa,u_xs;
	FLOATING   u_x,u_y, u_z, u_n[19] , u_squ;


	yc= (ly+1)/2 -1;//CHANGED! ORIGINALLY IT WAS yc= (ly +1)/2; AND zc= (ly +1)/2;
	zc= (ly+1)/2 -1;
	R_big=35;
	mass_flow=0.05;
	pi = acos(0.0);
	u_avg =reynolds*nu/(2*R_big);
	u_xa= (R_big*R_big) / (R_big*R_big-r_small*r_small)*u_avg/(1+mass_flow);
	A_out=pi*R_big*R_big;
	A_inn=pi*r_small*r_small;
	A_anu=A_out-A_inn;
	u_xs=A_anu*u_xa*mass_flow/A_inn;



	//.....compute weighting factors (depending on lattice geometry) for
	//     increasing/decreasing inlet densities


	//8etei se olo to domain thn idia taxuthta, thn opoia 8a allaksei meta
	//gia ton eswteriko swlhna

	//todo: vale sto katw for, to x na paizei metaksu timwn pou prosdiorizontai apo to rank!
	//px. gia x=0...1/rank... 1/rank+margin... 2/rank.... etc!


	for( z = 0; z< lz ; ++z){
		for( y = 0; y< ly ; ++y){
			for( x = 0; x< baffle_position+1 ; ++x){

				zr=z-zc;
				yr=y-yc;

				if(  yr*yr+zr*zr < r_small*r_small  ){

					rho=0.0;
					rho+=D3.Q0[index(z,y,x)]+D3.Q1[index(z,y,x)]+D3.Q2[index(z,y,x)]+D3.Q3[index(z,y,x)];
					rho+=D3.Q4[index(z,y,x)]+D3.Q5[index(z,y,x)]+D3.Q6[index(z,y,x)]+D3.Q7[index(z,y,x)];
					rho+=D3.Q8[index(z,y,x)]+D3.Q9[index(z,y,x)]+D3.Q10[index(z,y,x)]+D3.Q11[index(z,y,x)];
					rho+=D3.Q12[index(z,y,x)]+D3.Q13[index(z,y,x)]+D3.Q14[index(z,y,x)]+D3.Q15[index(z,y,x)];
					rho+=D3.Q16[index(z,y,x)]+D3.Q17[index(z,y,x)]+D3.Q18[index(z,y,x)];


					u_x = u_xs;
					u_y =  0.0;
					u_z =  0.0;

					u_n[0]= 0.0; //SHOULD NEVER USED!
					u_n[1] =   u_x; //u_xa
					u_n[2] =         u_y;
					u_n[3] = - u_x;
					u_n[4] =       - u_y;
					u_n[5] =   u_z;
					u_n[6] =       - u_z;
					u_n[7] =   u_x + u_y;
					u_n[8] = - u_x + u_y;
					u_n[9] = - u_x - u_y;
					u_n[10] =   u_x - u_y;
					u_n[11] =   u_x - u_z;
					u_n[12] = - u_x - u_z;
					u_n[13] = - u_x + u_z;
					u_n[14] =   u_x + u_z;
					u_n[15] =   u_z + u_y;
					u_n[16] = - u_z + u_y;
					u_n[17] = - u_z - u_y;
					u_n[18] =   u_z - u_y;

					u_squ = u_x*u_x + u_y*u_y + u_z*u_z;


					D3.Q0[index(z,y,x)]=(FLOATING) (t_0  * rho *(1.0  - u_squ / (2.0  * c_squ)));

					//...........axis speeds (factor: t_1)

					D3.Q1[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[1]/c_squ ) +  0.5* ( (u_n[1]*u_n[1])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));;
					D3.Q2[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[2]/c_squ ) +  0.5* ( (u_n[2]*u_n[2])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));;
					D3.Q3[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[3]/c_squ ) +  0.5* ( (u_n[3]*u_n[3])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));;
					D3.Q4[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[4]/c_squ ) +  0.5* ( (u_n[4]*u_n[4])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));;
					D3.Q5[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[5]/c_squ ) +  0.5* ( (u_n[5]*u_n[5])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));;
					D3.Q6[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[6]/c_squ ) +  0.5* ( (u_n[6]*u_n[6])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));;

					//...........diagonal speeds (factor: t_2)
					D3.Q7[index(z,y,x)]=(FLOATING) (t_2 * rho *	 (1.0+ ( u_n[7]/c_squ ) +  0.5* ( (u_n[7]*u_n[7])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));;
					D3.Q8[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[8]/c_squ ) +  0.5* ( (u_n[8]*u_n[8])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q9[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[9]/c_squ ) +  0.5* ( (u_n[9]*u_n[9])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q10[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[10]/c_squ ) +  0.5* ( (u_n[10]*u_n[10])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q11[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[11]/c_squ ) +  0.5* ( (u_n[11]*u_n[11])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q12[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[12]/c_squ ) +  0.5* ( (u_n[12]*u_n[12])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q13[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[13]/c_squ ) +  0.5* ( (u_n[13]*u_n[13])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q14[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[14]/c_squ ) +  0.5* ( (u_n[14]*u_n[14])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q15[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[15]/c_squ ) +  0.5* ( (u_n[15]*u_n[15])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q16[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[16]/c_squ ) +  0.5* ( (u_n[16]*u_n[16])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q17[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[17]/c_squ ) +  0.5* ( (u_n[17]*u_n[17])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q18[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[18]/c_squ ) +  0.5* ( (u_n[18]*u_n[18])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));

				}
			}
		}
	}

#ifdef DEBUG
	cout << "#LBM redistribute OK!" << endl;
#endif
}

void LBM::fortran_redistribute(int time){
	int R_big;
//	int baffle_position=59;
	FLOATING    mass_flow;

	//.....local variables

	int  y, z, x,  yc, zc, yr, zr;


	FLOATING   rho/*local density*/, u_avg,A_out,A_inn,A_anu,pi,u_xa,u_xs;
	FLOATING   u_x,u_y, u_z, u_n[19] , u_squ;


	yc= (ly+1)/2 ;//CHANGED! ORIGINALLY IT WAS yc= (ly +1)/2; AND zc= (ly +1)/2;
	zc= (ly+1)/2 ;
	R_big=35;
	mass_flow=0.05;
	pi = acos(0.0);
	u_avg =reynolds*nu/(2*R_big);
	u_xa= (R_big*R_big) / (R_big*R_big-r_small*r_small)*u_avg/(1+mass_flow);
	A_out=pi*R_big*R_big;
	A_inn=pi*r_small*r_small;
	A_anu=A_out-A_inn;
	u_xs=A_anu*u_xa*mass_flow/A_inn;



	//.....compute weighting factors (depending on lattice geometry) for
	//     increasing/decreasing inlet densities


	//8etei se olo to domain thn idia taxuthta, thn opoia 8a allaksei meta
	//gia ton eswteriko swlhna

	//todo: vale sto katw for, to x na paizei metaksu timwn pou prosdiorizontai apo to rank!
	//px. gia x=0...1/rank... 1/rank+margin... 2/rank.... etc!

	int end_x=0;
	if (time==0)
		end_x=lx;
	for( z = 0; z< lz ; ++z){
		for( y = 0; y< ly ; ++y){
			for( x = 0; x< end_x ; ++x){



				rho=0.0;
				rho+=D3.Q0[index(z,y,x)]+D3.Q1[index(z,y,x)]+D3.Q2[index(z,y,x)]+D3.Q3[index(z,y,x)];
				rho+=D3.Q4[index(z,y,x)]+D3.Q5[index(z,y,x)]+D3.Q6[index(z,y,x)]+D3.Q7[index(z,y,x)];
				rho+=D3.Q8[index(z,y,x)]+D3.Q9[index(z,y,x)]+D3.Q10[index(z,y,x)]+D3.Q11[index(z,y,x)];
				rho+=D3.Q12[index(z,y,x)]+D3.Q13[index(z,y,x)]+D3.Q14[index(z,y,x)]+D3.Q15[index(z,y,x)];
				rho+=D3.Q16[index(z,y,x)]+D3.Q17[index(z,y,x)]+D3.Q18[index(z,y,x)];




				u_x = u_xa;
				u_y =  0.0;
				u_z =  0.0;

				u_n[0]= 0.0; //SHOULD NEVER USED!
				u_n[1] =   u_x; //u_xa
				u_n[2] =         u_y;
				u_n[3] = - u_x;
				u_n[4] =       - u_y;
				u_n[5] =   u_z;
				u_n[6] =       - u_z;
				u_n[7] =   u_x + u_y;
				u_n[8] = - u_x + u_y;
				u_n[9] = - u_x - u_y;
				u_n[10] =   u_x - u_y;
				u_n[11] =   u_x - u_z;
				u_n[12] = - u_x - u_z;
				u_n[13] = - u_x + u_z;
				u_n[14] =   u_x + u_z;
				u_n[15] =   u_z + u_y;
				u_n[16] = - u_z + u_y;
				u_n[17] = - u_z - u_y;
				u_n[18] =   u_z - u_y;

				u_squ = u_x*u_x + u_y*u_y + u_z*u_z;


				D3.Q0[index(z,y,x)]=(FLOATING)           (t_0  * rho *(1.0  - u_squ / (2.0  * c_squ)));


				//...........axis speeds (factor: t_1)
				//				for( i = 1 ; i < 7 ; ++i)
				//					node_grid[index4D(z,y,x,i)] = (FLOATING) (t_1 * rho *	 (1.0+ ( u_n[i]/c_squ ) +  0.5* ( (u_n[i]*u_n[i])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));

				D3.Q1[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[1]/c_squ ) +  0.5* ( (u_n[1]*u_n[1])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
				D3.Q2[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[2]/c_squ ) +  0.5* ( (u_n[2]*u_n[2])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
				D3.Q3[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[3]/c_squ ) +  0.5* ( (u_n[3]*u_n[3])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
				D3.Q4[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[4]/c_squ ) +  0.5* ( (u_n[4]*u_n[4])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
				D3.Q5[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[5]/c_squ ) +  0.5* ( (u_n[5]*u_n[5])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
				D3.Q6[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[6]/c_squ ) +  0.5* ( (u_n[6]*u_n[6])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));


				//...........diagonal speeds (factor: t_2)
				//				for( i = 7 ; i < 19 ; ++i)
				//					node_grid[index4D(z,y,x,i)] =  (FLOATING) (t_2  * rho * (1.0+ ( u_n[i]/c_squ ) +  0.5* ( (u_n[i]*u_n[i])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));

				D3.Q7[index(z,y,x)]=(FLOATING)  (t_2 * rho *  (1.0+ ( u_n[7]/c_squ ) +  0.5* ( (u_n[7]*u_n[7])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));;
				D3.Q8[index(z,y,x)]=(FLOATING)  (t_2  * rho * (1.0+ ( u_n[8]/c_squ ) +  0.5* ( (u_n[8]*u_n[8])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
				D3.Q9[index(z,y,x)]=(FLOATING)  (t_2  * rho * (1.0+ ( u_n[9]/c_squ ) +  0.5* ( (u_n[9]*u_n[9])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
				D3.Q10[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[10]/c_squ ) +  0.5* ( (u_n[10]*u_n[10])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
				D3.Q11[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[11]/c_squ ) +  0.5* ( (u_n[11]*u_n[11])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
				D3.Q12[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[12]/c_squ ) +  0.5* ( (u_n[12]*u_n[12])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
				D3.Q13[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[13]/c_squ ) +  0.5* ( (u_n[13]*u_n[13])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
				D3.Q14[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[14]/c_squ ) +  0.5* ( (u_n[14]*u_n[14])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
				D3.Q15[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[15]/c_squ ) +  0.5* ( (u_n[15]*u_n[15])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
				D3.Q16[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[16]/c_squ ) +  0.5* ( (u_n[16]*u_n[16])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
				D3.Q17[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[17]/c_squ ) +  0.5* ( (u_n[17]*u_n[17])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
				D3.Q18[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[18]/c_squ ) +  0.5* ( (u_n[18]*u_n[18])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));








			}
		}
	}
	if (time==0)
		end_x=0;
	else
		end_x=59;
	for( z = 0; z< lz ; ++z){
		for( y = 0; y< ly ; ++y){
			for( x = 0; x< end_x ; ++x){

				zr=z+1-zc;
				yr=y+1-yc;

				if(  yr*yr+zr*zr < r_small*r_small  ){

					rho=0.0;
					rho+=D3.Q0[index(z,y,x)]+D3.Q1[index(z,y,x)]+D3.Q2[index(z,y,x)]+D3.Q3[index(z,y,x)];
					rho+=D3.Q4[index(z,y,x)]+D3.Q5[index(z,y,x)]+D3.Q6[index(z,y,x)]+D3.Q7[index(z,y,x)];
					rho+=D3.Q8[index(z,y,x)]+D3.Q9[index(z,y,x)]+D3.Q10[index(z,y,x)]+D3.Q11[index(z,y,x)];
					rho+=D3.Q12[index(z,y,x)]+D3.Q13[index(z,y,x)]+D3.Q14[index(z,y,x)]+D3.Q15[index(z,y,x)];
					rho+=D3.Q16[index(z,y,x)]+D3.Q17[index(z,y,x)]+D3.Q18[index(z,y,x)];




					u_x = u_xs;
					u_y =  0.0;
					u_z =  0.0;

					u_n[0]= 0.0; //SHOULD NEVER USED!
					u_n[1] =   u_x; //u_xa
					u_n[2] =         u_y;
					u_n[3] = - u_x;
					u_n[4] =       - u_y;
					u_n[5] =   u_z;
					u_n[6] =       - u_z;
					u_n[7] =   u_x + u_y;
					u_n[8] = - u_x + u_y;
					u_n[9] = - u_x - u_y;
					u_n[10] =   u_x - u_y;
					u_n[11] =   u_x - u_z;
					u_n[12] = - u_x - u_z;
					u_n[13] = - u_x + u_z;
					u_n[14] =   u_x + u_z;
					u_n[15] =   u_z + u_y;
					u_n[16] = - u_z + u_y;
					u_n[17] = - u_z - u_y;
					u_n[18] =   u_z - u_y;

					u_squ = u_x*u_x + u_y*u_y + u_z*u_z;


					D3.Q0[index(z,y,x)]=(FLOATING)           (t_0  * rho *(1.0  - u_squ / (2.0  * c_squ)));


					//...........axis speeds (factor: t_1)
					//				for( i = 1 ; i < 7 ; ++i)
					//					node_grid[index4D(z,y,x,i)] = (FLOATING) (t_1 * rho *	 (1.0+ ( u_n[i]/c_squ ) +  0.5* ( (u_n[i]*u_n[i])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));

					D3.Q1[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[1]/c_squ ) +  0.5* ( (u_n[1]*u_n[1])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q2[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[2]/c_squ ) +  0.5* ( (u_n[2]*u_n[2])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q3[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[3]/c_squ ) +  0.5* ( (u_n[3]*u_n[3])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q4[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[4]/c_squ ) +  0.5* ( (u_n[4]*u_n[4])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q5[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[5]/c_squ ) +  0.5* ( (u_n[5]*u_n[5])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q6[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[6]/c_squ ) +  0.5* ( (u_n[6]*u_n[6])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));


					//...........diagonal speeds (factor: t_2)
					//				for( i = 7 ; i < 19 ; ++i)
					//					node_grid[index4D(z,y,x,i)] =  (FLOATING) (t_2  * rho * (1.0+ ( u_n[i]/c_squ ) +  0.5* ( (u_n[i]*u_n[i])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));

					D3.Q7[index(z,y,x)]=(FLOATING)  (t_2 * rho *  (1.0+ ( u_n[7]/c_squ ) +  0.5* ( (u_n[7]*u_n[7])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));;
					D3.Q8[index(z,y,x)]=(FLOATING)  (t_2  * rho * (1.0+ ( u_n[8]/c_squ ) +  0.5* ( (u_n[8]*u_n[8])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q9[index(z,y,x)]=(FLOATING)  (t_2  * rho * (1.0+ ( u_n[9]/c_squ ) +  0.5* ( (u_n[9]*u_n[9])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q10[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[10]/c_squ ) +  0.5* ( (u_n[10]*u_n[10])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q11[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[11]/c_squ ) +  0.5* ( (u_n[11]*u_n[11])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q12[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[12]/c_squ ) +  0.5* ( (u_n[12]*u_n[12])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q13[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[13]/c_squ ) +  0.5* ( (u_n[13]*u_n[13])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q14[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[14]/c_squ ) +  0.5* ( (u_n[14]*u_n[14])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q15[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[15]/c_squ ) +  0.5* ( (u_n[15]*u_n[15])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q16[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[16]/c_squ ) +  0.5* ( (u_n[16]*u_n[16])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q17[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[17]/c_squ ) +  0.5* ( (u_n[17]*u_n[17])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
					D3.Q18[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[18]/c_squ ) +  0.5* ( (u_n[18]*u_n[18])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));






				}

			}
		}
	}






}





