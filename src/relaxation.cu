#include "global_defines.cuh"
#include <numeric>




void LBM::relaxation(){





	/*One-step density relaxation process

				c.......density relaxation: a single time relaxation with relaxation
				c       parameter omega is applied here. This step is only "local",
				c       nothing is propagated through the lattice.
				c*/

	if(data_location==GPU)
		copy_data_from_device_to_host();

	int  x,y,z;
	FLOATING  u_x=0.0, u_y=0.0, u_z=0.0, u_squ=0.0, rho=0.0, reciprocal_rho=0.0;




	FLOATING u_n[19];
	//	FLOATING n_equ[19];
	FLOATING buff[19];

	//	FLOATING u_n_squared[19];
	//FLOATING two_x_c_squ_sqared;
	const FLOATING omega_x_t_0=omega*t_0, omega_x_t_1=omega*t_1, omega_x_t_2=omega*t_2;
	FLOATING omega_x_rho_x_t_0, omega_x_rho_x_t_1, omega_x_rho_x_t_2;
	FLOATING temp_factor;
	FLOATING u_n__over__c_squ[19];
	FLOATING u_n__over__c_squ__squared_and_halved[19];

	//....square speed of sound
	/*      compute the out let velocity with a convevtive boundary condition
					c.....loop over all nodes
					c.....attention: actual densities are stored after the propagation
					c                step in the help-array n_hlp !*/


#pragma unroll
	for (z = 0 ; z< lz ; ++z){
#pragma unroll
		for (y = 0 ; y< ly ; ++y){
#pragma unroll
			for (x = 0 ; x< lx; ++x){

				/*c.........only free nodes are considered here
					!if (.not. obstacles[z][y][x]) then
					c...........integral local density
					c...........initialize variable ro*/
				//memory optimised implementation
				buff[0]=D3_hlp.Q0[index(z,y,x)];
				buff[1]=D3_hlp.Q1[index(z,y,x)];
				buff[2]=D3_hlp.Q2[index(z,y,x)];
				buff[3]=D3_hlp.Q3[index(z,y,x)];
				buff[4]=D3_hlp.Q4[index(z,y,x)];
				buff[5]=D3_hlp.Q5[index(z,y,x)];
				buff[6]=D3_hlp.Q6[index(z,y,x)];
				buff[7]=D3_hlp.Q7[index(z,y,x)];
				buff[8]=D3_hlp.Q8[index(z,y,x)];
				buff[9]=D3_hlp.Q9[index(z,y,x)];
				buff[10]=D3_hlp.Q10[index(z,y,x)];
				buff[11]=D3_hlp.Q11[index(z,y,x)];
				buff[12]=D3_hlp.Q12[index(z,y,x)];
				buff[13]=D3_hlp.Q13[index(z,y,x)];
				buff[14]=D3_hlp.Q14[index(z,y,x)];
				buff[15]=D3_hlp.Q15[index(z,y,x)];
				buff[16]=D3_hlp.Q16[index(z,y,x)];
				buff[17]=D3_hlp.Q17[index(z,y,x)];
				buff[18]=D3_hlp.Q18[index(z,y,x)];

				rho=accumulate(buff, buff+DENSITIES, 0.0);

				reciprocal_rho=1.0/rho;




				switch(obstacles[index(z,y,x)]){
				case 1:
					u_x = 0.0;
					u_y = 0.0;
					u_z = 0.0;
					break;
				default:
					u_x = 0.0;
					u_x =  reciprocal_rho*(buff[1] + buff[7] + buff[10] +buff[11] + buff[14]-
							(buff[3] + buff[8] + buff[9] +buff[12] + buff[13]));

					u_y = 0.0;
					u_y =  reciprocal_rho*(buff[2]+buff[8]+buff[7]+buff[16] + buff[15] -
							(buff[4] + buff[9] + buff[10] +buff[17] + buff[18]));

					u_z = 0.0;
					u_z =  reciprocal_rho*(buff[5]+buff[13]+buff[14]+buff[15]+buff[18]-
							(buff[6]+buff[12]+buff[11]+buff[16]+buff[17]));
					break;
				}//switch(obstacles[index(z,y,x)])

				//original implementation
				//				rho=0.0;
				//				rho+=D3_hlp.Q0[index(z,y,x)]+D3_hlp.Q1[index(z,y,x)]+D3_hlp.Q2[index(z,y,x)]+D3_hlp.Q3[index(z,y,x)];
				//				rho+=D3_hlp.Q4[index(z,y,x)]+D3_hlp.Q5[index(z,y,x)]+D3_hlp.Q6[index(z,y,x)]+D3_hlp.Q7[index(z,y,x)];
				//				rho+=D3_hlp.Q8[index(z,y,x)]+D3_hlp.Q9[index(z,y,x)]+D3_hlp.Q10[index(z,y,x)]+D3_hlp.Q11[index(z,y,x)];
				//				rho+=D3_hlp.Q12[index(z,y,x)]+D3_hlp.Q13[index(z,y,x)]+D3_hlp.Q14[index(z,y,x)]+D3_hlp.Q15[index(z,y,x)];
				//				rho+=D3_hlp.Q16[index(z,y,x)]+D3_hlp.Q17[index(z,y,x)]+D3_hlp.Q18[index(z,y,x)];
				//				reciprocal_rho=1.0/rho;

				//...........x-, and y- velocity components

				//				switch(obstacles[index(z,y,x)]){
				//				case 1:
				//					u_x = 0.0;
				//					u_y = 0.0;
				//					u_z = 0.0;
				//					break;
				//				default:
				//					u_x = (FLOATING) reciprocal_rho*(D3_hlp.Q1[index(z,y,x)] + D3_hlp.Q7[index(z,y,x)] + D3_hlp.Q10[index(z,y,x)] +
				//							D3_hlp.Q11[index(z,y,x)] + D3_hlp.Q14[index(z,y,x)] -
				//							(D3_hlp.Q3[index(z,y,x)] + D3_hlp.Q8[index(z,y,x)] + D3_hlp.Q9[index(z,y,x)] +
				//									D3_hlp.Q12[index(z,y,x)] + D3_hlp.Q13[index(z,y,x)]));
				//
				//					u_y = (FLOATING) reciprocal_rho*(D3_hlp.Q2[index(z,y,x)] + D3_hlp.Q8[index(z,y,x)] + D3_hlp.Q7[index(z,y,x)] +
				//							D3_hlp.Q16[index(z,y,x)] + D3_hlp.Q15[index(z,y,x)] -
				//							(D3_hlp.Q4[index(z,y,x)] + D3_hlp.Q9[index(z,y,x)] + D3_hlp.Q10[index(z,y,x)] +
				//									D3_hlp.Q17[index(z,y,x)] + D3_hlp.Q18[index(z,y,x)]));
				//
				//					u_z = (FLOATING) reciprocal_rho*(D3_hlp.Q5[index(z,y,x)] + D3_hlp.Q13[index(z,y,x)] + D3_hlp.Q14[index(z,y,x)] +
				//							D3_hlp.Q15[index(z,y,x)] + D3_hlp.Q18[index(z,y,x)] -
				//							(D3_hlp.Q6[index(z,y,x)] + D3_hlp.Q12[index(z,y,x)] + D3_hlp.Q11[index(z,y,x)] +
				//									D3_hlp.Q16[index(z,y,x)] + D3_hlp.Q17[index(z,y,x)]));
				//					break;
				//				}//switch(obstacles[index(z,y,x)])




				u_squ = (FLOATING)  u_x*u_x + u_y*u_y + u_z*u_z;
				temp_factor= 0.5*(2.0* c_squ - u_squ)/c_squ;
				//u_squ = (FLOATING)  pow(u_x,2) + pow(u_y,2) + pow(u_z,2);


				/*...........n- velocity compnents (n = lattice node connection vectors)
					c...........this is only necessary for clearence, and only 3 speeds would
					c...........be necessary*/


				//WARNING!!!! o pinakas autos exei tropopoihmena indices!!!!
				u_n[0]= 0.0; //SHOULD NEVER USED!
				u_n[1] =   u_x;
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

#pragma unroll
				for(int i=0; i<DENSITIES; ++i){
					u_n__over__c_squ[i]=reciprocal_c_squ*u_n[i];
					u_n__over__c_squ__squared_and_halved[i]=0.5*u_n__over__c_squ[i]*u_n__over__c_squ[i];
				}

				/*c...........equilibrium densities
					c...........this can be rewritten to improve computational performance
					c...........considerabely !
					c
					c...........zero velocity density
					c*/
				//memory optimised implementation! WARNING!!! different from the original case!

				//two_x_c_squ_sqared=2.0*c_squ*c_squ;
				omega_x_rho_x_t_0=omega_x_t_0*rho;
				omega_x_rho_x_t_1=omega_x_t_1*rho;
				omega_x_rho_x_t_2=omega_x_t_2*rho;


				//				//...........relaxation step


				//omega_x_rho_x_t_0*(1.0 - 0.5*u_squ/c_squ);
				D3.Q0[index(z,y,x)]=buff[0]*one_minus_omega+omega_x_rho_x_t_0*(u_n__over__c_squ__squared_and_halved[0]+u_n__over__c_squ[0]+temp_factor);




				D3.Q1[index(z,y,x)]=buff[1]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[1]+u_n__over__c_squ[1]+temp_factor);
				D3.Q2[index(z,y,x)]=buff[2]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[2]+u_n__over__c_squ[2]+temp_factor);
				D3.Q3[index(z,y,x)]=buff[3]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[3]+u_n__over__c_squ[3]+temp_factor);
				D3.Q4[index(z,y,x)]=buff[4]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[4]+u_n__over__c_squ[4]+temp_factor);
				D3.Q5[index(z,y,x)]=buff[5]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[5]+u_n__over__c_squ[5]+temp_factor);
				D3.Q6[index(z,y,x)]=buff[6]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[6]+u_n__over__c_squ[6]+temp_factor);



				D3.Q7[index(z,y,x)]= buff[ 7]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[ 7]+u_n__over__c_squ[ 7]+temp_factor);
				D3.Q8[index(z,y,x)]= buff[ 8]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[ 8]+u_n__over__c_squ[ 8]+temp_factor);
				D3.Q9[index(z,y,x)]= buff[ 9]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[ 9]+u_n__over__c_squ[ 9]+temp_factor);
				D3.Q10[index(z,y,x)]=buff[10]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[10]+u_n__over__c_squ[10]+temp_factor);
				D3.Q11[index(z,y,x)]=buff[11]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[11]+u_n__over__c_squ[11]+temp_factor);
				D3.Q12[index(z,y,x)]=buff[12]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[12]+u_n__over__c_squ[12]+temp_factor);
				D3.Q13[index(z,y,x)]=buff[13]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[13]+u_n__over__c_squ[13]+temp_factor);
				D3.Q14[index(z,y,x)]=buff[14]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[14]+u_n__over__c_squ[14]+temp_factor);
				D3.Q15[index(z,y,x)]=buff[15]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[15]+u_n__over__c_squ[15]+temp_factor);
				D3.Q16[index(z,y,x)]=buff[16]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[16]+u_n__over__c_squ[16]+temp_factor);
				D3.Q17[index(z,y,x)]=buff[17]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[17]+u_n__over__c_squ[17]+temp_factor);
				D3.Q18[index(z,y,x)]=buff[18]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[18]+u_n__over__c_squ[18]+temp_factor);

				//original implementation
				//				n_equ[0] = t_0  * rho*(1.0 - u_squ / (2.0 * c_squ));
				//
				//								//...........axis speeds (factor: t_1)
				//				#pragma unroll
				//								for (int i = 1 ; i< 7; ++i){
				//									n_equ[i] = t_1 * rho*(1.0 + u_n[i] / c_squ
				//											+ (u_n[i]*u_n[i]) / (2.0 * c_squ *c_squ)
				//											- u_squ / (2.0  * c_squ));
				//								}
				//
				//								//...........diagonal speeds (factor: t_2)
				//				#pragma unroll
				//								for (int i = 7 ; i< 19; ++i){
				//									n_equ[i] = t_2  * rho*(1.0 + u_n[i] / c_squ
				//											+ (u_n[i]*u_n[i]) / (2.0 * c_squ *c_squ)
				//											- u_squ / (2.0  * c_squ));
				//								}

				//				D3.Q0[index(z,y,x)]=D3_hlp.Q0[index(z,y,x)]+omega*(n_equ[0] - D3_hlp.Q0[index(z,y,x)]);
				//								D3.Q1[index(z,y,x)]=D3_hlp.Q1[index(z,y,x)]+omega*(n_equ[1] - D3_hlp.Q1[index(z,y,x)]);
				//								D3.Q2[index(z,y,x)]=D3_hlp.Q2[index(z,y,x)]+omega*(n_equ[2] - D3_hlp.Q2[index(z,y,x)]);
				//								D3.Q3[index(z,y,x)]=D3_hlp.Q3[index(z,y,x)]+omega*(n_equ[3] - D3_hlp.Q3[index(z,y,x)]);
				//								D3.Q4[index(z,y,x)]=D3_hlp.Q4[index(z,y,x)]+omega*(n_equ[4] - D3_hlp.Q4[index(z,y,x)]);
				//								D3.Q5[index(z,y,x)]=D3_hlp.Q5[index(z,y,x)]+omega*(n_equ[5] - D3_hlp.Q5[index(z,y,x)]);
				//								D3.Q6[index(z,y,x)]=D3_hlp.Q6[index(z,y,x)]+omega*(n_equ[6] - D3_hlp.Q6[index(z,y,x)]);
				//								D3.Q7[index(z,y,x)]=D3_hlp.Q7[index(z,y,x)]+omega*(n_equ[7] - D3_hlp.Q7[index(z,y,x)]);
				//								D3.Q8[index(z,y,x)]=D3_hlp.Q8[index(z,y,x)]+omega*(n_equ[8] - D3_hlp.Q8[index(z,y,x)]);
				//								D3.Q9[index(z,y,x)]=D3_hlp.Q9[index(z,y,x)]+omega*(n_equ[9] - D3_hlp.Q9[index(z,y,x)]);
				//								D3.Q10[index(z,y,x)]=D3_hlp.Q10[index(z,y,x)]+omega*(n_equ[10] - D3_hlp.Q10[index(z,y,x)]);
				//								D3.Q11[index(z,y,x)]=D3_hlp.Q11[index(z,y,x)]+omega*(n_equ[11] - D3_hlp.Q11[index(z,y,x)]);
				//								D3.Q12[index(z,y,x)]=D3_hlp.Q12[index(z,y,x)]+omega*(n_equ[12] - D3_hlp.Q12[index(z,y,x)]);
				//								D3.Q13[index(z,y,x)]=D3_hlp.Q13[index(z,y,x)]+omega*(n_equ[13] - D3_hlp.Q13[index(z,y,x)]);
				//								D3.Q14[index(z,y,x)]=D3_hlp.Q14[index(z,y,x)]+omega*(n_equ[14] - D3_hlp.Q14[index(z,y,x)]);
				//								D3.Q15[index(z,y,x)]=D3_hlp.Q15[index(z,y,x)]+omega*(n_equ[15] - D3_hlp.Q15[index(z,y,x)]);
				//								D3.Q16[index(z,y,x)]=D3_hlp.Q16[index(z,y,x)]+omega*(n_equ[16] - D3_hlp.Q16[index(z,y,x)]);
				//								D3.Q17[index(z,y,x)]=D3_hlp.Q17[index(z,y,x)]+omega*(n_equ[17] - D3_hlp.Q17[index(z,y,x)]);
				//								D3.Q18[index(z,y,x)]=D3_hlp.Q18[index(z,y,x)]+omega*(n_equ[18] - D3_hlp.Q18[index(z,y,x)]);

				if (x == lx-2) {
					u_previous_spatial_boundary[index2D(z,y)] = u_x;
					v_previous_spatial_boundary[index2D(z,y)] = u_y;
					w_previous_spatial_boundary[index2D(z,y)] = u_z;

					u_current[index2D(z,y)] = u_x;
					v_current[index2D(z,y)] = u_y;
					w_current[index2D(z,y)] = u_z;

				}//if (x == lx-2)
			}//for (x = 0 ; x< lx; ++x)
		}//for (y = 0 ; y< ly ; ++y)
	}//for (z = 0 ; z< lz ; ++z)
#ifdef DEBUG
	cout << " #LBM relaxation OK!" << endl;
#endif
}

void LBM::initial_relaxation(){





	/*One-step density relaxation process

				c.......density relaxation: a single time relaxation with relaxation
				c       parameter omega is applied here. This step is only "local",
				c       nothing is propagated through the lattice.
				c*/

	int  x,y,z;
	FLOATING  u_x=0.0, u_y=0.0, u_z=0.0, u_squ=0.0, rho=0.0;
	const FLOATING  tau=3.0*nu + 0.5, omega = 1.0 /tau; //	omega=1.0/(3.0*nu+0.5);



	FLOATING u_n[19], n_equ[19];

	//....square speed of sound
	/*      compute the out let velocity with a convevtive boundary condition
					c.....loop over all nodes
					c.....attention: actual densities are stored after the propagation
					c                step in the help-array n_hlp !*/



	for (z = 0 ; z< lz ; ++z){
		for (y = 0 ; y< ly ; ++y){
			for (x = 0 ; x< lx; ++x){

				/*c.........only free nodes are considered here
					!if (.not. obstacles[z][y][x]) then
					c...........integral local density
					c...........initialize variable ro*/


				rho=0.0;
				rho+=D3_hlp.Q0[index(z,y,x)]+D3_hlp.Q1[index(z,y,x)]+D3_hlp.Q2[index(z,y,x)]+D3_hlp.Q3[index(z,y,x)];
				rho+=D3_hlp.Q4[index(z,y,x)]+D3_hlp.Q5[index(z,y,x)]+D3_hlp.Q6[index(z,y,x)]+D3_hlp.Q7[index(z,y,x)];
				rho+=D3_hlp.Q8[index(z,y,x)]+D3_hlp.Q9[index(z,y,x)]+D3_hlp.Q10[index(z,y,x)]+D3_hlp.Q11[index(z,y,x)];
				rho+=D3_hlp.Q12[index(z,y,x)]+D3_hlp.Q13[index(z,y,x)]+D3_hlp.Q14[index(z,y,x)]+D3_hlp.Q15[index(z,y,x)];
				rho+=D3_hlp.Q16[index(z,y,x)]+D3_hlp.Q17[index(z,y,x)]+D3_hlp.Q18[index(z,y,x)];


				//...........x-, and y- velocity components


				if ( obstacles[index(z,y,x)]==1 ) {
					u_x = 0.0;
					u_y = 0.0;
					u_z = 0.0;
				}else{





					u_x = (FLOATING) (D3_hlp.Q1[index(z,y,x)] + D3_hlp.Q7[index(z,y,x)] + D3_hlp.Q10[index(z,y,x)] +
							D3_hlp.Q11[index(z,y,x)] + D3_hlp.Q14[index(z,y,x)] -
							(D3_hlp.Q3[index(z,y,x)] + D3_hlp.Q8[index(z,y,x)] + D3_hlp.Q9[index(z,y,x)] +
									D3_hlp.Q12[index(z,y,x)] + D3_hlp.Q13[index(z,y,x)])) / rho;

					u_y = (FLOATING) (D3_hlp.Q2[index(z,y,x)] + D3_hlp.Q8[index(z,y,x)] + D3_hlp.Q7[index(z,y,x)] +
							D3_hlp.Q16[index(z,y,x)] + D3_hlp.Q15[index(z,y,x)] -
							(D3_hlp.Q4[index(z,y,x)] + D3_hlp.Q9[index(z,y,x)] + D3_hlp.Q10[index(z,y,x)] +
									D3_hlp.Q17[index(z,y,x)] + D3_hlp.Q18[index(z,y,x)])) / rho;

					u_z = (FLOATING) (D3_hlp.Q5[index(z,y,x)] + D3_hlp.Q13[index(z,y,x)] + D3_hlp.Q14[index(z,y,x)] +
							D3_hlp.Q15[index(z,y,x)] + D3_hlp.Q18[index(z,y,x)] -
							(D3_hlp.Q6[index(z,y,x)] + D3_hlp.Q12[index(z,y,x)] + D3_hlp.Q11[index(z,y,x)] +
									D3_hlp.Q16[index(z,y,x)] + D3_hlp.Q17[index(z,y,x)])) / rho;

				}

				u_squ = (FLOATING)  u_x*u_x + u_y*u_y + u_z*u_z;
				//u_squ = (FLOATING)  pow(u_x,2) + pow(u_y,2) + pow(u_z,2);


				/*...........n- velocity compnents (n = lattice node connection vectors)
					c...........this is only necessary for clearence, and only 3 speeds would
					c...........be necessary*/

				//WARNING!!!! o pinakas autos exei tropopoihmena indices!!!!
				u_n[0]= 0.0; //SHOULD NEVER USED!
				u_n[1] =   u_x;
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

				/*c...........equilibrium densities
					c...........this can be rewritten to improve computational performance
					c...........considerabely !
					c
					c...........zero velocity density
					c*/
				n_equ[0] = t_0  * rho*(1.0 - u_squ / (2.0 * c_squ));

				//...........axis speeds (factor: t_1)
				for (int i = 1 ; i< 7; ++i){
					n_equ[i] = t_1 * rho*(1.0 + u_n[i] / c_squ
							+ (u_n[i]*u_n[i]) / (2.0 * c_squ *c_squ)
							- u_squ / (2.0  * c_squ));
				}

				//...........diagonal speeds (factor: t_2)
				for (int i = 7 ; i< 19; ++i){
					n_equ[i] = t_2  * rho*(1.0 + u_n[i] / c_squ
							+ (u_n[i]*u_n[i]) / (2.0 * c_squ *c_squ)
							- u_squ / (2.0  * c_squ));
				}


				//...........relaxation step




				D3.Q0[index(z,y,x)]=D3_hlp.Q0[index(z,y,x)]+omega*(n_equ[0] - D3_hlp.Q0[index(z,y,x)]);
				D3.Q1[index(z,y,x)]=D3_hlp.Q1[index(z,y,x)]+omega*(n_equ[1] - D3_hlp.Q1[index(z,y,x)]);
				D3.Q2[index(z,y,x)]=D3_hlp.Q2[index(z,y,x)]+omega*(n_equ[2] - D3_hlp.Q2[index(z,y,x)]);
				D3.Q3[index(z,y,x)]=D3_hlp.Q3[index(z,y,x)]+omega*(n_equ[3] - D3_hlp.Q3[index(z,y,x)]);
				D3.Q4[index(z,y,x)]=D3_hlp.Q4[index(z,y,x)]+omega*(n_equ[4] - D3_hlp.Q4[index(z,y,x)]);
				D3.Q5[index(z,y,x)]=D3_hlp.Q5[index(z,y,x)]+omega*(n_equ[5] - D3_hlp.Q5[index(z,y,x)]);
				D3.Q6[index(z,y,x)]=D3_hlp.Q6[index(z,y,x)]+omega*(n_equ[6] - D3_hlp.Q6[index(z,y,x)]);
				D3.Q7[index(z,y,x)]=D3_hlp.Q7[index(z,y,x)]+omega*(n_equ[7] - D3_hlp.Q7[index(z,y,x)]);
				D3.Q8[index(z,y,x)]=D3_hlp.Q8[index(z,y,x)]+omega*(n_equ[8] - D3_hlp.Q8[index(z,y,x)]);
				D3.Q9[index(z,y,x)]=D3_hlp.Q9[index(z,y,x)]+omega*(n_equ[9] - D3_hlp.Q9[index(z,y,x)]);
				D3.Q10[index(z,y,x)]=D3_hlp.Q10[index(z,y,x)]+omega*(n_equ[10] - D3_hlp.Q10[index(z,y,x)]);
				D3.Q11[index(z,y,x)]=D3_hlp.Q11[index(z,y,x)]+omega*(n_equ[11] - D3_hlp.Q11[index(z,y,x)]);
				D3.Q12[index(z,y,x)]=D3_hlp.Q12[index(z,y,x)]+omega*(n_equ[12] - D3_hlp.Q12[index(z,y,x)]);
				D3.Q13[index(z,y,x)]=D3_hlp.Q13[index(z,y,x)]+omega*(n_equ[13] - D3_hlp.Q13[index(z,y,x)]);
				D3.Q14[index(z,y,x)]=D3_hlp.Q14[index(z,y,x)]+omega*(n_equ[14] - D3_hlp.Q14[index(z,y,x)]);
				D3.Q15[index(z,y,x)]=D3_hlp.Q15[index(z,y,x)]+omega*(n_equ[15] - D3_hlp.Q15[index(z,y,x)]);
				D3.Q16[index(z,y,x)]=D3_hlp.Q16[index(z,y,x)]+omega*(n_equ[16] - D3_hlp.Q16[index(z,y,x)]);
				D3.Q17[index(z,y,x)]=D3_hlp.Q17[index(z,y,x)]+omega*(n_equ[17] - D3_hlp.Q17[index(z,y,x)]);
				D3.Q18[index(z,y,x)]=D3_hlp.Q18[index(z,y,x)]+omega*(n_equ[18] - D3_hlp.Q18[index(z,y,x)]);

				//at the penultimat slice, save previous and current slices
				if (x == lx-2) {
					u_previous_spatial_boundary[index2D(z,y)] = u_x;
					v_previous_spatial_boundary[index2D(z,y)] = u_y;
					w_previous_spatial_boundary[index2D(z,y)] = u_z;

					u_current[index2D(z,y)] = u_x;
					v_current[index2D(z,y)] = u_y;
					w_current[index2D(z,y)] = u_z;

					//the following 3 lines correspond to time_unit==0!!!
					u_previous_temporal_boundary[index2D(z,y)] = u_current[index2D(z,y)];
					v_previous_temporal_boundary[index2D(z,y)] = v_current[index2D(z,y)];
					w_previous_temporal_boundary[index2D(z,y)] = w_current[index2D(z,y)];

				}//if (x == lx-2)
			}//for (x = 0 ; x< lx; ++x)
		}//for (y = 0 ; y< ly ; ++y)
	}//for (z = 0 ; z< lz ; ++z)
#ifdef DEBUG
	cout << " #LBM relaxation OK!" << endl;
#endif
}

//__global__
//void relaxation_kernel(int lx, int ly, int lz, FLOATING reynolds, FLOATING nu, FLOATING r_small,
//		FLOATING t_0, FLOATING t_1, FLOATING t_2, FLOATING c_squ, FLOATING omega, FLOATING one_minus_omega,
//		FLOATING reciprocal_c_squ,lattice D3, lattice D3_hlp, int *obstacles_d,
//		FLOATING *u_previous_spatial_boundary, FLOATING *v_previous_spatial_boundary, FLOATING *w_previous_spatial_boundary,
//		FLOATING *u_current, FLOATING *v_current, FLOATING *w_current){
//
//
//
//
//
//	/*One-step density relaxation process
//
//				c.......density relaxation: a single time relaxation with relaxation
//				c       parameter omega is applied here. This step is only "local",
//				c       nothing is propagated through the lattice.
//				c*/
//
//	int  x,y,z;
//	FLOATING  u_x=0.0, u_y=0.0, u_z=0.0, u_squ=0.0, rho=0.0, reciprocal_rho=0.0;
//
//
//
//
//	FLOATING u_n[19];
//	//	FLOATING n_equ[19];
//	FLOATING buff[19];
//
//	//	FLOATING u_n_squared[19];
//	//FLOATING two_x_c_squ_sqared;
//	const FLOATING omega_x_t_0=omega*t_0, omega_x_t_1=omega*t_1, omega_x_t_2=omega*t_2;
//	FLOATING omega_x_rho_x_t_0, omega_x_rho_x_t_1, omega_x_rho_x_t_2;
//	FLOATING temp_factor;
//	FLOATING u_n__over__c_squ[19];
//	FLOATING u_n__over__c_squ__squared_and_halved[19];
//
//	//....square speed of sound
//	/*      compute the out let velocity with a convevtive boundary condition
//					c.....loop over all nodes
//					c.....attention: actual densities are stored after the propagation
//					c                step in the help-array n_hlp !*/
//
//
//#pragma unroll
//	for (z = 0 ; z< lz ; ++z){
//#pragma unroll
//		for (y = 0 ; y< ly ; ++y){
//#pragma unroll
//			for (x = 0 ; x< lx; ++x){
//
//				/*c.........only free nodes are considered here
//					!if (.not. obstacles[z][y][x]) then
//					c...........integral local density
//					c...........initialize variable ro*/
//				//memory optimised implementation
//				buff[0]=D3_hlp.Q0[index(z,y,x)];
//				buff[1]=D3_hlp.Q1[index(z,y,x)];
//				buff[2]=D3_hlp.Q2[index(z,y,x)];
//				buff[3]=D3_hlp.Q3[index(z,y,x)];
//				buff[4]=D3_hlp.Q4[index(z,y,x)];
//				buff[5]=D3_hlp.Q5[index(z,y,x)];
//				buff[6]=D3_hlp.Q6[index(z,y,x)];
//				buff[7]=D3_hlp.Q7[index(z,y,x)];
//				buff[8]=D3_hlp.Q8[index(z,y,x)];
//				buff[9]=D3_hlp.Q9[index(z,y,x)];
//				buff[10]=D3_hlp.Q10[index(z,y,x)];
//				buff[11]=D3_hlp.Q11[index(z,y,x)];
//				buff[12]=D3_hlp.Q12[index(z,y,x)];
//				buff[13]=D3_hlp.Q13[index(z,y,x)];
//				buff[14]=D3_hlp.Q14[index(z,y,x)];
//				buff[15]=D3_hlp.Q15[index(z,y,x)];
//				buff[16]=D3_hlp.Q16[index(z,y,x)];
//				buff[17]=D3_hlp.Q17[index(z,y,x)];
//				buff[18]=D3_hlp.Q18[index(z,y,x)];
//
//				rho=0.0;
//				for(int k=0; k<DENSITIES; ++k)
//					rho+=buff[k];
//
//				//	rho=accumulate(buff, buff+DENSITIES, 0.0);
//
//				reciprocal_rho=1.0/rho;
//
//
//
//
//				switch(obstacles_d[index(z,y,x)]){
//				case 1:
//					u_x = 0.0;
//					u_y = 0.0;
//					u_z = 0.0;
//					break;
//				default:
//					u_x = 0.0;
//					u_x =  reciprocal_rho*(buff[1] + buff[7] + buff[10] +buff[11] + buff[14]-
//							(buff[3] + buff[8] + buff[9] +buff[12] + buff[13]));
//
//					u_y = 0.0;
//					u_y =  reciprocal_rho*(buff[2]+buff[8]+buff[7]+buff[16] + buff[15] -
//							(buff[4] + buff[9] + buff[10] +buff[17] + buff[18]));
//
//					u_z = 0.0;
//					u_z =  reciprocal_rho*(buff[5]+buff[13]+buff[14]+buff[15]+buff[18]-
//							(buff[6]+buff[12]+buff[11]+buff[16]+buff[17]));
//					break;
//				}//switch(obstacles[index(z,y,x)])
//
//				//original implementation
//				//				rho=0.0;
//				//				rho+=D3_hlp.Q0[index(z,y,x)]+D3_hlp.Q1[index(z,y,x)]+D3_hlp.Q2[index(z,y,x)]+D3_hlp.Q3[index(z,y,x)];
//				//				rho+=D3_hlp.Q4[index(z,y,x)]+D3_hlp.Q5[index(z,y,x)]+D3_hlp.Q6[index(z,y,x)]+D3_hlp.Q7[index(z,y,x)];
//				//				rho+=D3_hlp.Q8[index(z,y,x)]+D3_hlp.Q9[index(z,y,x)]+D3_hlp.Q10[index(z,y,x)]+D3_hlp.Q11[index(z,y,x)];
//				//				rho+=D3_hlp.Q12[index(z,y,x)]+D3_hlp.Q13[index(z,y,x)]+D3_hlp.Q14[index(z,y,x)]+D3_hlp.Q15[index(z,y,x)];
//				//				rho+=D3_hlp.Q16[index(z,y,x)]+D3_hlp.Q17[index(z,y,x)]+D3_hlp.Q18[index(z,y,x)];
//				//				reciprocal_rho=1.0/rho;
//
//				//...........x-, and y- velocity components
//
//				//				switch(obstacles[index(z,y,x)]){
//				//				case 1:
//				//					u_x = 0.0;
//				//					u_y = 0.0;
//				//					u_z = 0.0;
//				//					break;
//				//				default:
//				//					u_x = (FLOATING) reciprocal_rho*(D3_hlp.Q1[index(z,y,x)] + D3_hlp.Q7[index(z,y,x)] + D3_hlp.Q10[index(z,y,x)] +
//				//							D3_hlp.Q11[index(z,y,x)] + D3_hlp.Q14[index(z,y,x)] -
//				//							(D3_hlp.Q3[index(z,y,x)] + D3_hlp.Q8[index(z,y,x)] + D3_hlp.Q9[index(z,y,x)] +
//				//									D3_hlp.Q12[index(z,y,x)] + D3_hlp.Q13[index(z,y,x)]));
//				//
//				//					u_y = (FLOATING) reciprocal_rho*(D3_hlp.Q2[index(z,y,x)] + D3_hlp.Q8[index(z,y,x)] + D3_hlp.Q7[index(z,y,x)] +
//				//							D3_hlp.Q16[index(z,y,x)] + D3_hlp.Q15[index(z,y,x)] -
//				//							(D3_hlp.Q4[index(z,y,x)] + D3_hlp.Q9[index(z,y,x)] + D3_hlp.Q10[index(z,y,x)] +
//				//									D3_hlp.Q17[index(z,y,x)] + D3_hlp.Q18[index(z,y,x)]));
//				//
//				//					u_z = (FLOATING) reciprocal_rho*(D3_hlp.Q5[index(z,y,x)] + D3_hlp.Q13[index(z,y,x)] + D3_hlp.Q14[index(z,y,x)] +
//				//							D3_hlp.Q15[index(z,y,x)] + D3_hlp.Q18[index(z,y,x)] -
//				//							(D3_hlp.Q6[index(z,y,x)] + D3_hlp.Q12[index(z,y,x)] + D3_hlp.Q11[index(z,y,x)] +
//				//									D3_hlp.Q16[index(z,y,x)] + D3_hlp.Q17[index(z,y,x)]));
//				//					break;
//				//				}//switch(obstacles[index(z,y,x)])
//
//
//
//
//				u_squ = (FLOATING)  u_x*u_x + u_y*u_y + u_z*u_z;
//				temp_factor= 0.5*(2.0* c_squ - u_squ)/c_squ;
//				//u_squ = (FLOATING)  pow(u_x,2) + pow(u_y,2) + pow(u_z,2);
//
//
//				/*...........n- velocity compnents (n = lattice node connection vectors)
//					c...........this is only necessary for clearence, and only 3 speeds would
//					c...........be necessary*/
//
//
//				//WARNING!!!! o pinakas autos exei tropopoihmena indices!!!!
//				u_n[0]= 0.0; //SHOULD NEVER USED!
//				u_n[1] =   u_x;
//				u_n[2] =         u_y;
//				u_n[3] = - u_x;
//				u_n[4] =       - u_y;
//				u_n[5] =   u_z;
//				u_n[6] =       - u_z;
//				u_n[7] =   u_x + u_y;
//				u_n[8] = - u_x + u_y;
//				u_n[9] = - u_x - u_y;
//				u_n[10] =   u_x - u_y;
//				u_n[11] =   u_x - u_z;
//				u_n[12] = - u_x - u_z;
//				u_n[13] = - u_x + u_z;
//				u_n[14] =   u_x + u_z;
//				u_n[15] =   u_z + u_y;
//				u_n[16] = - u_z + u_y;
//				u_n[17] = - u_z - u_y;
//				u_n[18] =   u_z - u_y;
//
//#pragma unroll
//				for(int i=0; i<DENSITIES; ++i){
//					u_n__over__c_squ[i]=reciprocal_c_squ*u_n[i];
//					u_n__over__c_squ__squared_and_halved[i]=0.5*u_n__over__c_squ[i]*u_n__over__c_squ[i];
//				}
//
//				/*c...........equilibrium densities
//					c...........this can be rewritten to improve computational performance
//					c...........considerabely !
//					c
//					c...........zero velocity density
//					c*/
//				//memory optimised implementation! WARNING!!! different from the original case!
//
//				//two_x_c_squ_sqared=2.0*c_squ*c_squ;
//				omega_x_rho_x_t_0=omega_x_t_0*rho;
//				omega_x_rho_x_t_1=omega_x_t_1*rho;
//				omega_x_rho_x_t_2=omega_x_t_2*rho;
//
//
//				//				//...........relaxation step
//
//
//				//omega_x_rho_x_t_0*(1.0 - 0.5*u_squ/c_squ);
//				D3.Q0[index(z,y,x)]=buff[0]*one_minus_omega+omega_x_rho_x_t_0*(u_n__over__c_squ__squared_and_halved[0]+u_n__over__c_squ[0]+temp_factor);
//
//
//
//
//				D3.Q1[index(z,y,x)]=buff[1]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[1]+u_n__over__c_squ[1]+temp_factor);
//				D3.Q2[index(z,y,x)]=buff[2]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[2]+u_n__over__c_squ[2]+temp_factor);
//				D3.Q3[index(z,y,x)]=buff[3]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[3]+u_n__over__c_squ[3]+temp_factor);
//				D3.Q4[index(z,y,x)]=buff[4]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[4]+u_n__over__c_squ[4]+temp_factor);
//				D3.Q5[index(z,y,x)]=buff[5]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[5]+u_n__over__c_squ[5]+temp_factor);
//				D3.Q6[index(z,y,x)]=buff[6]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[6]+u_n__over__c_squ[6]+temp_factor);
//
//
//
//				D3.Q7[index(z,y,x)]= buff[ 7]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[ 7]+u_n__over__c_squ[ 7]+temp_factor);
//				D3.Q8[index(z,y,x)]= buff[ 8]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[ 8]+u_n__over__c_squ[ 8]+temp_factor);
//				D3.Q9[index(z,y,x)]= buff[ 9]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[ 9]+u_n__over__c_squ[ 9]+temp_factor);
//				D3.Q10[index(z,y,x)]=buff[10]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[10]+u_n__over__c_squ[10]+temp_factor);
//				D3.Q11[index(z,y,x)]=buff[11]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[11]+u_n__over__c_squ[11]+temp_factor);
//				D3.Q12[index(z,y,x)]=buff[12]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[12]+u_n__over__c_squ[12]+temp_factor);
//				D3.Q13[index(z,y,x)]=buff[13]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[13]+u_n__over__c_squ[13]+temp_factor);
//				D3.Q14[index(z,y,x)]=buff[14]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[14]+u_n__over__c_squ[14]+temp_factor);
//				D3.Q15[index(z,y,x)]=buff[15]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[15]+u_n__over__c_squ[15]+temp_factor);
//				D3.Q16[index(z,y,x)]=buff[16]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[16]+u_n__over__c_squ[16]+temp_factor);
//				D3.Q17[index(z,y,x)]=buff[17]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[17]+u_n__over__c_squ[17]+temp_factor);
//				D3.Q18[index(z,y,x)]=buff[18]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[18]+u_n__over__c_squ[18]+temp_factor);
//
//				//original implementation
//				//				n_equ[0] = t_0  * rho*(1.0 - u_squ / (2.0 * c_squ));
//				//
//				//								//...........axis speeds (factor: t_1)
//				//				#pragma unroll
//				//								for (int i = 1 ; i< 7; ++i){
//				//									n_equ[i] = t_1 * rho*(1.0 + u_n[i] / c_squ
//				//											+ (u_n[i]*u_n[i]) / (2.0 * c_squ *c_squ)
//				//											- u_squ / (2.0  * c_squ));
//				//								}
//				//
//				//								//...........diagonal speeds (factor: t_2)
//				//				#pragma unroll
//				//								for (int i = 7 ; i< 19; ++i){
//				//									n_equ[i] = t_2  * rho*(1.0 + u_n[i] / c_squ
//				//											+ (u_n[i]*u_n[i]) / (2.0 * c_squ *c_squ)
//				//											- u_squ / (2.0  * c_squ));
//				//								}
//
//				//				D3.Q0[index(z,y,x)]=D3_hlp.Q0[index(z,y,x)]+omega*(n_equ[0] - D3_hlp.Q0[index(z,y,x)]);
//				//								D3.Q1[index(z,y,x)]=D3_hlp.Q1[index(z,y,x)]+omega*(n_equ[1] - D3_hlp.Q1[index(z,y,x)]);
//				//								D3.Q2[index(z,y,x)]=D3_hlp.Q2[index(z,y,x)]+omega*(n_equ[2] - D3_hlp.Q2[index(z,y,x)]);
//				//								D3.Q3[index(z,y,x)]=D3_hlp.Q3[index(z,y,x)]+omega*(n_equ[3] - D3_hlp.Q3[index(z,y,x)]);
//				//								D3.Q4[index(z,y,x)]=D3_hlp.Q4[index(z,y,x)]+omega*(n_equ[4] - D3_hlp.Q4[index(z,y,x)]);
//				//								D3.Q5[index(z,y,x)]=D3_hlp.Q5[index(z,y,x)]+omega*(n_equ[5] - D3_hlp.Q5[index(z,y,x)]);
//				//								D3.Q6[index(z,y,x)]=D3_hlp.Q6[index(z,y,x)]+omega*(n_equ[6] - D3_hlp.Q6[index(z,y,x)]);
//				//								D3.Q7[index(z,y,x)]=D3_hlp.Q7[index(z,y,x)]+omega*(n_equ[7] - D3_hlp.Q7[index(z,y,x)]);
//				//								D3.Q8[index(z,y,x)]=D3_hlp.Q8[index(z,y,x)]+omega*(n_equ[8] - D3_hlp.Q8[index(z,y,x)]);
//				//								D3.Q9[index(z,y,x)]=D3_hlp.Q9[index(z,y,x)]+omega*(n_equ[9] - D3_hlp.Q9[index(z,y,x)]);
//				//								D3.Q10[index(z,y,x)]=D3_hlp.Q10[index(z,y,x)]+omega*(n_equ[10] - D3_hlp.Q10[index(z,y,x)]);
//				//								D3.Q11[index(z,y,x)]=D3_hlp.Q11[index(z,y,x)]+omega*(n_equ[11] - D3_hlp.Q11[index(z,y,x)]);
//				//								D3.Q12[index(z,y,x)]=D3_hlp.Q12[index(z,y,x)]+omega*(n_equ[12] - D3_hlp.Q12[index(z,y,x)]);
//				//								D3.Q13[index(z,y,x)]=D3_hlp.Q13[index(z,y,x)]+omega*(n_equ[13] - D3_hlp.Q13[index(z,y,x)]);
//				//								D3.Q14[index(z,y,x)]=D3_hlp.Q14[index(z,y,x)]+omega*(n_equ[14] - D3_hlp.Q14[index(z,y,x)]);
//				//								D3.Q15[index(z,y,x)]=D3_hlp.Q15[index(z,y,x)]+omega*(n_equ[15] - D3_hlp.Q15[index(z,y,x)]);
//				//								D3.Q16[index(z,y,x)]=D3_hlp.Q16[index(z,y,x)]+omega*(n_equ[16] - D3_hlp.Q16[index(z,y,x)]);
//				//								D3.Q17[index(z,y,x)]=D3_hlp.Q17[index(z,y,x)]+omega*(n_equ[17] - D3_hlp.Q17[index(z,y,x)]);
//				//								D3.Q18[index(z,y,x)]=D3_hlp.Q18[index(z,y,x)]+omega*(n_equ[18] - D3_hlp.Q18[index(z,y,x)]);
//
//				if (x == lx-2) {
//					u_previous_spatial_boundary[index2D(z,y)] = u_x;
//					v_previous_spatial_boundary[index2D(z,y)] = u_y;
//					w_previous_spatial_boundary[index2D(z,y)] = u_z;
//
//					u_current[index2D(z,y)] = u_x;
//					v_current[index2D(z,y)] = u_y;
//					w_current[index2D(z,y)] = u_z;
//
//				}//if (x == lx-2)
//			}//for (x = 0 ; x< lx; ++x)
//		}//for (y = 0 ; y< ly ; ++y)
//	}//for (z = 0 ; z< lz ; ++z)
//#ifdef DEBUG
//	cout << " #LBM relaxation OK!" << endl;
//#endif
//}
//
//__global__
//void relaxation_kernel_v2(int lx, int ly, int lz, FLOATING reynolds, FLOATING nu, FLOATING r_small,
//		FLOATING t_0, FLOATING t_1, FLOATING t_2, FLOATING c_squ, FLOATING omega, FLOATING one_minus_omega,
//		FLOATING reciprocal_c_squ,
//		FLOATING *hlp_Q0, FLOATING *hlp_Q1, FLOATING *hlp_Q2, FLOATING *hlp_Q3,
//		FLOATING *hlp_Q4, FLOATING *hlp_Q5, FLOATING *hlp_Q6, FLOATING *hlp_Q7,
//		FLOATING *hlp_Q8, FLOATING *hlp_Q9, FLOATING *hlp_Q10, FLOATING *hlp_Q11,
//		FLOATING *hlp_Q12, FLOATING *hlp_Q13, FLOATING *hlp_Q14, FLOATING *hlp_Q15,
//		FLOATING *hlp_Q16, FLOATING *hlp_Q17, FLOATING *hlp_Q18,
//		FLOATING *Q0, FLOATING *Q1, FLOATING *Q2, FLOATING *Q3,
//		FLOATING *Q4, FLOATING *Q5, FLOATING *Q6, FLOATING *Q7,
//		FLOATING *Q8, FLOATING *Q9, FLOATING *Q10, FLOATING *Q11,
//		FLOATING *Q12, FLOATING *Q13, FLOATING *Q14, FLOATING *Q15,
//		FLOATING *Q16, FLOATING *Q17, FLOATING *Q18,
//		int *obstacles_d,
//		FLOATING *u_previous_spatial_boundary, FLOATING *v_previous_spatial_boundary, FLOATING *w_previous_spatial_boundary,
//		FLOATING *u_current, FLOATING *v_current, FLOATING *w_current){
//
//
//
//
//
//	/*One-step density relaxation process
//
//				c.......density relaxation: a single time relaxation with relaxation
//				c       parameter omega is applied here. This step is only "local",
//				c       nothing is propagated through the lattice.
//				c*/
//
//	int  x,y,z;
//	FLOATING  u_x=0.0, u_y=0.0, u_z=0.0, u_squ=0.0, rho=0.0, reciprocal_rho=0.0;
//
//
//
//
//	FLOATING u_n[19];
//	//	FLOATING n_equ[19];
//	FLOATING buff[19];
//
//	//	FLOATING u_n_squared[19];
//	//FLOATING two_x_c_squ_sqared;
//	const FLOATING omega_x_t_0=omega*t_0, omega_x_t_1=omega*t_1, omega_x_t_2=omega*t_2;
//	FLOATING omega_x_rho_x_t_0, omega_x_rho_x_t_1, omega_x_rho_x_t_2;
//	FLOATING temp_factor;
//	FLOATING u_n__over__c_squ[19];
//	FLOATING u_n__over__c_squ__squared_and_halved[19];
//
//	//....square speed of sound
//	/*      compute the out let velocity with a convevtive boundary condition
//					c.....loop over all nodes
//					c.....attention: actual densities are stored after the propagation
//					c                step in the help-array n_hlp !*/
//
//	//
//	//#pragma unroll
//	//	for (z = 0 ; z< lz ; ++z){
//	//#pragma unroll
//	//		for (y = 0 ; y< ly ; ++y){
//	//#pragma unroll
//	//			for (x = 0 ; x< lx; ++x){
//
//	const int tid=blockIdx.x*blockDim.x+threadIdx.x;
//	int rest;
//	int end_of_memory=lz*ly*(lx);
//
//	z=(int) (tid/(ly*lx));
//	rest=tid-z;
//	y=(int)(rest/lx);
//	x=rest-y;
//
//	if (tid<end_of_memory){
//		/*c.........only free nodes are considered here
//					!if (.not. obstacles[z][y][x]) then
//					c...........integral local density
//					c...........initialize variable ro*/
//		//memory optimised implementation
//		buff[0]=hlp_Q0[index(z,y,x)];
//		buff[1]=hlp_Q1[index(z,y,x)];
//		buff[2]=hlp_Q2[index(z,y,x)];
//		buff[3]=hlp_Q3[index(z,y,x)];
//		buff[4]=hlp_Q4[index(z,y,x)];
//		buff[5]=hlp_Q5[index(z,y,x)];
//		buff[6]=hlp_Q6[index(z,y,x)];
//		buff[7]=hlp_Q7[index(z,y,x)];
//		buff[8]=hlp_Q8[index(z,y,x)];
//		buff[9]=hlp_Q9[index(z,y,x)];
//		buff[10]=hlp_Q10[index(z,y,x)];
//		buff[11]=hlp_Q11[index(z,y,x)];
//		buff[12]=hlp_Q12[index(z,y,x)];
//		buff[13]=hlp_Q13[index(z,y,x)];
//		buff[14]=hlp_Q14[index(z,y,x)];
//		buff[15]=hlp_Q15[index(z,y,x)];
//		buff[16]=hlp_Q16[index(z,y,x)];
//		buff[17]=hlp_Q17[index(z,y,x)];
//		buff[18]=hlp_Q18[index(z,y,x)];
//
//		rho=0.0;
//		for(int k=0; k<DENSITIES; ++k)
//			rho+=buff[k];
//
//		//	rho=accumulate(buff, buff+DENSITIES, 0.0);
//
//		reciprocal_rho=1.0/rho;
//
//
//
//
//		switch(obstacles_d[index(z,y,x)]){
//		case 1:
//			u_x = 0.0;
//			u_y = 0.0;
//			u_z = 0.0;
//			break;
//		default:
//			u_x = 0.0;
//			u_x =  reciprocal_rho*(buff[1] + buff[7] + buff[10] +buff[11] + buff[14]-
//					(buff[3] + buff[8] + buff[9] +buff[12] + buff[13]));
//
//			u_y = 0.0;
//			u_y =  reciprocal_rho*(buff[2]+buff[8]+buff[7]+buff[16] + buff[15] -
//					(buff[4] + buff[9] + buff[10] +buff[17] + buff[18]));
//
//			u_z = 0.0;
//			u_z =  reciprocal_rho*(buff[5]+buff[13]+buff[14]+buff[15]+buff[18]-
//					(buff[6]+buff[12]+buff[11]+buff[16]+buff[17]));
//			break;
//		}//switch(obstacles[index(z,y,x)])
//
//		//original implementation
//		//				rho=0.0;
//		//				rho+=D3_hlp.Q0[index(z,y,x)]+D3_hlp.Q1[index(z,y,x)]+D3_hlp.Q2[index(z,y,x)]+D3_hlp.Q3[index(z,y,x)];
//		//				rho+=D3_hlp.Q4[index(z,y,x)]+D3_hlp.Q5[index(z,y,x)]+D3_hlp.Q6[index(z,y,x)]+D3_hlp.Q7[index(z,y,x)];
//		//				rho+=D3_hlp.Q8[index(z,y,x)]+D3_hlp.Q9[index(z,y,x)]+D3_hlp.Q10[index(z,y,x)]+D3_hlp.Q11[index(z,y,x)];
//		//				rho+=D3_hlp.Q12[index(z,y,x)]+D3_hlp.Q13[index(z,y,x)]+D3_hlp.Q14[index(z,y,x)]+D3_hlp.Q15[index(z,y,x)];
//		//				rho+=D3_hlp.Q16[index(z,y,x)]+D3_hlp.Q17[index(z,y,x)]+D3_hlp.Q18[index(z,y,x)];
//		//				reciprocal_rho=1.0/rho;
//
//		//...........x-, and y- velocity components
//
//		//				switch(obstacles[index(z,y,x)]){
//		//				case 1:
//		//					u_x = 0.0;
//		//					u_y = 0.0;
//		//					u_z = 0.0;
//		//					break;
//		//				default:
//		//					u_x = (FLOATING) reciprocal_rho*(D3_hlp.Q1[index(z,y,x)] + D3_hlp.Q7[index(z,y,x)] + D3_hlp.Q10[index(z,y,x)] +
//		//							D3_hlp.Q11[index(z,y,x)] + D3_hlp.Q14[index(z,y,x)] -
//		//							(D3_hlp.Q3[index(z,y,x)] + D3_hlp.Q8[index(z,y,x)] + D3_hlp.Q9[index(z,y,x)] +
//		//									D3_hlp.Q12[index(z,y,x)] + D3_hlp.Q13[index(z,y,x)]));
//		//
//		//					u_y = (FLOATING) reciprocal_rho*(D3_hlp.Q2[index(z,y,x)] + D3_hlp.Q8[index(z,y,x)] + D3_hlp.Q7[index(z,y,x)] +
//		//							D3_hlp.Q16[index(z,y,x)] + D3_hlp.Q15[index(z,y,x)] -
//		//							(D3_hlp.Q4[index(z,y,x)] + D3_hlp.Q9[index(z,y,x)] + D3_hlp.Q10[index(z,y,x)] +
//		//									D3_hlp.Q17[index(z,y,x)] + D3_hlp.Q18[index(z,y,x)]));
//		//
//		//					u_z = (FLOATING) reciprocal_rho*(D3_hlp.Q5[index(z,y,x)] + D3_hlp.Q13[index(z,y,x)] + D3_hlp.Q14[index(z,y,x)] +
//		//							D3_hlp.Q15[index(z,y,x)] + D3_hlp.Q18[index(z,y,x)] -
//		//							(D3_hlp.Q6[index(z,y,x)] + D3_hlp.Q12[index(z,y,x)] + D3_hlp.Q11[index(z,y,x)] +
//		//									D3_hlp.Q16[index(z,y,x)] + D3_hlp.Q17[index(z,y,x)]));
//		//					break;
//		//				}//switch(obstacles[index(z,y,x)])
//
//
//
//
//		u_squ = (FLOATING)  u_x*u_x + u_y*u_y + u_z*u_z;
//		temp_factor= 0.5*(2.0* c_squ - u_squ)/c_squ;
//		//u_squ = (FLOATING)  pow(u_x,2) + pow(u_y,2) + pow(u_z,2);
//
//
//		/*...........n- velocity compnents (n = lattice node connection vectors)
//					c...........this is only necessary for clearence, and only 3 speeds would
//					c...........be necessary*/
//
//
//		//WARNING!!!! o pinakas autos exei tropopoihmena indices!!!!
//		u_n[0]= 0.0; //SHOULD NEVER USED!
//		u_n[1] =   u_x;
//		u_n[2] =         u_y;
//		u_n[3] = - u_x;
//		u_n[4] =       - u_y;
//		u_n[5] =   u_z;
//		u_n[6] =       - u_z;
//		u_n[7] =   u_x + u_y;
//		u_n[8] = - u_x + u_y;
//		u_n[9] = - u_x - u_y;
//		u_n[10] =   u_x - u_y;
//		u_n[11] =   u_x - u_z;
//		u_n[12] = - u_x - u_z;
//		u_n[13] = - u_x + u_z;
//		u_n[14] =   u_x + u_z;
//		u_n[15] =   u_z + u_y;
//		u_n[16] = - u_z + u_y;
//		u_n[17] = - u_z - u_y;
//		u_n[18] =   u_z - u_y;
//
//#pragma unroll
//		for(int i=0; i<DENSITIES; ++i){
//			u_n__over__c_squ[i]=reciprocal_c_squ*u_n[i];
//			u_n__over__c_squ__squared_and_halved[i]=0.5*u_n__over__c_squ[i]*u_n__over__c_squ[i];
//		}
//
//		/*c...........equilibrium densities
//					c...........this can be rewritten to improve computational performance
//					c...........considerabely !
//					c
//					c...........zero velocity density
//					c*/
//		//memory optimised implementation! WARNING!!! different from the original case!
//
//		//two_x_c_squ_sqared=2.0*c_squ*c_squ;
//		omega_x_rho_x_t_0=omega_x_t_0*rho;
//		omega_x_rho_x_t_1=omega_x_t_1*rho;
//		omega_x_rho_x_t_2=omega_x_t_2*rho;
//
//
//		//				//...........relaxation step
//
//
//		//omega_x_rho_x_t_0*(1.0 - 0.5*u_squ/c_squ);
//		Q0[index(z,y,x)]=buff[0]*one_minus_omega+omega_x_rho_x_t_0*(u_n__over__c_squ__squared_and_halved[0]+u_n__over__c_squ[0]+temp_factor);
//
//
//
//
//		Q1[index(z,y,x)]=buff[1]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[1]+u_n__over__c_squ[1]+temp_factor);
//		Q2[index(z,y,x)]=buff[2]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[2]+u_n__over__c_squ[2]+temp_factor);
//		Q3[index(z,y,x)]=buff[3]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[3]+u_n__over__c_squ[3]+temp_factor);
//		Q4[index(z,y,x)]=buff[4]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[4]+u_n__over__c_squ[4]+temp_factor);
//		Q5[index(z,y,x)]=buff[5]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[5]+u_n__over__c_squ[5]+temp_factor);
//		Q6[index(z,y,x)]=buff[6]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[6]+u_n__over__c_squ[6]+temp_factor);
//
//
//
//		Q7[index(z,y,x)]= buff[ 7]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[ 7]+u_n__over__c_squ[ 7]+temp_factor);
//		Q8[index(z,y,x)]= buff[ 8]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[ 8]+u_n__over__c_squ[ 8]+temp_factor);
//		Q9[index(z,y,x)]= buff[ 9]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[ 9]+u_n__over__c_squ[ 9]+temp_factor);
//		Q10[index(z,y,x)]=buff[10]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[10]+u_n__over__c_squ[10]+temp_factor);
//		Q11[index(z,y,x)]=buff[11]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[11]+u_n__over__c_squ[11]+temp_factor);
//		Q12[index(z,y,x)]=buff[12]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[12]+u_n__over__c_squ[12]+temp_factor);
//		Q13[index(z,y,x)]=buff[13]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[13]+u_n__over__c_squ[13]+temp_factor);
//		Q14[index(z,y,x)]=buff[14]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[14]+u_n__over__c_squ[14]+temp_factor);
//		Q15[index(z,y,x)]=buff[15]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[15]+u_n__over__c_squ[15]+temp_factor);
//		Q16[index(z,y,x)]=buff[16]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[16]+u_n__over__c_squ[16]+temp_factor);
//		Q17[index(z,y,x)]=buff[17]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[17]+u_n__over__c_squ[17]+temp_factor);
//		Q18[index(z,y,x)]=buff[18]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[18]+u_n__over__c_squ[18]+temp_factor);
//
//		//original implementation
//		//				n_equ[0] = t_0  * rho*(1.0 - u_squ / (2.0 * c_squ));
//		//
//		//								//...........axis speeds (factor: t_1)
//		//				#pragma unroll
//		//								for (int i = 1 ; i< 7; ++i){
//		//									n_equ[i] = t_1 * rho*(1.0 + u_n[i] / c_squ
//		//											+ (u_n[i]*u_n[i]) / (2.0 * c_squ *c_squ)
//		//											- u_squ / (2.0  * c_squ));
//		//								}
//		//
//		//								//...........diagonal speeds (factor: t_2)
//		//				#pragma unroll
//		//								for (int i = 7 ; i< 19; ++i){
//		//									n_equ[i] = t_2  * rho*(1.0 + u_n[i] / c_squ
//		//											+ (u_n[i]*u_n[i]) / (2.0 * c_squ *c_squ)
//		//											- u_squ / (2.0  * c_squ));
//		//								}
//
//		//				D3.Q0[index(z,y,x)]=D3_hlp.Q0[index(z,y,x)]+omega*(n_equ[0] - D3_hlp.Q0[index(z,y,x)]);
//		//								D3.Q1[index(z,y,x)]=D3_hlp.Q1[index(z,y,x)]+omega*(n_equ[1] - D3_hlp.Q1[index(z,y,x)]);
//		//								D3.Q2[index(z,y,x)]=D3_hlp.Q2[index(z,y,x)]+omega*(n_equ[2] - D3_hlp.Q2[index(z,y,x)]);
//		//								D3.Q3[index(z,y,x)]=D3_hlp.Q3[index(z,y,x)]+omega*(n_equ[3] - D3_hlp.Q3[index(z,y,x)]);
//		//								D3.Q4[index(z,y,x)]=D3_hlp.Q4[index(z,y,x)]+omega*(n_equ[4] - D3_hlp.Q4[index(z,y,x)]);
//		//								D3.Q5[index(z,y,x)]=D3_hlp.Q5[index(z,y,x)]+omega*(n_equ[5] - D3_hlp.Q5[index(z,y,x)]);
//		//								D3.Q6[index(z,y,x)]=D3_hlp.Q6[index(z,y,x)]+omega*(n_equ[6] - D3_hlp.Q6[index(z,y,x)]);
//		//								D3.Q7[index(z,y,x)]=D3_hlp.Q7[index(z,y,x)]+omega*(n_equ[7] - D3_hlp.Q7[index(z,y,x)]);
//		//								D3.Q8[index(z,y,x)]=D3_hlp.Q8[index(z,y,x)]+omega*(n_equ[8] - D3_hlp.Q8[index(z,y,x)]);
//		//								D3.Q9[index(z,y,x)]=D3_hlp.Q9[index(z,y,x)]+omega*(n_equ[9] - D3_hlp.Q9[index(z,y,x)]);
//		//								D3.Q10[index(z,y,x)]=D3_hlp.Q10[index(z,y,x)]+omega*(n_equ[10] - D3_hlp.Q10[index(z,y,x)]);
//		//								D3.Q11[index(z,y,x)]=D3_hlp.Q11[index(z,y,x)]+omega*(n_equ[11] - D3_hlp.Q11[index(z,y,x)]);
//		//								D3.Q12[index(z,y,x)]=D3_hlp.Q12[index(z,y,x)]+omega*(n_equ[12] - D3_hlp.Q12[index(z,y,x)]);
//		//								D3.Q13[index(z,y,x)]=D3_hlp.Q13[index(z,y,x)]+omega*(n_equ[13] - D3_hlp.Q13[index(z,y,x)]);
//		//								D3.Q14[index(z,y,x)]=D3_hlp.Q14[index(z,y,x)]+omega*(n_equ[14] - D3_hlp.Q14[index(z,y,x)]);
//		//								D3.Q15[index(z,y,x)]=D3_hlp.Q15[index(z,y,x)]+omega*(n_equ[15] - D3_hlp.Q15[index(z,y,x)]);
//		//								D3.Q16[index(z,y,x)]=D3_hlp.Q16[index(z,y,x)]+omega*(n_equ[16] - D3_hlp.Q16[index(z,y,x)]);
//		//								D3.Q17[index(z,y,x)]=D3_hlp.Q17[index(z,y,x)]+omega*(n_equ[17] - D3_hlp.Q17[index(z,y,x)]);
//		//								D3.Q18[index(z,y,x)]=D3_hlp.Q18[index(z,y,x)]+omega*(n_equ[18] - D3_hlp.Q18[index(z,y,x)]);
//
//		if (x == lx-2) {
//			u_previous_spatial_boundary[index2D(z,y)] = u_x;
//			v_previous_spatial_boundary[index2D(z,y)] = u_y;
//			w_previous_spatial_boundary[index2D(z,y)] = u_z;
//
//			u_current[index2D(z,y)] = u_x;
//			v_current[index2D(z,y)] = u_y;
//			w_current[index2D(z,y)] = u_z;
//
//		}//if (x == lx-2)
//	}//if memory!
//	//			}//for (x = 0 ; x< lx; ++x)
//	//		}//for (y = 0 ; y< ly ; ++y)
//	//	}//for (z = 0 ; z< lz ; ++z)
//
//}
//
//__global__
//void relaxation_kernel_v3(int lx, int ly, int lz, FLOATING reynolds, FLOATING nu, FLOATING r_small,
//		FLOATING t_0, FLOATING t_1, FLOATING t_2, FLOATING c_squ, FLOATING omega, FLOATING one_minus_omega,
//		FLOATING reciprocal_c_squ,
//		FLOATING *hlp_Q0, FLOATING *hlp_Q1, FLOATING *hlp_Q2, FLOATING *hlp_Q3,
//		FLOATING *hlp_Q4, FLOATING *hlp_Q5, FLOATING *hlp_Q6, FLOATING *hlp_Q7,
//		FLOATING *hlp_Q8, FLOATING *hlp_Q9, FLOATING *hlp_Q10, FLOATING *hlp_Q11,
//		FLOATING *hlp_Q12, FLOATING *hlp_Q13, FLOATING *hlp_Q14, FLOATING *hlp_Q15,
//		FLOATING *hlp_Q16, FLOATING *hlp_Q17, FLOATING *hlp_Q18,
//		FLOATING *Q0, FLOATING *Q1, FLOATING *Q2, FLOATING *Q3,
//		FLOATING *Q4, FLOATING *Q5, FLOATING *Q6, FLOATING *Q7,
//		FLOATING *Q8, FLOATING *Q9, FLOATING *Q10, FLOATING *Q11,
//		FLOATING *Q12, FLOATING *Q13, FLOATING *Q14, FLOATING *Q15,
//		FLOATING *Q16, FLOATING *Q17, FLOATING *Q18,
//		int *obstacles,
//		FLOATING *u_previous_spatial_boundary, FLOATING *v_previous_spatial_boundary, FLOATING *w_previous_spatial_boundary,
//		FLOATING *u_current, FLOATING *v_current, FLOATING *w_current){
//
//
//
//
//
//	/*One-step density relaxation process
//
//				c.......density relaxation: a single time relaxation with relaxation
//				c       parameter omega is applied here. This step is only "local",
//				c       nothing is propagated through the lattice.
//				c*/
//
//	int  x,y,z;
//	FLOATING  u_x=0.0, u_y=0.0, u_z=0.0, u_squ=0.0, rho=0.0, reciprocal_rho=0.0;
//
//
//
//
//	FLOATING u_n[19];
//	//	FLOATING n_equ[19];
//	FLOATING buff[19];
//
//	//	FLOATING u_n_squared[19];
//	//FLOATING two_x_c_squ_sqared;
//	const FLOATING omega_x_t_0=omega*t_0, omega_x_t_1=omega*t_1, omega_x_t_2=omega*t_2;
//	FLOATING omega_x_rho_x_t_0, omega_x_rho_x_t_1, omega_x_rho_x_t_2;
//	FLOATING temp_factor;
//	FLOATING u_n__over__c_squ[19];
//	FLOATING u_n__over__c_squ__squared_and_halved[19];
//
//	//....square speed of sound
//	/*      compute the out let velocity with a convevtive boundary condition
//					c.....loop over all nodes
//					c.....attention: actual densities are stored after the propagation
//					c                step in the help-array n_hlp !*/
//
//	//
//	//#pragma unroll
//	//	for (z = 0 ; z< lz ; ++z){
//	//#pragma unroll
//	//		for (y = 0 ; y< ly ; ++y){
//	//#pragma unroll
//	//			for (x = 0 ; x< lx; ++x){
//
//	const int tid=blockIdx.x*blockDim.x+threadIdx.x;
//	int rest;
//	int end_of_memory=lz*ly*(lx);
//
//	z=(int) (tid/(ly*lx));
//	rest=tid-z;
//	y=(int)(rest/lx);
//	x=rest-y;
//
//	if (tid<end_of_memory){
//		/*c.........only free nodes are considered here
//					!if (.not. obstacles[z][y][x]) then
//					c...........integral local density
//					c...........initialize variable ro*/
//		//memory optimised implementation
//		buff[0]=hlp_Q0[tid];
//		buff[1]=hlp_Q1[tid];
//		buff[2]=hlp_Q2[tid];
//		buff[3]=hlp_Q3[tid];
//		buff[4]=hlp_Q4[tid];
//		buff[5]=hlp_Q5[tid];
//		buff[6]=hlp_Q6[tid];
//		buff[7]=hlp_Q7[tid];
//		buff[8]=hlp_Q8[tid];
//		buff[9]=hlp_Q9[tid];
//		buff[10]=hlp_Q10[tid];
//		buff[11]=hlp_Q11[tid];
//		buff[12]=hlp_Q12[tid];
//		buff[13]=hlp_Q13[tid];
//		buff[14]=hlp_Q14[tid];
//		buff[15]=hlp_Q15[tid];
//		buff[16]=hlp_Q16[tid];
//		buff[17]=hlp_Q17[tid];
//		buff[18]=hlp_Q18[tid];
//
//		rho=0.0;
//		for(int k=0; k<DENSITIES; ++k)
//			rho+=buff[k];
//
//		//	rho=accumulate(buff, buff+DENSITIES, 0.0);
//
//		reciprocal_rho=1.0/rho;
//
//
//
//
//		switch(obstacles[tid]){
//		case 1:
//			u_x = 0.0;
//			u_y = 0.0;
//			u_z = 0.0;
//			break;
//		default:
//			u_x = 0.0;
//			u_x =  reciprocal_rho*(buff[1] + buff[7] + buff[10] +buff[11] + buff[14]-
//					(buff[3] + buff[8] + buff[9] +buff[12] + buff[13]));
//
//			u_y = 0.0;
//			u_y =  reciprocal_rho*(buff[2]+buff[8]+buff[7]+buff[16] + buff[15] -
//					(buff[4] + buff[9] + buff[10] +buff[17] + buff[18]));
//
//			u_z = 0.0;
//			u_z =  reciprocal_rho*(buff[5]+buff[13]+buff[14]+buff[15]+buff[18]-
//					(buff[6]+buff[12]+buff[11]+buff[16]+buff[17]));
//			break;
//		}//switch(obstacles[index(z,y,x)])
//
//		//original implementation
//		//				rho=0.0;
//		//				rho+=D3_hlp.Q0[index(z,y,x)]+D3_hlp.Q1[index(z,y,x)]+D3_hlp.Q2[index(z,y,x)]+D3_hlp.Q3[index(z,y,x)];
//		//				rho+=D3_hlp.Q4[index(z,y,x)]+D3_hlp.Q5[index(z,y,x)]+D3_hlp.Q6[index(z,y,x)]+D3_hlp.Q7[index(z,y,x)];
//		//				rho+=D3_hlp.Q8[index(z,y,x)]+D3_hlp.Q9[index(z,y,x)]+D3_hlp.Q10[index(z,y,x)]+D3_hlp.Q11[index(z,y,x)];
//		//				rho+=D3_hlp.Q12[index(z,y,x)]+D3_hlp.Q13[index(z,y,x)]+D3_hlp.Q14[index(z,y,x)]+D3_hlp.Q15[index(z,y,x)];
//		//				rho+=D3_hlp.Q16[index(z,y,x)]+D3_hlp.Q17[index(z,y,x)]+D3_hlp.Q18[index(z,y,x)];
//		//				reciprocal_rho=1.0/rho;
//
//		//...........x-, and y- velocity components
//
//		//				switch(obstacles[index(z,y,x)]){
//		//				case 1:
//		//					u_x = 0.0;
//		//					u_y = 0.0;
//		//					u_z = 0.0;
//		//					break;
//		//				default:
//		//					u_x = (FLOATING) reciprocal_rho*(D3_hlp.Q1[index(z,y,x)] + D3_hlp.Q7[index(z,y,x)] + D3_hlp.Q10[index(z,y,x)] +
//		//							D3_hlp.Q11[index(z,y,x)] + D3_hlp.Q14[index(z,y,x)] -
//		//							(D3_hlp.Q3[index(z,y,x)] + D3_hlp.Q8[index(z,y,x)] + D3_hlp.Q9[index(z,y,x)] +
//		//									D3_hlp.Q12[index(z,y,x)] + D3_hlp.Q13[index(z,y,x)]));
//		//
//		//					u_y = (FLOATING) reciprocal_rho*(D3_hlp.Q2[index(z,y,x)] + D3_hlp.Q8[index(z,y,x)] + D3_hlp.Q7[index(z,y,x)] +
//		//							D3_hlp.Q16[index(z,y,x)] + D3_hlp.Q15[index(z,y,x)] -
//		//							(D3_hlp.Q4[index(z,y,x)] + D3_hlp.Q9[index(z,y,x)] + D3_hlp.Q10[index(z,y,x)] +
//		//									D3_hlp.Q17[index(z,y,x)] + D3_hlp.Q18[index(z,y,x)]));
//		//
//		//					u_z = (FLOATING) reciprocal_rho*(D3_hlp.Q5[index(z,y,x)] + D3_hlp.Q13[index(z,y,x)] + D3_hlp.Q14[index(z,y,x)] +
//		//							D3_hlp.Q15[index(z,y,x)] + D3_hlp.Q18[index(z,y,x)] -
//		//							(D3_hlp.Q6[index(z,y,x)] + D3_hlp.Q12[index(z,y,x)] + D3_hlp.Q11[index(z,y,x)] +
//		//									D3_hlp.Q16[index(z,y,x)] + D3_hlp.Q17[index(z,y,x)]));
//		//					break;
//		//				}//switch(obstacles[index(z,y,x)])
//
//
//
//
//		u_squ = (FLOATING)  u_x*u_x + u_y*u_y + u_z*u_z;
//		temp_factor= 0.5*(2.0* c_squ - u_squ)/c_squ;
//		//u_squ = (FLOATING)  pow(u_x,2) + pow(u_y,2) + pow(u_z,2);
//
//
//		/*...........n- velocity compnents (n = lattice node connection vectors)
//					c...........this is only necessary for clearence, and only 3 speeds would
//					c...........be necessary*/
//
//
//		//WARNING!!!! o pinakas autos exei tropopoihmena indices!!!!
//		u_n[0]= 0.0; //SHOULD NEVER USED!
//		u_n[1] =   u_x;
//		u_n[2] =         u_y;
//		u_n[3] = - u_x;
//		u_n[4] =       - u_y;
//		u_n[5] =   u_z;
//		u_n[6] =       - u_z;
//		u_n[7] =   u_x + u_y;
//		u_n[8] = - u_x + u_y;
//		u_n[9] = - u_x - u_y;
//		u_n[10] =   u_x - u_y;
//		u_n[11] =   u_x - u_z;
//		u_n[12] = - u_x - u_z;
//		u_n[13] = - u_x + u_z;
//		u_n[14] =   u_x + u_z;
//		u_n[15] =   u_z + u_y;
//		u_n[16] = - u_z + u_y;
//		u_n[17] = - u_z - u_y;
//		u_n[18] =   u_z - u_y;
//
//#pragma unroll
//		for(int i=0; i<DENSITIES; ++i){
//			u_n__over__c_squ[i]=reciprocal_c_squ*u_n[i];
//			u_n__over__c_squ__squared_and_halved[i]=0.5*u_n__over__c_squ[i]*u_n__over__c_squ[i];
//		}
//
//		/*c...........equilibrium densities
//					c...........this can be rewritten to improve computational performance
//					c...........considerabely !
//					c
//					c...........zero velocity density
//					c*/
//		//memory optimised implementation! WARNING!!! different from the original case!
//
//		//two_x_c_squ_sqared=2.0*c_squ*c_squ;
//		omega_x_rho_x_t_0=omega_x_t_0*rho;
//		omega_x_rho_x_t_1=omega_x_t_1*rho;
//		omega_x_rho_x_t_2=omega_x_t_2*rho;
//
//
//		//				//...........relaxation step
//
//
//		//omega_x_rho_x_t_0*(1.0 - 0.5*u_squ/c_squ);
//		Q0[tid]=buff[0]*one_minus_omega+omega_x_rho_x_t_0*(u_n__over__c_squ__squared_and_halved[0]+u_n__over__c_squ[0]+temp_factor);
//
//
//
//
//		Q1[tid]=buff[1]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[1]+u_n__over__c_squ[1]+temp_factor);
//		Q2[tid]=buff[2]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[2]+u_n__over__c_squ[2]+temp_factor);
//		Q3[tid]=buff[3]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[3]+u_n__over__c_squ[3]+temp_factor);
//		Q4[tid]=buff[4]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[4]+u_n__over__c_squ[4]+temp_factor);
//		Q5[tid]=buff[5]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[5]+u_n__over__c_squ[5]+temp_factor);
//		Q6[tid]=buff[6]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[6]+u_n__over__c_squ[6]+temp_factor);
//
//
//
//		Q7[tid]= buff[ 7]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[ 7]+u_n__over__c_squ[ 7]+temp_factor);
//		Q8[tid]= buff[ 8]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[ 8]+u_n__over__c_squ[ 8]+temp_factor);
//		Q9[tid]= buff[ 9]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[ 9]+u_n__over__c_squ[ 9]+temp_factor);
//		Q10[tid]=buff[10]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[10]+u_n__over__c_squ[10]+temp_factor);
//		Q11[tid]=buff[11]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[11]+u_n__over__c_squ[11]+temp_factor);
//		Q12[tid]=buff[12]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[12]+u_n__over__c_squ[12]+temp_factor);
//		Q13[tid]=buff[13]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[13]+u_n__over__c_squ[13]+temp_factor);
//		Q14[tid]=buff[14]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[14]+u_n__over__c_squ[14]+temp_factor);
//		Q15[tid]=buff[15]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[15]+u_n__over__c_squ[15]+temp_factor);
//		Q16[tid]=buff[16]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[16]+u_n__over__c_squ[16]+temp_factor);
//		Q17[tid]=buff[17]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[17]+u_n__over__c_squ[17]+temp_factor);
//		Q18[tid]=buff[18]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[18]+u_n__over__c_squ[18]+temp_factor);
//
//		//original implementation
//		//				n_equ[0] = t_0  * rho*(1.0 - u_squ / (2.0 * c_squ));
//		//
//		//								//...........axis speeds (factor: t_1)
//		//				#pragma unroll
//		//								for (int i = 1 ; i< 7; ++i){
//		//									n_equ[i] = t_1 * rho*(1.0 + u_n[i] / c_squ
//		//											+ (u_n[i]*u_n[i]) / (2.0 * c_squ *c_squ)
//		//											- u_squ / (2.0  * c_squ));
//		//								}
//		//
//		//								//...........diagonal speeds (factor: t_2)
//		//				#pragma unroll
//		//								for (int i = 7 ; i< 19; ++i){
//		//									n_equ[i] = t_2  * rho*(1.0 + u_n[i] / c_squ
//		//											+ (u_n[i]*u_n[i]) / (2.0 * c_squ *c_squ)
//		//											- u_squ / (2.0  * c_squ));
//		//								}
//
//		//				D3.Q0[index(z,y,x)]=D3_hlp.Q0[index(z,y,x)]+omega*(n_equ[0] - D3_hlp.Q0[index(z,y,x)]);
//		//								D3.Q1[index(z,y,x)]=D3_hlp.Q1[index(z,y,x)]+omega*(n_equ[1] - D3_hlp.Q1[index(z,y,x)]);
//		//								D3.Q2[index(z,y,x)]=D3_hlp.Q2[index(z,y,x)]+omega*(n_equ[2] - D3_hlp.Q2[index(z,y,x)]);
//		//								D3.Q3[index(z,y,x)]=D3_hlp.Q3[index(z,y,x)]+omega*(n_equ[3] - D3_hlp.Q3[index(z,y,x)]);
//		//								D3.Q4[index(z,y,x)]=D3_hlp.Q4[index(z,y,x)]+omega*(n_equ[4] - D3_hlp.Q4[index(z,y,x)]);
//		//								D3.Q5[index(z,y,x)]=D3_hlp.Q5[index(z,y,x)]+omega*(n_equ[5] - D3_hlp.Q5[index(z,y,x)]);
//		//								D3.Q6[index(z,y,x)]=D3_hlp.Q6[index(z,y,x)]+omega*(n_equ[6] - D3_hlp.Q6[index(z,y,x)]);
//		//								D3.Q7[index(z,y,x)]=D3_hlp.Q7[index(z,y,x)]+omega*(n_equ[7] - D3_hlp.Q7[index(z,y,x)]);
//		//								D3.Q8[index(z,y,x)]=D3_hlp.Q8[index(z,y,x)]+omega*(n_equ[8] - D3_hlp.Q8[index(z,y,x)]);
//		//								D3.Q9[index(z,y,x)]=D3_hlp.Q9[index(z,y,x)]+omega*(n_equ[9] - D3_hlp.Q9[index(z,y,x)]);
//		//								D3.Q10[index(z,y,x)]=D3_hlp.Q10[index(z,y,x)]+omega*(n_equ[10] - D3_hlp.Q10[index(z,y,x)]);
//		//								D3.Q11[index(z,y,x)]=D3_hlp.Q11[index(z,y,x)]+omega*(n_equ[11] - D3_hlp.Q11[index(z,y,x)]);
//		//								D3.Q12[index(z,y,x)]=D3_hlp.Q12[index(z,y,x)]+omega*(n_equ[12] - D3_hlp.Q12[index(z,y,x)]);
//		//								D3.Q13[index(z,y,x)]=D3_hlp.Q13[index(z,y,x)]+omega*(n_equ[13] - D3_hlp.Q13[index(z,y,x)]);
//		//								D3.Q14[index(z,y,x)]=D3_hlp.Q14[index(z,y,x)]+omega*(n_equ[14] - D3_hlp.Q14[index(z,y,x)]);
//		//								D3.Q15[index(z,y,x)]=D3_hlp.Q15[index(z,y,x)]+omega*(n_equ[15] - D3_hlp.Q15[index(z,y,x)]);
//		//								D3.Q16[index(z,y,x)]=D3_hlp.Q16[index(z,y,x)]+omega*(n_equ[16] - D3_hlp.Q16[index(z,y,x)]);
//		//								D3.Q17[index(z,y,x)]=D3_hlp.Q17[index(z,y,x)]+omega*(n_equ[17] - D3_hlp.Q17[index(z,y,x)]);
//		//								D3.Q18[index(z,y,x)]=D3_hlp.Q18[index(z,y,x)]+omega*(n_equ[18] - D3_hlp.Q18[index(z,y,x)]);
//
//		if (x == lx-2) {
//			u_previous_spatial_boundary[index2D(z,y)] = u_x;
//			v_previous_spatial_boundary[index2D(z,y)] = u_y;
//			w_previous_spatial_boundary[index2D(z,y)] = u_z;
//
//			u_current[index2D(z,y)] = u_x;
//			v_current[index2D(z,y)] = u_y;
//			w_current[index2D(z,y)] = u_z;
//
//		}//if (x == lx-2)
//	}//if memory!
//	//			}//for (x = 0 ; x< lx; ++x)
//	//		}//for (y = 0 ; y< ly ; ++y)
//	//	}//for (z = 0 ; z< lz ; ++z)
//
//}

__global__
void relaxation_kernel_v4(const int lx, const int ly, const int lz, const FLOATING reynolds, FLOATING nu, FLOATING r_small,
		FLOATING t_0, FLOATING t_1, FLOATING t_2, FLOATING c_squ, FLOATING omega, FLOATING one_minus_omega,
		FLOATING reciprocal_c_squ,
		const FLOATING *hlp_Q0, const FLOATING *hlp_Q1, const FLOATING *hlp_Q2, const FLOATING *hlp_Q3,
		const FLOATING *hlp_Q4, const FLOATING *hlp_Q5, const FLOATING *hlp_Q6, const FLOATING *hlp_Q7,
		const FLOATING *hlp_Q8, const FLOATING *hlp_Q9, const FLOATING *hlp_Q10, const FLOATING *hlp_Q11,
		const FLOATING *hlp_Q12, const FLOATING *hlp_Q13, const FLOATING *hlp_Q14, const FLOATING *hlp_Q15,
		const FLOATING *hlp_Q16, const FLOATING *hlp_Q17, const FLOATING *hlp_Q18,
		FLOATING *Q0, FLOATING *Q1, FLOATING *Q2, FLOATING *Q3,
		FLOATING *Q4, FLOATING *Q5, FLOATING *Q6, FLOATING *Q7,
		FLOATING *Q8, FLOATING *Q9, FLOATING *Q10, FLOATING *Q11,
		FLOATING *Q12, FLOATING *Q13, FLOATING *Q14, FLOATING *Q15,
		FLOATING *Q16, FLOATING *Q17, FLOATING *Q18,
		const int *obstacles,
		FLOATING *u_previous_spatial_boundary, FLOATING *v_previous_spatial_boundary, FLOATING *w_previous_spatial_boundary,
		FLOATING *u_current, FLOATING *v_current, FLOATING *w_current, FLOATING *u_current_temp){





	/*One-step density relaxation process

				c.......density relaxation: a single time relaxation with relaxation
				c       parameter omega is applied here. This step is only "local",
				c       nothing is propagated through the lattice.
				c*/

	int  x,y,z;
	FLOATING  u_x=0.0, u_y=0.0, u_z=0.0, u_squ=0.0, rho=0.0, reciprocal_rho=0.0;




	__shared__ FLOATING u_n[19];

	__shared__ FLOATING buff[19];

	//	extern __shared__ FLOATING shared_buffer[];

	const FLOATING omega_x_t_0=omega*t_0, omega_x_t_1=omega*t_1, omega_x_t_2=omega*t_2;
	FLOATING omega_x_rho_x_t_0, omega_x_rho_x_t_1, omega_x_rho_x_t_2;
	FLOATING temp_factor;
	FLOATING u_n__over__c_squ[19];
	FLOATING u_n__over__c_squ__squared_and_halved[19];

	//....square speed of sound
	/*      compute the out let velocity with a convevtive boundary condition
					c.....loop over all nodes
					c.....attention: actual densities are stored after the propagation
					c                step in the help-array n_hlp !*/



	const int tid=blockIdx.x*blockDim.x+threadIdx.x;
	int rest;
	int end_of_memory=lz*ly*(lx);

	z=(int) (tid/(ly*lx));
	rest=tid-z;
	y=(int)(rest/lx);
	x=rest-y;

	if (tid<end_of_memory){
		/*c.........only free nodes are considered here
					!if (.not. obstacles[z][y][x]) then
					c...........integral local density
					c...........initialize variable ro*/
		//memory optimised implementation
		buff[0]=hlp_Q0[tid];
		buff[1]=hlp_Q1[tid];
		buff[2]=hlp_Q2[tid];
		buff[3]=hlp_Q3[tid];
		buff[4]=hlp_Q4[tid];
		buff[5]=hlp_Q5[tid];
		buff[6]=hlp_Q6[tid];
		buff[7]=hlp_Q7[tid];
		buff[8]=hlp_Q8[tid];
		buff[9]=hlp_Q9[tid];
		buff[10]=hlp_Q10[tid];
		buff[11]=hlp_Q11[tid];
		buff[12]=hlp_Q12[tid];
		buff[13]=hlp_Q13[tid];
		buff[14]=hlp_Q14[tid];
		buff[15]=hlp_Q15[tid];
		buff[16]=hlp_Q16[tid];
		buff[17]=hlp_Q17[tid];
		buff[18]=hlp_Q18[tid];

		rho=0.0;
#pragma unroll
		for(int k=0; k<DENSITIES; ++k)
			rho+=buff[k];



		reciprocal_rho=1.0/rho;




		u_x =  (1-obstacles[tid])*reciprocal_rho*(buff[1] + buff[7] + buff[10] +buff[11] + buff[14]-
				(buff[3] + buff[8] + buff[9] +buff[12] + buff[13]));


		u_y =   (1-obstacles[tid])*reciprocal_rho*(buff[2]+buff[8]+buff[7]+buff[16] + buff[15] -
				(buff[4] + buff[9] + buff[10] +buff[17] + buff[18]));


		u_z =   (1-obstacles[tid])*reciprocal_rho*(buff[5]+buff[13]+buff[14]+buff[15]+buff[18]-
				(buff[6]+buff[12]+buff[11]+buff[16]+buff[17]));






		u_squ = (FLOATING)  u_x*u_x + u_y*u_y + u_z*u_z;
		temp_factor= 0.5*(2.0* c_squ - u_squ)/c_squ;
		//u_squ = (FLOATING)  pow(u_x,2) + pow(u_y,2) + pow(u_z,2);


		/*...........n- velocity compnents (n = lattice node connection vectors)
					c...........this is only necessary for clearence, and only 3 speeds would
					c...........be necessary*/


		//WARNING!!!! o pinakas autos exei tropopoihmena indices!!!!
		u_n[0]= 0.0; //SHOULD NEVER USED!
		u_n[1] =   u_x;
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

#pragma unroll
		for(int i=0; i<DENSITIES; ++i){
			u_n__over__c_squ[i]=reciprocal_c_squ*u_n[i];
			u_n__over__c_squ__squared_and_halved[i]=0.5*u_n__over__c_squ[i]*u_n__over__c_squ[i];
		}

		/*c...........equilibrium densities
					c...........this can be rewritten to improve computational performance
					c...........considerabely !
					c
					c...........zero velocity density
					c*/
		//memory optimised implementation! WARNING!!! different from the original case!

		//two_x_c_squ_sqared=2.0*c_squ*c_squ;
		omega_x_rho_x_t_0=omega_x_t_0*rho;
		omega_x_rho_x_t_1=omega_x_t_1*rho;
		omega_x_rho_x_t_2=omega_x_t_2*rho;


		//				//...........relaxation step


		//omega_x_rho_x_t_0*(1.0 - 0.5*u_squ/c_squ);
		Q0[tid]=buff[0]*one_minus_omega+omega_x_rho_x_t_0*(u_n__over__c_squ__squared_and_halved[0]+u_n__over__c_squ[0]+temp_factor);




		Q1[tid]=buff[1]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[1]+u_n__over__c_squ[1]+temp_factor);
		Q2[tid]=buff[2]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[2]+u_n__over__c_squ[2]+temp_factor);
		Q3[tid]=buff[3]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[3]+u_n__over__c_squ[3]+temp_factor);
		Q4[tid]=buff[4]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[4]+u_n__over__c_squ[4]+temp_factor);
		Q5[tid]=buff[5]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[5]+u_n__over__c_squ[5]+temp_factor);
		Q6[tid]=buff[6]*one_minus_omega+omega_x_rho_x_t_1*(u_n__over__c_squ__squared_and_halved[6]+u_n__over__c_squ[6]+temp_factor);



		Q7[tid]= buff[ 7]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[ 7]+u_n__over__c_squ[ 7]+temp_factor);
		Q8[tid]= buff[ 8]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[ 8]+u_n__over__c_squ[ 8]+temp_factor);
		Q9[tid]= buff[ 9]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[ 9]+u_n__over__c_squ[ 9]+temp_factor);
		Q10[tid]=buff[10]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[10]+u_n__over__c_squ[10]+temp_factor);
		Q11[tid]=buff[11]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[11]+u_n__over__c_squ[11]+temp_factor);
		Q12[tid]=buff[12]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[12]+u_n__over__c_squ[12]+temp_factor);
		Q13[tid]=buff[13]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[13]+u_n__over__c_squ[13]+temp_factor);
		Q14[tid]=buff[14]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[14]+u_n__over__c_squ[14]+temp_factor);
		Q15[tid]=buff[15]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[15]+u_n__over__c_squ[15]+temp_factor);
		Q16[tid]=buff[16]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[16]+u_n__over__c_squ[16]+temp_factor);
		Q17[tid]=buff[17]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[17]+u_n__over__c_squ[17]+temp_factor);
		Q18[tid]=buff[18]*one_minus_omega+omega_x_rho_x_t_2*(u_n__over__c_squ__squared_and_halved[18]+u_n__over__c_squ[18]+temp_factor);


		//todo : improve the following if by having a kernel to collect all the necessary data
		if (x == lx-2 and obstacles[index(z,y,(lx-1))]==0) {
			u_previous_spatial_boundary[index2D(z,y)] = u_x;
			v_previous_spatial_boundary[index2D(z,y)] = u_y;
			w_previous_spatial_boundary[index2D(z,y)] = u_z;

			u_current[index2D(z,y)] = u_x;
			v_current[index2D(z,y)] = u_y;
			w_current[index2D(z,y)] = u_z;

			u_current_temp[index2D(z,y)] = u_x;


		}//if (x == lx-2)
	}//if memory!


}


void LBM::cuda_relaxation(){

	if(data_location==CPU)
		copy_data_from_host_to_device();

	dim3 threads_type2(threads_for_streaming_collision_and_relaxation,1,1);
	dim3 grid_type2(blocks_for_streaming_collision_and_relaxation,1,1);

	relaxation_kernel_v4<<<grid_type2, threads_type2>>>( lx,  ly,  lz,  reynolds,  nu,  r_small,
			t_0,  t_1,  t_2,  c_squ,  omega,  one_minus_omega,
			reciprocal_c_squ,
			D3_hlp_d.Q0, D3_hlp_d.Q1, D3_hlp_d.Q2, D3_hlp_d.Q3,
			D3_hlp_d.Q4, D3_hlp_d.Q5, D3_hlp_d.Q6, D3_hlp_d.Q7,
			D3_hlp_d.Q8, D3_hlp_d.Q9, D3_hlp_d.Q10, D3_hlp_d.Q11,
			D3_hlp_d.Q12, D3_hlp_d.Q13, D3_hlp_d.Q14, D3_hlp_d.Q15,
			D3_hlp_d.Q16, D3_hlp_d.Q17, D3_hlp_d.Q18,
			D3_d.Q0, D3_d.Q1, D3_d.Q2, D3_d.Q3,
			D3_d.Q4, D3_d.Q5, D3_d.Q6, D3_d.Q7,
			D3_d.Q8, D3_d.Q9, D3_d.Q10, D3_d.Q11,
			D3_d.Q12, D3_d.Q13, D3_d.Q14, D3_d.Q15,
			D3_d.Q16, D3_d.Q17, D3_d.Q18,
			obstacles_d,
			u_previous_spatial_boundary_d, v_previous_spatial_boundary_d, w_previous_spatial_boundary_d,
			u_current_d, v_current_d, w_current_d, u_current_temp_d);




	cudaDeviceSynchronize();
}
