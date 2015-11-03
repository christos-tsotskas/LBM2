#include "global_defines.cuh"
#include "kernels3.cuh"


void	LBM::convective_BC(){

	if(data_location==GPU)
		copy_data_from_device_to_host();

	int  y,z,num;

	FLOATING  u_x,u_y,u_z,u_n[19];//,n_equ[19];
	FLOATING	u_squ,rho,Uc;



	//.....first compute the mean outflow velocity, Uc
	Uc = 0.0;
	num = 0;
	for (z = 0 ; z< lz ; ++z){
		for (y = 0 ; y< ly ; ++y){
			if (obstacles[index(z,y,(lx-1))]==0) {
				Uc +=   u_current[index2D(z,y)];
				++num ;
			}
		}
	}
	Uc /= num; //!if (num>0) check not needed


	cout << " CCCCPU(convective_BC, U_C_avg) Uc:"  << Uc << endl;

	for (z = 0 ; z< lz ; ++z){
		for (y = 0 ; y< ly ; ++y){
			if (!obstacles[index(z,y,(lx-1))]) {

				//.....compute the new velocities (based on convective BC)
				//originally proposed by Djenidi
//				u_current[index2D(z,y)] = (u_previous_temporal_boundary[index2D(z,y)] + Uc*u_previous_spatial_boundary[index2D(z,y)])/(1.0+Uc);
//				u_previous_temporal_boundary[index2D(z,y)] = u_current[index2D(z,y)];
//				v_current[index2D(z,y)] = (v_previous_temporal_boundary[index2D(z,y)] + Uc*v_previous_spatial_boundary[index2D(z,y)])/(1.0+Uc);
//				v_previous_temporal_boundary[index2D(z,y)] = v_current[index2D(z,y)];
//				w_current[index2D(z,y)] = (w_previous_temporal_boundary[index2D(z,y)] + Uc*w_previous_spatial_boundary[index2D(z,y)])/(1.0+Uc);
//				w_previous_temporal_boundary[index2D(z,y)] = w_current[index2D(z,y)];

				////.....compute the new velocities (based on NON convective BC)
				//suggested by Timos on 04102015
				u_current[index2D(z,y)] = Uc;
				v_current[index2D(z,y)] = 0;
				w_current[index2D(z,y)] = 0;


				rho=0.0;
				rho+=D3_hlp.Q0[index(z,y,lx-1)]+D3_hlp.Q1[index(z,y,lx-1)]+D3_hlp.Q2[index(z,y,lx-1)]+D3_hlp.Q3[index(z,y,lx-1)];
				rho+=D3_hlp.Q4[index(z,y,lx-1)]+D3_hlp.Q5[index(z,y,lx-1)]+D3_hlp.Q6[index(z,y,lx-1)]+D3_hlp.Q7[index(z,y,lx-1)];
				rho+=D3_hlp.Q8[index(z,y,lx-1)]+D3_hlp.Q9[index(z,y,lx-1)]+D3_hlp.Q10[index(z,y,lx-1)]+D3_hlp.Q11[index(z,y,lx-1)];
				rho+=D3_hlp.Q12[index(z,y,lx-1)]+D3_hlp.Q13[index(z,y,lx-1)]+D3_hlp.Q14[index(z,y,lx-1)]+D3_hlp.Q15[index(z,y,lx-1)];
				rho+=D3_hlp.Q16[index(z,y,lx-1)]+D3_hlp.Q17[index(z,y,lx-1)]+D3_hlp.Q18[index(z,y,lx-1)];

				u_x = u_current[index2D(z,y)];
				u_y = v_current[index2D(z,y)];
				u_z = w_current[index2D(z,y)];

				//...........square velocity
				u_squ = u_x * u_x + u_y * u_y + u_z * u_z;
				/*
								c...........n- velocity compnents (n = lattice node connection vectors)
								c...........this is only necessary for clearence, and only 3 speeds would
								c...........be necessary
				 */
				u_n[0]= 0.0; //SHOULD NEVER USED!
				u_n[1] =   u_x; //u_xs
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
					c...........zero velocity density
					c*/
				//original part!
				//				n_equ[0] = (FLOATING) (t_0  * rho*(1.0 - u_squ / (2.0 * c_squ)));
				//
				//
				//				//...........axis speeds (factor: t_1)
				//				//TODO: NA GINEI SUGXWNEUSH SE ENA CASE ME FOR
				//				for (i = 1 ; i< 7 ; ++i)
				//					n_equ[i] = (FLOATING) (t_1 * rho*(1.0 + u_n[i] / c_squ
				//							+ ( u_n[i] * u_n[i]) / (2.0 * (c_squ * c_squ))
				//							- u_squ / (2.0 * c_squ)));
				//
				//				//...........diagonal speeds (factor: t_2)
				//				for (i = 7 ; i< 19 ; ++i)
				//					n_equ[i] =  (FLOATING) (t_2 * rho*(1.0 + u_n[i] / c_squ
				//							+ ( u_n[i] * u_n[i]) / (2.0 * (c_squ * c_squ))
				//							- u_squ / (2.0 * c_squ)));
				//
				//
				//
				//				D3.Q0[index(z,y,lx-1)]=(FLOATING) n_equ[0];
				//
				//				//...........axis speeds (factor: t_1)
				//
				//				D3.Q1[index(z,y,lx-1)]=(FLOATING) n_equ[1];
				//				D3.Q2[index(z,y,lx-1)]=(FLOATING) n_equ[2];
				//				D3.Q3[index(z,y,lx-1)]=(FLOATING) n_equ[3];
				//				D3.Q4[index(z,y,lx-1)]=(FLOATING) n_equ[4];
				//				D3.Q5[index(z,y,lx-1)]=(FLOATING) n_equ[5];
				//				D3.Q6[index(z,y,lx-1)]=(FLOATING) n_equ[6];
				//
				//				//...........diagonal speeds (factor: t_2)
				//				D3.Q7[index(z,y,lx-1)]=(FLOATING) n_equ[7];
				//				D3.Q8[index(z,y,lx-1)]=(FLOATING) n_equ[8];
				//				D3.Q9[index(z,y,lx-1)]=(FLOATING) n_equ[9];
				//				D3.Q10[index(z,y,lx-1)]=(FLOATING) n_equ[10];
				//				D3.Q11[index(z,y,lx-1)]=(FLOATING) n_equ[11];
				//				D3.Q12[index(z,y,lx-1)]=(FLOATING) n_equ[12];
				//				D3.Q13[index(z,y,lx-1)]=(FLOATING) n_equ[13];
				//				D3.Q14[index(z,y,lx-1)]=(FLOATING) n_equ[14];
				//				D3.Q15[index(z,y,lx-1)]=(FLOATING) n_equ[15];
				//				D3.Q16[index(z,y,lx-1)]=(FLOATING) n_equ[16];
				//				D3.Q17[index(z,y,lx-1)]=(FLOATING) n_equ[17];
				//				D3.Q18[index(z,y,lx-1)]=(FLOATING) n_equ[18];

				//optimised!
				D3.Q0[index(z,y,lx-1)]=(t_0  * rho*(1.0 - u_squ / (2.0 * c_squ)));
				//...........axis speeds (factor: t_1)
				D3.Q1[index(z,y,lx-1)]=(t_1 * rho*(1.0 + u_n[1] / c_squ	+ ( u_n[1] * u_n[1]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
				D3.Q2[index(z,y,lx-1)]=(t_1 * rho*(1.0 + u_n[2] / c_squ	+ ( u_n[2] * u_n[2]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
				D3.Q3[index(z,y,lx-1)]=(t_1 * rho*(1.0 + u_n[3] / c_squ	+ ( u_n[3] * u_n[3]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
				D3.Q4[index(z,y,lx-1)]=(t_1 * rho*(1.0 + u_n[4] / c_squ	+ ( u_n[4] * u_n[4]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
				D3.Q5[index(z,y,lx-1)]=(t_1 * rho*(1.0 + u_n[5] / c_squ	+ ( u_n[5] * u_n[5]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
				D3.Q6[index(z,y,lx-1)]=(t_1 * rho*(1.0 + u_n[6] / c_squ	+ ( u_n[6] * u_n[6]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));

				//...........diagonal speeds (factor: t_2)
				D3.Q7[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[7] / c_squ 	+ ( u_n[7] * u_n[7]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
				D3.Q8[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[8] / c_squ 	+ ( u_n[8] * u_n[8]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
				D3.Q9[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[9] / c_squ 	+ ( u_n[9] * u_n[9]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
				D3.Q10[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[10] / c_squ 	+ ( u_n[10] * u_n[10]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
				D3.Q11[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[11] / c_squ 	+ ( u_n[11] * u_n[11]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
				D3.Q12[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[12] / c_squ 	+ ( u_n[12] * u_n[12]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
				D3.Q13[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[13] / c_squ 	+ ( u_n[13] * u_n[13]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
				D3.Q14[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[14] / c_squ 	+ ( u_n[14] * u_n[14]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
				D3.Q15[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[15] / c_squ 	+ ( u_n[15] * u_n[15]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
				D3.Q16[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[16] / c_squ 	+ ( u_n[16] * u_n[16]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
				D3.Q17[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[17] / c_squ 	+ ( u_n[17] * u_n[17]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
				D3.Q18[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[18] / c_squ 	+ ( u_n[18] * u_n[18]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));

			}
		}
	}

#ifdef DEBUG
	cout << " #LBM convective_bc OK!" << endl;
#endif

}






__global__
void	new_convective_BC_kernel_v1(FLOATING Uc, const int free_lattices_at_U_direction,int end_of_memory, int lx, int ly, int lz, FLOATING reynolds, FLOATING nu, FLOATING r_small,
		FLOATING t_0, FLOATING t_1, FLOATING t_2, FLOATING c_squ, FLOATING omega, FLOATING one_minus_omega,
		FLOATING reciprocal_c_squ,
		FLOATING *hlp_Q0, FLOATING *hlp_Q1, FLOATING *hlp_Q2, FLOATING *hlp_Q3,
		FLOATING *hlp_Q4, FLOATING *hlp_Q5, FLOATING *hlp_Q6, FLOATING *hlp_Q7,
		FLOATING *hlp_Q8, FLOATING *hlp_Q9, FLOATING *hlp_Q10, FLOATING *hlp_Q11,
		FLOATING *hlp_Q12, FLOATING *hlp_Q13, FLOATING *hlp_Q14, FLOATING *hlp_Q15,
		FLOATING *hlp_Q16, FLOATING *hlp_Q17, FLOATING *hlp_Q18,
		FLOATING *Q0, FLOATING *Q1, FLOATING *Q2, FLOATING *Q3,
		FLOATING *Q4, FLOATING *Q5, FLOATING *Q6, FLOATING *Q7,
		FLOATING *Q8, FLOATING *Q9, FLOATING *Q10, FLOATING *Q11,
		FLOATING *Q12, FLOATING *Q13, FLOATING *Q14, FLOATING *Q15,
		FLOATING *Q16, FLOATING *Q17, FLOATING *Q18,
		int *obstacles,
		FLOATING *u_previous_spatial_boundary, FLOATING *v_previous_spatial_boundary, FLOATING *w_previous_spatial_boundary,
		FLOATING *u_current, FLOATING *v_current, FLOATING *w_current,
		FLOATING *u_previous_temporal_boundary, FLOATING *v_previous_temporal_boundary, FLOATING *w_previous_temporal_boundary){

	if(blockIdx.x*blockDim.x+threadIdx.x==0){
		int  y,z  ;

		FLOATING  u_x,u_y,u_z,u_n[19];//,n_equ[19];
		FLOATING	u_squ,rho ;



		//.....first compute the mean outflow velocity, Uc
		//		Uc = 0.0;
		//
		//		for (z = 0 ; z< lz ; ++z){
		//			for (y = 0 ; y< ly ; ++y){
		////				if (obstacles[index(z,y,(lx-1))]==0) {
		//					Uc +=   u_current[index2D(z,y)];
		//
		////				}
		//			}
		//		}
		//		Uc /= free_lattices_at_U_direction; //!if (num>0) check not needed
		printf( "u-wise free lattices:%d \n",free_lattices_at_U_direction);

		printf( "within convective BC, Uc: %f\n", Uc);



		for (z = 0 ; z< lz ; ++z){
			for (y = 0 ; y< ly ; ++y){
				if (!obstacles[index(z,y,(lx-1))]) {

					//.....compute the new velocities (based on convective BC)
					u_current[index2D(z,y)] = (u_previous_temporal_boundary[index2D(z,y)] + Uc*u_previous_spatial_boundary[index2D(z,y)])/(1.0+Uc);
					u_previous_temporal_boundary[index2D(z,y)] = u_current[index2D(z,y)];
					v_current[index2D(z,y)] = (v_previous_temporal_boundary[index2D(z,y)] + Uc*v_previous_spatial_boundary[index2D(z,y)])/(1.0+Uc);
					v_previous_temporal_boundary[index2D(z,y)] = v_current[index2D(z,y)];
					w_current[index2D(z,y)] = (w_previous_temporal_boundary[index2D(z,y)] + Uc*w_previous_spatial_boundary[index2D(z,y)])/(1.0+Uc);
					w_previous_temporal_boundary[index2D(z,y)] = w_current[index2D(z,y)];

					rho=hlp_Q0[index(z,y,lx-1)]+hlp_Q1[index(z,y,lx-1)]+hlp_Q2[index(z,y,lx-1)]+hlp_Q3[index(z,y,lx-1)]+
							hlp_Q4[index(z,y,lx-1)]+hlp_Q5[index(z,y,lx-1)]+hlp_Q6[index(z,y,lx-1)]+hlp_Q7[index(z,y,lx-1)]+
							hlp_Q8[index(z,y,lx-1)]+hlp_Q9[index(z,y,lx-1)]+hlp_Q10[index(z,y,lx-1)]+hlp_Q11[index(z,y,lx-1)]+
							hlp_Q12[index(z,y,lx-1)]+hlp_Q13[index(z,y,lx-1)]+hlp_Q14[index(z,y,lx-1)]+hlp_Q15[index(z,y,lx-1)]+
							hlp_Q16[index(z,y,lx-1)]+hlp_Q17[index(z,y,lx-1)]+hlp_Q18[index(z,y,lx-1)];


					u_x = u_current[index2D(z,y)];
					u_y = v_current[index2D(z,y)];
					u_z = w_current[index2D(z,y)];

					//...........square velocity
					u_squ = u_x * u_x + u_y * u_y + u_z * u_z;
					/*
									c...........n- velocity compnents (n = lattice node connection vectors)
									c...........this is only necessary for clearence, and only 3 speeds would
									c...........be necessary
					 */
					u_n[0]= 0.0; //SHOULD NEVER USED!
					u_n[1] =   u_x; //u_xs
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
						c...........zero velocity density
						c*/


					//optimised!
					Q0[index(z,y,lx-1)]=(t_0  * rho*(1.0 - u_squ / (2.0 * c_squ)));
					//...........axis speeds (factor: t_1)
					Q1[index(z,y,lx-1)]=(t_1 * rho*(1.0 + u_n[1] / c_squ	+ ( u_n[1] * u_n[1]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
					Q2[index(z,y,lx-1)]=(t_1 * rho*(1.0 + u_n[2] / c_squ	+ ( u_n[2] * u_n[2]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
					Q3[index(z,y,lx-1)]=(t_1 * rho*(1.0 + u_n[3] / c_squ	+ ( u_n[3] * u_n[3]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
					Q4[index(z,y,lx-1)]=(t_1 * rho*(1.0 + u_n[4] / c_squ	+ ( u_n[4] * u_n[4]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
					Q5[index(z,y,lx-1)]=(t_1 * rho*(1.0 + u_n[5] / c_squ	+ ( u_n[5] * u_n[5]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
					Q6[index(z,y,lx-1)]=(t_1 * rho*(1.0 + u_n[6] / c_squ	+ ( u_n[6] * u_n[6]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));

					//...........diagonal speeds (factor: t_2)
					Q7[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[7] / c_squ 	+ ( u_n[7] * u_n[7]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
					Q8[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[8] / c_squ 	+ ( u_n[8] * u_n[8]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
					Q9[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[9] / c_squ 	+ ( u_n[9] * u_n[9]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
					Q10[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[10] / c_squ 	+ ( u_n[10] * u_n[10]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
					Q11[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[11] / c_squ 	+ ( u_n[11] * u_n[11]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
					Q12[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[12] / c_squ 	+ ( u_n[12] * u_n[12]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
					Q13[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[13] / c_squ 	+ ( u_n[13] * u_n[13]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
					Q14[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[14] / c_squ 	+ ( u_n[14] * u_n[14]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
					Q15[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[15] / c_squ 	+ ( u_n[15] * u_n[15]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
					Q16[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[16] / c_squ 	+ ( u_n[16] * u_n[16]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
					Q17[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[17] / c_squ 	+ ( u_n[17] * u_n[17]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
					Q18[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[18] / c_squ 	+ ( u_n[18] * u_n[18]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));

				}
			}
		}

	}
}

void compute_mean_outflow_velocity(const FLOATING *u_current_d, const int *obstacles_d){


}

__global__
void	new_convective_BC_kernel_v2(const FLOATING Uc, int end_of_memory, int lx, int ly, int lz, FLOATING reynolds, FLOATING nu, FLOATING r_small,
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
		FLOATING *u_current, FLOATING *v_current, FLOATING *w_current,
		FLOATING *u_previous_temporal_boundary, FLOATING *v_previous_temporal_boundary, FLOATING *w_previous_temporal_boundary){

	int tid=blockIdx.x*blockDim.x+threadIdx.x;

	int z=(int) tid / ly;
	int y=(int) tid % ly;

	if(tid<ly*lz){


		FLOATING  u_x,u_y,u_z,u_n[19];//,n_equ[19];
		FLOATING	u_squ,rho ;






		//		for (z = 0 ; z< lz ; ++z){
		//			for (y = 0 ; y< ly ; ++y){
		//				if (!obstacles[index(z,y,(lx-1))]) {

		//.....compute the new velocities (based on convective BC)
		u_current[index2D(z,y)] = (u_previous_temporal_boundary[index2D(z,y)] + Uc*u_previous_spatial_boundary[index2D(z,y)])/(1.0+Uc);
		u_previous_temporal_boundary[index2D(z,y)] = u_current[index2D(z,y)];
		v_current[index2D(z,y)] = (v_previous_temporal_boundary[index2D(z,y)] + Uc*v_previous_spatial_boundary[index2D(z,y)])/(1.0+Uc);
		v_previous_temporal_boundary[index2D(z,y)] = v_current[index2D(z,y)];
		w_current[index2D(z,y)] = (w_previous_temporal_boundary[index2D(z,y)] + Uc*w_previous_spatial_boundary[index2D(z,y)])/(1.0+Uc);
		w_previous_temporal_boundary[index2D(z,y)] = w_current[index2D(z,y)];

		rho=hlp_Q0[index(z,y,lx-1)]+hlp_Q1[index(z,y,lx-1)]+hlp_Q2[index(z,y,lx-1)]+hlp_Q3[index(z,y,lx-1)]+
				hlp_Q4[index(z,y,lx-1)]+hlp_Q5[index(z,y,lx-1)]+hlp_Q6[index(z,y,lx-1)]+hlp_Q7[index(z,y,lx-1)]+
				hlp_Q8[index(z,y,lx-1)]+hlp_Q9[index(z,y,lx-1)]+hlp_Q10[index(z,y,lx-1)]+hlp_Q11[index(z,y,lx-1)]+
				hlp_Q12[index(z,y,lx-1)]+hlp_Q13[index(z,y,lx-1)]+hlp_Q14[index(z,y,lx-1)]+hlp_Q15[index(z,y,lx-1)]+
				hlp_Q16[index(z,y,lx-1)]+hlp_Q17[index(z,y,lx-1)]+hlp_Q18[index(z,y,lx-1)];


		u_x = u_current[index2D(z,y)];
		u_y = v_current[index2D(z,y)];
		u_z = w_current[index2D(z,y)];

		//...........square velocity
		u_squ = u_x * u_x + u_y * u_y + u_z * u_z;
		/*
									c...........n- velocity compnents (n = lattice node connection vectors)
									c...........this is only necessary for clearence, and only 3 speeds would
									c...........be necessary
		 */
		u_n[0]= 0.0; //SHOULD NEVER USED!
		u_n[1] =   u_x; //u_xs
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
						c...........zero velocity density
						c*/


		//optimised!
		Q0[index(z,y,lx-1)]=(t_0  * rho*(1.0 - u_squ / (2.0 * c_squ)));
		//...........axis speeds (factor: t_1)
		Q1[index(z,y,lx-1)]=(t_1 * rho*(1.0 + u_n[1] / c_squ	+ ( u_n[1] * u_n[1]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
		Q2[index(z,y,lx-1)]=(t_1 * rho*(1.0 + u_n[2] / c_squ	+ ( u_n[2] * u_n[2]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
		Q3[index(z,y,lx-1)]=(t_1 * rho*(1.0 + u_n[3] / c_squ	+ ( u_n[3] * u_n[3]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
		Q4[index(z,y,lx-1)]=(t_1 * rho*(1.0 + u_n[4] / c_squ	+ ( u_n[4] * u_n[4]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
		Q5[index(z,y,lx-1)]=(t_1 * rho*(1.0 + u_n[5] / c_squ	+ ( u_n[5] * u_n[5]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
		Q6[index(z,y,lx-1)]=(t_1 * rho*(1.0 + u_n[6] / c_squ	+ ( u_n[6] * u_n[6]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));

		//...........diagonal speeds (factor: t_2)
		Q7[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[7] / c_squ 	+ ( u_n[7] * u_n[7]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
		Q8[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[8] / c_squ 	+ ( u_n[8] * u_n[8]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
		Q9[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[9] / c_squ 	+ ( u_n[9] * u_n[9]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
		Q10[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[10] / c_squ 	+ ( u_n[10] * u_n[10]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
		Q11[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[11] / c_squ 	+ ( u_n[11] * u_n[11]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
		Q12[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[12] / c_squ 	+ ( u_n[12] * u_n[12]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
		Q13[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[13] / c_squ 	+ ( u_n[13] * u_n[13]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
		Q14[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[14] / c_squ 	+ ( u_n[14] * u_n[14]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
		Q15[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[15] / c_squ 	+ ( u_n[15] * u_n[15]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
		Q16[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[16] / c_squ 	+ ( u_n[16] * u_n[16]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
		Q17[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[17] / c_squ 	+ ( u_n[17] * u_n[17]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));
		Q18[index(z,y,lx-1)]=(t_2 * rho*(1.0 + u_n[18] / c_squ 	+ ( u_n[18] * u_n[18]) / (2.0 * (c_squ * c_squ)) - u_squ / (2.0 * c_squ)));

		//} //if (!obstacles[index(z,y,(lx-1))])
		//			} // for y
		//		} // for z

	}
}


void LBM::cuda_convective_BC(){
	if(data_location==CPU)
		copy_data_from_host_to_device();


	dim3 threads_type2(threads_per_kernel,1,1);
	dim3 grid_type2(convective_boundary_conditions_blocks,1,1);



	FLOATING temp_Uc=reduce_sum(u_current_temp_d, lz*ly)/no_obstacle_lattices_at_penultimate_x_slice;



	//simplest working implementation
	//		new_convective_BC_kernel_v1<<<grid_type1, threads_type1>>>(temp_Uc,no_obstacle_lattices_at_penultimate_x_slice, ly*lz,  lx,  ly,  lz,
	//				reynolds,  nu,  r_small, t_0,  t_1,  t_2,
	//				c_squ,  omega,  one_minus_omega, reciprocal_c_squ,
	//				D3_hlp_d.Q0, D3_hlp_d.Q1, D3_hlp_d.Q2, D3_hlp_d.Q3,
	//				D3_hlp_d.Q4, D3_hlp_d.Q5, D3_hlp_d.Q6, D3_hlp_d.Q7,
	//				D3_hlp_d.Q8, D3_hlp_d.Q9, D3_hlp_d.Q10, D3_hlp_d.Q11,
	//				D3_hlp_d.Q12, D3_hlp_d.Q13, D3_hlp_d.Q14, D3_hlp_d.Q15,
	//				D3_hlp_d.Q16, D3_hlp_d.Q17, D3_hlp_d.Q18,
	//				D3_d.Q0, D3_d.Q1, D3_d.Q2, D3_d.Q3,
	//				D3_d.Q4, D3_d.Q5, D3_d.Q6, D3_d.Q7,
	//				D3_d.Q8, D3_d.Q9, D3_d.Q10, D3_d.Q11,
	//				D3_d.Q12, D3_d.Q13, D3_d.Q14, D3_d.Q15,
	//				D3_d.Q16, D3_d.Q17, D3_d.Q18,
	//				obstacles_d,
	//				u_previous_spatial_boundary_d,  v_previous_spatial_boundary_d,  w_previous_spatial_boundary_d,
	//				u_current_d,  v_current_d,  w_current_d,
	//				u_previous_temporal_boundary_d,  v_previous_temporal_boundary_d,  w_previous_temporal_boundary_d);

	//average of u_current_d with no obstacles
	//
	//
	new_convective_BC_kernel_v2<<<grid_type2, threads_type2>>>(temp_Uc, ly*lz,  lx,  ly,  lz,
			reynolds,  nu,  r_small, t_0,  t_1,  t_2,
			c_squ,  omega,  one_minus_omega, reciprocal_c_squ,
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
			u_previous_spatial_boundary_d,  v_previous_spatial_boundary_d,  w_previous_spatial_boundary_d,
			u_current_d,  v_current_d,  w_current_d,
			u_previous_temporal_boundary_d,  v_previous_temporal_boundary_d,  w_previous_temporal_boundary_d);


	cudaDeviceSynchronize();

}
