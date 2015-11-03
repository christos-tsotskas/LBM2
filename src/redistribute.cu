#include "global_defines.cuh"

__global__
void redistribute_kernel(const int lx, const int ly, const int lz,
		FLOATING reynolds, FLOATING nu, FLOATING r_small,
		FLOATING t_0, FLOATING t_1, FLOATING t_2, FLOATING c_squ,

		FLOATING *Q0, FLOATING *Q1, FLOATING *Q2, FLOATING *Q3,
		FLOATING *Q4, FLOATING *Q5, FLOATING *Q6, FLOATING *Q7,
		FLOATING *Q8, FLOATING *Q9, FLOATING *Q10, FLOATING *Q11,
		FLOATING *Q12, FLOATING *Q13, FLOATING *Q14, FLOATING *Q15,
		FLOATING *Q16, FLOATING *Q17, FLOATING *Q18){

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

	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	//int end_of_memory=lz*ly*(lx);
	int z=(int) (tid/(ly*lx));
	int	y=(tid-z*(ly*lx))/lx;
	int	x=tid-z*(ly*lx)-y*lx;

	int R_big;
	int baffle_position=59;

	FLOATING    mass_flow;

	//.....local variables

	int   yc, zc, yr, zr;


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

	//	int end_of_memory=lz*ly*(baffle_position+1);





	zr=z-zc;
	yr=y-yc;

	if(x< (baffle_position+1) and yr*yr+zr*zr < r_small*r_small and tid<lx*ly*lz){
		// id = z*+y+x


		//		for( z = 0; z< lz ; ++z){
		//			for( y = 0; y< ly ; ++y){
		//				for( x = 0; x< baffle_position+1 ; ++x){



		//		rho=0.0;
		rho=Q0[index(z,y,x)]+Q1[index(z,y,x)]+Q2[index(z,y,x)]+Q3[index(z,y,x)]+
				Q4[index(z,y,x)]+Q5[index(z,y,x)]+Q6[index(z,y,x)]+Q7[index(z,y,x)]+
				Q8[index(z,y,x)]+Q9[index(z,y,x)]+Q10[index(z,y,x)]+Q11[index(z,y,x)]+
				Q12[index(z,y,x)]+Q13[index(z,y,x)]+Q14[index(z,y,x)]+Q15[index(z,y,x)]+
				Q16[index(z,y,x)]+Q17[index(z,y,x)]+Q18[index(z,y,x)];


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


		Q0[index(z,y,x)]=(FLOATING) (t_0  * rho *(1.0  - u_squ / (2.0  * c_squ)));

		//...........axis speeds (factor: t_1)

		Q1[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[1]/c_squ ) +  0.5* ( (u_n[1]*u_n[1])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));;
		Q2[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[2]/c_squ ) +  0.5* ( (u_n[2]*u_n[2])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));;
		Q3[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[3]/c_squ ) +  0.5* ( (u_n[3]*u_n[3])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));;
		Q4[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[4]/c_squ ) +  0.5* ( (u_n[4]*u_n[4])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));;
		Q5[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[5]/c_squ ) +  0.5* ( (u_n[5]*u_n[5])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));;
		Q6[index(z,y,x)]=(FLOATING) (t_1 * rho *	 (1.0+ ( u_n[6]/c_squ ) +  0.5* ( (u_n[6]*u_n[6])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));;

		//...........diagonal speeds (factor: t_2)
		Q7[index(z,y,x)]=(FLOATING) (t_2 * rho *	 (1.0+ ( u_n[7]/c_squ ) +  0.5* ( (u_n[7]*u_n[7])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));;
		Q8[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[8]/c_squ ) +  0.5* ( (u_n[8]*u_n[8])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
		Q9[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[9]/c_squ ) +  0.5* ( (u_n[9]*u_n[9])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
		Q10[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[10]/c_squ ) +  0.5* ( (u_n[10]*u_n[10])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
		Q11[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[11]/c_squ ) +  0.5* ( (u_n[11]*u_n[11])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
		Q12[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[12]/c_squ ) +  0.5* ( (u_n[12]*u_n[12])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
		Q13[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[13]/c_squ ) +  0.5* ( (u_n[13]*u_n[13])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
		Q14[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[14]/c_squ ) +  0.5* ( (u_n[14]*u_n[14])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
		Q15[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[15]/c_squ ) +  0.5* ( (u_n[15]*u_n[15])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
		Q16[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[16]/c_squ ) +  0.5* ( (u_n[16]*u_n[16])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
		Q17[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[17]/c_squ ) +  0.5* ( (u_n[17]*u_n[17])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));
		Q18[index(z,y,x)]=(FLOATING) (t_2  * rho * (1.0+ ( u_n[18]/c_squ ) +  0.5* ( (u_n[18]*u_n[18])/(c_squ*c_squ) ) - 0.5 * ( u_squ/c_squ ) ));


		//				}
		//			}
		//		}
	}


}


void LBM::cuda_redistribute(){

	if(data_location==CPU)
		copy_data_from_host_to_device();

//	int lattice_nodes=lz*ly*lx;
//
//	int n_of_threads=128;
//	int n_of_blocks=ceil((lattice_nodes*1.0)/n_of_threads);
//	dim3 threads_type2(n_of_threads,1,1);
//	dim3 grid_type2(n_of_blocks,1,1);
//#ifdef REPORT
//	cout << "redistribute kernel with:" << lattice_nodes << " lattice nodes" << endl;
//	cout << "\tthreads:" << n_of_threads << endl;
//	cout << "\tblocks:" << n_of_blocks << endl;
//#endif

	dim3 threads_type2(threads_for_streaming_collision_and_relaxation,1,1);
		dim3 grid_type2(blocks_for_streaming_collision_and_relaxation,1,1);


	redistribute_kernel<<<grid_type2, threads_type2>>>(lx, ly, lz, reynolds, nu, r_small,
			t_0, t_1, t_2, c_squ,
			D3_d.Q0, D3_d.Q1, D3_d.Q2, D3_d.Q3,
			D3_d.Q4, D3_d.Q5, D3_d.Q6, D3_d.Q7,
			D3_d.Q8, D3_d.Q9, D3_d.Q10, D3_d.Q11,
			D3_d.Q12, D3_d.Q13, D3_d.Q14, D3_d.Q15,
			D3_d.Q16, D3_d.Q17, D3_d.Q18);

	cudaDeviceSynchronize();

}
