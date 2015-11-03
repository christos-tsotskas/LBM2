#include "global_defines.cuh"


void LBM::streaming(){



	/*Propagate fluid densities to their next neighbour nodes */
	/*c
	c.......density propagation: all fluid densities are propagated from
	c       non-occupied nodes along the lattice connection lines
	c       to their next neighbours.
	c*/

	if(data_location==GPU)
		copy_data_from_device_to_host();


	int  x,y,z;
	int x_e/*east*/,x_w/*west*/;
	int y_n/*north*/,y_s/*south*/;
	int z_l/*left*/,z_r/*right*/;





	//todo rwta ton Dirk pws to kanei auto to vhma
	for (z = 0 ; z< lz ; ++z){
		z_l =  (z+1)%lz; //1 8esh meta to trexon x
		z_r =   (z+lz-1) %(  lz) ;

		for (y = 0; y< ly; ++y){
			y_n =  (y+1)%ly; //1 8esh meta to trexon y
			y_s =   (y+ly-1) %(  ly) ;

			for (x = 0; x< lx; ++x){
				x_e =  (x+1)%lx; //1 8esh meta to trexon z
				x_w =   (x+lx-1) %(  lx) ;

				//regular streaming process

				/*
					.........density propagation

					.........zero: just copy
				 */
				D3_hlp.Q0[index(z,y,x)]=D3.Q0[index(z,y,x)];
				//node_hlp_grid[index4D(z,y,x,0)]=node_grid[index4D(z,y,x,0)];

				//.........in the x,y and z directions
				D3_hlp.Q1[index(z  ,y  ,x_e)] = D3.Q1[index(z,y,x)];
				D3_hlp.Q2[index(z  ,y_n,x  )] = D3.Q2[index(z,y,x)];
				D3_hlp.Q3[index(z  ,y  ,x_w)] = D3.Q3[index(z,y,x)];
				D3_hlp.Q4[index(z  ,y_s,x  )] = D3.Q4[index(z,y,x)];
				D3_hlp.Q5[index(z_l,y  ,x  )] = D3.Q5[index(z,y,x)];
				D3_hlp.Q6[index(z_r,y  ,x  )] = D3.Q6[index(z,y,x)];


				//......... in the x,y diagonals


				D3_hlp.Q7[ index(z  ,y_n,x_e)] = D3.Q7[ index(z,y,x)];
				D3_hlp.Q8[ index(z  ,y_n,x_w)] = D3.Q8[ index(z,y,x)];
				D3_hlp.Q9[ index(z  ,y_s,x_w)] = D3.Q9[ index(z,y,x)];
				D3_hlp.Q10[index(z  ,y_s,x_e)] = D3.Q10[index(z,y,x)];

				//......... in the x,z diagonals

				D3_hlp.Q11[index(z_r,y  ,x_e)] = D3.Q11[index(z,y,x)];
				D3_hlp.Q12[index(z_r,y  ,x_w)] = D3.Q12[index(z,y,x)];
				D3_hlp.Q13[index(z_l,y  ,x_w)] = D3.Q13[index(z,y,x)];
				D3_hlp.Q14[index(z_l,y  ,x_e)] = D3.Q14[index(z,y,x)];

				//......... in the y,z diagonals

				D3_hlp.Q15[index(z_l,y_n,x  )] = D3.Q15[index(z,y,x)];
				D3_hlp.Q16[index(z_r,y_n,x  )] = D3.Q16[index(z,y,x)];
				D3_hlp.Q17[index(z_r,y_s,x  )] = D3.Q17[index(z,y,x)];
				D3_hlp.Q18[index(z_l,y_s,x  )] = D3.Q18[index(z,y,x)];



			}//z-loop
		}//y loop
	}//x loop

	for (z = 0 ; z< lz ; ++z){
		//loop for x=0 and x=lx-1!!!! (first and last slice)
		for (y = 0; y< ly; ++y){
			//toslice 0 antigrafetai sto 0
			D3_hlp.Q0[index(z,y,0)]=D3.Q0[index(z,y,0)];
			D3_hlp.Q1[index(z,y,0)]=D3.Q1[index(z,y,0)];
			D3_hlp.Q2[index(z,y,0)]=D3.Q2[index(z,y,0)];
			D3_hlp.Q3[index(z,y,0)]=D3.Q3[index(z,y,0)];
			D3_hlp.Q4[index(z,y,0)]=D3.Q4[index(z,y,0)];
			D3_hlp.Q5[index(z,y,0)]=D3.Q5[index(z,y,0)];
			D3_hlp.Q6[index(z,y,0)]=D3.Q6[index(z,y,0)];
			D3_hlp.Q7[index(z,y,0)]=D3.Q7[index(z,y,0)];
			D3_hlp.Q8[index(z,y,0)]=D3.Q8[index(z,y,0)];
			D3_hlp.Q9[index(z,y,0)]=D3.Q9[index(z,y,0)];
			D3_hlp.Q10[index(z,y,0)]=D3.Q10[index(z,y,0)];
			D3_hlp.Q11[index(z,y,0)]=D3.Q11[index(z,y,0)];
			D3_hlp.Q12[index(z,y,0)]=D3.Q12[index(z,y,0)];
			D3_hlp.Q13[index(z,y,0)]=D3.Q13[index(z,y,0)];
			D3_hlp.Q14[index(z,y,0)]=D3.Q14[index(z,y,0)];
			D3_hlp.Q15[index(z,y,0)]=D3.Q15[index(z,y,0)];
			D3_hlp.Q16[index(z,y,0)]=D3.Q16[index(z,y,0)];
			D3_hlp.Q17[index(z,y,0)]=D3.Q17[index(z,y,0)];
			D3_hlp.Q18[index(z,y,0)]=D3.Q18[index(z,y,0)];
			//at x= lx I set the incomming density as the one (equilibrum) calculated after the collision
			D3_hlp.Q3[index(z,y,lx-1)] = D3.Q3[index(z,y,lx-1)];
			D3_hlp.Q8[index(z,y,lx-1)] = D3.Q8[index(z,y,lx-1)] ;
			D3_hlp.Q9[index(z,y,lx-1)] = D3.Q9[index(z,y,lx-1)] ;
			D3_hlp.Q12[index(z,y,lx-1)] = D3.Q12[index(z,y,lx-1)] ;
			D3_hlp.Q13[index(z,y,lx-1)] = D3.Q13[index(z,y,lx-1)] ;
			//			gia x=lx-1: the densities 0,   2,4,5,6,15,16,17,18
			//				aplws 8a metadw8oune
			//				"ka8eta" sto slice kai de 8a ginoun propagate se
			//				alla slices (tuxainei na voleuei auto)
			//				auto to kommati exei HDH ginei pio panw!
			//
			//
			//				during streaming some of the indices at
			//				(fixed) x=lx-1, y=ly-1, do not get updated.
			//				it doesn't matter since the last slice on x is not useful.
		}
	}

	/*
	for (z = 0 ; z< lz ; ++z){
		z_l =  (z+1)%lz; //1 8esh meta to trexon x
		z_r =   (z+lz-1) %(  lz) ;

		//loop for x=0 and x=lx-1!!!! (first and last slice)
		for (y = 0; y< ly; ++y){
			y_n =  (y+1)%ly; //1 8esh meta to trexon x
			y_s =   (y+ly-1) %(  ly) ;



			//toslice 0 antigrafetai sto 0


				D3_hlp.Q0[index(z,y,0)]=D3.Q0[index(z,y,0)];
				D3_hlp.Q1[index(z,y,0)]=D3.Q1[index(z,y,0)];
				D3_hlp.Q2[index(z,y,0)]=D3.Q2[index(z,y,0)];
				D3_hlp.Q3[index(z,y,0)]=D3.Q3[index(z,y,0)];
				D3_hlp.Q4[index(z,y,0)]=D3.Q4[index(z,y,0)];
				D3_hlp.Q5[index(z,y,0)]=D3.Q5[index(z,y,0)];
				D3_hlp.Q6[index(z,y,0)]=D3.Q6[index(z,y,0)];
				D3_hlp.Q7[index(z,y,0)]=D3.Q7[index(z,y,0)];
				D3_hlp.Q8[index(z,y,0)]=D3.Q8[index(z,y,0)];
				D3_hlp.Q9[index(z,y,0)]=D3.Q9[index(z,y,0)];
				D3_hlp.Q10[index(z,y,0)]=D3.Q10[index(z,y,0)];
				D3_hlp.Q11[index(z,y,0)]=D3.Q11[index(z,y,0)];
				D3_hlp.Q12[index(z,y,0)]=D3.Q12[index(z,y,0)];
				D3_hlp.Q13[index(z,y,0)]=D3.Q13[index(z,y,0)];
				D3_hlp.Q14[index(z,y,0)]=D3.Q14[index(z,y,0)];
				D3_hlp.Q15[index(z,y,0)]=D3.Q15[index(z,y,0)];
				D3_hlp.Q16[index(z,y,0)]=D3.Q16[index(z,y,0)];
				D3_hlp.Q17[index(z,y,0)]=D3.Q17[index(z,y,0)];
				D3_hlp.Q18[index(z,y,0)]=D3.Q18[index(z,y,0)];


			//at x= lx I set the incomming density as the one (equilibrum) calculated after the collision
			D3_hlp.Q3[index(z,y,lx-1)] = D3.Q3[index(z,y,lx-1)];
			D3_hlp.Q8[index(z,y,lx-1)] = D3.Q8[index(z,y,lx-1)] ;
			D3_hlp.Q9[index(z,y,lx-1)] = D3.Q9[index(z,y,lx-1)] ;
			D3_hlp.Q12[index(z,y,lx-1)] = D3.Q12[index(z,y,lx-1)] ;
			D3_hlp.Q13[index(z,y,lx-1)] = D3.Q13[index(z,y,lx-1)] ;

			//tis aristeres densities apo to slice 0  tis metaferei sto telos

			D3_hlp.Q3[ index(z  ,y  ,lx-1)] = D3.Q3[ index(z,y,0)];
			D3_hlp.Q8[ index(z  ,y_n,lx-1)] = D3.Q8[ index(z,y,0)];
			D3_hlp.Q9[ index(z  ,y_s,lx-1)] = D3.Q9[ index(z,y,0)];
			D3_hlp.Q12[index(z_r,y  ,lx-1)] = D3.Q12[index(z,y,0)];
			D3_hlp.Q13[index(z_l,y  ,lx-1)] = D3.Q13[index(z,y,0)];

			//tis deksies densities tou teleutaiou slice tis metaferei sto slice 0
			D3_hlp.Q1[ index(z  ,y  ,0)] = D3.Q1[ index(z,y,lx-1)];
			D3_hlp.Q7[ index(z  ,y_n,0)] = D3.Q7[ index(z,y,lx-1)];
			D3_hlp.Q10[index(z  ,y_s,0)] = D3.Q10[index(z,y,lx-1)];
			D3_hlp.Q11[index(z_r,y  ,0)] = D3.Q11[index(z,y,lx-1)];
			D3_hlp.Q14[index(z_l,y  ,0)] = D3.Q14[index(z,y,lx-1)];

//			gia x=lx-1: the densities 0,   2,4,5,6,15,16,17,18
//				aplws 8a metadw8oune
//				"ka8eta" sto slice kai de 8a ginoun propagate se
//				alla slices (tuxainei na voleuei auto)
//				auto to kommati exei HDH ginei pio panw!
//
//
//				during streaming some of the indices at
//				(fixed) x=lx-1, y=ly-1, do not get updated.
//				it doesn't matter since the last slice on x is not useful.

		}
	}*/
#ifdef DEBUG
	cout << " #LBM streaming OK!" << endl;
#endif
}

void LBM::streaming_first_part(){



	/*Propagate fluid densities to their next neighbour nodes */
	/*c
	c.......density propagation: all fluid densities are propagated from
	c       non-occupied nodes along the lattice connection lines
	c       to their next neighbours.
	c*/
	int  x,y,z;
	int x_e/*east*/,x_w/*west*/;
	int y_n/*north*/,y_s/*south*/;
	int z_l/*left*/,z_r/*right*/;
	//todo rwta ton Dirk pws to kanei auto to vhma
	for (z = 0 ; z< lz ; ++z){
		z_l =  (z+1)%lz; //1 8esh meta to trexon x
		z_r =   (z+lz-1) %(  lz) ;

		for (y = 0; y< ly; ++y){
			y_n =  (y+1)%ly; //1 8esh meta to trexon y
			y_s =   (y+ly-1) %(  ly) ;

			for (x = 0; x< lx; ++x){
				x_e =  (x+1)%lx; //1 8esh meta to trexon z
				x_w =   (x+lx-1) %(  lx) ;

				//regular streaming process

				/*
					.........density propagation

					.........zero: just copy
				 */
				D3_hlp.Q0[index(z,y,x)]=D3.Q0[index(z,y,x)];
				//				node_hlp_grid[index4D(z,y,x,0)]=node_grid[index4D(z,y,x,0)];

				//.........in the x,y and z directions
				D3_hlp.Q1[index(z  ,y  ,x_e)] = D3.Q1[index(z,y,x)];
				D3_hlp.Q2[index(z  ,y_n,x  )] = D3.Q2[index(z,y,x)];
				D3_hlp.Q3[index(z  ,y  ,x_w)] = D3.Q3[index(z,y,x)];
				D3_hlp.Q4[index(z  ,y_s,x  )] = D3.Q4[index(z,y,x)];
				D3_hlp.Q5[index(z_l,y  ,x  )] = D3.Q5[index(z,y,x)];
				D3_hlp.Q6[index(z_r,y  ,x  )] = D3.Q6[index(z,y,x)];


				//......... in the x,y diagonals


				D3_hlp.Q7[ index(z  ,y_n,x_e)] = D3.Q7[ index(z,y,x)];
				D3_hlp.Q8[ index(z  ,y_n,x_w)] = D3.Q8[ index(z,y,x)];
				D3_hlp.Q9[ index(z  ,y_s,x_w)] = D3.Q9[ index(z,y,x)];
				D3_hlp.Q10[index(z  ,y_s,x_e)] = D3.Q10[index(z,y,x)];

				//......... in the x,z diagonals

				D3_hlp.Q11[index(z_r,y  ,x_e)] = D3.Q11[index(z,y,x)];
				D3_hlp.Q12[index(z_r,y  ,x_w)] = D3.Q12[index(z,y,x)];
				D3_hlp.Q13[index(z_l,y  ,x_w)] = D3.Q13[index(z,y,x)];
				D3_hlp.Q14[index(z_l,y  ,x_e)] = D3.Q14[index(z,y,x)];

				//......... in the y,z diagonals

				D3_hlp.Q15[index(z_l,y_n,x  )] = D3.Q15[index(z,y,x)];
				D3_hlp.Q16[index(z_r,y_n,x  )] = D3.Q16[index(z,y,x)];
				D3_hlp.Q17[index(z_r,y_s,x  )] = D3.Q17[index(z,y,x)];
				D3_hlp.Q18[index(z_l,y_s,x  )] = D3.Q18[index(z,y,x)];



			}//z-loop
		}//y loop
	}//x loop


#ifdef DEBUG
	cout << " #LBM streaming OK!" << endl;
#endif
}

void LBM::streaming_last_part(){
	/*Propagate fluid densities to their next neighbour nodes */
	/*c
	c.......density propagation: all fluid densities are propagated from
	c       non-occupied nodes along the lattice connection lines
	c       to their next neighbours.
	c*/

	int   y,z;

	for (z = 0 ; z< lz ; ++z){
		//loop for x=0 and x=lx-1!!!! (first and last slice)
		for (y = 0; y< ly; ++y){
			//toslice 0 antigrafetai sto 0
			D3_hlp.Q0[index(z,y,0)]=D3.Q0[index(z,y,0)];
			D3_hlp.Q1[index(z,y,0)]=D3.Q1[index(z,y,0)];
			D3_hlp.Q2[index(z,y,0)]=D3.Q2[index(z,y,0)];
			D3_hlp.Q3[index(z,y,0)]=D3.Q3[index(z,y,0)];
			D3_hlp.Q4[index(z,y,0)]=D3.Q4[index(z,y,0)];
			D3_hlp.Q5[index(z,y,0)]=D3.Q5[index(z,y,0)];
			D3_hlp.Q6[index(z,y,0)]=D3.Q6[index(z,y,0)];
			D3_hlp.Q7[index(z,y,0)]=D3.Q7[index(z,y,0)];
			D3_hlp.Q8[index(z,y,0)]=D3.Q8[index(z,y,0)];
			D3_hlp.Q9[index(z,y,0)]=D3.Q9[index(z,y,0)];
			D3_hlp.Q10[index(z,y,0)]=D3.Q10[index(z,y,0)];
			D3_hlp.Q11[index(z,y,0)]=D3.Q11[index(z,y,0)];
			D3_hlp.Q12[index(z,y,0)]=D3.Q12[index(z,y,0)];
			D3_hlp.Q13[index(z,y,0)]=D3.Q13[index(z,y,0)];
			D3_hlp.Q14[index(z,y,0)]=D3.Q14[index(z,y,0)];
			D3_hlp.Q15[index(z,y,0)]=D3.Q15[index(z,y,0)];
			D3_hlp.Q16[index(z,y,0)]=D3.Q16[index(z,y,0)];
			D3_hlp.Q17[index(z,y,0)]=D3.Q17[index(z,y,0)];
			D3_hlp.Q18[index(z,y,0)]=D3.Q18[index(z,y,0)];
			//at x= lx I set the incomming density as the one (equilibrum) calculated after the collision
			D3_hlp.Q3[index(z,y,lx-1)] = D3.Q3[index(z,y,lx-1)];
			D3_hlp.Q8[index(z,y,lx-1)] = D3.Q8[index(z,y,lx-1)] ;
			D3_hlp.Q9[index(z,y,lx-1)] = D3.Q9[index(z,y,lx-1)] ;
			D3_hlp.Q12[index(z,y,lx-1)] = D3.Q12[index(z,y,lx-1)] ;
			D3_hlp.Q13[index(z,y,lx-1)] = D3.Q13[index(z,y,lx-1)] ;
			//			gia x=lx-1: the densities 0,   2,4,5,6,15,16,17,18
			//				aplws 8a metadw8oune
			//				"ka8eta" sto slice kai de 8a ginoun propagate se
			//				alla slices (tuxainei na voleuei auto)
			//				auto to kommati exei HDH ginei pio panw!
			//
			//
			//				during streaming some of the indices at
			//				(fixed) x=lx-1, y=ly-1, do not get updated.
			//				it doesn't matter since the last slice on x is not useful.
		}
	}
}




__global__
void streaming_kernel(int lx, int ly, int lz, FLOATING reynolds, FLOATING nu, FLOATING r_small, FLOATING t_0, FLOATING t_1, FLOATING t_2, FLOATING c_squ, lattice D3, lattice D3_hlp){



	/*Propagate fluid densities to their next neighbour nodes */
	/*c
	c.......density propagation: all fluid densities are propagated from
	c       non-occupied nodes along the lattice connection lines
	c       to their next neighbours.
	c*/




	int  x,y,z;
	int x_e/*east*/,x_w/*west*/;
	int y_n/*north*/,y_s/*south*/;
	int z_l/*left*/,z_r/*right*/;


	int rest;
	int end_of_memory=lz*ly*(lx);
	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	z=(int) (tid/(ly*lx));
	rest=tid-z;
	y=(int)(rest/lx);
	x=rest-y;

	//if(tid<end_of_memory){
	if( z<lz and y<ly and x<lx){

		//todo rwta ton Dirk pws to kanei auto to vhma
		//		for (z = 0 ; z< lz ; ++z){
		z_l =  (z+1)%lz; //1 8esh meta to trexon x
		z_r =   (z+lz-1) %(  lz) ;

		//			for (y = 0; y< ly; ++y){
		y_n =  (y+1)%ly; //1 8esh meta to trexon y
		y_s =   (y+ly-1) %(  ly) ;

		//				for (x = 0; x< lx; ++x){
		x_e =  (x+1)%lx; //1 8esh meta to trexon z
		x_w =   (x+lx-1) %(  lx) ;

		//regular streaming process

		/*
					.........density propagation

					.........zero: just copy
		 */
		D3_hlp.Q0[index(z,y,x)]=D3.Q0[index(z,y,x)];
		//node_hlp_grid[index4D(z,y,x,0)]=node_grid[index4D(z,y,x,0)];

		//.........in the x,y and z directions
		D3_hlp.Q1[index(z  ,y  ,x_e)] = D3.Q1[index(z,y,x)];
		D3_hlp.Q2[index(z  ,y_n,x  )] = D3.Q2[index(z,y,x)];
		D3_hlp.Q3[index(z  ,y  ,x_w)] = D3.Q3[index(z,y,x)];
		D3_hlp.Q4[index(z  ,y_s,x  )] = D3.Q4[index(z,y,x)];
		D3_hlp.Q5[index(z_l,y  ,x  )] = D3.Q5[index(z,y,x)];
		D3_hlp.Q6[index(z_r,y  ,x  )] = D3.Q6[index(z,y,x)];


		//......... in the x,y diagonals


		D3_hlp.Q7[ index(z  ,y_n,x_e)] = D3.Q7[ index(z,y,x)];
		D3_hlp.Q8[ index(z  ,y_n,x_w)] = D3.Q8[ index(z,y,x)];
		D3_hlp.Q9[ index(z  ,y_s,x_w)] = D3.Q9[ index(z,y,x)];
		D3_hlp.Q10[index(z  ,y_s,x_e)] = D3.Q10[index(z,y,x)];

		//......... in the x,z diagonals

		D3_hlp.Q11[index(z_r,y  ,x_e)] = D3.Q11[index(z,y,x)];
		D3_hlp.Q12[index(z_r,y  ,x_w)] = D3.Q12[index(z,y,x)];
		D3_hlp.Q13[index(z_l,y  ,x_w)] = D3.Q13[index(z,y,x)];
		D3_hlp.Q14[index(z_l,y  ,x_e)] = D3.Q14[index(z,y,x)];

		//......... in the y,z diagonals

		D3_hlp.Q15[index(z_l,y_n,x  )] = D3.Q15[index(z,y,x)];
		D3_hlp.Q16[index(z_r,y_n,x  )] = D3.Q16[index(z,y,x)];
		D3_hlp.Q17[index(z_r,y_s,x  )] = D3.Q17[index(z,y,x)];
		D3_hlp.Q18[index(z_l,y_s,x  )] = D3.Q18[index(z,y,x)];



		//				}//z-loop
		//			}//y loop
		//		}//x loop
	}

	if(tid<end_of_memory and x==0 ){
		//		for (z = 0 ; z< lz ; ++z){
		//
		//
		//			//loop for x=0 and x=lx-1!!!! (first and last slice)
		//			for (y = 0; y< ly; ++y){
		//toslice 0 antigrafetai sto 0
		D3_hlp.Q0[index(z,y,0)]=D3.Q0[index(z,y,0)];
		D3_hlp.Q1[index(z,y,0)]=D3.Q1[index(z,y,0)];
		D3_hlp.Q2[index(z,y,0)]=D3.Q2[index(z,y,0)];
		D3_hlp.Q3[index(z,y,0)]=D3.Q3[index(z,y,0)];
		D3_hlp.Q4[index(z,y,0)]=D3.Q4[index(z,y,0)];
		D3_hlp.Q5[index(z,y,0)]=D3.Q5[index(z,y,0)];
		D3_hlp.Q6[index(z,y,0)]=D3.Q6[index(z,y,0)];
		D3_hlp.Q7[index(z,y,0)]=D3.Q7[index(z,y,0)];
		D3_hlp.Q8[index(z,y,0)]=D3.Q8[index(z,y,0)];
		D3_hlp.Q9[index(z,y,0)]=D3.Q9[index(z,y,0)];
		D3_hlp.Q10[index(z,y,0)]=D3.Q10[index(z,y,0)];
		D3_hlp.Q11[index(z,y,0)]=D3.Q11[index(z,y,0)];
		D3_hlp.Q12[index(z,y,0)]=D3.Q12[index(z,y,0)];
		D3_hlp.Q13[index(z,y,0)]=D3.Q13[index(z,y,0)];
		D3_hlp.Q14[index(z,y,0)]=D3.Q14[index(z,y,0)];
		D3_hlp.Q15[index(z,y,0)]=D3.Q15[index(z,y,0)];
		D3_hlp.Q16[index(z,y,0)]=D3.Q16[index(z,y,0)];
		D3_hlp.Q17[index(z,y,0)]=D3.Q17[index(z,y,0)];
		D3_hlp.Q18[index(z,y,0)]=D3.Q18[index(z,y,0)];



		//			gia x=lx-1: the densities 0,   2,4,5,6,15,16,17,18
		//				aplws 8a metadw8oune
		//				"ka8eta" sto slice kai de 8a ginoun propagate se
		//				alla slices (tuxainei na voleuei auto)
		//				auto to kommati exei HDH ginei pio panw!
		//
		//
		//				during streaming some of the indices at
		//				(fixed) x=lx-1, y=ly-1, do not get updated.
		//				it doesn't matter since the last slice on x is not useful.

		//			}
		//		}
	}
	if(tid<end_of_memory and   x==lx-1  ){
		//at x= lx I set the incomming density as the one (equilibrum) calculated after the collision
		D3_hlp.Q3[index(z,y,lx-1)] = D3.Q3[index(z,y,lx-1)];
		D3_hlp.Q8[index(z,y,lx-1)] = D3.Q8[index(z,y,lx-1)] ;
		D3_hlp.Q9[index(z,y,lx-1)] = D3.Q9[index(z,y,lx-1)] ;
		D3_hlp.Q12[index(z,y,lx-1)] = D3.Q12[index(z,y,lx-1)] ;
		D3_hlp.Q13[index(z,y,lx-1)] = D3.Q13[index(z,y,lx-1)] ;
	}


#ifdef DEBUG
	cout << " #LBM streaming OK!" << endl;
#endif
}






__global__
void streaming_kernel_single_threaded_p1(int lx, int ly, int lz,
		FLOATING *hlp_Q0, FLOATING *hlp_Q1, FLOATING *hlp_Q2, FLOATING *hlp_Q3,
		FLOATING *hlp_Q4, FLOATING *hlp_Q5, FLOATING *hlp_Q6, FLOATING *hlp_Q7,
		FLOATING *hlp_Q8, FLOATING *hlp_Q9, FLOATING *hlp_Q10, FLOATING *hlp_Q11,
		FLOATING *hlp_Q12, FLOATING *hlp_Q13, FLOATING *hlp_Q14, FLOATING *hlp_Q15,
		FLOATING *hlp_Q16, FLOATING *hlp_Q17, FLOATING *hlp_Q18,
		FLOATING *Q0, FLOATING *Q1, FLOATING *Q2, FLOATING *Q3,
		FLOATING *Q4, FLOATING *Q5, FLOATING *Q6, FLOATING *Q7,
		FLOATING *Q8, FLOATING *Q9, FLOATING *Q10, FLOATING *Q11,
		FLOATING *Q12, FLOATING *Q13, FLOATING *Q14, FLOATING *Q15,
		FLOATING *Q16, FLOATING *Q17, FLOATING *Q18){


	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	int end_of_memory=lz*ly*(lx);
	int z=(int) (tid/(ly*lx));
	int	y=(tid-z*(ly*lx))/lx;
	int	x=tid-z*(ly*lx)-y*lx;

	//	int x,y,z;
	int x_e/*east*/,x_w/*west*/;
	int y_n/*north*/,y_s/*south*/;
	int z_l/*left*/,z_r/*right*/;

	//	__shared__ FLOATING shared_buffer[64];


	if(tid<end_of_memory){

		z_l =  (z+1)%lz; //1 8esh meta to trexon x
		z_r =   (z+lz-1) %(  lz) ;

		//				for (y = 0; y< ly; ++y){
		y_n =  (y+1)%ly; //1 8esh meta to trexon y
		y_s =   (y+ly-1) %(  ly) ;

		//					for (x = 0; x< lx; ++x){
		x_e =  (x+1)%lx; //1 8esh meta to trexon z
		x_w =   (x+lx-1) %(  lx) ;

		//regular streaming process

		//					.........density propagation
		//					.........zero: just copy
		hlp_Q0[index(z,y,x)]=Q0[index(z,y,x)];
		//				node_hlp_grid[index4D(z,y,x,0)]=node_grid[index4D(z,y,x,0)];

		//.........in the x,y and z directions
		hlp_Q1[index(z  ,y  ,x_e)] = Q1[index(z,y,x)];
		hlp_Q2[index(z  ,y_n,x  )] = Q2[index(z,y,x)];
		hlp_Q3[index(z  ,y  ,x_w)] = Q3[index(z,y,x)];
		hlp_Q4[index(z  ,y_s,x  )] = Q4[index(z,y,x)];
		hlp_Q5[index(z_l,y  ,x  )] = Q5[index(z,y,x)];
		hlp_Q6[index(z_r,y  ,x  )] = Q6[index(z,y,x)];

		//......... in the x,y diagonals
		hlp_Q7[ index(z  ,y_n,x_e)] = Q7[ index(z,y,x)];
		hlp_Q8[ index(z  ,y_n,x_w)] = Q8[ index(z,y,x)];
		hlp_Q9[ index(z  ,y_s,x_w)] = Q9[ index(z,y,x)];
		hlp_Q10[index(z  ,y_s,x_e)] = Q10[index(z,y,x)];

		//......... in the x,z diagonals
		hlp_Q11[index(z_r,y  ,x_e)] = Q11[index(z,y,x)];
		hlp_Q12[index(z_r,y  ,x_w)] = Q12[index(z,y,x)];
		hlp_Q13[index(z_l,y  ,x_w)] = Q13[index(z,y,x)];
		hlp_Q14[index(z_l,y  ,x_e)] = Q14[index(z,y,x)];

		//......... in the y,z diagonals
		hlp_Q15[index(z_l,y_n,x  )] = Q15[index(z,y,x)];
		hlp_Q16[index(z_r,y_n,x  )] = Q16[index(z,y,x)];
		hlp_Q17[index(z_r,y_s,x  )] = Q17[index(z,y,x)];
		hlp_Q18[index(z_l,y_s,x  )] = Q18[index(z,y,x)];
	}
}

__global__
void streaming_kernel_single_threaded_p1_shared(int lx, int ly, int lz,
		FLOATING *hlp_Q0, FLOATING *hlp_Q1, FLOATING *hlp_Q2, FLOATING *hlp_Q3,
		FLOATING *hlp_Q4, FLOATING *hlp_Q5, FLOATING *hlp_Q6, FLOATING *hlp_Q7,
		FLOATING *hlp_Q8, FLOATING *hlp_Q9, FLOATING *hlp_Q10, FLOATING *hlp_Q11,
		FLOATING *hlp_Q12, FLOATING *hlp_Q13, FLOATING *hlp_Q14, FLOATING *hlp_Q15,
		FLOATING *hlp_Q16, FLOATING *hlp_Q17, FLOATING *hlp_Q18,
		FLOATING *Q0, FLOATING *Q1, FLOATING *Q2, FLOATING *Q3,
		FLOATING *Q4, FLOATING *Q5, FLOATING *Q6, FLOATING *Q7,
		FLOATING *Q8, FLOATING *Q9, FLOATING *Q10, FLOATING *Q11,
		FLOATING *Q12, FLOATING *Q13, FLOATING *Q14, FLOATING *Q15,
		FLOATING *Q16, FLOATING *Q17, FLOATING *Q18){


	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	int end_of_memory=lz*ly*(lx);
	int z=(int) (tid/(ly*lx));
	int	y=(tid-z*(ly*lx))/lx;
	int	x=tid-z*(ly*lx)-y*lx;

	//	int x,y,z;
	int x_e/*east*/,x_w/*west*/;
	int y_n/*north*/,y_s/*south*/;
	int z_l/*left*/,z_r/*right*/;

	extern __shared__ FLOATING shared_buffer[];
	FLOATING *shared_Q0=shared_buffer;
	FLOATING *shared_Q1 = &shared_Q0[blockDim.x];
	FLOATING *shared_Q2 = &shared_Q1[blockDim.x];
	FLOATING *shared_Q3 = &shared_Q2[blockDim.x];
	FLOATING *shared_Q4 = &shared_Q3[blockDim.x];
	FLOATING *shared_Q5 = &shared_Q4[blockDim.x];
	FLOATING *shared_Q6 = &shared_Q5[blockDim.x];
	FLOATING *shared_Q7 = &shared_Q6[blockDim.x];
	FLOATING *shared_Q8 = &shared_Q7[blockDim.x];
	FLOATING *shared_Q9 = &shared_Q8[blockDim.x];
	FLOATING *shared_Q10 = &shared_Q9[blockDim.x];
	FLOATING *shared_Q11 = &shared_Q10[blockDim.x];
	FLOATING *shared_Q12 = &shared_Q11[blockDim.x];
	FLOATING *shared_Q13 = &shared_Q12[blockDim.x];
	FLOATING *shared_Q14 = &shared_Q13[blockDim.x];
	FLOATING *shared_Q15 = &shared_Q14[blockDim.x];
	FLOATING *shared_Q16 = &shared_Q15[blockDim.x];
	FLOATING *shared_Q17 = &shared_Q16[blockDim.x];
	FLOATING *shared_Q18 = &shared_Q17[blockDim.x];



	shared_Q0[threadIdx.x]=Q0[index(z,y,x)];
	shared_Q1[threadIdx.x]=Q1[index(z,y,x)];
	shared_Q2[threadIdx.x]=Q2[index(z,y,x)];
	shared_Q3[threadIdx.x]=Q3[index(z,y,x)];
	shared_Q4[threadIdx.x]=Q4[index(z,y,x)];
	shared_Q5[threadIdx.x]=Q5[index(z,y,x)];
	shared_Q6[threadIdx.x]=Q6[index(z,y,x)];
	shared_Q7[threadIdx.x]=Q7[index(z,y,x)];
	shared_Q8[threadIdx.x]=Q8[index(z,y,x)];
	shared_Q9[threadIdx.x]=Q9[index(z,y,x)];
	shared_Q10[threadIdx.x]=Q10[index(z,y,x)];
	shared_Q11[threadIdx.x]=Q11[index(z,y,x)];
	shared_Q12[threadIdx.x]=Q12[index(z,y,x)];
	shared_Q13[threadIdx.x]=Q13[index(z,y,x)];
	shared_Q14[threadIdx.x]=Q14[index(z,y,x)];
	shared_Q15[threadIdx.x]=Q15[index(z,y,x)];
	shared_Q16[threadIdx.x]=Q16[index(z,y,x)];
	shared_Q17[threadIdx.x]=Q17[index(z,y,x)];
	shared_Q18[threadIdx.x]=Q18[index(z,y,x)];
	__syncthreads();
	if(tid<end_of_memory){

		z_l =  (z+1)%lz; //1 8esh meta to trexon x
		z_r =   (z+lz-1) %(  lz) ;

		//				for (y = 0; y< ly; ++y){
		y_n =  (y+1)%ly; //1 8esh meta to trexon y
		y_s =   (y+ly-1) %(  ly) ;

		//					for (x = 0; x< lx; ++x){
		x_e =  (x+1)%lx; //1 8esh meta to trexon z
		x_w =   (x+lx-1) %(  lx) ;

		//regular streaming process

		//					.........density propagation
		//					.........zero: just copy
		hlp_Q0[index(z,y,x)]=shared_Q0[threadIdx.x];
		//				node_hlp_grid[index4D(z,y,x,0)]=node_grid[index4D(z,y,x,0)];

		//.........in the x,y and z directions
		hlp_Q1[index(z  ,y  ,x_e)] =shared_Q1[threadIdx.x];
		hlp_Q2[index(z  ,y_n,x  )] =shared_Q2[threadIdx.x];
		hlp_Q3[index(z  ,y  ,x_w)] =shared_Q3[threadIdx.x];
		hlp_Q4[index(z  ,y_s,x  )] =shared_Q4[threadIdx.x];
		hlp_Q5[index(z_l,y  ,x  )] =shared_Q5[threadIdx.x];
		hlp_Q6[index(z_r,y  ,x  )] =shared_Q6[threadIdx.x];

		//......... in the x,y diagonals
		hlp_Q7[ index(z  ,y_n,x_e)] =shared_Q7[ threadIdx.x];
		hlp_Q8[ index(z  ,y_n,x_w)] =shared_Q8[ threadIdx.x];
		hlp_Q9[ index(z  ,y_s,x_w)] =shared_Q9[ threadIdx.x];
		hlp_Q10[index(z  ,y_s,x_e)] =shared_Q10[threadIdx.x];

		//......... in the x,z diagonals
		hlp_Q11[index(z_r,y  ,x_e)] =shared_Q11[threadIdx.x];
		hlp_Q12[index(z_r,y  ,x_w)] =shared_Q12[threadIdx.x];
		hlp_Q13[index(z_l,y  ,x_w)] =shared_Q13[threadIdx.x];
		hlp_Q14[index(z_l,y  ,x_e)] =shared_Q14[threadIdx.x];

		//......... in the y,z diagonals
		hlp_Q15[index(z_l,y_n,x  )] =shared_Q15[threadIdx.x];
		hlp_Q16[index(z_r,y_n,x  )] =shared_Q16[threadIdx.x];
		hlp_Q17[index(z_r,y_s,x  )] =shared_Q17[threadIdx.x];
		hlp_Q18[index(z_l,y_s,x  )] =shared_Q18[threadIdx.x];
	}
}

__global__
void streaming_kernel_first_part_Q0(int lx, int ly, int lz, const FLOATING *Q0,  FLOATING *hlp_Q0){
	/*Propagate fluid densities to their next neighbour nodes */
	/*c....density propagation: all fluid densities are propagated from
	c       non-occupied nodes along the lattice connection lines
	c       to their next neighbours.*/
	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	int end_of_memory=lz*ly*(lx);
	int z=(int) (tid/(ly*lx));
	int	y=(tid-z*(ly*lx))/lx;
	int	x=tid-z*(ly*lx)-y*lx;

	extern __shared__ FLOATING shared_buffer[];

	shared_buffer[threadIdx.x]=Q0[index(z,y,x)];
	__syncthreads();

	if( tid<end_of_memory)
		hlp_Q0[index(z,y,x)] = shared_buffer[threadIdx.x];
}

__global__
void streaming_kernel_first_part_Q1(int lx, int ly, int lz, const FLOATING *Q1,  FLOATING *hlp_Q1){
	/*Propagate fluid densities to their next neighbour nodes */
	/*c....density propagation: all fluid densities are propagated from
	c       non-occupied nodes along the lattice connection lines
	c       to their next neighbours.*/


	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	int end_of_memory=lz*ly*(lx);
	int z=(int) (tid/(ly*lx));
	int	y=(tid-z*(ly*lx))/lx;
	int	x=tid-z*(ly*lx)-y*lx;

	int x_e/*east*/ =  (x+1)%lx; //1 8esh meta to trexon x

	extern __shared__ FLOATING shared_buffer[];

	shared_buffer[threadIdx.x]=Q1[index(z,y,x)];
	__syncthreads();

	if( tid<end_of_memory)
		hlp_Q1[index(z  ,y  ,x_e)] = shared_buffer[threadIdx.x];

}

__global__
void streaming_kernel_first_part_Q2(int lx, int ly, int lz, const FLOATING *Q2,  FLOATING *hlp_Q2){
	/*Propagate fluid densities to their next neighbour nodes */
	/*c....density propagation: all fluid densities are propagated from
	c       non-occupied nodes along the lattice connection lines
	c       to their next neighbours.*/


	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	int end_of_memory=lz*ly*(lx);
	int z=(int) (tid/(ly*lx));
	int	y=(tid-z*(ly*lx))/lx;
	int	x=tid-z*(ly*lx)-y*lx;

	int y_n/*north*/ =  (y+1)%ly; //1 8esh meta to trexon y

	extern __shared__ FLOATING shared_buffer[];

	shared_buffer[threadIdx.x]=Q2[index(z,y,x)];
	__syncthreads();

	if( tid<end_of_memory)
		hlp_Q2[index(z  ,y_n,x  )] = shared_buffer[threadIdx.x];

}

__global__
void streaming_kernel_first_part_Q3(int lx, int ly, int lz, const FLOATING *Q3,  FLOATING *hlp_Q3){
	/*Propagate fluid densities to their next neighbour nodes */
	/*c....density propagation: all fluid densities are propagated from
	c       non-occupied nodes along the lattice connection lines
	c       to their next neighbours.*/

	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	int end_of_memory=lz*ly*(lx);
	int z=(int) (tid/(ly*lx));
	int	y=(tid-z*(ly*lx))/lx;
	int	x=tid-z*(ly*lx)-y*lx;

	int x_w/*west*/ =   (x+lx-1) %(  lx) ;

	extern __shared__ FLOATING shared_buffer[];

	shared_buffer[threadIdx.x]=Q3[index(z,y,x)];
	__syncthreads();

	if( tid<end_of_memory)
		hlp_Q3[index(z  ,y  ,x_w)] = shared_buffer[threadIdx.x];

}

__global__
void streaming_kernel_first_part_Q4(int lx, int ly, int lz, const FLOATING *Q4,  FLOATING *hlp_Q4){
	/*Propagate fluid densities to their next neighbour nodes */
	/*c....density propagation: all fluid densities are propagated from
	c       non-occupied nodes along the lattice connection lines
	c       to their next neighbours.*/


	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	int end_of_memory=lz*ly*(lx);
	int z=(int) (tid/(ly*lx));
	int	y=(tid-z*(ly*lx))/lx;
	int	x=tid-z*(ly*lx)-y*lx;

	int y_s/*south*/ =   (y+ly-1) %(  ly) ;

	extern __shared__ FLOATING shared_buffer[];

	shared_buffer[threadIdx.x]=Q4[index(z,y,x)];
	__syncthreads();

	if( tid<end_of_memory)
		hlp_Q4[index(z  ,y_s,x  )] = shared_buffer[threadIdx.x];

}

__global__
void streaming_kernel_first_part_Q5(int lx, int ly, int lz, const FLOATING *Q5,  FLOATING *hlp_Q5){
	/*Propagate fluid densities to their next neighbour nodes */
	/*c....density propagation: all fluid densities are propagated from
	c       non-occupied nodes along the lattice connection lines
	c       to their next neighbours.*/

	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	int end_of_memory=lz*ly*(lx);
	int z=(int) (tid/(ly*lx));
	int	y=(tid-z*(ly*lx))/lx;
	int	x=tid-z*(ly*lx)-y*lx;
	int z_l/*left*/ =  (z+1)%lz;

	extern __shared__ FLOATING shared_buffer[];

	shared_buffer[threadIdx.x]=Q5[index(z,y,x)];
	__syncthreads();

	if( tid<end_of_memory)
		hlp_Q5[index(z_l,y  ,x  )] = shared_buffer[threadIdx.x];

}

__global__
void streaming_kernel_first_part_Q6(int lx, int ly, int lz, const FLOATING *Q6,  FLOATING *hlp_Q6){
	/*Propagate fluid densities to their next neighbour nodes */
	/*c....density propagation: all fluid densities are propagated from
	c       non-occupied nodes along the lattice connection lines
	c       to their next neighbours.*/

	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	int end_of_memory=lz*ly*(lx);
	int z=(int) (tid/(ly*lx));
	int	y=(tid-z*(ly*lx))/lx;
	int	x=tid-z*(ly*lx)-y*lx;
	int z_r/*right*/ =   (z+lz-1) %(  lz) ;

	extern __shared__ FLOATING shared_buffer[];

	shared_buffer[threadIdx.x]=Q6[index(z,y,x)];
	__syncthreads();

	if( tid<end_of_memory)
		hlp_Q6[index(z_r,y  ,x  )] = shared_buffer[threadIdx.x];

}

__global__
void streaming_kernel_first_part_Q7(int lx, int ly, int lz, const FLOATING *Q7,  FLOATING *hlp_Q7){
	/*Propagate fluid densities to their next neighbour nodes */
	/*c....density propagation: all fluid densities are propagated from
	c       non-occupied nodes along the lattice connection lines
	c       to their next neighbours.*/
	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	int end_of_memory=lz*ly*(lx);
	int z=(int) (tid/(ly*lx));
	int	y=(tid-z*(ly*lx))/lx;
	int	x=tid-z*(ly*lx)-y*lx;

	int y_n/*north*/ =  (y+1)%ly;
	int x_e/*east*/ =  (x+1)%lx;

	extern __shared__ FLOATING shared_buffer[];

	shared_buffer[threadIdx.x]=Q7[index(z,y,x)];
	__syncthreads();

	if( tid<end_of_memory)
		hlp_Q7[ index(z  ,y_n,x_e)] = shared_buffer[threadIdx.x];

}

__global__
void streaming_kernel_first_part_Q8(int lx, int ly, int lz, const FLOATING *Q8,  FLOATING *hlp_Q8){
	/*Propagate fluid densities to their next neighbour nodes */
	/*c....density propagation: all fluid densities are propagated from
	c       non-occupied nodes along the lattice connection lines
	c       to their next neighbours.*/


	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	int end_of_memory=lz*ly*(lx);
	int z=(int) (tid/(ly*lx));
	int	y=(tid-z*(ly*lx))/lx;
	int	x=tid-z*(ly*lx)-y*lx;

	int y_n/*north*/ =  (y+1)%ly;
	int x_w/*west*/ =   (x+lx-1) %(  lx) ;

	extern __shared__ FLOATING shared_buffer[];

	shared_buffer[threadIdx.x]=Q8[index(z,y,x)];
	__syncthreads();

	if( tid<end_of_memory)
		hlp_Q8[ index(z  ,y_n,x_w)] = shared_buffer[threadIdx.x];

}

__global__
void streaming_kernel_first_part_Q9(int lx, int ly, int lz, const FLOATING *Q9,  FLOATING *hlp_Q9){
	/*Propagate fluid densities to their next neighbour nodes */
	/*c....density propagation: all fluid densities are propagated from
	c       non-occupied nodes along the lattice connection lines
	c       to their next neighbours.*/
	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	int end_of_memory=lz*ly*(lx);
	int z=(int) (tid/(ly*lx));
	int	y=(tid-z*(ly*lx))/lx;
	int	x=tid-z*(ly*lx)-y*lx;

	int y_s/*south*/ =   (y+ly-1) %(  ly) ;
	int x_w/*west*/ =   (x+lx-1) %(  lx) ;

	extern __shared__ FLOATING shared_buffer[];

	shared_buffer[threadIdx.x]=Q9[index(z,y,x)];
	__syncthreads();

	if( tid<end_of_memory)
		hlp_Q9[ index(z  ,y_s,x_w)] = shared_buffer[threadIdx.x];

}

__global__
void streaming_kernel_first_part_Q10(int lx, int ly, int lz, const FLOATING *Q10,  FLOATING *hlp_Q10){
	/*Propagate fluid densities to their next neighbour nodes */
	/*c....density propagation: all fluid densities are propagated from
	c       non-occupied nodes along the lattice connection lines
	c       to their next neighbours.*/
	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	int end_of_memory=lz*ly*(lx);
	int z=(int) (tid/(ly*lx));
	int	y=(tid-z*(ly*lx))/lx;
	int	x=tid-z*(ly*lx)-y*lx;

	int y_s/*south*/ =   (y+ly-1) %(  ly) ;
	int x_e/*east*/ =  (x+1)%lx;

	extern __shared__ FLOATING shared_buffer[];

	shared_buffer[threadIdx.x]=Q10[index(z,y,x)];
	__syncthreads();

	if( tid<end_of_memory)
		hlp_Q10[index(z  ,y_s,x_e)] = shared_buffer[threadIdx.x];

}

__global__
void streaming_kernel_first_part_Q11(int lx, int ly, int lz, const FLOATING *Q11,  FLOATING *hlp_Q11){
	/*Propagate fluid densities to their next neighbour nodes */
	/*c....density propagation: all fluid densities are propagated from
	c       non-occupied nodes along the lattice connection lines
	c       to their next neighbours.*/
	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	int end_of_memory=lz*ly*(lx);
	int z=(int) (tid/(ly*lx));
	int	y=(tid-z*(ly*lx))/lx;
	int	x=tid-z*(ly*lx)-y*lx;

	int z_r/*right*/ =   (z+lz-1) %(  lz) ;
	int x_e/*east*/ =  (x+1)%lx;

	extern __shared__ FLOATING shared_buffer[];

	shared_buffer[threadIdx.x]=Q11[index(z,y,x)];
	__syncthreads();

	if( tid<end_of_memory)
		hlp_Q11[index(z_r,y  ,x_e)] = shared_buffer[threadIdx.x];

}

__global__
void streaming_kernel_first_part_Q12(int lx, int ly, int lz, const FLOATING *Q12,  FLOATING *hlp_Q12){
	/*Propagate fluid densities to their next neighbour nodes */
	/*c....density propagation: all fluid densities are propagated from
	c       non-occupied nodes along the lattice connection lines
	c       to their next neighbours.*/

	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	int end_of_memory=lz*ly*(lx);
	int z=(int) (tid/(ly*lx));
	int	y=(tid-z*(ly*lx))/lx;
	int	x=tid-z*(ly*lx)-y*lx;

	int	z_r/*right*/ =   (z+lz-1) %(  lz) ;
	int x_w/*west*/ =   (x+lx-1) %(  lx) ;

	extern __shared__ FLOATING shared_buffer[];

	shared_buffer[threadIdx.x]=Q12[index(z,y,x)];
	__syncthreads();

	if( tid<end_of_memory)
		hlp_Q12[index(z_r,y  ,x_w)] = shared_buffer[threadIdx.x];

}

__global__
void streaming_kernel_first_part_Q13(int lx, int ly, int lz, const FLOATING *Q13,  FLOATING *hlp_Q13){
	/*Propagate fluid densities to their next neighbour nodes */
	/*c....density propagation: all fluid densities are propagated from
	c       non-occupied nodes along the lattice connection lines
	c       to their next neighbours.*/
	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	int end_of_memory=lz*ly*(lx);
	int z=(int) (tid/(ly*lx));
	int	y=(tid-z*(ly*lx))/lx;
	int	x=tid-z*(ly*lx)-y*lx;


	int z_l/*left*/ =  (z+1)%lz;
	int x_w/*west*/ =   (x+lx-1) %(  lx) ;

	extern __shared__ FLOATING shared_buffer[];

	shared_buffer[threadIdx.x]=Q13[index(z,y,x)];
	__syncthreads();

	if( tid<end_of_memory)
		hlp_Q13[index(z_l,y  ,x_w)] = shared_buffer[threadIdx.x];

}

__global__
void streaming_kernel_first_part_Q14(int lx, int ly, int lz, const FLOATING *Q14,  FLOATING *hlp_Q14){
	/*Propagate fluid densities to their next neighbour nodes */
	/*c....density propagation: all fluid densities are propagated from
	c       non-occupied nodes along the lattice connection lines
	c       to their next neighbours.*/
	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	int end_of_memory=lz*ly*(lx);
	int z=(int) (tid/(ly*lx));
	int	y=(tid-z*(ly*lx))/lx;
	int	x=tid-z*(ly*lx)-y*lx;



	int z_l/*left*/ =  (z+1)%lz;
	int x_e/*east*/ =  (x+1)%lx;

	extern __shared__ FLOATING shared_buffer[];

	shared_buffer[threadIdx.x]=Q14[index(z,y,x)];
	__syncthreads();

	if( tid<end_of_memory)
		hlp_Q14[index(z_l,y  ,x_e)] = shared_buffer[threadIdx.x];

}

__global__
void streaming_kernel_first_part_Q15(int lx, int ly, int lz, const FLOATING *Q15,  FLOATING *hlp_Q15){
	/*Propagate fluid densities to their next neighbour nodes */
	/*c....density propagation: all fluid densities are propagated from
	c       non-occupied nodes along the lattice connection lines
	c       to their next neighbours.*/
	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	int end_of_memory=lz*ly*(lx);
	int z=(int) (tid/(ly*lx));
	int	y=(tid-z*(ly*lx))/lx;
	int	x=tid-z*(ly*lx)-y*lx;

	int z_l/*left*/ =  (z+1)%lz;
	int y_n/*north*/ =  (y+1)%ly;

	extern __shared__ FLOATING shared_buffer[];

	shared_buffer[threadIdx.x]=Q15[index(z,y,x)];
	__syncthreads();

	if( tid<end_of_memory)
		hlp_Q15[index(z_l,y_n,x  )] = shared_buffer[threadIdx.x];

}

__global__
void streaming_kernel_first_part_Q16(int lx, int ly, int lz, const FLOATING *Q16,  FLOATING *hlp_Q16){
	/*Propagate fluid densities to their next neighbour nodes */
	/*c....density propagation: all fluid densities are propagated from
	c       non-occupied nodes along the lattice connection lines
	c       to their next neighbours.*/

	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	int end_of_memory=lz*ly*(lx);
	int z=(int) (tid/(ly*lx));
	int	y=(tid-z*(ly*lx))/lx;
	int	x=tid-z*(ly*lx)-y*lx;

	int z_r /*right*/ =   (z+lz-1) %(  lz) ;
	int y_n /*north*/=  (y+1)%ly;

	extern __shared__ FLOATING shared_buffer[];

	shared_buffer[threadIdx.x]=Q16[index(z,y,x)];
	__syncthreads();

	if( tid<end_of_memory)
		hlp_Q16[index(z_r,y_n,x  )] = shared_buffer[threadIdx.x];

}

__global__
void streaming_kernel_first_part_Q17(int lx, int ly, int lz, const FLOATING *Q17,  FLOATING *hlp_Q17){
	/*Propagate fluid densities to their next neighbour nodes */
	/*c....density propagation: all fluid densities are propagated from
	c       non-occupied nodes along the lattice connection lines
	c       to their next neighbours.*/

	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	int end_of_memory=lz*ly*(lx);
	int z=(int) (tid/(ly*lx));
	int	y=(tid-z*(ly*lx))/lx;
	int	x=tid-z*(ly*lx)-y*lx;

	int	z_r/*right*/ =   (z+lz-1) %(  lz) ;
	int y_s/*south*/ =   (y+ly-1) %(  ly) ;

	extern __shared__ FLOATING shared_buffer[];

	shared_buffer[threadIdx.x]=Q17[index(z,y,x)];
	__syncthreads();

	if( tid<end_of_memory)
		hlp_Q17[index(z_r,y_s,x  )] = shared_buffer[threadIdx.x];

}

__global__
void streaming_kernel_first_part_Q18(int lx, int ly, int lz, const FLOATING *Q18,  FLOATING *hlp_Q18){
	/*Propagate fluid densities to their next neighbour nodes */
	/*c....density propagation: all fluid densities are propagated from
	c       non-occupied nodes along the lattice connection lines
	c       to their next neighbours.*/
	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	int end_of_memory=lz*ly*(lx);
	int z=(int) (tid/(ly*lx));
	int	y=(tid-z*(ly*lx))/lx;
	int	x=tid-z*(ly*lx)-y*lx;

	int z_l/*left*/ =  (z+1)%lz;
	int y_s/*south*/ =   (y+ly-1) %(  ly) ;

	extern __shared__ FLOATING shared_buffer[];

	shared_buffer[threadIdx.x]=Q18[index(z,y,x)];
	__syncthreads();

	if( tid<end_of_memory)
		hlp_Q18[index(z_l,y_s,x  )] = shared_buffer[threadIdx.x];
}



__global__
void streaming_kernel_single_threaded_p2_v2(int lx, int ly, int lz,
		FLOATING *hlp_Q0, FLOATING *hlp_Q1, FLOATING *hlp_Q2, FLOATING *hlp_Q3,
		FLOATING *hlp_Q4, FLOATING *hlp_Q5, FLOATING *hlp_Q6, FLOATING *hlp_Q7,
		FLOATING *hlp_Q8, FLOATING *hlp_Q9, FLOATING *hlp_Q10, FLOATING *hlp_Q11,
		FLOATING *hlp_Q12, FLOATING *hlp_Q13, FLOATING *hlp_Q14, FLOATING *hlp_Q15,
		FLOATING *hlp_Q16, FLOATING *hlp_Q17, FLOATING *hlp_Q18,
		FLOATING *Q0, FLOATING *Q1, FLOATING *Q2, FLOATING *Q3,
		FLOATING *Q4, FLOATING *Q5, FLOATING *Q6, FLOATING *Q7,
		FLOATING *Q8, FLOATING *Q9, FLOATING *Q10, FLOATING *Q11,
		FLOATING *Q12, FLOATING *Q13, FLOATING *Q14, FLOATING *Q15,
		FLOATING *Q16, FLOATING *Q17, FLOATING *Q18){

	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	int z=tid/ly;
	int  y=tid-z*ly;

	if(tid<ly*lz){


		//		for (z = 0 ; z< lz ; ++z){
		//loop for x=0 and x=lx-1!!!! (first and last slice)
		//			for (y = 0; y< ly; ++y){
		//toslice 0 antigrafetai sto 0
		hlp_Q0[index(z,y,0)]=Q0[index(z,y,0)];
		hlp_Q1[index(z,y,0)]=Q1[index(z,y,0)];
		hlp_Q2[index(z,y,0)]=Q2[index(z,y,0)];
		hlp_Q3[index(z,y,0)]=Q3[index(z,y,0)];
		hlp_Q4[index(z,y,0)]=Q4[index(z,y,0)];
		hlp_Q5[index(z,y,0)]=Q5[index(z,y,0)];
		hlp_Q6[index(z,y,0)]=Q6[index(z,y,0)];
		hlp_Q7[index(z,y,0)]=Q7[index(z,y,0)];
		hlp_Q8[index(z,y,0)]=Q8[index(z,y,0)];
		hlp_Q9[index(z,y,0)]=Q9[index(z,y,0)];
		hlp_Q10[index(z,y,0)]=Q10[index(z,y,0)];
		hlp_Q11[index(z,y,0)]=Q11[index(z,y,0)];
		hlp_Q12[index(z,y,0)]=Q12[index(z,y,0)];
		hlp_Q13[index(z,y,0)]=Q13[index(z,y,0)];
		hlp_Q14[index(z,y,0)]=Q14[index(z,y,0)];
		hlp_Q15[index(z,y,0)]=Q15[index(z,y,0)];
		hlp_Q16[index(z,y,0)]=Q16[index(z,y,0)];
		hlp_Q17[index(z,y,0)]=Q17[index(z,y,0)];
		hlp_Q18[index(z,y,0)]=Q18[index(z,y,0)];
		//at x= lx I set the incomming density as the one (equilibrum) calculated after the collision
		hlp_Q3[index(z,y,lx-1)] = Q3[index(z,y,lx-1)];
		hlp_Q8[index(z,y,lx-1)] = Q8[index(z,y,lx-1)] ;
		hlp_Q9[index(z,y,lx-1)] = Q9[index(z,y,lx-1)] ;
		hlp_Q12[index(z,y,lx-1)] = Q12[index(z,y,lx-1)] ;
		hlp_Q13[index(z,y,lx-1)] = Q13[index(z,y,lx-1)] ;
		//			gia x=lx-1: the densities 0,   2,4,5,6,15,16,17,18
		//				aplws 8a metadw8oune
		//				"ka8eta" sto slice kai de 8a ginoun propagate se
		//				alla slices (tuxainei na voleuei auto)
		//				auto to kommati exei HDH ginei pio panw!
		//
		//
		//				during streaming some of the indices at
		//				(fixed) x=lx-1, y=ly-1, do not get updated.
		//				it doesn't matter since the last slice on x is not useful.


	}
}

void LBM::cuda_streaming(){

	if(data_location==CPU)
		copy_data_from_host_to_device();

	dim3 threads_type2(threads_for_streaming_collision_and_relaxation,1,1);
	dim3 grid_type2(blocks_for_streaming_collision_and_relaxation,1,1);

	streaming_kernel_first_part_Q0<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>( lx,  ly,  lz, D3_d.Q0,  D3_hlp_d.Q0);
	streaming_kernel_first_part_Q1<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>( lx,  ly,  lz, D3_d.Q1,  D3_hlp_d.Q1);
	streaming_kernel_first_part_Q2<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>( lx,  ly,  lz, D3_d.Q2,  D3_hlp_d.Q2);
	streaming_kernel_first_part_Q3<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>( lx,  ly,  lz, D3_d.Q3,  D3_hlp_d.Q3);
	streaming_kernel_first_part_Q4<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>( lx,  ly,  lz, D3_d.Q4,  D3_hlp_d.Q4);
	streaming_kernel_first_part_Q5<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>( lx,  ly,  lz, D3_d.Q5,  D3_hlp_d.Q5);
	streaming_kernel_first_part_Q6<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>( lx,  ly,  lz, D3_d.Q6,  D3_hlp_d.Q6);
	streaming_kernel_first_part_Q7<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>( lx,  ly,  lz, D3_d.Q7,  D3_hlp_d.Q7);
	streaming_kernel_first_part_Q8<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>( lx,  ly,  lz, D3_d.Q8,  D3_hlp_d.Q8);
	streaming_kernel_first_part_Q9<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>( lx,  ly,  lz, D3_d.Q9,  D3_hlp_d.Q9);
	streaming_kernel_first_part_Q10<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>( lx,  ly,  lz, D3_d.Q10,  D3_hlp_d.Q10);
	streaming_kernel_first_part_Q11<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>( lx,  ly,  lz, D3_d.Q11,  D3_hlp_d.Q11);
	streaming_kernel_first_part_Q12<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>( lx,  ly,  lz, D3_d.Q12,  D3_hlp_d.Q12);
	streaming_kernel_first_part_Q13<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>( lx,  ly,  lz, D3_d.Q13,  D3_hlp_d.Q13);
	streaming_kernel_first_part_Q14<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>( lx,  ly,  lz, D3_d.Q14,  D3_hlp_d.Q14);
	streaming_kernel_first_part_Q15<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>( lx,  ly,  lz, D3_d.Q15,  D3_hlp_d.Q15);
	streaming_kernel_first_part_Q16<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>( lx,  ly,  lz, D3_d.Q16,  D3_hlp_d.Q16);
	streaming_kernel_first_part_Q17<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>( lx,  ly,  lz, D3_d.Q17,  D3_hlp_d.Q17);
	streaming_kernel_first_part_Q18<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>( lx,  ly,  lz, D3_d.Q18,  D3_hlp_d.Q18);
	cudaDeviceSynchronize();

	if(data_location==CPU)
		copy_data_from_host_to_device();

	int n_of_threads=64;
	int n_of_blocks=(lz*ly)/n_of_threads;
	if ( (lattice_nodes%n_of_threads)!=0 )
		++n_of_blocks;

	dim3 threads_type3(n_of_threads,1,1);
	dim3 grid_type3(n_of_blocks,1,1);

	streaming_kernel_single_threaded_p2_v2<<<grid_type3, threads_type3>>>(lx,ly,lz,
			D3_hlp_d.Q0, D3_hlp_d.Q1, D3_hlp_d.Q2, D3_hlp_d.Q3,
			D3_hlp_d.Q4, D3_hlp_d.Q5, D3_hlp_d.Q6, D3_hlp_d.Q7,
			D3_hlp_d.Q8, D3_hlp_d.Q9, D3_hlp_d.Q10, D3_hlp_d.Q11,
			D3_hlp_d.Q12, D3_hlp_d.Q13, D3_hlp_d.Q14, D3_hlp_d.Q15,
			D3_hlp_d.Q16, D3_hlp_d.Q17, D3_hlp_d.Q18,
			D3_d.Q0, D3_d.Q1, D3_d.Q2, D3_d.Q3,
			D3_d.Q4, D3_d.Q5, D3_d.Q6, D3_d.Q7,
			D3_d.Q8, D3_d.Q9, D3_d.Q10, D3_d.Q11,
			D3_d.Q12, D3_d.Q13, D3_d.Q14, D3_d.Q15,
			D3_d.Q16, D3_d.Q17, D3_d.Q18);

	cudaDeviceSynchronize();
}

