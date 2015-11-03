#include "global_defines.cuh"

void init_density( const int &lx, const int &ly, const int &lz, FLOATING *node_grid, FLOATING *node_hlp_grid){

	const FLOATING t_0 =  1.0  / 3.0  ;
		const FLOATING t_1 =  1.0 /  18.0 ;
		const FLOATING t_2 =  1.0 / 36.0 ;

	//Initialize density distribution function n with equilibrium	 for zero velocity                                                *
	int  x,y,z,i;
	//.....compute weighting factors (depending on lattice geometry)

	//.....loop over computational domain
	for (z = 0 ; z< lz ; ++z){
		for (y = 0 ; y< ly ; ++y){
			for (x = 0 ; x< lx ; ++x){

				//.........zero velocity density

				node_grid[index4D(z,y,x,0)] = t_0;
				node_hlp_grid[index4D(z,y,x,0)] = t_0;

				//.........equilibrium densities for axis speeds
#pragma unroll
				for( i=1; i<7 ; ++i){
					node_grid[index4D(z,y,x,i)] = t_1;
					node_hlp_grid[index4D(z,y,x,i)] = t_1;
				}

				//.........equilibrium densities for diagonal speeds
#pragma unroll
				for( i=7; i<19 ; ++i){
					node_grid[index4D(z,y,x,i)] = t_2;
					node_hlp_grid[index4D(z,y,x,i)] = t_2;
				}

				//todo, replaces the following me to memcopy!!! h alliws memcpy!!!!

			}
		}
		/*
	// Use this only when reading a 2d file to initiate 3d simulation.
	if (lzd ==1 ){
		for ( x=0 ; x<lx ; ++x)
			for ( y=0 ; y<ly ; ++y)
				for ( z=1 ; z<lz ; ++z)
					for( i=0 ; i<19 ; ++i)
						node[x][y][z][i]= (FLOATING) node[x][y][lzd][i];
	}
		 */
	}
}

void LBM::initialise_microscopic_density_arrays_in_the_host(){
	//Initialize density distribution function n with equilibrium	 for zero velocity                                                *
	int  x,y,z;
	//.....compute weighting factors (depending on lattice geometry)

	//.....loop over computational domain
	for (z = 0 ; z< lz ; ++z){
		for (y = 0 ; y< ly ; ++y){
			for (x = 0 ; x< lx ; ++x){

				//zero velocity density
				D3.Q0[index(z,y,x)]=t_0;
				D3_hlp.Q0[index(z,y,x)]=t_0;


				//equilibrium densities for axis speeds
				D3.Q1[index(z,y,x)]=t_1;
				D3.Q2[index(z,y,x)]=t_1;
				D3.Q3[index(z,y,x)]=t_1;
				D3.Q4[index(z,y,x)]=t_1;
				D3.Q5[index(z,y,x)]=t_1;
				D3.Q6[index(z,y,x)]=t_1;

				D3_hlp.Q1[index(z,y,x)]=t_1;
				D3_hlp.Q2[index(z,y,x)]=t_1;
				D3_hlp.Q3[index(z,y,x)]=t_1;
				D3_hlp.Q4[index(z,y,x)]=t_1;
				D3_hlp.Q5[index(z,y,x)]=t_1;
				D3_hlp.Q6[index(z,y,x)]=t_1;


				//equilibrium densities for diagonal speeds
				D3.Q7[index(z,y,x)]=t_2;
				D3.Q8[index(z,y,x)]=t_2;
				D3.Q9[index(z,y,x)]=t_2;
				D3.Q10[index(z,y,x)]=t_2;
				D3.Q11[index(z,y,x)]=t_2;
				D3.Q12[index(z,y,x)]=t_2;
				D3.Q13[index(z,y,x)]=t_2;
				D3.Q14[index(z,y,x)]=t_2;
				D3.Q15[index(z,y,x)]=t_2;
				D3.Q16[index(z,y,x)]=t_2;
				D3.Q17[index(z,y,x)]=t_2;
				D3.Q18[index(z,y,x)]=t_2;

				D3_hlp.Q7[index(z,y,x)]=t_2;
				D3_hlp.Q8[index(z,y,x)]=t_2;
				D3_hlp.Q9[index(z,y,x)]=t_2;
				D3_hlp.Q10[index(z,y,x)]=t_2;
				D3_hlp.Q11[index(z,y,x)]=t_2;
				D3_hlp.Q12[index(z,y,x)]=t_2;
				D3_hlp.Q13[index(z,y,x)]=t_2;
				D3_hlp.Q14[index(z,y,x)]=t_2;
				D3_hlp.Q15[index(z,y,x)]=t_2;
				D3_hlp.Q16[index(z,y,x)]=t_2;
				D3_hlp.Q17[index(z,y,x)]=t_2;
				D3_hlp.Q18[index(z,y,x)]=t_2;
			}
		}
		/*
	// Use this only when reading a 2d file to initiate 3d simulation.
	if (lzd ==1 ){
		for ( x=0 ; x<lx ; ++x)
			for ( y=0 ; y<ly ; ++y)
				for ( z=1 ; z<lz ; ++z)
					for( i=0 ; i<19 ; ++i)
						node[x][y][z][i]= (FLOATING) node[x][y][lzd][i];
	}
		 */
	}

#ifdef DEBUG
	cout << " LBM init density OK!" << endl;
#endif

}
