#include "global_defines.cuh"
#include "kernels3.cuh"

FLOATING LBM::sum_microscopic_densities_for_a_single_node(const int x, const int y, const int z){
	FLOATING sum = 0.0;
	sum = sum + D3.Q0[index(z,y,x)];

	sum+=D3.Q1[index(z,y,x)];
	sum+=D3.Q2[index(z,y,x)];
	sum+=D3.Q3[index(z,y,x)];
	sum+=D3.Q4[index(z,y,x)];
	sum+=D3.Q5[index(z,y,x)];
	sum+=D3.Q6[index(z,y,x)];

	sum+=D3.Q7[index(z,y,x)];
	sum+=D3.Q8[index(z,y,x)];
	sum+=D3.Q9[index(z,y,x)];
	sum+=D3.Q10[index(z,y,x)];
	sum+=D3.Q11[index(z,y,x)];
	sum+=D3.Q12[index(z,y,x)];
	sum+=D3.Q13[index(z,y,x)];
	sum+=D3.Q14[index(z,y,x)];
	sum+=D3.Q15[index(z,y,x)];
	sum+=D3.Q16[index(z,y,x)];
	sum+=D3.Q17[index(z,y,x)];
	sum+=D3.Q18[index(z,y,x)];

	return sum;
}




int LBM::calculate_macroscopic_density_in_the_host(FLOATING &density){
	/*compute integral density*/

	if(data_location==GPU)
		copy_data_from_device_to_host();
	//.....local variables
	int  x,y,z;
	FLOATING n_sum=0.0;

	//.....loop over computational domain
	//...........loop over all densities

	//loop only for the last quarter of the domain
#pragma unroll
	for ( z =0 ; z < lz ; ++z){
#pragma unroll
		for (y = 0; y < ly ; ++y){
#pragma unroll
			for ( x = (lx*3/4)-1 ; x< lx ; ++x){

				n_sum+=sum_microscopic_densities_for_a_single_node(x,y,z);

			}
		}
	}





	cout.setf(ios::fixed,ios::floatfield);
	cout.precision(10);
	cout << "density check: Integral density=" << n_sum << " " << n_sum/(lz*ly*(lx/4))  << endl;
	//check for nan

	if( isnan(n_sum/(lx*lz*ly))==0 ){
		//it is NOT nan, it is a good number!
		density=n_sum/(lz*ly*(lx/4));
		return 1;
	}else{
		//NAN WAS FOUND:the rest of this branch is for debugging purposes (print the nan on an external file)
		cout << " nsum:" << n_sum << endl;
#pragma unroll
		for ( z =0 ; z < lz ; ++z){
#pragma unroll
			for (y = 0; y < ly ; ++y){
#pragma unroll
				for ( x = (lx*3/4)-1 ; x< lx ; ++x){
					n_sum+=sum_microscopic_densities_for_a_single_node(x,y,z);

					if(isnan(n_sum)==0)
						cout << " x,y,z: " << x << " " << y << " " << z << endl;
				}
			}
		}

		return 0;

	}

}


FLOATING sum_one_density(int big_array_length, int numThreads, int numBlocks, int whichKernel, FLOATING *D3_Q,FLOATING* h_odata, FLOATING* d_odata){
	FLOATING partial_sum=0.0;

	reduce<FLOATING>(big_array_length, numThreads, numBlocks, whichKernel, D3_Q, d_odata);
	cudaDeviceSynchronize();
	cudaMemcpy( h_odata, d_odata, numBlocks*sizeof(FLOATING), cudaMemcpyDeviceToHost);
	for(int i=0; i<numBlocks; i++)
		partial_sum+= h_odata[i];

	return partial_sum;
}


__global__
void simple_check_density(const FLOATING *input_array, const int lx, const int ly, const int lz, FLOATING *output_array){
	int tid=blockIdx.x*blockDim.x+threadIdx.x;

	int z=(int) (tid/(ly*lx));
	int	y=(int) (tid-z*(ly*lx))/lx;
	int	x=(int) tid-z*(ly*lx)-y*lx;

	if( (3/4*lx-1)<=x and x< lx)
		output_array[tid]=input_array[tid];
}

FLOATING LBM::calulate_partial_sum(const FLOATING *input_array_d){
	FLOATING partial_sum=0.0;
	//calculate cuda threads and blocks
	//	simple_check_density<<< 1,1 >>> (input_array_d, 680, 73, 73, temp_check_density_d_full);
	partial_sum=reduce_sum(temp_check_density_d_full, 680*73*73);
	return partial_sum;
}


__global__
void collect_data_convergence_interest(const CUDA_FLOATING *input_array, const int lx, const int ly, const int lz,
		FLOATING *output_array){

	//collects a range of elements from input_array and inserts them into the output
	//then the reduce will be applied upon these data
	int tid=blockIdx.x*blockDim.x+threadIdx.x;

	int z_short= (int)(4*tid)/(lx*ly);
	int y_short=(int) (4*tid-z_short*lx*ly) / lx ;
	int x_short=((4*tid+3)%lx);

	extern __shared__ FLOATING shared_buffer[];

	shared_buffer[threadIdx.x]=input_array[index(z_short,y_short,x_short)];
	__syncthreads();

	output_array[tid]=shared_buffer[threadIdx.x];

	//	output_array[tid]=input_array[index(z_short,y_short,x_short)];



	__syncthreads();

}

FLOATING LBM::reduce_a_fraction_of_the_domain(const CUDA_FLOATING *domain_of_interest_d){

	int lattice_nodes=lz*ly*lx/4;
	int n_of_threads=threads_per_kernel;
	int n_of_blocks= (lattice_nodes)/n_of_threads;

	if( (lattice_nodes%n_of_threads)!=0)
		++n_of_blocks;

	const int size_of_allocated_shared_memory=(n_of_threads)*sizeof(CUDA_FLOATING);

	dim3 threads_type2(n_of_threads,1,1);
	dim3 grid_type2(n_of_blocks,1,1);

	collect_data_convergence_interest<<<grid_type2, threads_type2,size_of_allocated_shared_memory>>>(domain_of_interest_d,
			lx, ly,  lz, temp_check_density_d);

	FLOATING partial_sum=reduce_sum(temp_check_density_d, lz*ly*lx/4);
	return partial_sum;


}

void LBM::cuda_check_density(const int iteration){

//	int const array_length=(lx*ly*lz);
//	int const FLOATING_array_size=array_length*sizeof(FLOATING);



	FLOATING n_sum=0.0;


	n_sum+=reduce_a_fraction_of_the_domain(D3_d.Q0);
	n_sum+=reduce_a_fraction_of_the_domain(D3_d.Q1);
	n_sum+=reduce_a_fraction_of_the_domain(D3_d.Q2);
	n_sum+=reduce_a_fraction_of_the_domain(D3_d.Q3);
	n_sum+=reduce_a_fraction_of_the_domain(D3_d.Q4);
	n_sum+=reduce_a_fraction_of_the_domain(D3_d.Q5);
	n_sum+=reduce_a_fraction_of_the_domain(D3_d.Q6);
	n_sum+=reduce_a_fraction_of_the_domain(D3_d.Q7);
	n_sum+=reduce_a_fraction_of_the_domain(D3_d.Q8);
	n_sum+=reduce_a_fraction_of_the_domain(D3_d.Q9);
	n_sum+=reduce_a_fraction_of_the_domain(D3_d.Q10);
	n_sum+=reduce_a_fraction_of_the_domain(D3_d.Q11);
	n_sum+=reduce_a_fraction_of_the_domain(D3_d.Q12);
	n_sum+=reduce_a_fraction_of_the_domain(D3_d.Q13);
	n_sum+=reduce_a_fraction_of_the_domain(D3_d.Q14);
	n_sum+=reduce_a_fraction_of_the_domain(D3_d.Q15);
	n_sum+=reduce_a_fraction_of_the_domain(D3_d.Q16);
	n_sum+=reduce_a_fraction_of_the_domain(D3_d.Q17);
	n_sum+=reduce_a_fraction_of_the_domain(D3_d.Q18);
	//	n_sum=1.1;

	FLOATING temp_cpu_density=0.0;
	//CPU WORKING
//	if(data_location==GPU)
//		copy_data_from_device_to_host();
//	if( check_density(temp_cpu_density)==0)
//		exit(10);

	//save results in file - start
	const int quarter_length=lx*ly*lz/4;

	char buffer1[256];
	snprintf(buffer1, sizeof(buffer1), "LBM2_%s_convergence.log", case_name.c_str());

#ifdef PRODUCE_OUTPUT_FILES
	ofstream convergence_file;
	convergence_file.open( buffer1);
	//convergence_file.open( buffer1 , ofstream::app);
	convergence_file.precision(10);
#endif //PRODUCE_OUTPUT_FILES

	if (n_sum!=0){
		cout <<" cuda Integral density:" << n_sum/quarter_length << endl;
#ifdef PRODUCE_OUTPUT_FILES
		convergence_file<< iteration << "\t" << n_sum/quarter_length << "\t" << temp_cpu_density << endl;
#endif //PRODUCE_OUTPUT_FILES
	}else{
		cout <<" n_sum=0! potentially wrong! check code!" << endl;
#ifdef PRODUCE_OUTPUT_FILES
		convergence_file<< iteration << "\t" << n_sum/quarter_length << " WARNING! "<< endl;
#endif //PRODUCE_OUTPUT_FILES
	}
#ifdef PRODUCE_OUTPUT_FILES
	convergence_file.close();
#endif //PRODUCE_OUTPUT_FILES
	//save results in file - end

}
