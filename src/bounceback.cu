#include "global_defines.cuh"


void LBM::bounceback(){
	/*Fluid densities are rotated. By the next propagation step, this  *
	 *     results in a bounce back from obstacle nodes.*/

	/*
			.......bounce back from obstacles: this is the no-slip boundary-
			condition.
			The velocity vector of all fluid densities is inverted, so all
			the fluid densities will be sent back to the node  where they
			were located before the last propagation step, but with opposite
			velocity vector
			... there exist lots of other possibilities.
		 */

	if(data_location==GPU)
		copy_data_from_device_to_host();

	int  x,y,z;

	//.....loop over all nodes
	for( z = 0; z< lz ; ++z){
		for( y = 0; y< ly ; ++y){
			for( x = 0; x< lx ; ++x){

				//.........consider only obstacle nodes
				if (obstacles[index(z,y,x)]==1){

					//...........rotate all ensities and write back to node
					D3.Q1[index(z,y,x)] = D3_hlp.Q3[index(z,y,x)];
					D3.Q2[index(z,y,x)] = D3_hlp.Q4[index(z,y,x)];
					D3.Q3[index(z,y,x)] = D3_hlp.Q1[index(z,y,x)];
					D3.Q4[index(z,y,x)] = D3_hlp.Q2[index(z,y,x)];
					D3.Q5[index(z,y,x)] = D3_hlp.Q6[index(z,y,x)];
					D3.Q6[index(z,y,x)] = D3_hlp.Q5[index(z,y,x)];
					D3.Q7[index(z,y,x)] = D3_hlp.Q9[index(z,y,x)];
					D3.Q8[index(z,y,x)] = D3_hlp.Q10[index(z,y,x)];
					D3.Q9[index(z,y,x)] = D3_hlp.Q7[index(z,y,x)];
					D3.Q10[index(z,y,x)] = D3_hlp.Q8[index(z,y,x)];
					D3.Q11[index(z,y,x)] = D3_hlp.Q13[index(z,y,x)];
					D3.Q12[index(z,y,x)] = D3_hlp.Q14[index(z,y,x)];
					D3.Q13[index(z,y,x)] = D3_hlp.Q11[index(z,y,x)];
					D3.Q14[index(z,y,x)] = D3_hlp.Q12[index(z,y,x)];
					D3.Q15[index(z,y,x)] = D3_hlp.Q17[index(z,y,x)];
					D3.Q16[index(z,y,x)] = D3_hlp.Q18[index(z,y,x)];
					D3.Q17[index(z,y,x)] = D3_hlp.Q15[index(z,y,x)];
					D3.Q18[index(z,y,x)] = D3_hlp.Q16[index(z,y,x)];
				}
			}
		}
	}
#ifdef DEBUG
	cout << " #LBM bounceback OK!" << endl;
#endif
}









__global__
void bounceback_kernel_v4_shared(const int end_of_memory, const CUDA_FLOATING *source_data, CUDA_FLOATING *destination_data,
		const int *obstacles){
	/*Fluid densities are rotated. By the next propagation step, this  *
	 *     results in a bounce back from obstacle nodes.*/

	/*
		.......bounce back from obstacles: this is the no-slip boundary-
		condition.
		The velocity vector of all fluid densities is inverted, so all
		the fluid densities will be sent back to the node  where they
		were located before the last propagation step, but with opposite
		velocity vector
		... there exist lots of other possibilities.
	 */

	const int tid=blockIdx.x*blockDim.x+threadIdx.x;
	extern __shared__ CUDA_FLOATING shared_buffer[];
	shared_buffer[threadIdx.x]=source_data[tid];

	__syncthreads();
	if (tid<end_of_memory and obstacles[tid]){
		destination_data[tid]=shared_buffer[threadIdx.x];
	}
}



void LBM::cuda_bounceback(){

	if(data_location==CPU)
		copy_data_from_host_to_device();

	dim3 threads_type2(threads_for_streaming_collision_and_relaxation,1,1);
	dim3 grid_type2(blocks_for_streaming_collision_and_relaxation,1,1);



	bounceback_kernel_v4_shared<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>(lattice_nodes,
			D3_hlp_d.Q3, D3_d.Q1,   obstacles_d);
	bounceback_kernel_v4_shared<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>(lattice_nodes,
			D3_hlp_d.Q4, D3_d.Q2,   obstacles_d);
	bounceback_kernel_v4_shared<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>(lattice_nodes,
			D3_hlp_d.Q1, D3_d.Q3,   obstacles_d);
	bounceback_kernel_v4_shared<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>(lattice_nodes,
			D3_hlp_d.Q2, D3_d.Q4,   obstacles_d);
	bounceback_kernel_v4_shared<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>(lattice_nodes,
			D3_hlp_d.Q6, D3_d.Q5,   obstacles_d);
	bounceback_kernel_v4_shared<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>(lattice_nodes,
			D3_hlp_d.Q5, D3_d.Q6,   obstacles_d);
	bounceback_kernel_v4_shared<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>(lattice_nodes,
			D3_hlp_d.Q9, D3_d.Q7,   obstacles_d);
	bounceback_kernel_v4_shared<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>(lattice_nodes,
			D3_hlp_d.Q10, D3_d.Q8,   obstacles_d);
	bounceback_kernel_v4_shared<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>(lattice_nodes,
			D3_hlp_d.Q7, D3_d.Q9,   obstacles_d);
	bounceback_kernel_v4_shared<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>(lattice_nodes,
			D3_hlp_d.Q8, D3_d.Q10,   obstacles_d);
	bounceback_kernel_v4_shared<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>(lattice_nodes,
			D3_hlp_d.Q13, D3_d.Q11,   obstacles_d);
	bounceback_kernel_v4_shared<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>(lattice_nodes,
			D3_hlp_d.Q14, D3_d.Q12,   obstacles_d);
	bounceback_kernel_v4_shared<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>(lattice_nodes,
			D3_hlp_d.Q11, D3_d.Q13,   obstacles_d);
	bounceback_kernel_v4_shared<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>(lattice_nodes,
			D3_hlp_d.Q12, D3_d.Q14,   obstacles_d);
	bounceback_kernel_v4_shared<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>(lattice_nodes,
			D3_hlp_d.Q17, D3_d.Q15,   obstacles_d);
	bounceback_kernel_v4_shared<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>(lattice_nodes,
			D3_hlp_d.Q18, D3_d.Q16,   obstacles_d);
	bounceback_kernel_v4_shared<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>(lattice_nodes,
			D3_hlp_d.Q15, D3_d.Q17,   obstacles_d);
	bounceback_kernel_v4_shared<<<grid_type2, threads_type2, size_of_allocated_shared_memory_for_streaming_collision_and_relaxation>>>(lattice_nodes,
			D3_hlp_d.Q16, D3_d.Q18,   obstacles_d);

	cudaDeviceSynchronize();
#ifdef DEBUG
cout << " #LBM bounceback OK!" << endl;
#endif
}

