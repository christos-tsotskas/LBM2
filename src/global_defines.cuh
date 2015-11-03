/*
 * global_defines.h
 *
 *  Created on: Jun 21, 2012
 *      Author: christos
 */

#ifndef GLOBAL_DEFINES_H_
#define GLOBAL_DEFINES_H_

#include <iostream>
#include <cstdlib>
//#include "mpi.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <cstring>
#include <sys/resource.h> //sets stack size
#include <set>
#include <vector>


//#define DEBUG
//#define REPORT
#define MANY_THREADS 256

#define USING_DOUBLE
//#define USING_FLOAT

#define FLOATING float
#define CUDA_FLOATING float

//#define CPU_part
#define GPU_part
#define PRODUCE_OUTPUT_FILES //this should be deactivated when LBM is on interface mode

#define CPU 0
#define GPU 1


#define DENSITIES 19
#define XBAFFLE 59




#define index(z,y,x) (z)*ly*lx+(y)*lx+(x)

#define index2D(z,y) (z)*ly+(y)
#define index4D(z,y,x,i) (z)*ly*lx*DENSITIES+(y)*lx*DENSITIES+(x)*DENSITIES+(i)

#define SUM
using namespace std;

//extern const FLOATING t_0;
//extern const FLOATING t_1;
//extern const FLOATING t_2;
//extern const FLOATING c_squ;
//debug
void temp_compare(FLOATING *a, FLOATING *b);

void cuda_device_querry();
inline int _ConvertSMVer2Cores(int major, int minor);
void read_external_geometry_file_specification_for_LBM(int &lx, int &ly, int &lz, int &n_of_densities, const string filename);

class lattice{
	const int lx, ly, lz;
public:
	CUDA_FLOATING *Q0, *Q1, *Q2, *Q3, *Q4, *Q5, *Q6, *Q7, *Q8, *Q9;
	CUDA_FLOATING *Q10, *Q11, *Q12, *Q13, *Q14, *Q15, *Q16, *Q17, *Q18;
	lattice(int LX,int LY, int LZ);
	lattice(int LX,int LY, int LZ, int dump);
	void initialise(FLOATING *Q);
	~lattice();
};

class LBM{
public:
	time_t time_start,time_end;
	FLOATING time_elapsed;
	int max_iterations, check_step;
	const int lx, ly, lz;
	const int lattice_nodes;
	int no_obstacle_lattices_at_penultimate_x_slice;
	int threads_for_streaming_collision_and_relaxation;
	int blocks_for_streaming_collision_and_relaxation;
	int size_of_allocated_shared_memory_for_streaming_collision_and_relaxation;
	int convective_boundary_conditions_blocks;

private:
	FLOATING nu, r_small, reynolds, s, density;
	const FLOATING t_0, t_1, t_2, c_squ, reciprocal_c_squ;
	int baffle;//on X
	int threads_per_kernel;
	int time_unit;
	int two_dimensional_length;
	int three_dimensional_length;
	int const FLOATING_slice_size;
	int const int_array_size;
	string case_name;
	const FLOATING tau, omega, one_minus_omega;
	FLOATING pr_diff, pr_out, pr_in, vor;


	//  tau=3.0*nu + 0.5, omega = 1.0 /tau, ; //	omega=1.0/(3.0*nu+0.5);

	lattice D3, D3_hlp;
	int *obstacles;

	FLOATING *u_current;
	FLOATING *v_current;
	FLOATING *w_current;
	// u_previous_spatial_boundary: at boundary - 1 (in x)
	FLOATING *u_previous_spatial_boundary;
	FLOATING *v_previous_spatial_boundary;
	FLOATING *w_previous_spatial_boundary;
	// u_prev: at boundary - 1 (in time)
	FLOATING *u_previous_temporal_boundary;
	FLOATING *v_previous_temporal_boundary;
	FLOATING *w_previous_temporal_boundary;

	FLOATING *u_current_d;
	FLOATING *u_current_temp_d;
	FLOATING *v_current_d;
	FLOATING *w_current_d;
	// u_previous_spatial_boundary: at boundary - 1 (in x)
	FLOATING *u_previous_spatial_boundary_d;
	FLOATING *v_previous_spatial_boundary_d;
	FLOATING *w_previous_spatial_boundary_d;
	// u_prev: at boundary - 1 (in time)
	FLOATING *u_previous_temporal_boundary_d;
	FLOATING *v_previous_temporal_boundary_d;
	FLOATING *w_previous_temporal_boundary_d;


	FLOATING *temp_cpu_u_current_d;
	FLOATING *temp_cpu_v_current_d;
	FLOATING *temp_cpu_w_current_d;

	FLOATING *temp_cpu_u_previous_temporal_boundary_d;
	FLOATING *temp_cpu_v_previous_temporal_boundary_d;
	FLOATING *temp_cpu_w_previous_temporal_boundary_d;

	FLOATING *temp_cpu_u_previous_spatial_boundary_d;
	FLOATING *temp_cpu_v_previous_spatial_boundary_d;
	FLOATING *temp_cpu_w_previous_spatial_boundary_d;

	FLOATING *temp_check_density_d;
	FLOATING *temp_check_density_d_full;

	int data_location;//0-CPU, 1-GPU
	FLOATING *temp_Uc_d;
	int *obstacles_d;

	lattice D3_d, D3_hlp_d;

	//host memories - for saving the final results
	FLOATING *Ux, *Uy, *Uz, *Pressure, *Wx, *Wy, *Wz;

	void create_an_example_configuration_files(const string filename);

	void display_the_structure_of_an_example_configuration_file();

	void read_external_configuration_file_for_the_solver(const string filename);
	void reset_convergence_file();
	void calculate_CUDA_quantities();
	void display_CUDA_specifications();
	void display_LBM_specifications();

	void allocate_device_arrays();
	void initialise_host_data();
	void initialise_device_data();
	void initialise_all_data_arrays();


	void abstract_check_density();

	void abstract_debug_computations();
	void abstract_redistribute();
	void abstract_streaming();
	void abstract_bounce_back();
	void abstract_relaxation();
	void abstract_convective_boundary_conditions();
	void core_computations();


	void delete_device_data();
	void delete_host_memories();


	void abstract_initialise();


	template <class T>
	void initialise_array(T *array, const int length,const T init_value);
	template <typename T>
	void allocate_and_initialise(T *array, const int length);

	void count_no_obstacles_at_penultimate_x_slice();
	void copy_data_from_host_to_device();
	void small_copy_data_from_host_to_device();
	void copy_data_from_device_to_host();
	void small_copy_data_from_device_to_host();

	void cuda_func(int n_of_blocks, int n_of_threads, int a);

	FLOATING calulate_partial_sum(const FLOATING *input_array_d);
	FLOATING reduce_a_fraction_of_the_domain(const CUDA_FLOATING *domain_of_interest_d);
	void cuda_check_density(const int iteration);
	void cuda_redistribute();
	void cuda_streaming();
	void cuda_bounceback();
	void cuda_relaxation();
	void cuda_convective_BC();





	void read_reactor();
	void create_reactor_geometry_in_the_host();
	void initialise_microscopic_density_arrays_in_the_host();
	void initial_redistribute2();
	void initial_redistribute();
	FLOATING sum_microscopic_densities_for_a_single_node(const int x, const int y, const int z);
	int calculate_macroscopic_density_in_the_host(FLOATING &density);
	void redistribute();
	void fortran_redistribute(int time);
	void streaming();
	void streaming_last_part();
	void streaming_first_part();
	void streaming_A();
	void streaming_B();
	void streaming_C();
	void bounceback();
	void initial_relaxation();
	void relaxation();
	void convective_BC();

	FLOATING convective_BC_A();
	void convective_BC_B(FLOATING Uc);

	FLOATING calculate_vorticity();
	FLOATING calculate_pressure_differences();

	void write_log_file(const int iteration);
	void write_objectives();
	void calculate_macroscopic_quantities(const int& iteration);

	void write_VTK_SI(const int& iteration);


	void geometry_file_in_VTK();


	void write_output(const int &time_step, const FLOATING &vor, const FLOATING &pr_diff, const FLOATING &time_elapsed);
	//Debug
	void compare_obstacles(int *outter_obst);
	void compare_nodes(FLOATING *outter_node);
	void compare_nodes_hlp(FLOATING *outter_node_hlp);

public:
	LBM(const int &LX, const int &LY, const int &LZ, const FLOATING &DENSITY, const FLOATING &T_0,
			const FLOATING &T_1, const FLOATING &T_2, const FLOATING &C_SQU);
	~LBM();
	void compute_domain();
	void export_solution();


};

#endif /* GLOBAL_DEFINES_H_ */
