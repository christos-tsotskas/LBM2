#include "global_defines.cuh"
#include <cstdio>
#include <cstdlib>


void temp_compare(FLOATING *a, FLOATING *b){
	int x,y,z;
	int missed=0;
	int lx=680, ly=73, lz=73;

	for (z = 0 ; z< lz ; ++z){
		for (y = 0 ; y< ly ; ++y){
			for (x = 0 ; x< lx; ++x){
				if(abs(a[index(z,y,x)]-b[index(z,y,x)])>0.00001)
					++missed;

			}
		}
	}
	cout <<"totally missed:" << missed << endl;

	cout << " FULL MATCH!" << endl;
}

// Beginning of GPU Architecture definitions
inline int _ConvertSMVer2Cores(int major, int minor)
{
	// Defines for GPU Architecture types (using the SM version to determine the # of cores per SM
	typedef struct
	{
		int SM; // 0xMm (hexidecimal notation), M = SM Major version, and m = SM minor version
		int Cores;
	} sSMtoCores;

	sSMtoCores nGpuArchCoresPerSM[] =
	{
			{ 0x10,  8 }, // Tesla Generation (SM 1.0) G80 class
			{ 0x11,  8 }, // Tesla Generation (SM 1.1) G8x class
			{ 0x12,  8 }, // Tesla Generation (SM 1.2) G9x class
			{ 0x13,  8 }, // Tesla Generation (SM 1.3) GT200 class
			{ 0x20, 32 }, // Fermi Generation (SM 2.0) GF100 class
			{ 0x21, 48 }, // Fermi Generation (SM 2.1) GF10x class
			{ 0x30, 192}, // Kepler Generation (SM 3.0) GK10x class
			{ 0x35, 192}, // Kepler Generation (SM 3.5) GK11x class
			{   -1, -1 }
	};

	int index = 0;

	while (nGpuArchCoresPerSM[index].SM != -1)
	{
		if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor))
		{
			return nGpuArchCoresPerSM[index].Cores;
		}

		index++;
	}

	// If we don't find the values, we default use the previous one to run properly
	printf("MapSMtoCores for SM %d.%d is undefined.  Default to use %d Cores/SM\n", major, minor, nGpuArchCoresPerSM[7].Cores);
	return nGpuArchCoresPerSM[7].Cores;
}
// end of GPU Architecture definitions


void cuda_device_querry(){
	cout << "CUDA DEVISE TEST - START" << endl;

	int deviceCount = 0;
	cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

	if (error_id != cudaSuccess)
	{
		cout <<"\tcudaGetDeviceCount returned" <<(int)error_id  << " %d\n-> " <<  cudaGetErrorString(error_id) << endl;
		cout <<"\tResult = FAIL\n" << endl;
		exit(EXIT_FAILURE);
	}

	// This function call returns 0 if there are no CUDA capable devices.
	if (deviceCount == 0)
	{
		cout << "\tThere are no available device(s) that support CUDA\n" << endl;
	}
	else
	{
		cout << "\tDetected " << deviceCount << " CUDA Capable device(s)\n" << endl;
	}

	int dev, driverVersion = 0, runtimeVersion = 0;

	for (dev = 0; dev < deviceCount; ++dev){
		cudaSetDevice(dev);
		cudaDeviceProp deviceProp;
		cudaGetDeviceProperties(&deviceProp, dev);

		printf("\tDevice %d: \"%s\"\n", dev, deviceProp.name);

		cout <<"\tDevice " << dev << ": \" " << deviceProp.name <<    " \" "<<endl;


		// Console log
		cudaDriverGetVersion(&driverVersion);
		cudaRuntimeGetVersion(&runtimeVersion);
		printf("\tCUDA Driver Version / Runtime Version          %d.%d / %d.%d\n", driverVersion/1000, (driverVersion%100)/10, runtimeVersion/1000, (runtimeVersion%100)/10);
		printf("\tCUDA Capability Major/Minor version number:    %d.%d\n", deviceProp.major, deviceProp.minor);

		char msg[256];
		sprintf(msg, "\tTotal amount of global memory:                 %.0f MBytes (%llu bytes)\n",
				(float)deviceProp.totalGlobalMem/1048576.0f, (unsigned long long) deviceProp.totalGlobalMem);
		printf("\t%s", msg);

		printf("\t(%2d) Multiprocessors x (%3d) CUDA Cores/MP:    %d CUDA Cores\n",
				deviceProp.multiProcessorCount,
				_ConvertSMVer2Cores(deviceProp.major, deviceProp.minor),
				_ConvertSMVer2Cores(deviceProp.major, deviceProp.minor) * deviceProp.multiProcessorCount);
		printf("\tGPU Clock rate:                                %.0f MHz (%0.2f GHz)\n", deviceProp.clockRate * 1e-3f, deviceProp.clockRate * 1e-6f);


		printf("\tMax Texture Dimension Size (x,y,z)             1D=(%d), 2D=(%d,%d), 3D=(%d,%d,%d)\n",
				deviceProp.maxTexture1D   , deviceProp.maxTexture2D[0], deviceProp.maxTexture2D[1],
				deviceProp.maxTexture3D[0], deviceProp.maxTexture3D[1], deviceProp.maxTexture3D[2]);
		printf("\tMax Layered Texture Size (dim) x layers        1D=(%d) x %d, 2D=(%d,%d) x %d\n",
				deviceProp.maxTexture1DLayered[0], deviceProp.maxTexture1DLayered[1],
				deviceProp.maxTexture2DLayered[0], deviceProp.maxTexture2DLayered[1], deviceProp.maxTexture2DLayered[2]);

		printf("\tTotal amount of constant memory:               %lu bytes\n", deviceProp.totalConstMem);
		printf("\tTotal amount of shared memory per block:       %lu bytes\n", deviceProp.sharedMemPerBlock);
		printf("\tTotal number of registers available per block: %d\n", deviceProp.regsPerBlock);
		printf("\tWarp size:                                     %d\n", deviceProp.warpSize);
		printf("\tMaximum number of threads per multiprocessor:  %d\n", deviceProp.maxThreadsPerMultiProcessor);
		printf("\tMaximum number of threads per block:           %d\n", deviceProp.maxThreadsPerBlock);
		printf("\tMaximum sizes of each dimension of a block:    %d x %d x %d\n",
				deviceProp.maxThreadsDim[0],
				deviceProp.maxThreadsDim[1],
				deviceProp.maxThreadsDim[2]);
		printf("\tMaximum sizes of each dimension of a grid:     %d x %d x %d\n",
				deviceProp.maxGridSize[0],
				deviceProp.maxGridSize[1],
				deviceProp.maxGridSize[2]);
		printf("\tMaximum memory pitch:                          %lu bytes\n", deviceProp.memPitch);
		printf("\tTexture alignment:                             %lu bytes\n", deviceProp.textureAlignment);
		printf("\tConcurrent copy and kernel execution:          %s with %d copy engine(s)\n", (deviceProp.deviceOverlap ? "Yes" : "No"), deviceProp.asyncEngineCount);
		printf("\tRun time limit on kernels:                     %s\n", deviceProp.kernelExecTimeoutEnabled ? "Yes" : "No");
		printf("\tIntegrated GPU sharing Host Memory:            %s\n", deviceProp.integrated ? "Yes" : "No");
		printf("\tSupport host page-locked memory mapping:       %s\n", deviceProp.canMapHostMemory ? "Yes" : "No");
		printf("\tAlignment requirement for Surfaces:            %s\n", deviceProp.surfaceAlignment ? "Yes" : "No");
		printf("\tDevice has ECC support:                        %s\n", deviceProp.ECCEnabled ? "Enabled" : "Disabled");

		printf("\tDevice supports Unified Addressing (UVA):      %s\n", deviceProp.unifiedAddressing ? "Yes" : "No");
		printf("\tDevice PCI Bus ID / PCI location ID:           %d / %d\n", deviceProp.pciBusID, deviceProp.pciDeviceID);

		const char *sComputeMode[] =
		{
				"Default (multiple host threads can use ::cudaSetDevice() with device simultaneously)",
				"Exclusive (only one host thread in one process is able to use ::cudaSetDevice() with this device)",
				"Prohibited (no host thread can use ::cudaSetDevice() with this device)",
				"Exclusive Process (many threads in one process is able to use ::cudaSetDevice() with this device)",
				"Unknown",
				NULL
		};
		printf("\tCompute Mode:\n");
		printf("\t\t< %s >\n", sComputeMode[deviceProp.computeMode]);
	}


	// csv masterlog info
	// *****************************
	// exe and CUDA driver name
	printf("\n");
	std::string sProfileString = "deviceQuery, CUDA Driver = CUDART";
	char cTemp[16];

	// driver version
	sProfileString += ", CUDA Driver Version = ";
#ifdef WIN32
	sprintf_s(cTemp, 10, "%d.%d", driverVersion/1000, (driverVersion%100)/10);
#else
	sprintf(cTemp, "%d.%d", driverVersion/1000, (driverVersion%100)/10);
#endif
	sProfileString +=  cTemp;

	// Runtime version
	sProfileString += ", CUDA Runtime Version = ";
#ifdef WIN32
	sprintf_s(cTemp, 10, "%d.%d", runtimeVersion/1000, (runtimeVersion%100)/10);
#else
	sprintf(cTemp, "%d.%d", runtimeVersion/1000, (runtimeVersion%100)/10);
#endif
	sProfileString +=  cTemp;

	// Device count
	sProfileString += ", NumDevs = ";
#ifdef WIN32
	sprintf_s(cTemp, 10, "%d", deviceCount);
#else
	sprintf(cTemp, "%d", deviceCount);
#endif
	sProfileString += cTemp;

	// Print Out all device Names
	for (dev = 0; dev < deviceCount; ++dev)
	{
#ifdef _WIN32
		sprintf_s(cTemp, 13, ", Device%d = ", dev);
#else
		sprintf(cTemp, ", Device%d = ", dev);
#endif
		cudaDeviceProp deviceProp;
		cudaGetDeviceProperties(&deviceProp, dev);
		sProfileString += cTemp;
		sProfileString += deviceProp.name;
	}

	sProfileString += "\n";
	printf("%s", sProfileString.c_str());

	printf("\tResult = PASS\n");

	// finish

	cout << "CUDA DEVISE TEST - END" << endl;
}

void read_external_geometry_file_specification_for_LBM(int &lx, int &ly, int &lz, int &n_of_densities, const string filename){
	vector<string> geometry_parameters;
	ifstream conf_file(filename.c_str());
	string buff;
	if(conf_file.is_open()){
		while(conf_file>>buff){
			geometry_parameters.push_back(buff);
		}
		cout << "Geometry Parameters Read:" << endl;
		lx=atoi(geometry_parameters[0].c_str());
		cout << "\t domain length in X: " << lx << endl;

		ly=atoi(geometry_parameters[1].c_str());
		cout << "\t domain length in Y: " << ly << endl;

		lz=atoi(geometry_parameters[2].c_str());
		cout << "\t domain length in Z: " << lz << endl;

		n_of_densities=atoi(geometry_parameters[3].c_str());
		cout << "\t number of densities on each node: " << n_of_densities << endl;

		cout <<"total:" << geometry_parameters.size() << " parameters were read" << endl;
		conf_file.close();
	}else{
		cout << "The file "<< filename << " was not found" << endl;
		cout << "Create a new file at the root directory with 3 lines (one number on every line), each corresponding to the respective dimension of X,Y,Z" << endl;
		exit (-1);
	}
}


lattice::lattice(int LX,int LY, int LZ):
																																																		lx(LX), ly(LY), lz(LZ){
	Q0=new FLOATING[lz*ly*lx];
	Q1=new FLOATING[lz*ly*lx];
	Q2=new FLOATING[lz*ly*lx];
	Q3=new FLOATING[lz*ly*lx];
	Q4=new FLOATING[lz*ly*lx];
	Q5=new FLOATING[lz*ly*lx];
	Q6=new FLOATING[lz*ly*lx];
	Q7=new FLOATING[lz*ly*lx];
	Q8=new FLOATING[lz*ly*lx];
	Q9=new FLOATING[lz*ly*lx];
	Q10=new FLOATING[lz*ly*lx];
	Q11=new FLOATING[lz*ly*lx];
	Q12=new FLOATING[lz*ly*lx];
	Q13=new FLOATING[lz*ly*lx];
	Q14=new FLOATING[lz*ly*lx];
	Q15=new FLOATING[lz*ly*lx];
	Q16=new FLOATING[lz*ly*lx];
	Q17=new FLOATING[lz*ly*lx];
	Q18=new FLOATING[lz*ly*lx];

	initialise(Q0);
	initialise(Q1);
	initialise(Q2);
	initialise(Q3);
	initialise(Q4);
	initialise(Q5);
	initialise(Q6);
	initialise(Q7);
	initialise(Q8);
	initialise(Q9);
	initialise(Q10);
	initialise(Q11);
	initialise(Q12);
	initialise(Q13);
	initialise(Q14);
	initialise(Q15);
	initialise(Q16);
	initialise(Q17);
	initialise(Q18);
}

lattice::lattice(int LX,int LY, int LZ, int dump):
																																				lx(LX), ly(LY), lz(LZ){
	int FLOATING_array_size=lx*ly*lz*sizeof(FLOATING);
	cudaMalloc((void **)&Q0, FLOATING_array_size);
	cudaMalloc((void **)&Q1, FLOATING_array_size);
	cudaMalloc((void **)&Q2, FLOATING_array_size);
	cudaMalloc((void **)&Q3, FLOATING_array_size);
	cudaMalloc((void **)&Q4, FLOATING_array_size);
	cudaMalloc((void **)&Q5, FLOATING_array_size);
	cudaMalloc((void **)&Q6, FLOATING_array_size);
	cudaMalloc((void **)&Q7, FLOATING_array_size);
	cudaMalloc((void **)&Q8, FLOATING_array_size);
	cudaMalloc((void **)&Q9, FLOATING_array_size);
	cudaMalloc((void **)&Q10, FLOATING_array_size);
	cudaMalloc((void **)&Q11, FLOATING_array_size);
	cudaMalloc((void **)&Q12, FLOATING_array_size);
	cudaMalloc((void **)&Q13, FLOATING_array_size);
	cudaMalloc((void **)&Q14, FLOATING_array_size);
	cudaMalloc((void **)&Q15, FLOATING_array_size);
	cudaMalloc((void **)&Q16, FLOATING_array_size);
	cudaMalloc((void **)&Q17, FLOATING_array_size);
	cudaMalloc((void **)&Q18, FLOATING_array_size);
}
lattice::~lattice(){
	delete [] Q0;
	delete [] Q1;
	delete [] Q2;
	delete [] Q3;
	delete [] Q4;
	delete [] Q5;
	delete [] Q6;
	delete [] Q7;
	delete [] Q8;
	delete [] Q9;
	delete [] Q10;
	delete [] Q11;
	delete [] Q12;
	delete [] Q13;
	delete [] Q14;
	delete [] Q15;
	delete [] Q16;
	delete [] Q17;
	delete [] Q18;

	printf("host memories deleted!\n");
}

void lattice::initialise(FLOATING *Q){

	for(int z=0; z<lz; ++z)
		for(int y=0; y<ly; ++y)
			for(int x=0; x<lx; ++x)
				Q[ index(z,y,x)]=0.0;
}

void LBM::create_an_example_configuration_files(const string filename){
	ofstream example_file(filename.c_str());

	if ( example_file.is_open()){
		example_file << "10" << endl;
		example_file << "100" << endl;
		example_file << "0.0175" << endl;
		example_file << "7" << endl;
		example_file << "100" << endl;
		example_file << "26" << endl;
		example_file << "59" << endl;
		example_file << "512" << endl;
		example_file << "datum_design_case_name1" << endl;
		example_file.close();
	}
}

void LBM::display_the_structure_of_an_example_configuration_file(){
	cout<< "Create a new file (with 9 lines) at the root directory following the template below:"
			<< endl;
	cout << "10 <--line 1: number of iterations" << endl;
	cout << "100 <--line 2: check frequency" << endl;
	cout << "0.0175 <--line 3: nu" << endl;
	cout << "7 <--line 4: r_small (from baffle geometry)" << endl;
	cout << "100 <--line 5: Reynolds Number" << endl;
	cout << "26 <--line 6: S (from baffle geometry)" << endl;
	cout << "59 <--line 7: baffle possition" << endl;
	cout << "512 <--line 8: CUDA threads per kernel" << endl;
	cout << "datum_design <--line 9: case name (ONE WORD!)" << endl;
}

void LBM::read_external_configuration_file_for_the_solver(const string filename) {
	vector<string> configuration_parameters;
	ifstream conf_file(filename.c_str());
	string buff;
	if (conf_file.is_open()) {
		while (conf_file >> buff) {
			configuration_parameters.push_back(buff);
		}
		cout << "Configuration Parameters Read:" << endl;
		max_iterations = atoi(configuration_parameters[0].c_str());
		cout << "\titerations: " << max_iterations << endl;
		//check step:perform check_density and export

		check_step = atoi(configuration_parameters[1].c_str());
		cout << "\tcheck step: " << check_step << endl;

		nu = atof(configuration_parameters[2].c_str());
		cout << "\tnu: " << nu << endl;

		r_small = atof(configuration_parameters[3].c_str());
		cout << "\tr_small: " << r_small << endl;

		reynolds = atof(configuration_parameters[4].c_str());
		cout << "\treynolds: " << reynolds << endl;

		s = atof(configuration_parameters[5].c_str());
		cout << "\ts: " << s << endl;

		baffle = atoi(configuration_parameters[6].c_str());
		cout << "\tbaffle position on X=" << baffle << endl;

		threads_per_kernel = atoi(configuration_parameters[7].c_str());
		cout << "\tCUDA threads per kernel: " << threads_per_kernel << endl;

		case_name = configuration_parameters[8].c_str();
		cout << "Case: " << case_name << endl;

		cout << "total:" << configuration_parameters.size()
								<< " parameters were read" << endl;
		conf_file.close();
	} else {
		cout << "The file "<<filename <<" was not found" << endl;
		display_the_structure_of_an_example_configuration_file();
		create_an_example_configuration_files(filename);
		exit(-2);
	}
}

void LBM::delete_device_data(){

	cudaFree(D3_d.Q0);
	cudaFree(D3_d.Q1);
	cudaFree(D3_d.Q2);
	cudaFree(D3_d.Q3);
	cudaFree(D3_d.Q4);
	cudaFree(D3_d.Q5);
	cudaFree(D3_d.Q6);
	cudaFree(D3_d.Q7);
	cudaFree(D3_d.Q8);
	cudaFree(D3_d.Q9);
	cudaFree(D3_d.Q10);
	cudaFree(D3_d.Q11);
	cudaFree(D3_d.Q12);
	cudaFree(D3_d.Q13);
	cudaFree(D3_d.Q14);
	cudaFree(D3_d.Q15);
	cudaFree(D3_d.Q16);
	cudaFree(D3_d.Q17);
	cudaFree(D3_d.Q18);

	cudaFree(D3_hlp_d.Q0);
	cudaFree(D3_hlp_d.Q1);
	cudaFree(D3_hlp_d.Q2);
	cudaFree(D3_hlp_d.Q3);
	cudaFree(D3_hlp_d.Q4);
	cudaFree(D3_hlp_d.Q5);
	cudaFree(D3_hlp_d.Q6);
	cudaFree(D3_hlp_d.Q7);
	cudaFree(D3_hlp_d.Q8);
	cudaFree(D3_hlp_d.Q9);
	cudaFree(D3_hlp_d.Q10);
	cudaFree(D3_hlp_d.Q11);
	cudaFree(D3_hlp_d.Q12);
	cudaFree(D3_hlp_d.Q13);
	cudaFree(D3_hlp_d.Q14);
	cudaFree(D3_hlp_d.Q15);
	cudaFree(D3_hlp_d.Q16);
	cudaFree(D3_hlp_d.Q17);
	cudaFree(D3_hlp_d.Q18);

	cudaFree(u_current_d);
	cudaFree(u_current_temp_d);
	cudaFree(v_current_d);
	cudaFree(w_current_d);

	cudaFree(u_previous_spatial_boundary_d);
	cudaFree(v_previous_spatial_boundary_d);
	cudaFree(w_previous_spatial_boundary_d);

	cudaFree(u_previous_temporal_boundary_d);
	cudaFree(v_previous_temporal_boundary_d);
	cudaFree(w_previous_temporal_boundary_d);

	cudaFree(temp_check_density_d);
	cudaFree(temp_check_density_d_full);

	cudaFree(obstacles_d);

	cudaFree(temp_Uc_d);


	printf("cuda memories deleted!\n");
}


__global__
void cuda_initialise_array(FLOATING *input_array, const int length, const FLOATING value){
	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	input_array[tid]=value;
	__syncthreads();

}



template <class T>
void LBM::initialise_array(T *array, const int length,const T init_value){

	for(int i=0; i<length; ++i)
		array[i]=init_value;

}


void LBM::abstract_initialise(){
	//objectives


	FLOATING temp_density=0.0;

	initialise_microscopic_density_arrays_in_the_host();

	calculate_macroscopic_density_in_the_host(temp_density);

	create_reactor_geometry_in_the_host();

#ifdef PRODUCE_OUTPUT_FILES
	geometry_file_in_VTK();
#endif //PRODUCE_OUTPUT_FILES
	cout << "0th loop" <<endl;
	//myLBM.initial_redistribute();
	fortran_redistribute(0);
	//first loop!!!
	calculate_macroscopic_density_in_the_host(temp_density);
	streaming();
	bounceback();
	initial_relaxation();
	convective_BC();
	cout <<"# iteration 0" << endl;

	count_no_obstacles_at_penultimate_x_slice();
	copy_data_from_host_to_device();

}

void LBM::abstract_check_density(){
	//check density
	if( time_unit%check_step==0){
#ifdef CPU_part
		calculate_macroscopic_density_in_the_host(density);
#endif
#ifdef GPU_part
		cuda_check_density(time_unit);
#endif
	}
}


void LBM::abstract_debug_computations(){
#ifdef DEBUG
	compare_obstacles(obstacles);
	compare_nodes_hlp(n_hlp);
	compare_nodes(node);
#endif
}

void LBM::abstract_redistribute(){
	//redistribute
#ifdef CPU_part
	redistribute();
#endif
#ifdef GPU_part
	cuda_redistribute();
#endif
}

void LBM::abstract_streaming(){
	//streaming
#ifdef CPU_part
	streaming();
#endif
#ifdef GPU_part
	cuda_streaming();
#endif
}

void LBM::abstract_bounce_back(){
	//bounce back
#ifdef CPU_part
	bounceback();
#endif
#ifdef GPU_part
	cuda_bounceback();
#endif
}

void LBM::abstract_relaxation(){
	//relaxation
#ifdef CPU_part
	relaxation();
#endif
#ifdef GPU_part
	cuda_relaxation();
#endif
}

void LBM::abstract_convective_boundary_conditions(){
	//convective BC
#ifdef CPU_part
	convective_BC();
#endif
#ifdef GPU_part
	cuda_convective_BC();
#endif
}

void LBM::core_computations(){
	//LBM CORE
	abstract_redistribute();

	abstract_streaming();

	abstract_bounce_back();

	abstract_relaxation();

	abstract_convective_boundary_conditions();
}

void LBM::compute_domain(){
	//starting from second loop!
	for (time_unit = 1; time_unit<max_iterations ; ++time_unit){
		//		cout <<"# iteration " << time_unit << endl;
		cout <<time_unit <<". ";

		abstract_check_density();

		abstract_debug_computations();

		core_computations();

	}
	cout << endl;
}

void LBM::export_solution(){


	calculate_macroscopic_quantities(time_unit);
#ifdef PRODUCE_OUTPUT_FILES
	write_VTK_SI(time_unit);
#endif //PRODUCE_OUTPUT_FILES
}

template void
LBM::initialise_array<double>(double *array, const int length,const double init_value);

template void
LBM::initialise_array<float>(float *array, const int length,const float init_value);


template <typename T>
void LBM::allocate_and_initialise(T *array, const int length){

	cudaMalloc((void **)&array, length*sizeof(T));

	int n_of_threads=threads_per_kernel;
	int n_of_blocks=ceil((length*1.0)/n_of_threads);

	dim3 threads_type2(n_of_threads,1,1);
	dim3 grid_type2(n_of_blocks,1,1);
	//kane to template!
	//cuda_initialise_array<<<grid_type2,threads_type2>>>(array, length, 0.0);
}

void LBM::calculate_CUDA_quantities() {
	threads_for_streaming_collision_and_relaxation=threads_per_kernel;
	blocks_for_streaming_collision_and_relaxation= (three_dimensional_length)/threads_per_kernel;
	if ((three_dimensional_length%threads_per_kernel)!=0)
		++blocks_for_streaming_collision_and_relaxation;
	size_of_allocated_shared_memory_for_streaming_collision_and_relaxation=threads_per_kernel*sizeof(FLOATING);

	convective_boundary_conditions_blocks=two_dimensional_length/threads_per_kernel;
	if ( (two_dimensional_length%threads_per_kernel)!=0 )
		++convective_boundary_conditions_blocks;
}

void LBM::reset_convergence_file(){
	//delete previous convergence.txt and create a new one
	if( remove("LBM2_convergence.txt")!=0 )
		cout <<"couldn't delete LBM2_convergence.txt" << endl;
	else
		cout<< "creating LBM2_convergence.txt" << endl;

	ofstream convergence_file("LBM2_convergence.txt");
	convergence_file<<"#iteration ; converegence_value" << endl;
	convergence_file.close();
}

void LBM::display_CUDA_specifications(){
	cout <<"CUDA specifications:" <<endl;
	cout <<"\tstreaming/collision/relaxation:" << endl;
	cout <<"\t\tthreads: "<<threads_for_streaming_collision_and_relaxation<<endl;
	cout <<"\t\tblocks: "<<blocks_for_streaming_collision_and_relaxation<<endl;
	cout <<"\t\tshare memory size: "<< size_of_allocated_shared_memory_for_streaming_collision_and_relaxation <<endl;

	cout <<"\tconvective boundary conditions:" << endl;
	cout <<"\t\tthreads:" << threads_per_kernel << endl;
	cout <<"\t\tblocks: "<< convective_boundary_conditions_blocks << endl;
}

void LBM::initialise_host_data(){
	initialise_array(obstacles,   three_dimensional_length,0 );

	initialise_array<FLOATING>(u_current,   two_dimensional_length,0.0);
	initialise_array<FLOATING>(v_current,   two_dimensional_length,0.0);
	initialise_array<FLOATING>(w_current,   two_dimensional_length,0.0);
	// u_previous_spatial_boundary: at boundary - 1 (in x)
	initialise_array<FLOATING>(u_previous_spatial_boundary,   two_dimensional_length,0.0);
	initialise_array<FLOATING>(v_previous_spatial_boundary,   two_dimensional_length,0.0);
	initialise_array<FLOATING>(w_previous_spatial_boundary,   two_dimensional_length,0.0);
	// u_prev: at boundary - 1 (in time)
	initialise_array<FLOATING>(u_previous_temporal_boundary,   two_dimensional_length,0.0);
	initialise_array<FLOATING>(v_previous_temporal_boundary,   two_dimensional_length,0.0);
	initialise_array<FLOATING>(w_previous_temporal_boundary,   two_dimensional_length,0.0);

}

void LBM::allocate_device_arrays(){
	//allocate additional cuda memories
	cudaMalloc((void **)&u_current_d, FLOATING_slice_size);

	cudaMalloc((void **)&u_current_temp_d, FLOATING_slice_size);
	cudaMalloc((void **)&v_current_d, FLOATING_slice_size);
	cudaMalloc((void **)&w_current_d, FLOATING_slice_size);

	cudaMalloc((void **)&u_previous_spatial_boundary_d, FLOATING_slice_size);
	cudaMalloc((void **)&v_previous_spatial_boundary_d, FLOATING_slice_size);
	cudaMalloc((void **)&w_previous_spatial_boundary_d, FLOATING_slice_size);

	cudaMalloc((void **)&u_previous_temporal_boundary_d, FLOATING_slice_size);
	cudaMalloc((void **)&v_previous_temporal_boundary_d, FLOATING_slice_size);
	cudaMalloc((void **)&w_previous_temporal_boundary_d, FLOATING_slice_size);

	cudaMalloc((void **)&temp_check_density_d, lx*ly*lz/4 *sizeof(FLOATING));
	cudaMalloc((void **)&temp_check_density_d_full, lx*ly*lz*sizeof(FLOATING));

	cudaMalloc((void **)&obstacles_d, int_array_size);

	cudaMalloc((void **)&temp_Uc_d, 2*sizeof(FLOATING));
}

void LBM::initialise_device_data(){
	dim3 threads_type2(threads_for_streaming_collision_and_relaxation,1,1);
	dim3 grid_type2(blocks_for_streaming_collision_and_relaxation,1,1);
	cuda_initialise_array<<<grid_type2,threads_type2>>>(temp_check_density_d_full, lx*ly*lz, 0.0);
}


void LBM::initialise_all_data_arrays(){
	initialise_host_data();
	allocate_device_arrays();
	initialise_device_data();
}

void LBM::display_LBM_specifications(){
	cout << "constructing LBM(built-in quantities)" << endl;

	cout << "\tdensity" << density << endl;
	cout << "\tt_0=" << t_0 << endl;
	cout << "\tt_1=" << t_1 << endl;
	cout << "\tt_2=" << t_2 << endl;
	cout << "\tc_squ=" << c_squ << endl;
	cout << "\ttau=" << tau << endl;
	cout << "\tomega=" << omega << endl;
}



LBM::LBM(const int &LX, const int &LY, const int &LZ, const FLOATING &DENSITY, const FLOATING &T_0,
		const FLOATING &T_1, const FLOATING &T_2, const FLOATING &C_SQU):
		time_elapsed(0),
		max_iterations(1000),
		check_step(100),
		lx(LX), ly(LY), lz(LZ),
		lattice_nodes(lx*ly*lz), no_obstacle_lattices_at_penultimate_x_slice(0),
		threads_for_streaming_collision_and_relaxation(512),
		blocks_for_streaming_collision_and_relaxation(32),
		size_of_allocated_shared_memory_for_streaming_collision_and_relaxation(48*1024),
		convective_boundary_conditions_blocks(32),
		nu(0.0175), r_small(6.67897), reynolds(195.732), s(23.7849), density(DENSITY),
		t_0(density*T_0), t_1(density*T_1), t_2(density*T_2), c_squ(C_SQU), reciprocal_c_squ(1.0/c_squ),
		baffle(XBAFFLE), threads_per_kernel(MANY_THREADS), time_unit(0),
		two_dimensional_length(ly*lz),
		three_dimensional_length(lx*ly*lz),
		FLOATING_slice_size((two_dimensional_length)*sizeof(FLOATING)),
		int_array_size((three_dimensional_length)*sizeof(int)),
		tau(3.0*nu + 0.5), omega(1.0 /tau), one_minus_omega (1.0-omega),
		pr_diff(0.0), pr_out(0.0), pr_in(0.0), vor(0.0),
		D3(lx, ly, lz), D3_hlp(lx, ly, lz), obstacles(new int[lz*ly*lx]),
		u_current(new FLOATING[ly*lz]),  v_current(new FLOATING[ly*lz]),  w_current(new FLOATING[ly*lz]),
		u_previous_spatial_boundary(new FLOATING[ly*lz]), v_previous_spatial_boundary(new FLOATING[ly*lz]), w_previous_spatial_boundary(new FLOATING[ly*lz]),
		u_previous_temporal_boundary(new FLOATING[ly*lz]), v_previous_temporal_boundary(new FLOATING[ly*lz]), w_previous_temporal_boundary(new FLOATING[ly*lz]),
		u_current_d(NULL), u_current_temp_d(NULL), v_current_d(NULL), w_current_d(NULL),
		// u_previous_spatial_boundary: at boundary - 1 (in x)
		u_previous_spatial_boundary_d(NULL), v_previous_spatial_boundary_d(NULL), w_previous_spatial_boundary_d(NULL),
		// u_prev: at boundary - 1 (in time)
		u_previous_temporal_boundary_d(NULL), v_previous_temporal_boundary_d(NULL), w_previous_temporal_boundary_d(NULL),
		temp_cpu_u_current_d(NULL), temp_cpu_v_current_d(NULL), temp_cpu_w_current_d(NULL),
		temp_cpu_u_previous_temporal_boundary_d(NULL), temp_cpu_v_previous_temporal_boundary_d(NULL), temp_cpu_w_previous_temporal_boundary_d(NULL),
		temp_cpu_u_previous_spatial_boundary_d(NULL), temp_cpu_v_previous_spatial_boundary_d(NULL), temp_cpu_w_previous_spatial_boundary_d(NULL),
		temp_check_density_d(NULL), temp_check_density_d_full(NULL),
		data_location(CPU),
		temp_Uc_d(NULL), obstacles_d(NULL),
		D3_d(lx, ly, lz, 0), D3_hlp_d(lx, ly, lz, 0),
		Ux(new FLOATING[lx*ly*lz]),
		Uy(new FLOATING[lx*ly*lz]),
		Uz(new FLOATING[lx*ly*lz]),
		Pressure(new FLOATING[lx*ly*lz]),
		Wx(new FLOATING[lx*ly*lz]),
		Wy(new FLOATING[lx*ly*lz]),
		Wz(new FLOATING[lx*ly*lz]){
	cout << "***LBM Starting***" << endl;

	time (&time_start);

	read_external_configuration_file_for_the_solver("LBM2_configuration.txt");

	reset_convergence_file();

	calculate_CUDA_quantities();

	display_CUDA_specifications();

	initialise_all_data_arrays();

	display_LBM_specifications();

	abstract_initialise();
}

void LBM::delete_host_memories(){
	delete [] obstacles;

	delete [] u_current;
	delete [] v_current;
	delete [] w_current;

	delete [] u_previous_spatial_boundary;
	delete [] v_previous_spatial_boundary;
	delete [] w_previous_spatial_boundary;

	delete [] u_previous_temporal_boundary;
	delete [] v_previous_temporal_boundary;
	delete [] w_previous_temporal_boundary;

	delete [] Ux;
	delete [] Uy;
	delete [] Uz;
	delete [] Pressure;
	delete [] Wx;
	delete [] Wy;
	delete [] Wz;
}

LBM::~LBM(){

	delete_host_memories();
	delete_device_data();


	cout <<"all memories were deallocated!" <<endl;


	cout << endl << "LBM2 ended in "<< time_elapsed<< "secs !" << endl; // prints
}

void LBM::count_no_obstacles_at_penultimate_x_slice(){

	no_obstacle_lattices_at_penultimate_x_slice = 0;
#pragma unroll
	for (int z = 0 ; z< lz ; ++z){
#pragma unroll
		for (int y = 0 ; y< ly ; ++y){
			if (obstacles[index(z,y,(lx-1))]==0) {
				++no_obstacle_lattices_at_penultimate_x_slice ;
			}
		}
	}
	cout << "number of free lattices at U direction at the penultimate slice:" << no_obstacle_lattices_at_penultimate_x_slice <<endl;

}

void LBM::copy_data_from_host_to_device(){//copy data to CUDA variables
	int const array_length=(lx*ly*lz);
	int const slice_length=(ly*lz);
	int const FLOATING_array_size=array_length*sizeof(FLOATING);
	int const FLOATING_slice_size=slice_length*sizeof(FLOATING);
	int const int_array_size=array_length*sizeof(int);




	cudaMemcpy(D3_d.Q0 ,D3.Q0,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q1 ,D3.Q1,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q2 ,D3.Q2,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q3 ,D3.Q3,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q4 ,D3.Q4,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q5 ,D3.Q5,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q6 ,D3.Q6,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q7 ,D3.Q7,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q8 ,D3.Q8,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q9 ,D3.Q9,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q10 ,D3.Q10,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q11 ,D3.Q11,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q12 ,D3.Q12,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q13 ,D3.Q13,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q14 ,D3.Q14,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q15 ,D3.Q15,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q16 ,D3.Q16,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q17 ,D3.Q17,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q18 ,D3.Q18,FLOATING_array_size,cudaMemcpyHostToDevice);

	cudaMemcpy(D3_hlp_d.Q0 ,D3_hlp.Q0,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q1 ,D3_hlp.Q1,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q2 ,D3_hlp.Q2,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q3 ,D3_hlp.Q3,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q4 ,D3_hlp.Q4,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q5 ,D3_hlp.Q5,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q6 ,D3_hlp.Q6,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q7 ,D3_hlp.Q7,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q8 ,D3_hlp.Q8,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q9 ,D3_hlp.Q9,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q10 ,D3_hlp.Q10,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q11 ,D3_hlp.Q11,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q12 ,D3_hlp.Q12,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q13 ,D3_hlp.Q13,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q14 ,D3_hlp.Q14,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q15 ,D3_hlp.Q15,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q16 ,D3_hlp.Q16,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q17 ,D3_hlp.Q17,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q18 ,D3_hlp.Q18,FLOATING_array_size,cudaMemcpyHostToDevice);


	cudaMemcpy(obstacles_d ,obstacles,int_array_size,cudaMemcpyHostToDevice);

	cudaMemcpy(u_current_d ,u_current,FLOATING_slice_size,cudaMemcpyHostToDevice);
	cudaMemcpy(u_current_temp_d ,u_current,FLOATING_slice_size,cudaMemcpyHostToDevice);
	cudaMemcpy(v_current_d ,v_current,FLOATING_slice_size,cudaMemcpyHostToDevice);
	cudaMemcpy(w_current_d ,w_current,FLOATING_slice_size,cudaMemcpyHostToDevice);

	cudaMemcpy(u_previous_spatial_boundary_d ,u_previous_spatial_boundary,FLOATING_slice_size,cudaMemcpyHostToDevice);
	cudaMemcpy(v_previous_spatial_boundary_d ,v_previous_spatial_boundary,FLOATING_slice_size,cudaMemcpyHostToDevice);
	cudaMemcpy(w_previous_spatial_boundary_d ,w_previous_spatial_boundary,FLOATING_slice_size,cudaMemcpyHostToDevice);

	cudaMemcpy(u_previous_temporal_boundary_d ,u_previous_temporal_boundary,FLOATING_slice_size,cudaMemcpyHostToDevice);
	cudaMemcpy(v_previous_temporal_boundary_d ,v_previous_temporal_boundary,FLOATING_slice_size,cudaMemcpyHostToDevice);
	cudaMemcpy(w_previous_temporal_boundary_d ,w_previous_temporal_boundary,FLOATING_slice_size,cudaMemcpyHostToDevice);

	data_location=GPU;
	printf("all data were copied to device\n");
}

void LBM::small_copy_data_from_host_to_device(){//copy data to CUDA variables
	int const array_length=(lx*ly*lz);

	int const FLOATING_array_size=array_length*sizeof(FLOATING);



	cudaMemcpy(D3_d.Q0 ,D3.Q0,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q1 ,D3.Q1,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q2 ,D3.Q2,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q3 ,D3.Q3,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q4 ,D3.Q4,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q5 ,D3.Q5,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q6 ,D3.Q6,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q7 ,D3.Q7,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q8 ,D3.Q8,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q9 ,D3.Q9,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q10 ,D3.Q10,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q11 ,D3.Q11,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q12 ,D3.Q12,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q13 ,D3.Q13,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q14 ,D3.Q14,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q15 ,D3.Q15,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q16 ,D3.Q16,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q17 ,D3.Q17,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_d.Q18 ,D3.Q18,FLOATING_array_size,cudaMemcpyHostToDevice);

	cudaMemcpy(D3_hlp_d.Q0 ,D3_hlp.Q0,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q1 ,D3_hlp.Q1,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q2 ,D3_hlp.Q2,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q3 ,D3_hlp.Q3,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q4 ,D3_hlp.Q4,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q5 ,D3_hlp.Q5,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q6 ,D3_hlp.Q6,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q7 ,D3_hlp.Q7,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q8 ,D3_hlp.Q8,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q9 ,D3_hlp.Q9,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q10 ,D3_hlp.Q10,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q11 ,D3_hlp.Q11,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q12 ,D3_hlp.Q12,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q13 ,D3_hlp.Q13,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q14 ,D3_hlp.Q14,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q15 ,D3_hlp.Q15,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q16 ,D3_hlp.Q16,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q17 ,D3_hlp.Q17,FLOATING_array_size,cudaMemcpyHostToDevice);
	cudaMemcpy(D3_hlp_d.Q18 ,D3_hlp.Q18,FLOATING_array_size,cudaMemcpyHostToDevice);

	data_location=GPU;
	printf("all data were copied to device\n");
}

void LBM::copy_data_from_device_to_host(){
	int const array_length=(lx*ly*lz);
	int const slice_length=(ly*lz);
	int const FLOATING_array_size=array_length*sizeof(FLOATING);
	int const FLOATING_slice_size=slice_length*sizeof(FLOATING);
	int const int_array_size=array_length*sizeof(int);



	cudaMemcpy(D3.Q0 ,D3_d.Q0,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q1 ,D3_d.Q1,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q2 ,D3_d.Q2,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q3 ,D3_d.Q3,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q4 ,D3_d.Q4,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q5 ,D3_d.Q5,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q6 ,D3_d.Q6,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q7 ,D3_d.Q7,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q8 ,D3_d.Q8,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q9 ,D3_d.Q9,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q10 ,D3_d.Q10,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q11 ,D3_d.Q11,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q12 ,D3_d.Q12,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q13 ,D3_d.Q13,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q14 ,D3_d.Q14,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q15 ,D3_d.Q15,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q16 ,D3_d.Q16,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q17 ,D3_d.Q17,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q18 ,D3_d.Q18,FLOATING_array_size,cudaMemcpyDeviceToHost);

	cudaMemcpy(D3_hlp.Q0 ,D3_hlp_d.Q0,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q1 ,D3_hlp_d.Q1,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q2 ,D3_hlp_d.Q2,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q3 ,D3_hlp_d.Q3,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q4 ,D3_hlp_d.Q4,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q5 ,D3_hlp_d.Q5,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q6 ,D3_hlp_d.Q6,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q7 ,D3_hlp_d.Q7,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q8 ,D3_hlp_d.Q8,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q9 ,D3_hlp_d.Q9,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q10 ,D3_hlp_d.Q10,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q11 ,D3_hlp_d.Q11,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q12 ,D3_hlp_d.Q12,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q13 ,D3_hlp_d.Q13,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q14 ,D3_hlp_d.Q14,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q15 ,D3_hlp_d.Q15,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q16 ,D3_hlp_d.Q16,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q17 ,D3_hlp_d.Q17,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q18 ,D3_hlp_d.Q18,FLOATING_array_size,cudaMemcpyDeviceToHost);

	cudaMemcpy(obstacles ,obstacles_d,int_array_size,cudaMemcpyDeviceToHost);

	cudaMemcpy(u_current ,u_current_d,FLOATING_slice_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(v_current ,v_current_d,FLOATING_slice_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(w_current ,w_current_d,FLOATING_slice_size,cudaMemcpyDeviceToHost);

	cudaMemcpy(u_previous_spatial_boundary ,u_previous_spatial_boundary_d,FLOATING_slice_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(v_previous_spatial_boundary ,v_previous_spatial_boundary_d,FLOATING_slice_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(w_previous_spatial_boundary ,w_previous_spatial_boundary_d,FLOATING_slice_size,cudaMemcpyDeviceToHost);

	cudaMemcpy(u_previous_temporal_boundary ,u_previous_temporal_boundary_d,FLOATING_slice_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(v_previous_temporal_boundary ,v_previous_temporal_boundary_d,FLOATING_slice_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(w_previous_temporal_boundary ,w_previous_temporal_boundary_d,FLOATING_slice_size,cudaMemcpyDeviceToHost);

	data_location=CPU;
	printf("all data were copied to host\n");
}

void LBM::small_copy_data_from_device_to_host(){
	int const array_length=(lx*ly*lz);

	int const FLOATING_array_size=array_length*sizeof(FLOATING);


	cudaMemcpy(D3.Q0 ,D3_d.Q0,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q1 ,D3_d.Q1,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q2 ,D3_d.Q2,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q3 ,D3_d.Q3,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q4 ,D3_d.Q4,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q5 ,D3_d.Q5,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q6 ,D3_d.Q6,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q7 ,D3_d.Q7,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q8 ,D3_d.Q8,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q9 ,D3_d.Q9,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q10 ,D3_d.Q10,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q11 ,D3_d.Q11,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q12 ,D3_d.Q12,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q13 ,D3_d.Q13,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q14 ,D3_d.Q14,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q15 ,D3_d.Q15,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q16 ,D3_d.Q16,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q17 ,D3_d.Q17,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3.Q18 ,D3_d.Q18,FLOATING_array_size,cudaMemcpyDeviceToHost);

	cudaMemcpy(D3_hlp.Q0 ,D3_hlp_d.Q0,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q1 ,D3_hlp_d.Q1,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q2 ,D3_hlp_d.Q2,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q3 ,D3_hlp_d.Q3,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q4 ,D3_hlp_d.Q4,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q5 ,D3_hlp_d.Q5,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q6 ,D3_hlp_d.Q6,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q7 ,D3_hlp_d.Q7,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q8 ,D3_hlp_d.Q8,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q9 ,D3_hlp_d.Q9,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q10 ,D3_hlp_d.Q10,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q11 ,D3_hlp_d.Q11,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q12 ,D3_hlp_d.Q12,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q13 ,D3_hlp_d.Q13,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q14 ,D3_hlp_d.Q14,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q15 ,D3_hlp_d.Q15,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q16 ,D3_hlp_d.Q16,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q17 ,D3_hlp_d.Q17,FLOATING_array_size,cudaMemcpyDeviceToHost);
	cudaMemcpy(D3_hlp.Q18 ,D3_hlp_d.Q18,FLOATING_array_size,cudaMemcpyDeviceToHost);


	data_location=CPU;
	printf("all data were copied to host\n");
}

void LBM::compare_obstacles(int *outter_obst){

	int x,y,z;

	for ( z =0 ; z < lz ; ++z){
		for (y = 0; y < ly ; ++y){
			for ( x = 0 ; x< lx ; ++x){
				if( obstacles[index(z,y,x)]!=outter_obst[index(z,y,x)]){
					cout << "obstacle miss-match @" << x << " " << y << " " << z <<endl;
					exit(-2);
				}
			}
		}
	}

	cout << "obstacles ok" << endl;

}

void LBM::compare_nodes(FLOATING *outter_node){
	int x,y,z,i;

	for ( z =0 ; z < lz ; ++z){
		for (y = 0; y < ly ; ++y){
			for ( x = 0 ; x< lx ; ++x){

				i=0;
				if( abs(outter_node[index4D(z,y,x,i)]-D3.Q0[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" << i <<  " :" << outter_node[index4D(z,y,x,i)] << " vs "<< D3.Q0[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=1;
				if( abs(outter_node[index4D(z,y,x,i)]-D3.Q1[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node[index4D(z,y,x,i)] << " vs "<< D3.Q1[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=2;
				if( abs(outter_node[index4D(z,y,x,i)]-D3.Q2[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node[index4D(z,y,x,i)] << " vs "<< D3.Q2[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=3;
				if( abs(outter_node[index4D(z,y,x,i)]-D3.Q3[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node[index4D(z,y,x,i)] << " vs "<< D3.Q3[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=4;
				if( abs(outter_node[index4D(z,y,x,i)]-D3.Q4[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node[index4D(z,y,x,i)] << " vs "<< D3.Q4[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=5;
				if( abs(outter_node[index4D(z,y,x,i)]-D3.Q5[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node[index4D(z,y,x,i)] << " vs "<< D3.Q5[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=6;
				if( abs(outter_node[index4D(z,y,x,i)]-D3.Q6[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node[index4D(z,y,x,i)] << " vs "<< D3.Q6[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=7;
				if( abs(outter_node[index4D(z,y,x,i)]-D3.Q7[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node[index4D(z,y,x,i)] << " vs "<< D3.Q7[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=8;
				if( abs(outter_node[index4D(z,y,x,i)]-D3.Q8[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node[index4D(z,y,x,i)] << " vs "<< D3.Q8[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=9;
				if( abs(outter_node[index4D(z,y,x,i)]-D3.Q9[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node[index4D(z,y,x,i)] << " vs "<< D3.Q9[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=10;
				if( abs(outter_node[index4D(z,y,x,i)]-D3.Q10[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node[index4D(z,y,x,i)] << " vs "<< D3.Q10[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=11;
				if( abs(outter_node[index4D(z,y,x,i)]-D3.Q11[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node[index4D(z,y,x,i)] << " vs "<< D3.Q11[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=12;
				if( abs(outter_node[index4D(z,y,x,i)]-D3.Q12[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node[index4D(z,y,x,i)] << " vs "<< D3.Q12[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=13;
				if( abs(outter_node[index4D(z,y,x,i)]-D3.Q13[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node[index4D(z,y,x,i)] << " vs "<< D3.Q13[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=14;
				if( abs(outter_node[index4D(z,y,x,i)]-D3.Q14[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node[index4D(z,y,x,i)] << " vs "<< D3.Q14[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=15;
				if( abs(outter_node[index4D(z,y,x,i)]-D3.Q15[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node[index4D(z,y,x,i)] << " vs "<< D3.Q15[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=16;
				if( abs(outter_node[index4D(z,y,x,i)]-D3.Q16[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node[index4D(z,y,x,i)] << " vs "<< D3.Q16[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=17;
				if( abs(outter_node[index4D(z,y,x,i)]-D3.Q17[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node[index4D(z,y,x,i)] << " vs "<< D3.Q17[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=18;
				if( abs(outter_node[index4D(z,y,x,i)]-D3.Q18[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node[index4D(z,y,x,i)] << " vs "<< D3.Q18[index(z,y,x)] <<endl;
					exit(-1000);
				}

			}
		}
	}
	cout << "nodes ok" << endl;
}

void LBM::compare_nodes_hlp(FLOATING *outter_node_hlp){
	int x,y,z,i;

	for ( z =0 ; z < lz ; ++z){
		for (y = 0; y < ly ; ++y){
			for ( x = 0 ; x< lx ; ++x){

				i=0;
				if( abs(outter_node_hlp[index4D(z,y,x,i)]-D3_hlp.Q0[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node_hlp[index4D(z,y,x,i)] << " vs "<< D3_hlp.Q0[index(z,y,x)] <<endl;

					exit(-1000);
				}

				i=1;
				if( abs(outter_node_hlp[index4D(z,y,x,i)]-D3_hlp.Q1[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node_hlp[index4D(z,y,x,i)] << " vs "<< D3_hlp.Q1[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=2;
				if( abs(outter_node_hlp[index4D(z,y,x,i)]-D3_hlp.Q2[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node_hlp[index4D(z,y,x,i)] << " vs "<< D3_hlp.Q2[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=3;
				if( abs(outter_node_hlp[index4D(z,y,x,i)]-D3_hlp.Q3[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node_hlp[index4D(z,y,x,i)] << " vs "<< D3_hlp.Q3[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=4;
				if( abs(outter_node_hlp[index4D(z,y,x,i)]-D3_hlp.Q4[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node_hlp[index4D(z,y,x,i)] << " vs "<< D3_hlp.Q4[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=5;
				if( abs(outter_node_hlp[index4D(z,y,x,i)]-D3_hlp.Q5[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node_hlp[index4D(z,y,x,i)] << " vs "<< D3_hlp.Q5[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=6;
				if( abs(outter_node_hlp[index4D(z,y,x,i)]-D3_hlp.Q6[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node_hlp[index4D(z,y,x,i)] << " vs "<< D3_hlp.Q6[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=7;
				if( abs(outter_node_hlp[index4D(z,y,x,i)]-D3_hlp.Q7[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node_hlp[index4D(z,y,x,i)] << " vs "<< D3_hlp.Q7[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=8;
				if( abs(outter_node_hlp[index4D(z,y,x,i)]-D3_hlp.Q8[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node_hlp[index4D(z,y,x,i)] << " vs "<< D3_hlp.Q8[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=9;
				if( abs(outter_node_hlp[index4D(z,y,x,i)]-D3_hlp.Q9[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node_hlp[index4D(z,y,x,i)] << " vs "<< D3_hlp.Q9[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=10;
				if( abs(outter_node_hlp[index4D(z,y,x,i)]-D3_hlp.Q10[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node_hlp[index4D(z,y,x,i)] << " vs "<< D3_hlp.Q10[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=11;
				if( abs(outter_node_hlp[index4D(z,y,x,i)]-D3_hlp.Q11[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node_hlp[index4D(z,y,x,i)] << " vs "<< D3_hlp.Q11[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=12;
				if( abs(outter_node_hlp[index4D(z,y,x,i)]-D3_hlp.Q12[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node_hlp[index4D(z,y,x,i)] << " vs "<< D3_hlp.Q12[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=13;
				if( abs(outter_node_hlp[index4D(z,y,x,i)]-D3_hlp.Q13[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node_hlp[index4D(z,y,x,i)] << " vs "<< D3_hlp.Q13[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=14;
				if( abs(outter_node_hlp[index4D(z,y,x,i)]-D3_hlp.Q14[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node_hlp[index4D(z,y,x,i)] << " vs "<< D3_hlp.Q14[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=15;
				if( abs(outter_node_hlp[index4D(z,y,x,i)]-D3_hlp.Q15[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node_hlp[index4D(z,y,x,i)] << " vs "<< D3_hlp.Q15[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=16;
				if( abs(outter_node_hlp[index4D(z,y,x,i)]-D3_hlp.Q16[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node_hlp[index4D(z,y,x,i)] << " vs "<< D3_hlp.Q16[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=17;
				if( abs(outter_node_hlp[index4D(z,y,x,i)]-D3_hlp.Q17[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node_hlp[index4D(z,y,x,i)] << " vs "<< D3_hlp.Q17[index(z,y,x)] <<endl;
					exit(-1000);
				}

				i=18;
				if( abs(outter_node_hlp[index4D(z,y,x,i)]-D3_hlp.Q18[index(z,y,x)])>0.00001){
					cout << "node miss-match @ x:" << x << " y:" << y << " z:" << z << " i:" <<i <<  " :" << outter_node_hlp[index4D(z,y,x,i)] << " vs "<< D3.Q18[index(z,y,x)] <<endl;
					exit(-1000);
				}

			}
		}
	}
	cout << "n_hlp ok" << endl;
}
