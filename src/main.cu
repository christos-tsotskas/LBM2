//============================================================================
// Name        : LBM2.cpp
// Author      : Christos Tsotskas
// Version     :
// Copyright   : L-GPL
// Description : CUDA Lattice Boltzmann in C++
//				(code adapted from the original F90 version by Prof. L.Djenidi)
//
// Input(2 files): LBM2_configuration.txt, LBM2_geometry.txt
// Output(4 files): geometry, convergence, report, data.
// The output files start with the "LBM2_", followed by the case name and their type (one of the 4 above).
// In addition, the report and data include the iteration counter
//============================================================================

#include "global_defines.cuh"

using namespace std;

	//	const FLOATING nu=0.0175, r_small=6.67897, reynolds=195.732, s=23.7849; //original
	//	const FLOATING nu=0.0175, r_small=9.7, reynolds=25.73, s=21; //case A
	//	const FLOATING nu=0.0175, r_small=5.7, reynolds=195.7, s=24; //case B
	//	const FLOATING nu=0.0175, r_small=11, reynolds=195.7, s=23; //case C

int main() {
	cuda_device_querry();

	int lx=1000, ly=1000,  lz=1000, n_of_densities=100;

	read_external_geometry_file_specification_for_LBM(lx, ly, lz, n_of_densities, "LBM2_geometry.txt");

	LBM myLBM( lx, ly, lz, 1.0,  1.0/3.0 , 1.0/18.0, 1.0/36.0, 1.0/3.0);

	myLBM.compute_domain();

	myLBM.export_solution();

	return 0;
}
