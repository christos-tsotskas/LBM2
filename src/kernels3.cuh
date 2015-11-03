#include "global_defines.cuh"

void getNumBlocksAndThreads(int whichKernel, int n, int maxBlocks, int maxThreads, int &blocks, int &threads);
template <class T>
void
reduce(int size, int threads, int blocks,
		int whichKernel, T *d_idata, T *d_odata);
FLOATING cuda_reduce_max(FLOATING *big_array, int big_array_size);
FLOATING reduce_sum(FLOATING *d_idata, const int big_array_length);
