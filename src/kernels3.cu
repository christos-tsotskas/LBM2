/*
 * Copyright 2008, Karen Hains, UWA . All rights reserved.
 *
 * NOTICE TO USER:
 *
 * This source code is subject to NVIDIA ownership rights under U.S. and
 * international Copyright laws. Users and possessors of this source code
 * are hereby granted a nonexclusive, royalty-free license to use this code
 * in individual and commercial software.
 *
 * WE MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE
 * CODE FOR ANY PURPOSE. IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR
 * IMPLIED WARRANTY OF ANY KIND.
 */

#ifndef _REDUCE_KERNELS_
#define _REDUCE_KERNELS_

#include "global_defines.cuh"
#include <cstdio>
#include <cmath>

///////////////////////////////////////////////////////////
// Simple Hello World kernel
// @param gpu_odata output data in global memory
///////////////////////////////////////////////////////////

#define SUM
//#define MINIMUM
//#define MAXIMUM


#ifndef MIN
#define MIN(x,y) ((x < y) ? x : y)
#endif

#ifndef MAX
#define MAX(x,y) ((x > y) ? x : y)
#endif

// Utility class used to avoid linker errors with extern
// unsized shared memory arrays with templated type
template<class T>
struct SharedMemory
{
	__device__ inline operator       T*()
	{
		extern __shared__ int __smem[];
		return (T*)__smem;
	}

	__device__ inline operator const T*() const
    																		{
		extern __shared__ int __smem[];
		return (T*)__smem;
    																		}
};

// specialize for FLOATING to avoid unaligned memory
// access compile errors

#ifdef USING_DOUBLE
template<>
struct SharedMemory<double>
{
	__device__ inline operator       double*()
    														{
		extern __shared__ double __smem_d[];
		return (double*)__smem_d;
    														}

	__device__ inline operator const double*() const
    														{
		extern __shared__ double __smem_d[];
		return (double*)__smem_d;
    														}
};
#endif //USING_DOUBLE

#ifndef USING_DOUBLE
template<>
struct SharedMemory<float>
{
	__device__ inline operator       float*()
    														{
		extern __shared__ float __smem_d[];
		return (float*)__smem_d;
    														}

	__device__ inline operator const float*() const
    														{
		extern __shared__ float __smem_d[];
		return (float*)__smem_d;
    														}
};
#endif //USING_FLOAT

/*
    Parallel sum reduction using shared memory
    - takes log(n) steps for n input elements
    - uses n threads
    - only works for power-of-2 arrays
 */

/* This reduction interleaves which threads are active by using the modulo
   operator.  This operator is very expensive on GPUs, and the interleaved
   inactivity means that no whole warps are active, which is also very
   inefficient */

unsigned int nextPow2( unsigned int x ) {
	--x;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return ++x;
}

void getNumBlocksAndThreads(int whichKernel, int n, int maxBlocks, int maxThreads, int &blocks, int &threads)
{
	//get device capability, to avoid block/grid size excceed the upbound
	cudaDeviceProp prop;
	int device;
	cudaGetDevice(&device);
	cudaGetDeviceProperties(&prop, device);

	if (whichKernel < 3)
	{
		threads = (n < maxThreads) ? nextPow2(n) : maxThreads;
		blocks = (n + threads - 1) / threads;
	}
	else
	{
		threads = (n < maxThreads*2) ? nextPow2((n + 1)/ 2) : maxThreads;
		blocks = (n + (threads * 2 - 1)) / (threads * 2);
	}

	if ((float)threads*blocks > (float)prop.maxGridSize[0] * prop.maxThreadsPerBlock)
	{
		printf("n is too large, please choose a smaller number!\n");
	}

	if (blocks > prop.maxGridSize[0])
	{
		printf("Grid size <%d> excceeds the device capability <%d>, set block size as %d (original %d)\n",
				blocks, prop.maxGridSize[0], threads*2, threads);

		blocks /= 2;
		threads *= 2;
	}

	if (whichKernel == 6)
	{
		blocks = MIN(maxBlocks, blocks);
	}
#ifdef REPORT
	printf("CUDA Kernels will be launched with:\n");
	printf("\tnumber of blocks=%d\n", blocks);
	printf("\tnumber of threads=%d\n", threads);
#endif //REPORT

}

template <class T>
__global__ void
reduce0(T *g_idata, T *g_odata, unsigned int n)
{
	T *sdata = SharedMemory<T>();

	// load shared mem
	unsigned int tid = threadIdx.x;
	unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;

#ifdef SUM
	sdata[tid] = (i < n) ? g_idata[i] : 0;
#endif
#ifdef MINIMUM
	sdata[tid] = (i < n) ? g_idata[i] : 0;
#endif
#ifdef MAXIMUM
	sdata[tid] = (i < n) ? g_idata[i] : 0;
#endif

	__syncthreads();

	// do reduction in shared mem
	for(unsigned int s=1; s < blockDim.x; s *= 2) {
		// modulo arithmetic is slow!
		if ((tid % (2*s)) == 0) {
#ifdef SUM
			sdata[tid] += sdata[tid + s];
#endif
#ifdef MINIMUM
			sdata[tid] = min(sdata[tid] ,  sdata[tid + s]);
#endif
#ifdef MAXIMUM
			sdata[tid] = max(sdata[tid] ,  sdata[tid + s]);
#endif

		}
		__syncthreads();
	}

	// write result for this block to global mem
	if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

template <class T, unsigned int blockSize, bool nIsPow2>
__global__ void
reduce61(T *g_idata, T *g_odata, unsigned int n)
{
	//original
	T *sdata = SharedMemory<T>();

	// perform first level of reduction,
	// reading from global memory, writing to shared memory
	unsigned int tid = threadIdx.x;
	unsigned int i = blockIdx.x*blockSize*2 + threadIdx.x;
	unsigned int gridSize = blockSize*2*gridDim.x;

	T mySum = 0;

	// we reduce multiple elements per thread.  The number is determined by the
	// number of active thread blocks (via gridDim).  More blocks will result
	// in a larger gridSize and therefore fewer elements per thread
	while (i < n)
	{
		mySum += g_idata[i];

		// ensure we don't read out of bounds -- this is optimized away for powerOf2 sized arrays
		if (nIsPow2 || i + blockSize < n)
			mySum += g_idata[i+blockSize];

		i += gridSize;
	}

	// each thread puts its local sum into shared memory
	sdata[tid] = mySum;
	__syncthreads();


	// do reduction in shared mem
	if (blockSize >= 512)
	{
		if (tid < 256)
		{
			sdata[tid] = mySum = mySum + sdata[tid + 256];
		}

		__syncthreads();
	}

	if (blockSize >= 256)
	{
		if (tid < 128)
		{
			sdata[tid] = mySum = mySum + sdata[tid + 128];
		}

		__syncthreads();
	}

	if (blockSize >= 128)
	{
		if (tid <  64)
		{
			sdata[tid] = mySum = mySum + sdata[tid +  64];
		}

		__syncthreads();
	}

	if (tid < 32)
	{
		// now that we are using warp-synchronous programming (below)
		// we need to declare our shared memory volatile so that the compiler
		// doesn't reorder stores to it and induce incorrect behavior.
		volatile T *smem = sdata;

		if (blockSize >=  64)
		{
			smem[tid] = mySum = mySum + smem[tid + 32];
		}

		if (blockSize >=  32)
		{
			smem[tid] = mySum = mySum + smem[tid + 16];
		}

		if (blockSize >=  16)
		{
			smem[tid] = mySum = mySum + smem[tid +  8];
		}

		if (blockSize >=   8)
		{
			smem[tid] = mySum = mySum + smem[tid +  4];
		}

		if (blockSize >=   4)
		{
			smem[tid] = mySum = mySum + smem[tid +  2];
		}

		if (blockSize >=   2)
		{
			smem[tid] = mySum = mySum + smem[tid +  1];
		}
	}

	// write result for this block to global mem
	if (tid == 0)
		g_odata[blockIdx.x] = sdata[0];
}



template <class T, unsigned int blockSize, bool nIsPow2>
__global__ void
reduce6(T *g_idata, T *g_odata, unsigned int n)
{
	T *sdata = SharedMemory<T>();

	// perform first level of reduction,
	// reading from global memory, writing to shared memory
	unsigned int tid = threadIdx.x;
	unsigned int i = blockIdx.x*blockSize*2 + threadIdx.x;
	unsigned int gridSize = blockSize*2*gridDim.x;

#ifdef SUM
	T mySum = 0;
#endif
#ifdef MINIMUM
	T myMin = g_idata[i];//g_idata[0];
#endif
#ifdef MAXIMUM
	T myMax = g_idata[i];//g_idata[0];
#endif
	// we reduce multiple elements per thread.  The number is determined by the
	// number of active thread blocks (via gridDim).  More blocks will result
	// in a larger gridSize and therefore fewer elements per thread
	while (i < n)
	{
#ifdef SUM
		mySum = mySum+ g_idata[i];
#endif
#ifdef MINIMUM
		myMin=min(myMin, g_idata[i]);
#endif
#ifdef MAXIMUM
		myMax=max(myMax, g_idata[i]);
#endif

		// ensure we don't read out of bounds -- this is optimized away for powerOf2 sized arrays
		if (nIsPow2 || i + blockSize < n){
#ifdef SUM
			mySum = mySum + g_idata[i+blockSize];
#endif
#ifdef MINIMUM
			myMin=min(myMin, g_idata[i+blockSize]);
#endif
#ifdef MAXIMUM
			myMax=max(myMax, g_idata[i+blockSize]);
#endif
		}
		i += gridSize;
	}

	// each thread puts its local sum into shared memory
#ifdef SUM
	sdata[tid] = mySum;
#endif
#ifdef MINIMUM
	sdata[tid] = myMin;
#endif
#ifdef MAXIMUM
	sdata[tid] = myMax;
#endif
	__syncthreads();


	// do reduction in shared mem
#ifdef SUM
	if (blockSize >= 512) { if (tid < 256) { sdata[tid] = mySum = mySum + sdata[tid + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (tid < 128) { sdata[tid] = mySum = mySum + sdata[tid + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (tid <  64) { sdata[tid] = mySum = mySum + sdata[tid +  64]; } __syncthreads(); }
#endif
#ifdef MINIMUM
	if (blockSize >= 512) { if (tid < 256) { sdata[tid] = myMin = min(myMin , sdata[tid + 256]); } __syncthreads(); }
	if (blockSize >= 256) { if (tid < 128) { sdata[tid] = myMin = min(myMin , sdata[tid + 128]); } __syncthreads(); }
	if (blockSize >= 128) { if (tid <  64) { sdata[tid] = myMin = min(myMin , sdata[tid +  64]); } __syncthreads(); }
#endif
#ifdef MAXIMUM
	if (blockSize >= 512) { if (tid < 256) { sdata[tid] = myMax = max(myMax , sdata[tid + 256]); } __syncthreads(); }
	if (blockSize >= 256) { if (tid < 128) { sdata[tid] = myMax = max(myMax , sdata[tid + 128]); } __syncthreads(); }
	if (blockSize >= 128) { if (tid <  64) { sdata[tid] = myMax = max(myMax , sdata[tid +  64]); } __syncthreads(); }
#endif
	if (tid < 32)
	{
		// now that we are using warp-synchronous programming (below)
		// we need to declare our shared memory volatile so that the compiler
		// doesn't reorder stores to it and induce incorrect behavior.
		volatile T* smem = sdata;
#ifdef SUM
		if (blockSize >=  64) { smem[tid] = mySum = mySum + smem[tid + 32]; }
		if (blockSize >=  32) { smem[tid] = mySum = mySum + smem[tid + 16]; }
		if (blockSize >=  16) { smem[tid] = mySum = mySum + smem[tid +  8]; }
		if (blockSize >=   8) { smem[tid] = mySum = mySum + smem[tid +  4]; }
		if (blockSize >=   4) { smem[tid] = mySum = mySum + smem[tid +  2]; }
		if (blockSize >=   2) { smem[tid] = mySum = mySum + smem[tid +  1]; }
#endif
#ifdef MINIMUM
		if (blockSize >=  64) { smem[tid] = myMin = min(myMin , smem[tid + 32]); }
		if (blockSize >=  32) { smem[tid] = myMin = min(myMin , smem[tid + 16]); }
		if (blockSize >=  16) { smem[tid] = myMin = min(myMin , smem[tid +  8]); }
		if (blockSize >=   8) { smem[tid] = myMin = min(myMin , smem[tid +  4]); }
		if (blockSize >=   4) { smem[tid] = myMin = min(myMin , smem[tid +  2]); }
		if (blockSize >=   2) { smem[tid] = myMin = min(myMin , smem[tid +  1]); }
#endif
#ifdef MAXIMUM
		if (blockSize >=  64) { smem[tid] = myMax = max(myMax , smem[tid + 32]); }
		if (blockSize >=  32) { smem[tid] = myMax = max(myMax , smem[tid + 16]); }
		if (blockSize >=  16) { smem[tid] = myMax = max(myMax , smem[tid +  8]); }
		if (blockSize >=   8) { smem[tid] = myMax = max(myMax , smem[tid +  4]); }
		if (blockSize >=   4) { smem[tid] = myMax = max(myMax , smem[tid +  2]); }
		if (blockSize >=   2) { smem[tid] = myMax = max(myMax , smem[tid +  1]); }
#endif
	}

	// write result for this block to global mem
	if (tid == 0)
		g_odata[blockIdx.x] = sdata[0];
}






bool isPow2(unsigned int x)
{
	return ((x&(x-1))==0);
}



////////////////////////////////////////////////////////////////////////////////
// Wrapper function for kernel launch
////////////////////////////////////////////////////////////////////////////////
template <class T>
void
reduce(int size, int threads, int blocks,
		int whichKernel, T *d_idata, T *d_odata)
{
	dim3 dimBlock(threads, 1, 1);
	dim3 dimGrid(blocks, 1, 1);

	// when there is only one warp per block, we need to allocate two warps
	// worth of shared memory so that we don't index shared memory out of bounds
	int smemSize = (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);

	// choose which of the optimized versions of reduction to launch
	switch (whichKernel)
	{
	case 0:
		reduce0<T><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
		break;
	case 6:
	default:
		if (isPow2(size))
		{
			switch (threads)
			{
			case 512:
				reduce6<T, 512, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
			case 256:
				reduce6<T, 256, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
			case 128:
				reduce6<T, 128, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
			case 64:
				reduce6<T,  64, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
			case 32:
				reduce6<T,  32, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
			case 16:
				reduce6<T,  16, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
			case  8:
				reduce6<T,   8, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
			case  4:
				reduce6<T,   4, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
			case  2:
				reduce6<T,   2, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
			case  1:
				reduce6<T,   1, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
			}
		}
		else
		{
			switch (threads)
			{
			case 512:
				reduce6<T, 512, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
			case 256:
				reduce6<T, 256, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
			case 128:
				reduce6<T, 128, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
			case 64:
				reduce6<T,  64, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
			case 32:
				reduce6<T,  32, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
			case 16:
				reduce6<T,  16, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
			case  8:
				reduce6<T,   8, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
			case  4:
				reduce6<T,   4, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
			case  2:
				reduce6<T,   2, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
			case  1:
				reduce6<T,   1, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
			}
		}
		break;
	}
}




// Instantiate the reduction function for 3 types
template void
reduce<int>(int size, int threads, int blocks,
		int whichKernel, int *d_idata, int *d_odata);

template void
reduce<float>(int size, int threads, int blocks,
		int whichKernel, float *d_idata, float *d_odata);

template void
reduce<double>(int size, int threads, int blocks,
		int whichKernel, double *d_idata, double *d_odata);


//void getNumBlocksAndThreads(int whichKernel, int n, int maxBlocks, int maxThreads, int &blocks, int &threads);



FLOATING reduce_sum(FLOATING *d_idata, const int big_array_length){
	int whichKernel = 6;
	int maxThreads = 256;  // number of threads per block
	int maxBlocks = min(33554432/maxThreads , 65535);
	int numBlocks, numThreads;
	int cpuFinalThreshold=1;


	getNumBlocksAndThreads(whichKernel, big_array_length, maxBlocks, maxThreads, numBlocks, numThreads);

	if (numBlocks == 1) cpuFinalThreshold = 1;
	FLOATING* d_odata = NULL;
	cudaMalloc((void**) &d_odata, numBlocks*sizeof(FLOATING));
	FLOATING* h_odata = (FLOATING*) malloc(numBlocks*sizeof(FLOATING));

	reduce<FLOATING>(big_array_length, numThreads, numBlocks, whichKernel, d_idata, d_odata);
	cudaDeviceSynchronize();

	FLOATING gpu_result_sum=0.0;
	bool needReadBack = true;

	int s=numBlocks;
	int kernel = whichKernel;

	while (s > cpuFinalThreshold)
	{
		int threads = 0, blocks = 0;
		getNumBlocksAndThreads(kernel, s, maxBlocks, maxThreads, blocks, threads);

		reduce<FLOATING>(s, threads, blocks, kernel, d_odata, d_odata);

		if (kernel < 3)
		{
			s = (s + threads - 1) / threads;
		}
		else
		{
			s = (s + (threads*2-1)) / (threads*2);
		}
	}

	if (s > 1)
	{
		// copy result from device to host
		cudaMemcpy(h_odata, d_odata, s * sizeof(FLOATING), cudaMemcpyDeviceToHost);

		for (int i=0; i < s; i++)
		{
			gpu_result_sum += h_odata[i];
		}

		needReadBack = false;
	}

	if (needReadBack)
	{
		// copy final sum from device to host
		cudaMemcpy(&gpu_result_sum, d_odata, sizeof(FLOATING), cudaMemcpyDeviceToHost);
	}

	cudaDeviceSynchronize();


#ifdef REPORT
	printf ( " [universal]gpu_result_sum: %30.10f \n", gpu_result_sum);
#endif //REPORT
	cudaFree(d_odata);
	return gpu_result_sum;
}


#endif // #ifndef _HELLOWORLD_KERNEL_H_
