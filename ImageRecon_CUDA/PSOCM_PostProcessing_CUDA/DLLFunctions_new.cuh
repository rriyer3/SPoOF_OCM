/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
* Institute:   University of Illinois Urbana Champaign                            *
* School:      Engineering                                                        *
* Authors:     Rishee Iyer											     		  *
* Emails:      rriyer3                                                            *
* Project:     Photron FastCAM                                                    *
* Advisor:     Prof. Stephen Boppart                                              *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "CUDAKernels_new.cuh"

uint8_t* gpuRawIm, * gpuRawBack;
cufftComplex* gpuIn, * gpuBack;
cufftComplex* gpuComplexIm, * gpuFFTIm, * gpuOUTIm1, * gpuOUTIm2;
float* gpuWin1, * gpuWin2, * gpuOUTabs1, * gpuOUTabs2, * gRad1, * gRad2;
long* gnumX, * gnumY, * gnumx, * gnumy, * gCent1X, * gCent1Y, * gCent2X, * gCent2Y;
cufftHandle cuPlan2DIN, cuPlan2DOUT1, cuPlan2DOUT2;
uint8_t* gDFEh;

DLLEXPORT cudaError_t FFPSOCM_InitializeEverything(long numX, long numY, long numx, long numy, uint8_t *Back, float* Rad1, float* Rad2)
{
	if (PRINTSTUFFEH) printf("Initializing GPU and memory\n");
	
	long *nnumX = (long*)malloc(sizeof(long));
	long* nnumY = (long*)malloc(sizeof(long));
	*nnumX = numX;
	*nnumY = numY;

	long* nnumx = (long*)malloc(sizeof(long));
	long* nnumy = (long*)malloc(sizeof(long));
	*nnumx = numx;
	*nnumy = numy;

	long sizeXY = numX * numY;
	long sizexy = numx * numy;

	cudaSetDevice(0);
	cudaMalloc((void**)&gpuRawIm, sizeof(uint8_t) * sizeXY * 3 / 2);
	cudaMalloc((void**)&gpuRawBack, sizeof(uint8_t) * sizeXY * 3 / 2);
	cudaMalloc((void**)&gpuIn, sizeof(cufftComplex) * sizeXY);
	cudaMalloc((void**)&gpuBack, sizeof(cufftComplex) * sizeXY);
	cudaMalloc((void**)&gpuComplexIm, sizeof(cufftComplex)*sizeXY);
	cudaMalloc((void**)&gpuFFTIm, sizeof(cufftComplex)*sizeXY);
	
	if (PRINTSTUFFEH) printf("Allocating memory set 1: %s\n", cudaGetErrorString(cudaPeekAtLastError()));

	cudaMalloc((void**)&gpuOUTIm1, sizeof(cufftComplex)*sizexy);
	cudaMalloc((void**)&gpuOUTIm2, sizeof(cufftComplex)*sizexy);
	cudaMalloc((void**)&gpuWin1, sizeof(float)*sizexy);
	cudaMalloc((void**)&gpuWin2, sizeof(float)*sizexy);
	cudaMalloc((void**)&gpuOUTabs1, sizeof(float)*sizexy);
	cudaMalloc((void**)&gpuOUTabs2, sizeof(float)*sizexy);
	if (PRINTSTUFFEH) printf("Allocating memory set 2: %s\n", cudaGetErrorString(cudaPeekAtLastError()));

	cudaMalloc((void**)&gRad1, sizeof(float));
	cudaMalloc((void**)&gRad2, sizeof(float));
	cudaMalloc((void**)&gnumX, sizeof(long));
	cudaMalloc((void**)&gnumY, sizeof(long));
	cudaMalloc((void**)&gnumx, sizeof(long));
	cudaMalloc((void**)&gnumy, sizeof(long));
	cudaMalloc((void**)&gCent1X, sizeof(long));
	cudaMalloc((void**)&gCent1Y, sizeof(long));
	cudaMalloc((void**)&gCent2X, sizeof(long));
	cudaMalloc((void**)&gCent2Y, sizeof(long));
	cudaMalloc((void**)&gDFEh, sizeof(uint8_t));
	if (PRINTSTUFFEH) printf("Allocating memory set 3: %s\n", cudaGetErrorString(cudaPeekAtLastError()));

	cufftPlan2d(&cuPlan2DIN, numY, numX, CUFFT_C2C);
	cufftPlan2d(&cuPlan2DOUT1, numy, numx, CUFFT_C2C);
	cufftPlan2d(&cuPlan2DOUT2, numy, numx, CUFFT_C2C);
	if (PRINTSTUFFEH) printf("Planning FFT: %s\n", cudaGetErrorString(cudaPeekAtLastError()));

	cudaMemcpy(gnumX, nnumX, sizeof(long), cudaMemcpyHostToDevice);
	cudaMemcpy(gnumY, nnumY, sizeof(long), cudaMemcpyHostToDevice);
	if (PRINTSTUFFEH) printf("Copying init stuff XY: %s\n", cudaGetErrorString(cudaPeekAtLastError()));

	cudaMemcpy(gnumx, nnumx, sizeof(long), cudaMemcpyHostToDevice);
	cudaMemcpy(gnumy, nnumy, sizeof(long), cudaMemcpyHostToDevice);
	if (PRINTSTUFFEH) printf("Copying init stuff xy: %s\n", cudaGetErrorString(cudaPeekAtLastError()));

	dim3 dimBlock_XY1p5(BLOCK_SIZE * BLOCK_SIZE);
	dim3 dimGrid_XY1p5((numX*numY/2) / dimBlock_XY1p5.x);
	
	cudaMemcpy(gpuRawBack, Back, sizeXY * sizeof(uint8_t) * 3 / 2, cudaMemcpyHostToDevice);
	GetInput_12bit << < dimGrid_XY1p5, dimBlock_XY1p5 >> > (gpuRawBack, gpuBack);

	dim3 dimBlock_xy(BLOCK_SIZE, BLOCK_SIZE);
	dim3 dimGrid_xy((numx) / dimBlock_xy.x, (numy) / dimBlock_xy.y);

	cudaMemcpy(gRad1, Rad1, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(gRad2, Rad2, sizeof(float), cudaMemcpyHostToDevice);

	GenerateWin << <dimGrid_xy, dimBlock_xy >> > (gpuWin1, gRad1, gnumx, gnumy);
	if (PRINTSTUFFEH && cudaPeekAtLastError() != cudaSuccess) printf("Generating Window 1 failed: %s\n", cudaGetErrorString(cudaPeekAtLastError()));

	GenerateWin << <dimGrid_xy, dimBlock_xy >> > (gpuWin2, gRad2, gnumx, gnumy);
	if (PRINTSTUFFEH && cudaPeekAtLastError() != cudaSuccess) printf("Generating Window 2 failed: %s\n", cudaGetErrorString(cudaPeekAtLastError()));

	if (PRINTSTUFFEH) printf("Copying init stuff: %s\n", cudaGetErrorString(cudaPeekAtLastError()));
	free(nnumX);
	free(nnumY);
	free(nnumx);
	free(nnumy);
	if (PRINTSTUFFEH) printf("Anything that happens henceforth happens in the processing step\n");
	return cudaGetLastError();
}

DLLEXPORT cudaError_t FFPSOCM_ProcessFrame(uint8_t *Input, cufftComplex *ImRaw, cufftComplex *Im1, cufftComplex* Im2,
	long numX, long numY, long numx, long numy, long *Cent1X, long *Cent1Y, long *Cent2X, long *Cent2Y,
	float *Rad1, float *Rad2, uint8_t *DFEh)
{
	dim3 dimBlock_XY(BLOCK_SIZE, BLOCK_SIZE);
	dim3 dimGrid_XY((numX) / dimBlock_XY.x, (numY) / dimBlock_XY.y);

	dim3 dimBlock_xy(BLOCK_SIZE, BLOCK_SIZE);
	dim3 dimGrid_xy((numx) / dimBlock_xy.x, (numy) / dimBlock_xy.y);

	dim3 dimBlock_linXY(BLOCK_SIZE*BLOCK_SIZE);
	dim3 dimGrid_linXY((numX * numY) / dimBlock_linXY.x);

	dim3 dimBlock_linxy(BLOCK_SIZE * BLOCK_SIZE);
	dim3 dimGrid_linxy((numx * numy) / dimBlock_linxy.x);

	cudaMemcpy(gCent1X, Cent1X, sizeof(long), cudaMemcpyHostToDevice);
	cudaMemcpy(gCent2X, Cent2X, sizeof(long), cudaMemcpyHostToDevice);
	cudaMemcpy(gCent1Y, Cent1Y, sizeof(long), cudaMemcpyHostToDevice);
	cudaMemcpy(gCent2Y, Cent2Y, sizeof(long), cudaMemcpyHostToDevice);
	cudaMemcpy(gRad1, Rad1, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(gRad2, Rad2, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(gDFEh, DFEh, sizeof(uint8_t), cudaMemcpyHostToDevice);
	if (PRINTSTUFFEH && cudaPeekAtLastError() != cudaSuccess) printf("Copying parameters failed: %s\n", cudaGetErrorString(cudaPeekAtLastError()));

	dim3 dimBlock_XY1p5(BLOCK_SIZE * BLOCK_SIZE);
	dim3 dimGrid_XY1p5((numX * numY / 2) / dimBlock_XY1p5.x);

	cudaMemcpy(gpuRawIm, Input, numX * numY * sizeof(uint8_t) * 3 / 2, cudaMemcpyHostToDevice);
	if (PRINTSTUFFEH && cudaPeekAtLastError() != cudaSuccess) printf("Copying input failed: %s\n", cudaGetErrorString(cudaPeekAtLastError()));

	GetInput_12bit << < dimGrid_XY1p5, dimBlock_XY1p5 >> > (gpuRawIm, gpuComplexIm);
	if (PRINTSTUFFEH && cudaPeekAtLastError() != cudaSuccess) printf("Restructuring input failed: %s\n", cudaGetErrorString(cudaPeekAtLastError()));
	
	BGSub << < dimGrid_linXY, dimBlock_linXY >> > (gpuComplexIm, gpuBack);
	if (PRINTSTUFFEH && cudaPeekAtLastError() != cudaSuccess) printf("Background subtraction failed: %s\n", cudaGetErrorString(cudaPeekAtLastError()));
	
	cudaMemcpy(ImRaw, gpuComplexIm, sizeof(cufftComplex) * numX * numY, cudaMemcpyDeviceToHost);
	if (PRINTSTUFFEH && cudaPeekAtLastError() != cudaSuccess) printf("Copying output 2 failed: %s\n", cudaGetErrorString(cudaPeekAtLastError()));
	
	cufftExecC2C(cuPlan2DIN, gpuComplexIm, gpuFFTIm, CUFFT_FORWARD);
	if (PRINTSTUFFEH && cudaPeekAtLastError() != cudaSuccess) printf("Calculating FFT of raw data failed: %s\n", cudaGetErrorString(cudaPeekAtLastError()));

	ApplyWin_shift << <dimGrid_xy, dimBlock_xy >> > (gpuFFTIm, gpuWin1, gpuOUTIm1, gnumX, gnumY,
		gnumx, gnumy, gCent1X, gCent1Y, gDFEh);
	if (PRINTSTUFFEH && cudaPeekAtLastError() != cudaSuccess) printf("Applying Window 1 failed: %s\n", cudaGetErrorString(cudaPeekAtLastError()));

	ApplyWin_shift << <dimGrid_xy, dimBlock_xy >> > (gpuFFTIm, gpuWin2, gpuOUTIm2, gnumX, gnumY,
		gnumx, gnumy, gCent2X, gCent2Y, gDFEh);
	if (PRINTSTUFFEH && cudaPeekAtLastError() != cudaSuccess) printf("Generating Window 2 failed: %s\n", cudaGetErrorString(cudaPeekAtLastError()));

	cufftExecC2C(cuPlan2DOUT1, gpuOUTIm1, gpuOUTIm1, CUFFT_INVERSE);
	if (PRINTSTUFFEH && cudaPeekAtLastError() != cudaSuccess) printf("Calculating iFFT of P1 failed: %s\n", cudaGetErrorString(cudaPeekAtLastError()));

	cufftExecC2C(cuPlan2DOUT2, gpuOUTIm2, gpuOUTIm2, CUFFT_INVERSE);
	if (PRINTSTUFFEH && cudaPeekAtLastError() != cudaSuccess) printf("Calculating iFFT of P1 failed: %s\n", cudaGetErrorString(cudaPeekAtLastError()));

	cudaMemcpy(Im1, gpuOUTIm1, sizeof(cufftComplex) * numx * numy, cudaMemcpyDeviceToHost);
	if (PRINTSTUFFEH && cudaPeekAtLastError() != cudaSuccess) printf("Copying output 1 failed: %s\n", cudaGetErrorString(cudaPeekAtLastError()));

	cudaMemcpy(Im2, gpuOUTIm2, sizeof(cufftComplex) * numx * numy, cudaMemcpyDeviceToHost);
	if (PRINTSTUFFEH && cudaPeekAtLastError() != cudaSuccess) printf("Copying output 2 failed: %s\n", cudaGetErrorString(cudaPeekAtLastError()));

	return cudaGetLastError();
}

DLLEXPORT cudaError_t FFPSOCM_DestroyEverything()
{

	if (PRINTSTUFFEH) printf("Destroying Everything\n");
	cudaFree(gpuBack);
	cudaFree(gpuRawIm);
	cudaFree(gpuComplexIm);
	cudaFree(gpuFFTIm);
	cudaFree(gpuOUTIm1);
	cudaFree(gpuOUTIm2);
	cudaFree(gpuWin1);
	cudaFree(gpuWin2);
	cudaFree(gpuOUTabs1);
	cudaFree(gpuOUTabs2);
	cudaFree(gRad1);
	cudaFree(gRad2);
	cudaFree(gnumX);
	cudaFree(gnumY);
	cudaFree(gCent1X);
	cudaFree(gCent1Y);
	cudaFree(gCent2X);
	cudaFree(gCent2Y);
	cudaFree(gpuRawBack);
	cudaFree(gpuIn);
	cudaFree(gDFEh);
	cufftDestroy(cuPlan2DIN);
	cufftDestroy(cuPlan2DOUT1);
	cufftDestroy(cuPlan2DOUT2);
	if (PRINTSTUFFEH) printf("Destroyed\n");
	cudaDeviceReset();
	return cudaGetLastError();
}
