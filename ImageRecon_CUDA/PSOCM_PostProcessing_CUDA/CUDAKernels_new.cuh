/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
* Institute:   University of Illinois Urbana Champaign                            *
* School:      Engineering                                                        *
* Authors:     Rishee Iyer														  *
* Emails:      rriyer3                                                            *
* Project:     Photron FastCAM                                                    *
* Advisor:     Prof. Stephen Boppart                                              *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <windows.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <cufft.h>
#include <cublas_v2.h>
#include <cuda_profiler_api.h>
#include <cassert>

#define PI 3.14159265f

#define DLLEXPORT extern "C" __declspec(dllexport)
#define BLOCK_SIZE 16
#define PRINTSTUFFEH TRUE
#define LOGFILENAME "C:\\PSOCM Control Software\\GPULogFile.txt"
#define IMFACT 4096.0f
#define IMFACT_DB 512.0f
#define COSRATIO 0.92f
#define WRITERAWOUT FALSE

extern uint8_t *gpuRawIm, *gpuRawBack;
extern cufftComplex* gpuIn, * gpuBack;
extern cufftComplex *gpuComplexIm, *gpuFFTIm, *gpuOUTIm1, *gpuOUTIm2;
extern float *gpuWin1, *gpuWin2, *gpuOUTabs1, *gpuOUTabs2, *gRad1, *gRad2;
extern long *gnumX, *gnumY, *gnumx, *gnumy, *gCent1X, *gCent1Y, *gCent2X, *gCent2Y;
extern cufftHandle cuPlan2DIN, cuPlan2DOUT1, cuPlan2DOUT2;
extern uint8_t* gDFEh;

/* Addressing x*numY + y */

__global__ void GetInput(uint16_t *In8, uint16_t *gpuBack, cufftComplex *InC, long *numX, long *numY)
{
	long xid = blockIdx.x * blockDim.x + threadIdx.x;
	long yid = blockIdx.y * blockDim.y + threadIdx.y;
	InC[xid + yid* *numX] = { (float)(In8[xid**numY + yid]) - (float)(gpuBack[xid**numY + yid]), 0.0f };
}

__global__ void GetInput_12bit(uint8_t* In8, cufftComplex* InC)
{
	long xid = blockIdx.x * blockDim.x + threadIdx.x;
	InC[2 * (xid - 1)] = { (float)((((uint16_t)(In8[3 * xid] & 0xFFu)) << 4u) + (((uint16_t)(In8[3 * xid + 1] & 0xF0u)) >> 4u)), 0.0f };
	InC[2 * (xid - 1) + 1] = { (float)((((uint16_t)(In8[3 * xid + 1] & 0x0Fu)) << 8u) + (((uint16_t)(In8[3 * xid + 2] & 0xFFu)) << 0u)), 0.0f };
}

__global__ void BGSub(cufftComplex* InC, cufftComplex* gpuBack)
{
	long xid = blockIdx.x * blockDim.x + threadIdx.x;
	InC[xid] = { (InC[xid].x - gpuBack[xid].x)/ IMFACT, 0.0f };
}

__global__ void GenerateWin(float *Window, float *Radius, long *numx, long *numy)
{
	long xid = blockIdx.x * blockDim.x + threadIdx.x;
	long yid = blockIdx.y * blockDim.y + threadIdx.y;
	float Condish = sqrtf((float)((xid - *numx / 2) * (xid - *numx / 2) +
		(yid - *numy / 2) * (yid - *numy / 2)));
	if (Condish <= COSRATIO * *Radius)
	{
		Window[xid + yid * *numx] = 1.0f;
	}
	else if (Condish > COSRATIO** Radius && Condish <= *Radius)
	{
		Window[xid + yid * *numx] = cosf(( Condish - COSRATIO ** Radius) / ((1 - COSRATIO) ** Radius) * PI * 0.5f);
	}
	else
	{
		Window[xid + yid * *numx] = 0.0f;
	}
}

__global__ void ApplyWin_shift(cufftComplex *Input, float *Win, cufftComplex *Output, 
	long *numX, long *numY,	long *numx, long *numy, long *CenterX, long *CenterY, uint8_t *DFEh)
{
	long XOG = (blockIdx.x * blockDim.x + threadIdx.x);
	long YOG = (blockIdx.y * blockDim.y + threadIdx.y);
	
	long xin = (XOG + *numX + (*CenterX)) % *numX;
	long yin = (YOG + *numY + (*CenterY)) % *numY;
	//long xin = (XOG + *numX - (*CenterX - 1)) % *numX;
	//long yin = (YOG + *numY - (*CenterY - 1)) % *numY;

	long xwin = (XOG) % *numx;
	long ywin = (YOG) % *numy;

	long xout = (XOG < *numx / 2) ? XOG + *numx / 2 : XOG - *numx / 2;
	long yout = (YOG < *numy / 2) ? YOG + *numy / 2 : YOG - *numy / 2;

	Output[xout + yout * *numx] = { Input[xin + yin * *numX].x * Win[xwin + ywin * *numx],
		Input[xin + yin * *numX].y * Win[xwin + ywin * *numx] };

	if (*DFEh && ((xout == 0 || xout == 1 || xout == 2 || xout == *numx - 1 || xout == *numx - 2 || xout == *numx - 3) 
		&& (yout == 0 || yout == 1 || yout == 2 || yout == *numy -1 || yout == *numy - 2 || yout == *numy - 3)))
		Output[xout + yout * *numx] = { 0.0f , 0.0f };
}

__global__ void AbsOutput_Raw(cufftComplex *Input, float *Output)
{
	long idx = (blockIdx.x * blockDim.x + threadIdx.x);

	Output[idx] = cuCabsf(Input[idx]);
}

__global__ void AbsOutput(cufftComplex *Input, float *Output)
{
	long idx = (blockIdx.x * blockDim.x + threadIdx.x);

	Output[idx] = cuCabsf(Input[idx]) / IMFACT;
}

__global__ void PhaseOutput(cufftComplex* Input, float* Output, float* Thresh)
{
	long idx = (blockIdx.x * blockDim.x + threadIdx.x);

	Output[idx] = ((cuCabsf(Input[idx]) / IMFACT) >= *Thresh) ? atan2f(Input[idx].y, Input[idx].x) : -10.0f;
}

__global__ void AbsOutput_dB(cufftComplex *Input, float *Output, long *numX, long *numY)
{
	long xin = (blockIdx.x * blockDim.x + threadIdx.x);
	long yin = (blockIdx.y * blockDim.y + threadIdx.y);
	long xout = (xin < *numX / 2) ? xin + *numX / 2 : xin - *numX / 2;
	long yout = (yin < *numY / 2) ? yin + *numY / 2 : yin - *numY / 2;
	
	Output[xout + yout * *numX] = 20.0f * log10f(cuCabsf(Input[xin + yin* *numX]) / IMFACT_DB);
}

__global__ void PhaseOutput_dB(cufftComplex *Input, float *Output, long *numX, long *numY)
{
	long xin = (blockIdx.x * blockDim.x + threadIdx.x);
	long yin = (blockIdx.y * blockDim.y + threadIdx.y);
	long xout = (xin < *numX / 2) ? xin + *numX / 2 : xin - *numX / 2;
	long yout = (yin < *numY / 2) ? yin + *numY / 2 : yin - *numY / 2;

	Output[xout + yout * *numX] = atan2f(Input[xin + yin * *numX].y, Input[xin + yin * *numX].x);
}
