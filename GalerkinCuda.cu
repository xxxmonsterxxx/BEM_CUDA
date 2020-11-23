#include "GalerkinCuda.cuh"
#include "GalerkinData.h"

__global__ void GalerkinCuda::ElementInfluenceMatrix(float* infMatr, float* beinfo, float* fright)
{
	uint beNum = blockIdx.x / 3; // each be generate 3 equation
	uint funcNum = blockIdx.x % 3 + 1; // equation number (1..3) for fixed boundary element

	uint beNumLocal = blockIdx.y / 3; // be num which influent on fixed global be
	uint funcNumLocal = blockIdx.y % 3 + 1; // func num of be which influent on fixed global be

	uint termNum = threadIdx.x; // numeric integral term (1..numIntDiscr)

	// TO-DO: CHECK if coeff number is correct???!!!
	uint coeffNumGlobal = blockIdx.x * gridDim.y + blockIdx.y; // global number of coefficient in full influence matrix

	uint beNumInfoK = 8; // shift multiplier = size of beinfo struct

	// global boundary element info
	float xBE = beinfo[beNum * beNumInfoK+2];
	float yBE = beinfo[beNum * beNumInfoK+3];
	float alphaBE = beinfo[beNum * beNumInfoK+6];
	float lngBE = beinfo[beNum * beNumInfoK+7];

	//local boundary element info
	float xBELoc = beinfo[beNumLocal * beNumInfoK+2];
	float yBELoc = beinfo[beNumLocal * beNumInfoK+3];
	float alphaBELoc = beinfo[beNumLocal * beNumInfoK+6];
	float lngBELoc = beinfo[beNumLocal * beNumInfoK+7];

	// info for discret integral
	float discrStep = 2 * lngBE / blockDim.x;
	float xSub = -lngBE + discrStep * termNum;
	float ySub = 0;

	float xSubTransofrmed = 0, ySubTransformed = 0;

	Transform2D(xBE,yBE,alphaBE,
				0,0,0,
				xSub,ySub,
				xSubTransofrmed,ySubTransformed);

	float increment = 0;
	float frightIncrement = 0;

	switch (funcNum) {
		case 1:
			increment = discrStep * f1(xSub, lngBE) * IG(xSubTransofrmed, ySubTransformed, xBELoc, yBELoc, lngBELoc, alphaBELoc, funcNumLocal);
			if (blockIdx.x == blockIdx.y)
				frightIncrement = discrStep * f1(xSub, lngBE) * Problem::InitCondition(xSubTransofrmed, ySubTransformed);
			break;
		case 2:
			increment = discrStep * f2(xSub, lngBE) * IG(xSubTransofrmed, ySubTransformed, xBELoc, yBELoc, lngBELoc, alphaBELoc, funcNumLocal);
			if (blockIdx.x == blockIdx.y)
				frightIncrement = discrStep * f2(xSub, lngBE) * Problem::InitCondition(xSubTransofrmed, ySubTransformed);
			break;
		case 3:
			increment = discrStep * f3(xSub, lngBE) * IG(xSubTransofrmed, ySubTransformed, xBELoc, yBELoc, lngBELoc, alphaBELoc, funcNumLocal);
			if (blockIdx.x == blockIdx.y)
				frightIncrement = discrStep * f3(xSub, lngBE) * Problem::InitCondition(xSubTransofrmed, ySubTransformed);
			break;
		default:
			break;
	}

	atomicAdd(&infMatr[coeffNumGlobal], increment);
	atomicAdd(&fright[blockIdx.x], frightIncrement);
}

using namespace GalerkinMethod;

void GalerkinCuda::CalculateInfMatrix()
{ 
	if (!initialisedData) {
		printf("\nFalse while reading input data");
		return;
	}

	ResetData();

	cudaSetDevice(0);

	// try to parallel maximal effective
	// we have 3N [N - number of boundary elements] equations
	// each equation is a result of numeric integral and is a sum of p_j*k [j=1..N]
	// so we need to calculate each k

	//cudaDeviceProp prop;
	//cudaGetDeviceProperties(&prop, 0);
	//printf("Device is %s\nnumber of blocks %dx%dx%d (each %dx%dx%d) = number of threads %d\n", prop.name,
	//	prop.maxGridSize[0],
	//	prop.maxGridSize[1],
	//	prop.maxGridSize[2],
	//	prop.maxThreadsDim[0],
	//	prop.maxThreadsDim[1],
	//	prop.maxThreadsDim[2],
	//	prop.maxThreadsPerBlock);

	dim3 blockSize = dim3(numIntDiscr, 1, 1); // each cofficient is a summ of numIntDiscr terms
	dim3 gridSize = dim3(beNum*3, beNum*3, 1); // each boundary element have 3 equation which consist of (beNum * 3) coefficients

	// data pointers for a kernel
	float* cudaInfMatr;
	float* cudaBeInfo;
	float* cudaFright;
	cudaMalloc((void**)&cudaInfMatr, infMatrSize * sizeof(float));
	cudaMalloc((void**)&cudaBeInfo, beInfoSize * sizeof(float));
	cudaMalloc((void**)&cudaFright, fRightSize * sizeof(float));
	cudaMemcpy(cudaInfMatr, infMatr, infMatrSize * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(cudaBeInfo, beInfo, beInfoSize * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(cudaFright, fRight, fRightSize * sizeof(float), cudaMemcpyHostToDevice);

	ElementInfluenceMatrix <<< gridSize, blockSize >>> (cudaInfMatr, cudaBeInfo, cudaFright);


	cudaError_t cudaerr = cudaDeviceSynchronize();
	if (cudaerr != cudaSuccess)
		printf("kernel launch failed with error \"%s\".\n",
			cudaGetErrorString(cudaerr));

	cudaEvent_t syncEvent;
	cudaEventCreate(&syncEvent);
	cudaEventRecord(syncEvent, 0);
	cudaEventSynchronize(syncEvent);

	cudaMemcpy(infMatr, cudaInfMatr, infMatrSize * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(fRight, cudaFright, fRightSize * sizeof(float), cudaMemcpyDeviceToHost);
	cudaDeviceReset();
	printf("\nSolved success!");
}

__global__ void GalerkinCuda::NodePotential(float* nodes, float* beinfo, float* coeffs)
{
	uint nodeInd = (blockIdx.x * gridDim.y + blockIdx.y) * 3;
	float x = nodes[nodeInd + 0];
	float y = nodes[nodeInd + 1];


	// local boundary element info
	uint beNumInfoK = 8; // shift multiplier = size of beinfo struct
	uint localBE = threadIdx.x;
	float xBELoc = beinfo[localBE * beNumInfoK + 2];
	float yBELoc = beinfo[localBE * beNumInfoK + 3];
	float alphaBELoc = beinfo[localBE * beNumInfoK + 6];
	float lngBELoc = beinfo[localBE * beNumInfoK + 7];

	uint func = threadIdx.y + 1;

	uint coeffInd = threadIdx.x * 3;

	float increment = 0;

	switch (func) {
		case 1:
			increment = coeffs[coeffInd + 0] * IG(x, y, xBELoc, yBELoc, lngBELoc, alphaBELoc, 1);
			break;
		case 2:
			increment = coeffs[coeffInd + 1] * IG(x, y, xBELoc, yBELoc, lngBELoc, alphaBELoc, 2);
			break;
		case 3:
			increment = coeffs[coeffInd + 2] * IG(x, y, xBELoc, yBELoc, lngBELoc, alphaBELoc, 3);
			break;
	}

	atomicAdd(&nodes[nodeInd + 2], increment);
}

void GalerkinCuda::CalculatePotentialField()
{
	if (!initialisedData || !initialisedCoeffs) {
		printf("\nFalse while reading input data");
		return;
	}

	ResetData();

	cudaSetDevice(0);

	// try to parallel maximal effective
	// we have 3N [N - number of boundary elements] equations
	// each equation is a result of numeric integral and is a sum of p_j*k [j=1..N]
	// so we need to calculate each k

	//cudaDeviceProp prop;
	//cudaGetDeviceProperties(&prop, 0);
	//printf("Device is %s\nnumber of blocks %dx%dx%d (each %dx%dx%d) = number of threads %d\n", prop.name,
	//	prop.maxGridSize[0],
	//	prop.maxGridSize[1],
	//	prop.maxGridSize[2],
	//	prop.maxThreadsDim[0],
	//	prop.maxThreadsDim[1],
	//	prop.maxThreadsDim[2],
	//	prop.maxThreadsPerBlock);

	dim3 blockSize = dim3(beNum, 3, 1); // each cofficient is a summ of numIntDiscr terms
	dim3 gridSize = dim3(fdSizeX, fdSizeY, 1); // each boundary element have 3 equation which consist of (beNum * 3) coefficients

	// data pointers for a kernel
	float* cudaPotField;
	float* cudaBeInfo;
	float* cudaCoeffs;
	cudaMalloc((void**)&cudaPotField, potFieldSize * sizeof(float));
	cudaMalloc((void**)&cudaBeInfo, beInfoSize * sizeof(float));
	cudaMalloc((void**)&cudaCoeffs, coeffsSize * sizeof(float));
	cudaMemcpy(cudaPotField, potField, potFieldSize * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(cudaBeInfo, beInfo, beInfoSize * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(cudaCoeffs, coeffs, coeffsSize * sizeof(float), cudaMemcpyHostToDevice);

	NodePotential <<< gridSize, blockSize >>> (cudaPotField, cudaBeInfo, cudaCoeffs);


	cudaError_t cudaerr = cudaDeviceSynchronize();
	if (cudaerr != cudaSuccess)
		printf("kernel launch failed with error \"%s\".\n",
			cudaGetErrorString(cudaerr));

	cudaEvent_t syncEvent;
	cudaEventCreate(&syncEvent);
	cudaEventRecord(syncEvent, 0);
	cudaEventSynchronize(syncEvent);

	cudaMemcpy(potField, cudaPotField, potFieldSize * sizeof(float), cudaMemcpyDeviceToHost);
	cudaDeviceReset();
	printf("\nSolved success!");
}