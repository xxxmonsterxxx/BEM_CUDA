#include "GalerkinCudaSmooth.cuh"

#include "GalerkinData.h"

__global__ void GalerkinCudaSmooth::ElementInfluenceMatrix(float* infMatr, float* beinfo, float* fright)
{
	uint beNum = blockIdx.x; // each be generate equation
	uint beNumLocal = blockIdx.y; // be num which influent on fixed global be
	uint termNum = threadIdx.x; // numeric integral term (1..numIntDiscr)

	uint coeffNumGlobal = blockIdx.x * gridDim.y + blockIdx.y; // global number of coefficient in full influence matrix

	uint beNumInfoK1 = 19; // shift multiplier = size of beinfo struct odd
	uint beNumInfoK2 = 7;  // shift multiplier = size of beinfo struct even
	uint beNumId = (beNum % 2) ? ((beNum-1)/2) * (beNumInfoK1 + beNumInfoK2) + beNumInfoK1 : (beNum/2) * (beNumInfoK1 + beNumInfoK2);
	uint beNumLocId = (beNumLocal % 2) ? ((beNumLocal-1)/2) * (beNumInfoK1 + beNumInfoK2) + beNumInfoK1 : (beNumLocal/2) * (beNumInfoK1 + beNumInfoK2);

	//local boundary element info
	float xBELocL, yBELocL, lngBELocL, alphaBELocL; // local be info for LEFT side
	float xBELoc, yBELoc, lngBELoc, alphaBELoc;		// local be info for CENTER
	float xBELocR, yBELocR, lngBELocR, alphaBELocR; // local be info for RIGHT side
	uint BELocType = beinfo[beNumLocId + 2];

	if (BELocType == 1) {
		xBELoc = beinfo[beNumLocId + 3];
		yBELoc = beinfo[beNumLocId + 4];
		lngBELoc = beinfo[beNumLocId + 5];
		alphaBELoc = beinfo[beNumLocId + 6];
	}
	else if (BELocType == 2) {
		xBELocL = beinfo[beNumLocId + 5];
		yBELocL = beinfo[beNumLocId + 6];
		alphaBELocL = beinfo[beNumLocId + 9];
		lngBELocL = beinfo[beNumLocId + 10];

		xBELocR = beinfo[beNumLocId + 13];
		yBELocR = beinfo[beNumLocId + 14];
		alphaBELocR = beinfo[beNumLocId + 17];
		lngBELocR = beinfo[beNumLocId + 18];
	}

	// global boundary element info
	float xBE, yBE, lngBE, alphaBE;	// global be info
	uint  BEType = beinfo[beNumId + 2];

	// info for discret integral
	bool side = (termNum < (blockDim.x / 2));
	float discrStep;
	float xSub, ySub = 0;

	if (BEType == 1) {
		xBE = beinfo[beNumId + 0];
		yBE = beinfo[beNumId + 1];
		lngBE = beinfo[beNumId + 5];
		alphaBE = beinfo[beNumId + 6];

		discrStep = 2 * lngBE / blockDim.x;

		xSub = -lngBE + termNum * discrStep;
	}
	else if (BEType == 2) {
		if (side) { // left semilength
			xBE = beinfo[beNumId + 5];
			yBE = beinfo[beNumId + 6];
			alphaBE = beinfo[beNumId + 9];
			lngBE = beinfo[beNumId + 10];
			discrStep = lngBE / (blockDim.x / 2);

			xSub =  termNum * discrStep;
		}
		else { // right
			xBE = beinfo[beNumId + 13];
			yBE = beinfo[beNumId + 14];
			alphaBE = beinfo[beNumId + 17];
			lngBE = beinfo[beNumId + 18];
			discrStep = lngBE / (blockDim.x / 2);

			xSub = -lngBE + (termNum - blockDim.x/2) * discrStep;
		}
	}

	float xSubTransofrmed = 0, ySubTransformed = 0;

	Transform2D(xBE, yBE, alphaBE,
		0, 0, 0,
		xSub, ySub,
		xSubTransofrmed, ySubTransformed);

	float increment = 0;
	float frightIncrement = 0;
	float localInf;

	if (BELocType == 1)
		localInf = IG(xSubTransofrmed, ySubTransformed, xBELoc, yBELoc, lngBELoc, alphaBELoc, 2);
	else if (BELocType == 2)
		localInf = (IG(xSubTransofrmed, ySubTransformed, xBELocL, yBELocL, lngBELocL, alphaBELocL, 3) +
				    IG(xSubTransofrmed, ySubTransformed, xBELocR, yBELocR, lngBELocR, alphaBELocR, 1)   );

	if (BEType == 1) {
		increment = discrStep * f2(xSub, lngBE) * localInf;

		if (blockIdx.x == blockIdx.y)
			frightIncrement = discrStep * f2(xSub, lngBE) * Problem::InitCondition(xSubTransofrmed, ySubTransformed);
	}
	else if (BEType == 2) {
		if (side) {
			increment = discrStep * f3(xSub, lngBE) * localInf;
			if (blockIdx.x == blockIdx.y)
				frightIncrement = discrStep * f3(xSub, lngBE) * Problem::InitCondition(xSubTransofrmed, ySubTransformed);
		}
		else {
			increment = discrStep * f1(xSub, lngBE) * localInf;
			if (blockIdx.x == blockIdx.y)
				frightIncrement = discrStep * f1(xSub, lngBE) * Problem::InitCondition(xSubTransofrmed, ySubTransformed);
		}
	}

	atomicAdd(&infMatr[coeffNumGlobal], increment);
	atomicAdd(&fright[blockIdx.x], frightIncrement);
}

using namespace GalerkinMethod;

void GalerkinCudaSmooth::CalculateInfMatrix()
{
	if (!initialisedSmoothData) {
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
	dim3 gridSize = dim3(beNum, beNum, 1); // each boundary element generate a equation which consist of beNum coefficients

	// data pointers for a kernel
	float* cudaInfMatr;
	float* cudaBeInfo;
	float* cudaFright;
	cudaMalloc((void**)&cudaInfMatr, infMatrSize * sizeof(float));
	cudaMalloc((void**)&cudaBeInfo, beInfoSize * sizeof(float));
	cudaMalloc((void**)&cudaFright, beNum * sizeof(float));
	cudaMemcpy(cudaInfMatr, infMatr, infMatrSize * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(cudaBeInfo, beInfo, beInfoSize * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(cudaFright, fRight, beNum * sizeof(float), cudaMemcpyHostToDevice);

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
	cudaMemcpy(fRight, cudaFright, beNum * sizeof(float), cudaMemcpyDeviceToHost);
	cudaDeviceReset();
	printf("\nSolved success!");
}