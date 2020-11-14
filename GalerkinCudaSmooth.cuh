#pragma once

#include "Problem.h"
#include "Additional.h"

namespace GalerkinCudaSmooth {
	// kernel function to calculate one element of influence matrix
	__global__ void ElementInfluenceMatrix(float* infMatr, float* beinfo, float* fright);
	void CalculateInfMatrix();
}