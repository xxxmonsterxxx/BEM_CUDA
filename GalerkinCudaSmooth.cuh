#pragma once

#include "Problem.h"
#include "Additional.h"

namespace GalerkinCudaSmooth {
	// kernel function to calculate one element of influence matrix
	__global__ void ElementInfluenceMatrix(float* infMatr, float* beinfo, float* fright);
	void CalculateInfMatrix();

	// kernel function to calculate node potential
	__global__ void NodePotential(float* nodes, float* beinfo, float* coeffs);
	void CalculatePotentialField();
}