#pragma once
/*
* Galerkins boundary elements method to solve potential rectangle 2D problem
*/
#include "Additional.h"
#include "Problem.h"

namespace GalerkinCuda {

	// kernel function to calculate one element of influence matrix
	__global__ void ElementInfluenceMatrix(float* infMatr, float* beinfo, float* fright);
	void CalculateInfMatrix();
}