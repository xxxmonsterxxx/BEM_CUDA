//#define _USE_MATH_DEFINES
//
//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"
//
//#include <math.h>
//#include <stdio.h>
//#include <conio.h>
//#include <iostream>
//#include <fstream>
//#include <ctime>
//#include <string>
//
//void GetInputAndCalcInfluence();
//void GetInputAndCalcDistr();
//void CreateMatrix(double* Matrixij, double* BECoords, double* BE, int BeNumber, int BEInfoSize);
//void CreateNodes(double* Matrixij, int size, double* Coeffs, int BeNumber, double* BE, int BEInfoSize);
//
//// cuda function create influence matrix with pointer 'Matrixij' for boundary condition in points with coords 'BECoords',
//// but boundary elements are in coords 'BE'
//__global__ void MatrixCreation(double* Matrixij, double* BECoords, double* BE)
//{
//	int i = blockIdx.x;
//	int j0 = threadIdx.x;
//	int index = i * blockDim.x + j0;
//
//	// coords of calculation point which is be coords with shift to avoid undeterminated state
//	double x = BECoords[i];
//	double y = BECoords[i + 1];
//
//	if (j0 % 2)
//	{
//		int j = (j0 - 1) / 2 * (19 + 7) + 19;
//		Matrixij[index] = (-6 * atanf((-BE[j + 5] + (-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6])) / ((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6]))) * ((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6])) *
//			(3 * powf(BE[j + 5], 2) + powf((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6]), 2) - 3 * powf((-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6]), 2)) +
//			(BE[j + 5] - (-BE[j + 3] + x) * cosf(BE[j + 6]) - (-BE[j + 4] + y) * sinf(BE[j + 6])) * (-16 * powf(BE[j + 5], 2) - 6 * powf((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6]), 2) +
//				5 * BE[j + 5] * ((-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6])) + 11 * powf((-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6]), 2)) -
//			3 * logf(powf((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6]), 2) + powf(-BE[j + 5] + (-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6]), 2)) *
//			(powf(BE[j + 5], 3) + 3 * powf((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6]), 2) * ((-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6])) -
//				powf((-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6]), 3) + 3 * powf(BE[j + 5], 2) * (-BE[j + 5] + (-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6])))) / (36. * powf(BE[j + 5], 2) * M_PI) -
//			(-6 * atanf((BE[j + 5] + (-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6])) / ((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6]))) * ((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6])) *
//				(3 * powf(BE[j + 5], 2) + powf((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6]), 2) - 3 * powf((-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6]), 2)) +
//				(-BE[j + 5] - (-BE[j + 3] + x) * cosf(BE[j + 6]) - (-BE[j + 4] + y) * sinf(BE[j + 6])) * (-16 * powf(BE[j + 5], 2) - 6 * powf((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6]), 2) -
//					5 * BE[j + 5] * ((-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6])) + 11 * powf((-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6]), 2)) -
//				3 * logf(powf((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6]), 2) + powf(BE[j + 5] + (-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6]), 2)) *
//				(-powf(BE[j + 5], 3) + 3 * powf((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6]), 2) * ((-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6])) -
//					powf((-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6]), 3) + 3 * powf(BE[j + 5], 2) * (BE[j + 5] + (-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6])))) / (36. * powf(BE[j + 5], 2) * M_PI);
//	}
//	else
//	{
//		int j = j0 / 2 * (19 + 7);
//		Matrixij[index] = -(12 * atanf((BE[j + 10] + (-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9])) / ((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]))) * ((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9])) *
//			(powf((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]), 2) - 3 * BE[j + 10] * ((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9])) -
//				3 * powf((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]), 2)) +
//			(BE[j + 10] + (-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9])) * (9 * BE[j + 10] * (-BE[j + 10] + 3 * ((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]))) +
//				2 * (2 * powf(BE[j + 10], 2) - 6 * powf((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]), 2) - 5 * BE[j + 10] * ((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9])) +
//					11 * powf((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]), 2))) +
//			logf(powf((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]), 2) + powf(BE[j + 10] + (-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]), 2)) *
//			(9 * BE[j + 10] * (powf(BE[j + 10], 2) + powf((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]), 2) - powf((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]), 2)) +
//				6 * (-powf(BE[j + 10], 3) + 3 * powf((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]), 2) * ((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9])) -
//					powf((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]), 3)))) / (144. * powf(BE[j + 10], 2) * M_PI) +
//			(12 * atanf((-BE[j + 10] + (-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9])) / ((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]))) * ((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9])) *
//				(powf((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]), 2) - 3 * BE[j + 10] * ((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9])) -
//					3 * powf((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]), 2)) +
//				(-BE[j + 10] + (-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9])) * (9 * BE[j + 10] * (BE[j + 10] + 3 * ((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]))) +
//					2 * (2 * powf(BE[j + 10], 2) - 6 * powf((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]), 2) + 5 * BE[j + 10] * ((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9])) +
//						11 * powf((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]), 2))) +
//				logf(powf((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]), 2) + powf(-BE[j + 10] + (-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]), 2)) *
//				(9 * BE[j + 10] * (powf(BE[j + 10], 2) + powf((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]), 2) - powf((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]), 2)) +
//					6 * (powf(BE[j + 10], 3) + 3 * powf((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]), 2) * ((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9])) -
//						powf((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]), 3)))) / (144. * powf(BE[j + 10], 2) * M_PI) -
//			(-12 * atanf((-BE[j + 18] + (-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17])) / ((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]))) *
//				((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17])) * (powf((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]), 2) +
//					3 * BE[j + 18] * ((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17])) - 3 * powf((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]), 2)) +
//				(-BE[j + 18] + (-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17])) * (9 * BE[j + 18] * (BE[j + 18] + 3 * ((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]))) -
//					2 * (2 * powf(BE[j + 18], 2) - 6 * powf((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]), 2) + 5 * BE[j + 18] * ((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17])) +
//						11 * powf((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]), 2))) +
//				logf(powf((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]), 2) + powf(-BE[j + 18] + (-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]), 2)) *
//				(9 * BE[j + 18] * (powf(BE[j + 18], 2) + powf((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]), 2) - powf((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]), 2)) +
//					6 * (-powf(BE[j + 18], 3) - 3 * powf((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]), 2) * ((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17])) +
//						powf((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]), 3)))) / (144. * powf(BE[j + 18], 2) * M_PI) +
//			(-12 * atanf((BE[j + 18] + (-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17])) / ((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]))) *
//				((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17])) * (powf((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]), 2) +
//					3 * BE[j + 18] * ((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17])) - 3 * powf((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]), 2)) +
//				(BE[j + 18] + (-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17])) * (9 * BE[j + 18] * (-BE[j + 18] + 3 * ((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]))) -
//					2 * (2 * powf(BE[j + 18], 2) - 6 * powf((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]), 2) - 5 * BE[j + 18] * ((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17])) +
//						11 * powf((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]), 2))) +
//				logf(powf((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]), 2) + powf(BE[j + 18] + (-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]), 2)) *
//				(9 * BE[j + 18] * (powf(BE[j + 18], 2) + powf((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]), 2) - powf((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]), 2)) +
//					6 * (powf(BE[j + 18], 3) - 3 * powf((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]), 2) * ((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17])) +
//						powf((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]), 3)))) / (144. * powf(BE[j + 18], 2) * M_PI);
//	}
//}
//
//// cuda function calculate one term for one be coord and add it to value belong to be coord
//__global__ void CalculateNodes(double* Matrixij, double* Coeff, double* BE)
//{
//	int i = blockIdx.x; // index by x axis
//	int j0 = blockIdx.y; // index by be element
//	int k = threadIdx.x; // index by y axis
//
//	int index = (i * blockDim.x + k) * 3;
//
//	double x = Matrixij[index];
//	double y = Matrixij[index + 1];
//
//	if (j0 % 2) {
//		int j = (j0 - 1) / 2 * (19 + 7) + 19;
//		double increment = Coeff[j0] * ((-6 * atanf((-BE[j + 5] + (-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6])) / ((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6]))) * ((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6])) *
//			(3 * powf(BE[j + 5], 2) + powf((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6]), 2) - 3 * powf((-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6]), 2)) +
//			(BE[j + 5] - (-BE[j + 3] + x) * cosf(BE[j + 6]) - (-BE[j + 4] + y) * sinf(BE[j + 6])) * (-16 * powf(BE[j + 5], 2) - 6 * powf((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6]), 2) +
//				5 * BE[j + 5] * ((-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6])) + 11 * powf((-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6]), 2)) -
//			3 * logf(powf((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6]), 2) + powf(-BE[j + 5] + (-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6]), 2)) *
//			(powf(BE[j + 5], 3) + 3 * powf((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6]), 2) * ((-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6])) -
//				powf((-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6]), 3) + 3 * powf(BE[j + 5], 2) * (-BE[j + 5] + (-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6])))) / (36. * powf(BE[j + 5], 2) * M_PI) -
//			(-6 * atanf((BE[j + 5] + (-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6])) / ((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6]))) * ((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6])) *
//				(3 * powf(BE[j + 5], 2) + powf((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6]), 2) - 3 * powf((-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6]), 2)) +
//				(-BE[j + 5] - (-BE[j + 3] + x) * cosf(BE[j + 6]) - (-BE[j + 4] + y) * sinf(BE[j + 6])) * (-16 * powf(BE[j + 5], 2) - 6 * powf((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6]), 2) -
//					5 * BE[j + 5] * ((-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6])) + 11 * powf((-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6]), 2)) -
//				3 * logf(powf((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6]), 2) + powf(BE[j + 5] + (-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6]), 2)) *
//				(-powf(BE[j + 5], 3) + 3 * powf((-BE[j + 4] + y) * cosf(BE[j + 6]) + (BE[j + 3] - x) * sinf(BE[j + 6]), 2) * ((-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6])) -
//					powf((-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6]), 3) + 3 * powf(BE[j + 5], 2) * (BE[j + 5] + (-BE[j + 3] + x) * cosf(BE[j + 6]) + (-BE[j + 4] + y) * sinf(BE[j + 6])))) / (36. * powf(BE[j + 5], 2) * M_PI));
//		atomicAdd(&(Matrixij[index + 2]), increment);
//	}
//	else {
//		int j = j0 / 2 * (19 + 7);
//		double increment = Coeff[j0] * (-(12 * atanf((BE[j + 10] + (-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9])) / ((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]))) * ((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9])) *
//			(powf((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]), 2) - 3 * BE[j + 10] * ((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9])) -
//				3 * powf((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]), 2)) +
//			(BE[j + 10] + (-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9])) * (9 * BE[j + 10] * (-BE[j + 10] + 3 * ((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]))) +
//				2 * (2 * powf(BE[j + 10], 2) - 6 * powf((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]), 2) - 5 * BE[j + 10] * ((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9])) +
//					11 * powf((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]), 2))) +
//			logf(powf((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]), 2) + powf(BE[j + 10] + (-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]), 2)) *
//			(9 * BE[j + 10] * (powf(BE[j + 10], 2) + powf((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]), 2) - powf((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]), 2)) +
//				6 * (-powf(BE[j + 10], 3) + 3 * powf((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]), 2) * ((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9])) -
//					powf((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]), 3)))) / (144. * powf(BE[j + 10], 2) * M_PI) +
//			(12 * atanf((-BE[j + 10] + (-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9])) / ((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]))) * ((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9])) *
//				(powf((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]), 2) - 3 * BE[j + 10] * ((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9])) -
//					3 * powf((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]), 2)) +
//				(-BE[j + 10] + (-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9])) * (9 * BE[j + 10] * (BE[j + 10] + 3 * ((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]))) +
//					2 * (2 * powf(BE[j + 10], 2) - 6 * powf((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]), 2) + 5 * BE[j + 10] * ((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9])) +
//						11 * powf((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]), 2))) +
//				logf(powf((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]), 2) + powf(-BE[j + 10] + (-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]), 2)) *
//				(9 * BE[j + 10] * (powf(BE[j + 10], 2) + powf((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]), 2) - powf((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]), 2)) +
//					6 * (powf(BE[j + 10], 3) + 3 * powf((-BE[j + 6] + y) * cosf(BE[j + 9]) + (BE[j + 5] - x) * sinf(BE[j + 9]), 2) * ((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9])) -
//						powf((-BE[j + 5] + x) * cosf(BE[j + 9]) + (-BE[j + 6] + y) * sinf(BE[j + 9]), 3)))) / (144. * powf(BE[j + 10], 2) * M_PI) -
//			(-12 * atanf((-BE[j + 18] + (-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17])) / ((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]))) *
//				((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17])) * (powf((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]), 2) +
//					3 * BE[j + 18] * ((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17])) - 3 * powf((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]), 2)) +
//				(-BE[j + 18] + (-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17])) * (9 * BE[j + 18] * (BE[j + 18] + 3 * ((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]))) -
//					2 * (2 * powf(BE[j + 18], 2) - 6 * powf((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]), 2) + 5 * BE[j + 18] * ((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17])) +
//						11 * powf((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]), 2))) +
//				logf(powf((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]), 2) + powf(-BE[j + 18] + (-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]), 2)) *
//				(9 * BE[j + 18] * (powf(BE[j + 18], 2) + powf((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]), 2) - powf((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]), 2)) +
//					6 * (-powf(BE[j + 18], 3) - 3 * powf((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]), 2) * ((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17])) +
//						powf((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]), 3)))) / (144. * powf(BE[j + 18], 2) * M_PI) +
//			(-12 * atanf((BE[j + 18] + (-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17])) / ((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]))) *
//				((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17])) * (powf((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]), 2) +
//					3 * BE[j + 18] * ((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17])) - 3 * powf((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]), 2)) +
//				(BE[j + 18] + (-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17])) * (9 * BE[j + 18] * (-BE[j + 18] + 3 * ((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]))) -
//					2 * (2 * powf(BE[j + 18], 2) - 6 * powf((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]), 2) - 5 * BE[j + 18] * ((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17])) +
//						11 * powf((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]), 2))) +
//				logf(powf((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]), 2) + powf(BE[j + 18] + (-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]), 2)) *
//				(9 * BE[j + 18] * (powf(BE[j + 18], 2) + powf((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]), 2) - powf((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]), 2)) +
//					6 * (powf(BE[j + 18], 3) - 3 * powf((-BE[j + 14] + y) * cosf(BE[j + 17]) + (BE[j + 13] - x) * sinf(BE[j + 17]), 2) * ((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17])) +
//						powf((-BE[j + 13] + x) * cosf(BE[j + 17]) + (-BE[j + 14] + y) * sinf(BE[j + 17]), 3)))) / (144. * powf(BE[j + 18], 2) * M_PI));
//		atomicAdd(&(Matrixij[index + 2]), increment);
//	}
//}
//
//// function using CUDA function "MatrixCreation" create influence matrix by pointer Matrixij[size*size]
//void CreateMatrix(double* Matrixij, double* BECoords, double* BE, int BeNumber, int BEInfoSize)
//{
//	double* dev_a, * dev_b, * dev_c;
//
//	cudaSetDevice(0);
//
//	cudaMalloc((void**)&dev_a, BeNumber * BeNumber * sizeof(double));
//	cudaMalloc((void**)&dev_b, 2 * BeNumber * sizeof(double));
//	cudaMalloc((void**)&dev_c, BEInfoSize * sizeof(double));
//	cudaMemcpy(dev_a, Matrixij, BeNumber * BeNumber * sizeof(double), cudaMemcpyHostToDevice);
//	cudaMemcpy(dev_b, BECoords, BeNumber * 2 * sizeof(double), cudaMemcpyHostToDevice);
//	cudaMemcpy(dev_c, BE, BEInfoSize * sizeof(double), cudaMemcpyHostToDevice);
//
//	dim3 blockSize = dim3(BeNumber, 1, 1);
//	dim3 gridSize = dim3(BeNumber, 1, 1);
//
//	MatrixCreation << <gridSize, blockSize >> > (dev_a, dev_b, dev_c);
//
//	cudaEvent_t syncEvent;
//	cudaEventCreate(&syncEvent);
//	cudaEventRecord(syncEvent, 0);
//	cudaEventSynchronize(syncEvent);
//
//	cudaMemcpy(Matrixij, dev_a, BeNumber * BeNumber * sizeof(double), cudaMemcpyDeviceToHost);
//}
//
//void CreateNodes(double* Matrixij, int size, double* Coeffs, int BeNumber, double* BE, int BEInfoSize)
//{
//	double* dev_a, * dev_b, * dev_c;
//
//	cudaSetDevice(0);
//
//	cudaMalloc((void**)&dev_a, 3 * size * size * sizeof(double));
//	cudaMalloc((void**)&dev_b, BeNumber * sizeof(double));
//	cudaMalloc((void**)&dev_c, BEInfoSize * sizeof(double));
//	cudaMemcpy(dev_a, Matrixij, 3 * size * size * sizeof(double), cudaMemcpyHostToDevice);
//	cudaMemcpy(dev_b, Coeffs, BeNumber * sizeof(double), cudaMemcpyHostToDevice);
//	cudaMemcpy(dev_c, BE, BEInfoSize * sizeof(double), cudaMemcpyHostToDevice);
//
//	dim3 blockSize = dim3(size, 1, 1);
//	dim3 gridSize = dim3(size, BeNumber, 1);
//
//	CalculateNodes << < gridSize, blockSize >> > (dev_a, dev_b, dev_c);
//
//	cudaEvent_t syncEvent;
//	cudaEventCreate(&syncEvent);
//	cudaEventRecord(syncEvent, 0);
//	cudaEventSynchronize(syncEvent);
//
//	cudaMemcpy(Matrixij, dev_a, 3 * size * size * sizeof(double), cudaMemcpyDeviceToHost);
//}
//
//void GetInputAndCalcInfluence()
//{
//	double* Coords, * BE;
//	double* Matrixij;
//
//	int bediscr_array[] = { 5,10,15,20,25,30,35 };
//	string input_file_name;
//	string output_file_name;
//
//	for (int i = 0; i < 7; i++)
//	{
//		int bediscr = bediscr_array[i];
//		int CoordsNumber, BEInfoSize;
//
//		input_file_name = "D:/Docs/article_03_20/data/shifted_coords" + to_string(bediscr) + ".txt";
//		ifstream in;
//		in.open(input_file_name);
//
//		in >> CoordsNumber;
//		Coords = new double[CoordsNumber];
//
//		for (int i = 0; i < CoordsNumber; i++)
//			in >> Coords[i];
//
//		in.close();
//
//		input_file_name = "D:/Docs/article_03_20/data/beinfo" + to_string(bediscr) + ".txt";
//		in.open(input_file_name);
//
//		in >> BEInfoSize;
//		BE = new double[BEInfoSize];
//
//		for (int i = 0; i < BEInfoSize; i++)
//			in >> BE[i];
//
//		in.close();
//
//		int benumb = CoordsNumber / 2;
//
//		Matrixij = new double[benumb * benumb];
//
//		unsigned int start_time;
//		unsigned int end_time;
//		unsigned int search_time = 0;
//
//		/////////////// CUDA method //////////////////////////
//		CreateMatrix(Matrixij, Coords, BE, benumb, BEInfoSize);
//		start_time = clock();
//		CreateMatrix(Matrixij, Coords, BE, benumb, BEInfoSize);
//		end_time = clock();
//		search_time = end_time - start_time;
//
//		output_file_name = "D:/Docs/article_03_20/CUDA_output/matr" + to_string(bediscr) + ".txt";
//		ofstream out;
//		out.open(output_file_name);
//
//		for (int i = 0; i < benumb; i++)
//			for (int j = 0; j < benumb; j++)
//				out << Matrixij[i * benumb + j] << '\n';
//		out << search_time;
//		out.close();
//
//		/////////////// sequential method ////////////////////
//		start_time = clock();
//		MatrixCreationSeq(Matrixij, Coords, BE, benumb);
//		end_time = clock();
//		search_time = end_time - start_time;
//
//		output_file_name = "D:/Docs/article_03_20/SEQ_output/matr" + to_string(bediscr) + ".txt";
//		out.open(output_file_name);
//
//		for (int i = 0; i < benumb; i++)
//			for (int j = 0; j < benumb; j++)
//				out << Matrixij[i * benumb + j] << '\n';
//		out << search_time;
//		out.close();
//
//		///////////////////////////////////////////////////////
//
//		delete(Matrixij);
//		delete(BE);
//		delete(Coords);
//	}
//}
//
//void GetInputAndCalcDistr()
//{
//	int bediscr_array[] = { 5,10,15,20,25,30,35 };
//	int areadiscr_array[] = { 10, 20, 30, 40, 50 };
//	double* BE, * Coefs;
//	int CoefsNumb, BEInfoSize;
//	string input_file_name;
//	string output_file_name;
//	int bediscr, areadiscr;
//
//	for (int k = 0; k < 7; k++)
//	{
//		bediscr = bediscr_array[k];
//		input_file_name = "D:/Docs/article_03_20/data/coefs" + to_string(bediscr) + ".txt";
//		ifstream in;
//		in.open(input_file_name);
//
//		in >> CoefsNumb;
//		Coefs = new double[CoefsNumb];
//
//		for (int i = 0; i < CoefsNumb; i++)
//			in >> Coefs[i];
//
//		in.close();
//
//		input_file_name = "D:/Docs/article_03_20/data/beinfo" + to_string(bediscr) + ".txt";
//		in.open(input_file_name);
//
//		in >> BEInfoSize;
//		BE = new double[BEInfoSize];
//
//		for (int i = 0; i < BEInfoSize; i++)
//			in >> BE[i];
//
//		in.close();
//
//		for (int s = 0; s < 5; s++)
//		{
//			areadiscr = areadiscr_array[s];
//			float discrx = (WIDTH - 2 * EPS) / (areadiscr - 1);
//			float discry = (HEIGHT - 2 * EPS) / (areadiscr - 1);
//
//			int matrixsize = areadiscr * areadiscr * 3;
//			double* Matrixij = new double[matrixsize];
//
//			for (int i = 0; i < areadiscr; i++)
//				for (int j = 0; j < areadiscr; j++)
//				{
//					int idx = (i * areadiscr + j) * 3;
//					Matrixij[idx] = -WIDTH / 2 + EPS + discrx * i;
//					Matrixij[idx + 1] = -EPS - discry * j;
//					Matrixij[idx + 2] = 0;
//				}
//
//			unsigned int start_time;
//			unsigned int end_time;
//			unsigned int search_time = 0;
//
//			CreateNodes(Matrixij, areadiscr, Coefs, CoefsNumb, BE, BEInfoSize); // startup calculation to activate cuda memory
//			start_time = clock();
//			CreateNodes(Matrixij, areadiscr, Coefs, CoefsNumb, BE, BEInfoSize);
//			end_time = clock();
//			search_time = end_time - start_time;
//
//			output_file_name = "D:/Docs/article_03_20/CUDA_Output/node" + to_string(bediscr) + "-" + to_string(areadiscr) + ".txt";
//			ofstream out;
//			out.open(output_file_name);
//
//			for (int i = 0; i < matrixsize; i = i + 3)
//				out << Matrixij[i] << " " << Matrixij[i + 1] << " " << Matrixij[i + 2] << '\n';
//			out << search_time;
//
//			out.close();
//
//			/////////////// sequential method ////////////////////
//			start_time = clock();
//			CalculateNodesSeq(Matrixij, areadiscr, Coefs, CoefsNumb, BE);
//			end_time = clock();
//			search_time = end_time - start_time;
//
//			output_file_name = "D:/Docs/article_03_20/SEQ_output/node" + to_string(bediscr) + "-" + to_string(areadiscr) + ".txt";
//			out.open(output_file_name);
//
//			for (int i = 0; i < matrixsize; i = i + 3)
//				out << Matrixij[i] << " " << Matrixij[i + 1] << " " << Matrixij[i + 2] << '\n';
//			out << search_time;
//			out.close();
//
//			///////////////////////////////////////////////////////
//
//			delete(Matrixij);
//		}
//
//		delete(Coefs);
//		delete(BE);
//	}
//}