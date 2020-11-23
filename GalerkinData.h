#pragma once

#include "Additional.h"

namespace GalerkinMethod {
	extern bool initialisedCoeffs;
	extern uint coeffsSize;
	extern float* coeffs;

	extern bool initialisedData;
	extern uint numIntDiscr; // discretization of be for numeric integrate

	extern uint beNum; // number of boundary elements
	extern uint beInfoSize;
	extern float* beInfo; // boundary elements info data

	extern float* infMatr;
	extern uint infMatrSize;

	extern float* fRight; // right part of matrix equation
	extern uint fRightSize;

	extern uint fdSizeX;
	extern uint fdSizeY;
	extern uint potFieldSize; // cuz each point should to have 3 number
	extern float* potField; // potential field data === (x,y,P) P - potential

	void ResetData();

// ---------------------------------Galerkin classic algorithm-------------------------------------

	void InitInputData(); // read files with initial data
	void Export();

// ---------------------------------Galerkin smooth algorithm--------------------------------------

	extern bool initialisedSmoothData;

	void InitInputSmoothData(); // read files with initial data for using smooth Galerkin algorithm
}