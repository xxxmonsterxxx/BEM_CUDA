#include "GalerkinData.h"

bool GalerkinMethod::initialisedData = false;
uint GalerkinMethod::numIntDiscr = 100;
uint GalerkinMethod::beNum = 0;
uint GalerkinMethod::beInfoSize = 0;
float* GalerkinMethod::beInfo = nullptr;
uint GalerkinMethod::infMatrSize = 0;
float* GalerkinMethod::infMatr = nullptr;
uint GalerkinMethod::fdSizeX = 0;
uint GalerkinMethod::fdSizeY = 0;
uint GalerkinMethod::potFieldSize = GalerkinMethod::fdSizeX * GalerkinMethod::fdSizeY * 3;
float* GalerkinMethod::potField = nullptr;
float* GalerkinMethod::fRight = nullptr;
uint GalerkinMethod::fRightSize = 0;
bool GalerkinMethod::initialisedCoeffs = false;
float* GalerkinMethod::coeffs = nullptr;
uint GalerkinMethod::coeffsSize = 0;

bool GalerkinMethod::initialisedSmoothData = false;

void GalerkinMethod::ResetData()
{
	if (!initialisedData && !initialisedSmoothData)
		return;

	for (uint i = 0; i < infMatrSize; i++)
		infMatr[i] = 0;

	for (uint i = 0; i < fRightSize; i++)
		fRight[i] = 0;

	if (!initialisedCoeffs)
		return;

	for (uint i = 2; i < potFieldSize; i+=3)
		potField[i] = 0;
}

void GalerkinMethod::InitInputData()
{
	std::string input_file_name = "D:/repos/BEM_CUDA/be_info.txt";

	std::ifstream input_stream(input_file_name);
	if (!input_stream.is_open()) {
		printf("\nFile not found");
		return;
	}

	input_stream >> beNum;

	if (beNum == 0) {
		printf("\nBoundary elements information not found.");
		input_stream.close();
		return;
	}

	// each beinfo element consist of 8 numbers {{x-,y-,x0,y0,x+,y+,a,lng}}
	beInfoSize = beNum * 8;
	beInfo = new float[beInfoSize];
	for (uint i = 0; i < beInfoSize; i++) {
		input_stream >> beInfo[i];
	}

	// infMatr is a matrix (3*beNum,3*beNum) cuz each be has 3 coefficient
	infMatrSize = 3 * beNum * 3 * beNum;
	infMatr = new float[infMatrSize];
	for (uint i = 0; i < infMatrSize; i++)
		infMatr[i] = 0;

	fRightSize = 3 * beNum;
	fRight = new float[fRightSize];
	for (uint i = 0; i < fRightSize; i++)
		fRight[i] = 0;

	initialisedData = true;
	printf("\nInput OK! BE num = %d Coeff num = %d", beNum, beNum * 3);
	input_stream.close();

	input_file_name = "D:/repos/BEM_CUDA/bem_coeffs.txt";
	input_stream.open(input_file_name);
	input_stream >> coeffsSize;
	if (coeffsSize == 0 || coeffsSize != beNum * 3) {
		printf("\nERROR! Incorrect coefficients!");
		input_stream.close();
	} else {
		coeffs = new float[coeffsSize];

		for (int i = 0; i < coeffsSize; i++)
			input_stream >> coeffs[i];

		input_stream.close();

		initialisedCoeffs = true;
		printf("\nInput coeffs OK!");
	}

	fdSizeX = 10 + 1;
	fdSizeY = 10 + 1;
	potFieldSize = 3 * fdSizeX * fdSizeY;
	float stepX = 10. / (fdSizeX-1);
	float stepY =  8. / (fdSizeY-1);
	potField = new float[potFieldSize];
	for (uint i = 0; i < fdSizeX; i++) {
		for (uint j = 0; j < fdSizeY; j++) {
			uint index = (i * fdSizeY + j) * 3;
			potField[index + 0] = -5 + i * stepX;
			potField[index + 1] = -8 + j * stepY;
			potField[index + 2] = 0;
		}
	}
}

void GalerkinMethod::Export()
{
	std::string output_file_name = "D:/repos/BEM_CUDA/outresult.txt";
	std::ofstream out;
	out.open(output_file_name);

	for (int i = 0; i < infMatrSize; i++)
		out << infMatr << '\n';
	out.close();
}

void GalerkinMethod::InitInputSmoothData()
{
	std::string input_file_name = "D:/repos/BEM_CUDA/be_info_smooth.txt";

	std::ifstream input_stream(input_file_name);
	if (!input_stream.is_open()) {
		printf("\nFile not found");
		return;
	}

	input_stream >> beNum;

	if (beNum == 0) {
		printf("\nBoundary elements information not found.");
		return;
	}

	// each odd beinfo element consist of 19 numbers, each even 7
	beInfoSize = beNum / 2 * (19+7);
	beInfo = new float[beInfoSize];
	for (uint i = 0; i < beInfoSize; i++) {
		input_stream >> beInfo[i];
	}

	infMatrSize = beNum * beNum;
	infMatr = new float[infMatrSize];
	for (uint i = 0; i < infMatrSize; i++)
		infMatr[i] = 0;

	fRightSize = beNum;
	fRight = new float[fRightSize];
	for (int i = 0; i < beNum; i++)
		fRight[i] = 0;

	initialisedSmoothData = true;
	printf("\nSmooth data input OK! BE num = %d Coeff num = %d", beNum, beNum);
}