#include "GalerkinData.h"
#include "GalerkinCuda.cuh"
#include "GalerkinSerial.h"
#include "GalerkinCudaSmooth.cuh"
#include "GalerkinSmoothSerial.h"
#include <ctime>

using namespace GalerkinMethod;

int main()
{
    //uint start_time;
    //uint end_time;
    //uint spand_time;

    ////----------------------------------------------------
    ////----------------------------------------------------
    //GalerkinMethod::InitInputData();
    ////----------------------------------------------------
    //start_time = clock();
    //GalerkinCuda::CalculateInfMatrix();
    //end_time = clock();
    //spand_time = end_time - start_time;
    //printf("\nCuda inf matr computation time = %d\n", spand_time);

    //start_time = clock();
    //GalerkinSerial::CalculateInfMatrix();
    //end_time = clock();
    //spand_time = end_time - start_time;
    //printf("\nSeq. inf matr computation time = %d\n", spand_time);
    ////----------------------------------------------------
    //GalerkinMethod::InitInputSmoothData();
    ////----------------------------------------------------
    //start_time = clock();
    //GalerkinCudaSmooth::CalculateInfMatrix();
    //end_time = clock();
    //spand_time = end_time - start_time;
    //printf("\nCuda smooth inf matr computation time = %d\n", spand_time);

    //start_time = clock();
    //GalerkinSmoothSerial::CalculateInfMatrix();
    //end_time = clock();
    //spand_time = end_time - start_time;
    //printf("\nSeq. smooth inf matr computation time = %d\n", spand_time);
    ////----------------------------------------------------
    ////----------------------------------------------------
    //start_time = clock();
    //GalerkinCuda::CalculatePotentialField();
    //end_time = clock();
    //spand_time = end_time - start_time;
    //printf("\nCuda calc potential field computation time = %d\n", spand_time);

    //start_time = clock();
    //GalerkinSerial::CalculatePotentialField();
    //end_time = clock();
    //spand_time = end_time - start_time;
    //printf("\Seq. calc potential field computation time = %d\n", spand_time);
    //----------------------------------------------------
    //start_time = clock();
    GalerkinMethod::InitInputSmoothData();
    GalerkinSmoothSerial::CalculatePotentialField();
    for (int i = 0; i < potFieldSize; i++) {
        printf("\n%f", potField[i]);
    }

    //end_time = clock();
    //spand_time = end_time - start_time;
    //printf("\nCuda smooth calc potential field computation time = %d\n", spand_time);

    /*start_time = clock();
    GalerkinSmoothSerial::CalculatePotentialField();
    end_time = clock();
    spand_time = end_time - start_time;
    printf("\Seq. smooth calc potential field computation time = %d\n", spand_time);*/
    //----------------------------------------------------
    

    return 0;
}