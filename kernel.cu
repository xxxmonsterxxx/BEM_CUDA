#include "GalerkinData.h"
#include "GalerkinCuda.cuh"
#include "GalerkinSerial.h"
#include "GalerkinCudaSmooth.cuh"
#include "GalerkinSmoothSerial.h"
#include <ctime>

int main()
{
    uint start_time;
    uint end_time;
    uint spand_time;

    /*GalerkinMethod::InitInputData();
    for (int i = 0; i < 4; i++) {
        start_time = clock();
        GalerkinCuda::Solve();
        end_time = clock();
        spand_time = end_time - start_time;
        printf("\nCuda computation time = %d\n", spand_time);
    }

    start_time = clock();
    GalerkinSerial::Solve();
    end_time = clock();
    spand_time = end_time - start_time;
    printf("\nSeq. computation time = %d\n", spand_time);*/

    GalerkinMethod::InitInputSmoothData();
    for (int i = 0; i < 4; i++) {
        start_time = clock();
        GalerkinCudaSmooth::CalculateInfMatrix();
        end_time = clock();
        spand_time = end_time - start_time;
        printf("\nCuda computation time = %d\n", spand_time);
    }

    start_time = clock();
    GalerkinSmoothSerial::CalculateInfMatrix();
    end_time = clock();
    spand_time = end_time - start_time;
    printf("\nSeq. computation time = %d\n", spand_time);

    return 0;
}