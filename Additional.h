#pragma once

#define uint unsigned int

#define _USE_MATH_DEFINES

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <math.h>
#include <string>
#include <fstream>

#define EPSILON 1e-6

/*
Transform2DClassic - transform coordinate (x,y)
from original system (0, 0, 0 rads)
to destination system (x_dst, y_dst, alpha_dst)
result stored to (x_,y_)
*/
static __device__ __host__ void Transform2DClassic(float x_dst, float y_dst, float alpha_dst,
												   float x, float y, float& x_, float& y_)
{
	x_ = (x - x_dst) * cosf(alpha_dst) + (y - y_dst) * sinf(alpha_dst);
	y_ = (x - x_dst) * -sinf(alpha_dst) + (y - y_dst) * cosf(alpha_dst);
}

/*
Transform2DClassic - transform coordinate (x,y)
from system (x_src, y_src, alpha_src)
to original system (0, 0, 0 rads)
result stored to (x_,y_)
*/
static __device__ __host__ void Transform2DClassicRevers(float x_src, float y_src, float alpha_src,
														   float x, float y, float& x_, float& y_)
{
	float x_origin_original, y_origin_original; // coords of original system in source system
	Transform2DClassic(x_src, y_src, alpha_src, 0, 0, x_origin_original, y_origin_original);
	Transform2DClassic(x_origin_original, y_origin_original, -alpha_src, x, y, x_, y_);
}

/*
Transform2D - transform coordinate (x,y) 
from source system (x_src, y_src, alpha_src)
to destination system (x_dst, y_dst, alpha_dst)
result stored to (x_,y_)
*/
static __device__ __host__ void Transform2D(float x_src, float y_src, float alpha_src,
							float x_dst, float y_dst, float alpha_dst,
							float x, float y, float &x_, float &y_)
{
	float x_original, y_original; // find searched coords in original system
	Transform2DClassicRevers(x_src, y_src, alpha_src, x, y, x_original, y_original);
	// now transform this coord in destination system
	Transform2DClassic(x_dst, y_dst, alpha_dst, x_original, y_original, x_, y_);
}

static __device__ __host__ bool IsEqualf(float f1, float f2)
{
	return fabsf(f1 - f2) < EPSILON;
}

static __device__ __host__ float f1(float x, float lng)
{
	return -x / (2 * lng) + powf(x, 2) / (2 * powf(lng, 2));
}

static __device__ __host__ float f2(float x, float lng)
{
	return 1 - powf(x, 2) / powf(lng, 2);
}

static __device__ __host__ float f3(float x, float lng)
{
	return x / (2 * lng) + powf(x, 2) / (2 * powf(lng, 2));
}

/*
* Integral of canonic solve with first form function
*/
static __device__ __host__ float G1(float x, float y, float lng, float ksi)
{
	return -((-ksi + x) * (9 * lng * (ksi + 3 * x) - 2 * (2 * powf(ksi, 2) + 5 * ksi * x + 11 * powf(x, 2) - 6 * powf(y, 2))) -
		12 * y * (3 * lng * x - 3 * powf(x, 2) + powf(y, 2)) * atanf((-ksi + x) / y) +
		(9 * lng * (powf(ksi, 2) - powf(x, 2) + powf(y, 2)) + 6 * (-powf(ksi, 3) + powf(x, 3) - 3 * x * powf(y, 2))) * logf(powf(-ksi + x, 2) + powf(y, 2))) /
		(144. * powf(lng, 2) * M_PI);
}

/*
* Integral of canonic solve with second form function
*/
static __device__ __host__ float G2(float x, float y, float lng, float ksi)
{
	return ((ksi - x) * (2 * powf(ksi, 2) - 18 * powf(lng, 2) + 5 * ksi * x + 11 * powf(x, 2) - 6 * powf(y, 2)) -
		6 * y * (3 * powf(lng, 2) - 3 * powf(x, 2) + powf(y, 2)) * atanf((-ksi + x) / y) -
		3 * (powf(ksi, 3) - powf(x, 3) + 3 * powf(lng, 2) * (-ksi + x) + 3 * x * powf(y, 2)) * logf(powf(-ksi + x, 2) + powf(y, 2))) / (36. * powf(lng, 2) * M_PI);
}

/*
* * Integral of canonic solve with third form function
*/
static __device__ __host__ float G3(float x, float y, float lng, float ksi)
{
	return ((-ksi + x) * (9 * lng * (ksi + 3 * x) + 2 * (2 * powf(ksi, 2) + 5 * ksi * x + 11 * powf(x, 2) - 6 * powf(y, 2))) +
		12 * y * (-3 * lng * x - 3 * powf(x, 2) + powf(y, 2)) * atanf((-ksi + x) / y) +
		(9 * lng * (powf(ksi, 2) - powf(x, 2) + powf(y, 2)) + 6 * (powf(ksi, 3) - powf(x, 3) + 3 * x * powf(y, 2))) * logf(powf(-ksi + x, 2) + powf(y, 2))) /
		(144. * powf(lng, 2) * M_PI);
}

/*
* Function to avoid undefined value by x==ksi
*/
static __device__ __host__ float Ggood(float x, float y, float lng, float ksi, uint funcNumber)
{
	if (IsEqualf(x, ksi)) {
		switch (funcNumber)
		{
			case 1:
				return G1(x - EPSILON, y, lng, ksi);
			case 2:
				return G2(x - EPSILON, y, lng, ksi);
			case 3:
				return G3(x - EPSILON, y, lng, ksi);
			default:
				return 0;
		}
	}
	else {
		switch (funcNumber)
		{
		case 1:
			return G1(x, y, lng, ksi);
		case 2:
			return G2(x, y, lng, ksi);
		case 3:
			return G3(x, y, lng, ksi);
		default:
			return 0;
		}
	}
}

/*
*/
static __device__ __host__ float IG(float x, float y, float xC, float yC, float lng, float alpha, uint funcNumber)
{
	float xTransformed, yTransformed;
	Transform2D(0,0,0,
				xC,yC,alpha,
				x,y,
				xTransformed, yTransformed);

	return Ggood(xTransformed, yTransformed, lng, lng, funcNumber) - Ggood(xTransformed, yTransformed, lng, -lng, funcNumber);
}
