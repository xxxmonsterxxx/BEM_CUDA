#pragma once

#include "Additional.h"

namespace Problem {

	static __device__ float InitCondition(float x, float y) {
		if (IsEqualf(y, 0) &&
			(IsEqualf(x, -2.5) || x > -2.5) &&
			(IsEqualf(x, 2.5) || x < 2.5)) {
			return 5 * sqrtf(1 - x * x / (float)(2.5 * 2.5));
		}
		else {
			return 0;
		}
	}
}