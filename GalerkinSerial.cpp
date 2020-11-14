#include "GalerkinSerial.h"
#include "Problem.h"

using namespace GalerkinMethod;

void GalerkinSerial::CalculateInfMatrix()
{
	if (!initialisedData) {
		printf("\nFalse while reading input data");
		return;
	}

	ResetData();

	// equation system consist of (beNum *3) equations --> each BE generate 3 equation
	for (uint be = 0; be < beNum; be++) {

		// global boundary element info
		uint beNumInfoK = 8; // shift multiplier = size of beinfo struct
		float xBE = beInfo[be * beNumInfoK + 2];
		float yBE = beInfo[be * beNumInfoK + 3];
		float alphaBE = beInfo[be * beNumInfoK + 6];
		float lngBE = beInfo[be * beNumInfoK + 7];
		float discrStep = 2 * lngBE / numIntDiscr;

		// each equation of one BE correspond to one form function and consist of (beNum * 3) terms
		for (uint beFunc = 1; beFunc <= 3; beFunc++) {

			// each term is a sum of corresponds elements of each local BE and correspond local BE function
			for (uint beLoc = 0; beLoc < beNum; beLoc++) {

				//local boundary element info
				float xBELoc = beInfo[beLoc * beNumInfoK + 2];
				float yBELoc = beInfo[beLoc * beNumInfoK + 3];
				float alphaBELoc = beInfo[beLoc * beNumInfoK + 6];
				float lngBELoc = beInfo[beLoc * beNumInfoK + 7];

				for (uint beLocFunc = 1; beLocFunc <= 3; beLocFunc++) {

					// each this term is a numeric integral which is a sum of numIntDiscr terms
					for (uint i = 0; i < numIntDiscr; i++) {
						// info for discret integral
						float xSub = -lngBE + discrStep * numIntDiscr;
						float ySub = 0;

						float xSubTransofrmed = 0, ySubTransformed = 0;

						Transform2D(xBE, yBE, alphaBE,
							0, 0, 0,
							xSub, ySub,
							xSubTransofrmed, ySubTransformed);

						float increment = 0;
						float frightIncrement = 0;

						switch (beFunc) {
						case 1:
							increment = discrStep * f1(xSub, lngBE) * IG(xSubTransofrmed, ySubTransformed, xBELoc, yBELoc, lngBELoc, alphaBELoc, beLocFunc);
							if (be == beLoc)
								frightIncrement = discrStep * f1(xSub, lngBE) * Problem::InitCondition(xSubTransofrmed, ySubTransformed);
							break;
						case 2:
							increment = discrStep * f2(xSub, lngBE) * IG(xSubTransofrmed, ySubTransformed, xBELoc, yBELoc, lngBELoc, alphaBELoc, beLocFunc);
							if (be == beLoc)
								frightIncrement = discrStep * f2(xSub, lngBE) * Problem::InitCondition(xSubTransofrmed, ySubTransformed);
							break;
						case 3:
							increment = discrStep * f3(xSub, lngBE) * IG(xSubTransofrmed, ySubTransformed, xBELoc, yBELoc, lngBELoc, alphaBELoc, beLocFunc);
							if (be == beLoc)
								frightIncrement = discrStep * f3(xSub, lngBE) * Problem::InitCondition(xSubTransofrmed, ySubTransformed);
							break;
						default:
							break;
						}

						//uint index = be*beNum + (beFunc-1)*3 + beLoc*beNum + (beLocFunc-1);
						uint index = (beLocFunc-1) + beLoc*3 + (beFunc-1)*beNum*3 + be*3*beNum*3;
						uint fRightIndex = (beFunc-1) + be*3;
						infMatr[index] += increment;
						fRight[fRightIndex] += frightIncrement;
					}
				}
			}
		}
	}
}