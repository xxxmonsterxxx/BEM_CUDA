#include "GalerkinSmoothSerial.h"
#include "Problem.h"

using namespace GalerkinMethod;

void GalerkinSmoothSerial::CalculateInfMatrix()
{
	if (!initialisedSmoothData) {
		printf("\nFalse while reading input data");
		return;
	}

	ResetData();

	uint beNumInfoK1 = 19; // shift multiplier = size of beinfo struct odd
	uint beNumInfoK2 = 7;  // shift multiplier = size of beinfo struct even
	// equation system consist of (beNum *3) equations --> each BE generate 3 equation
	for (uint be = 0; be < beNum; be++) {
		uint beNumId = (be % 2) ? ((be - 1) / 2) * (beNumInfoK1 + beNumInfoK2) + beNumInfoK1 : (be / 2) * (beNumInfoK1 + beNumInfoK2);

		// global boundary element info
		float xBE, yBE, lngBE, alphaBE;	// global be info
		uint  BEType = beInfo[beNumId + 2];

		// each term is a sum of corresponds elements of each local BE and correspond local BE function
		for (uint beLoc = 0; beLoc < beNum; beLoc++) {
			uint beNumLocId = (beLoc % 2) ? ((beLoc - 1) / 2) * (beNumInfoK1 + beNumInfoK2) + beNumInfoK1 : (beLoc / 2) * (beNumInfoK1 + beNumInfoK2);

			//local boundary element info
			float xBELocL, yBELocL, lngBELocL, alphaBELocL; // local be info for LEFT side
			float xBELoc, yBELoc, lngBELoc, alphaBELoc;		// local be info for CENTER
			float xBELocR, yBELocR, lngBELocR, alphaBELocR; // local be info for RIGHT side
			uint BELocType = beInfo[beNumLocId + 2];

			if (BELocType == 1) {
				xBELoc = beInfo[beNumLocId + 3];
				yBELoc = beInfo[beNumLocId + 4];
				lngBELoc = beInfo[beNumLocId + 5];
				alphaBELoc = beInfo[beNumLocId + 6];
			}
			else {
				xBELocL = beInfo[beNumLocId + 5];
				yBELocL = beInfo[beNumLocId + 6];
				alphaBELocL = beInfo[beNumLocId + 9];
				lngBELocL = beInfo[beNumLocId + 10];

				xBELocR = beInfo[beNumLocId + 13];
				yBELocR = beInfo[beNumLocId + 14];
				alphaBELocR = beInfo[beNumLocId + 17];
				lngBELocR = beInfo[beNumLocId + 18];
			}
			
			// each this term is a numeric integral which is a sum of numIntDiscr terms
			for (uint i = 0; i < numIntDiscr; i++) {
				// info for discret integral
				bool side = (i < (numIntDiscr / 2));
				float discrStep;
				float xSub, ySub = 0;

				if (BEType == 1) {
					xBE = beInfo[beNumId + 0];
					yBE = beInfo[beNumId + 1];
					lngBE = beInfo[beNumId + 5];
					alphaBE = beInfo[beNumId + 6];
					discrStep = 2 * lngBE / numIntDiscr;
					xSub = -lngBE + i * discrStep;
				}
				else {
					if (side) { // left semilength
						xBE = beInfo[beNumId + 5];
						yBE = beInfo[beNumId + 6];
						alphaBE = beInfo[beNumId + 9];
						lngBE = beInfo[beNumId + 10];
						discrStep = lngBE / (numIntDiscr / 2);

						xSub = i * discrStep;
					}
					else { // right
						xBE = beInfo[beNumId + 13];
						yBE = beInfo[beNumId + 14];
						alphaBE = beInfo[beNumId + 17];
						lngBE = beInfo[beNumId + 18];
						discrStep = lngBE / (numIntDiscr / 2);

						xSub = -lngBE + (i - numIntDiscr / 2) * discrStep;
					}
				}


				float xSubTransofrmed = 0, ySubTransformed = 0;

				Transform2D(xBE, yBE, alphaBE,
					0, 0, 0,
					xSub, ySub,
					xSubTransofrmed, ySubTransformed);

				float increment = 0;
				float frightIncrement = 0;
				float localInf = 0;

				if (BELocType == 1)
					localInf = IG(xSubTransofrmed, ySubTransformed, xBELoc, yBELoc, lngBELoc, alphaBELoc, 2);
				else if (BELocType == 2)
					localInf = (IG(xSubTransofrmed, ySubTransformed, xBELocL, yBELocL, lngBELocL, alphaBELocL, 3) +
						IG(xSubTransofrmed, ySubTransformed, xBELocR, yBELocR, lngBELocR, alphaBELocR, 1));

				if (BEType == 1) {
					increment = discrStep * f2(xSub, lngBE) * localInf;

					if (be == beLoc)
						frightIncrement = discrStep * f2(xSub, lngBE) * Problem::InitCondition(xSubTransofrmed, ySubTransformed);
				}
				else if (BEType == 2) {
					if (side) {
						increment = discrStep * f3(xSub, lngBE) * localInf;
						if (be == beLoc)
							frightIncrement = discrStep * f3(xSub, lngBE) * Problem::InitCondition(xSubTransofrmed, ySubTransformed);
					}
					else {
						increment = discrStep * f1(xSub, lngBE) * localInf;
						if (be == beLoc)
							frightIncrement = discrStep * f1(xSub, lngBE) * Problem::InitCondition(xSubTransofrmed, ySubTransformed);
					}
				}

				//uint index = be*beNum + (beFunc-1)*3 + beLoc*beNum + (beLocFunc-1);
				uint index = be * beNum + beLoc; // global number of coefficient in full influence matrix
				infMatr[index] += increment;
				fRight[be] += frightIncrement;
			}
		}
	}
}