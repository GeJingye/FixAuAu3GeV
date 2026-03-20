#ifndef StAnaCuts_H
#define StAnaCuts_H

/* **************************************************
 *
 *  Authors: Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************
 */
/* ******************************************************************************************
 * read PicoDst document about FixedTarget AuAu 3.85GeV for produciton within TOF acceptance*
 * ******************************************************************************************
 */
#include "Rtypes.h"
#include <string>
#include <array>

namespace anaCuts
{
	const std::array<UInt_t, 4> trigNumber = {
		820017
	};// MinimumBias with etof 3.85GeVAuAu

	// event cuts 
	Float_t const Vz_up = 201;// < cm.
	Float_t const Vz_low = 199.5;// < cm.
	Float_t const Vr = 2; // cm
	// tracks cuts
	Int_t   const NHitsFit = 20;//12
	Int_t   const NHitsDedx = 14;//6
	Float_t const NHitsFitRatio = 0.52;
	Float_t const Dca = 1.0;//3.0
	Float_t const GPt = 0.0;
	Float_t const Eta = -3.0;
    Float_t const PhiVCutMRange = 0.2;
	Int_t const nCenBins = 9;
}
#endif
