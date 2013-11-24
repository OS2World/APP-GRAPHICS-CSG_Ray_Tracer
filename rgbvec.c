/*

rgbvec.c - RGB Colour vector datatype

*/

/*...sincludes:0:*/
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <malloc.h>
#include <memory.h>
#include <math.h>
#define	_RGBVEC_
#include "rgbvec.h"

/*...vrgbvec\46\h:0:*/
/*...e*/

/*...sscale_rgbvec:0:*/
RGBVEC scale_rgbvec(RGBVEC rgbvec, double scalar)
	{
	rgbvec.r *= scalar;
	rgbvec.g *= scalar;
	rgbvec.b *= scalar;

	return rgbvec;
	}
/*...e*/

/*...srgb_interp_1d:0:*/
RGBVEC	rgbvec_interp_1d(RGBVEC rgbvec_1d[2], double w1)
	{
	double	w0 = 1.0 - w1;
	RGBVEC	rgbvec;

	rgbvec.r = w0 * rgbvec_1d[0].r + w1 * rgbvec_1d[1].r;
	rgbvec.g = w0 * rgbvec_1d[0].g + w1 * rgbvec_1d[1].g;
	rgbvec.b = w0 * rgbvec_1d[0].b + w1 * rgbvec_1d[1].b;

	return rgbvec;
	}
/*...e*/
/*...srgb_interp_2d:0:*/
RGBVEC	rgbvec_interp_2d(
	RGBVEC rgbvec_2d[2][2],
	double w1major,
	double w1minor
	)
	{
	RGBVEC	rgbvec_1d[2];

	rgbvec_1d[0] = rgbvec_interp_1d(rgbvec_2d[0], w1minor);
	rgbvec_1d[1] = rgbvec_interp_1d(rgbvec_2d[1], w1minor);

	return rgbvec_interp_1d(rgbvec_1d, w1major);
	}
/*...e*/
