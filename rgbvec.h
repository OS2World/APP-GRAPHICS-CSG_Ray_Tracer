/*

rgbvec.h - RGB Colour vector datatype and operation functions

*/

typedef struct { double	r, g, b; } RGBVEC;

#ifndef _RGBVEC_

extern RGBVEC scale_rgbvec(RGBVEC rgbvec, double scalar);

extern RGBVEC rgbvec_interp_1d(RGBVEC rgbvec_1d[2], double w1);

extern RGBVEC rgbvec_interp_2d(
	RGBVEC rgbvec_2d[2][2],
	double w1major,
	double w1minor
	);

#endif
