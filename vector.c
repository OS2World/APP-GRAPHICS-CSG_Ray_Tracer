/*

vector.c - Vector datatype

*/

/*...sincludes:0:*/
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#define	_VECTOR_
#include "vector.h"
/*...e*/

/*...smagnitude:0:*/
double	magnitude(VECTOR p)
	{
	return sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
	}
/*...e*/
/*...sunit_vector:0:*/
VECTOR	unit_vector(VECTOR p)
	{
	double	mag = magnitude(p);

	if ( mag == 0.0 )
		return p;	/* Zero vector cannot be unit-ised */

	p.x /= mag;
	p.y /= mag;
	p.z /= mag;

	return p;
	}
/*...e*/
/*...svector_product:0:*/
VECTOR	vector_product(VECTOR p, VECTOR q)
	{
	VECTOR	product;

	product.x = p.y * q.z - p.z * q.y;
	product.y = p.z * q.x - p.x * q.z;
	product.z = p.x * q.y - p.y * q.x;

	return product;
	}
/*...e*/
/*...sscalar_product:0:*/
double	scalar_product(VECTOR p, VECTOR q)
	{
	return p.x * q.x + p.y * q.y + p.z * q.z;
	}
/*...e*/
/*...svector_sum:0:*/
VECTOR	vector_sum(VECTOR p, VECTOR q)
	{
	p.x += q.x;
	p.y += q.y;
	p.z += q.z;

	return p;
	}
/*...e*/
/*...snegate_vector:0:*/
VECTOR	negate_vector(VECTOR p)
	{
	p.x = -p.x;
	p.y = -p.y;
	p.z = -p.z;

	return p;
	}
/*...e*/
/*...svector_difference:0:*/
VECTOR	vector_difference(VECTOR p, VECTOR q)
	{
	p.x -= q.x;
	p.y -= q.y;
	p.z -= q.z;

	return p;
	}
/*...e*/
/*...sscale_vector:0:*/
VECTOR	scale_vector(VECTOR p, double scalar)
	{
	p.x *= scalar;
	p.y *= scalar;
	p.z *= scalar;

	return p;
	}
/*...e*/
/*...sinv_scale_vector:0:*/
VECTOR	inv_scale_vector(VECTOR p, double inv_scalar)
	{
	p.x /= inv_scalar;
	p.y /= inv_scalar;
	p.z /= inv_scalar;

	return p;
	}
/*...e*/
/*...st_along_pq:0:*/
VECTOR	t_along_pq(VECTOR p, VECTOR q, double t)
	{
	p.x += q.x * t;
	p.y += q.y * t;
	p.z += q.z * t;

	return p;
	}
/*...e*/
/*...srot_x_vector:0:*/
VECTOR	rot_x_vector(VECTOR v, double angle)
	{
	double	ca = cos(angle);
	double	sa = sin(angle);
	VECTOR	rotated;

	rotated.x = v.x;
	rotated.y = ca * v.y - sa * v.z;
	rotated.z = sa * v.y + ca * v.z;
	return rotated;
	}
/*...e*/
/*...srot_y_vector:0:*/
VECTOR	rot_y_vector(VECTOR v, double angle)
	{
	double	ca = cos(angle);
	double	sa = sin(angle);
	VECTOR	rotated;

	rotated.x = sa * v.z + ca * v.x;
	rotated.y = v.y;
	rotated.z = ca * v.z - sa * v.x;
	return rotated;
	}
/*...e*/
/*...srot_z_vector:0:*/
VECTOR	rot_z_vector(VECTOR v, double angle)
	{
	double	ca = cos(angle);
	double	sa = sin(angle);
	VECTOR	rotated;

	rotated.x = ca * v.x - sa * v.y;
	rotated.y = sa * v.x + ca * v.y;
	rotated.z = v.z;
	return rotated;
	}
/*...e*/
