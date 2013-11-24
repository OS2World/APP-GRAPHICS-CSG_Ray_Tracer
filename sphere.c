/*

sphere.c - Sphere logic
			  2      2      2  2
Defined solid for :- (x-a) +(y-b) +(z-c) -d <=0

*/

/*...sincludes:0:*/
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <malloc.h>
#include <memory.h>
#include <math.h>
#include "rt.h"
#include "vector.h"
#include "sil.h"
#define	_SPHERE_
#include "sphere.h"

/*...vrt\46\h:0:*/
/*...vvector\46\h:0:*/
/*...vsil\46\h:0:*/
/*...vsphere\46\h:0:*/
/*...e*/

/*...screate_sphere      \45\ create any general sphere:0:*/
SPHERE	*create_sphere(double r)
	{
	SPHERE	*sphere;

	if ( (sphere = malloc(sizeof(SPHERE))) == NULL )
		return NULL;

	sphere->a = 0.0;
	sphere->b = 0.0;
	sphere->c = 0.0;
	sphere->d = r;
	return sphere;
	}
/*...e*/
/*...scopy_sphere        \45\ make a copy of a sphere:0:*/
SPHERE	*copy_sphere(SPHERE *sphere)
	{
	SPHERE	*copy;

	if ( (copy = malloc(sizeof(SPHERE))) == NULL )
		return NULL;

	memcpy(copy, sphere, sizeof(SPHERE));
	return copy;
	}
/*...e*/
/*...sdestroy_sphere     \45\ destroy a sphere:0:*/
void	destroy_sphere(SPHERE *sphere)
	{
	free(sphere);
	}
/*...e*/

/*...strans_x_sphere     \45\ translate by amount in x direction:0:*/
void	trans_x_sphere(SPHERE *sphere, double t)
	{
	sphere->a += t;
	}
/*...e*/
/*...strans_y_sphere     \45\ translate by amount in y direction:0:*/
void	trans_y_sphere(SPHERE *sphere, double t)
	{
	sphere->b += t;
	}
/*...e*/
/*...strans_z_sphere     \45\ translate by amount in z direction:0:*/
void	trans_z_sphere(SPHERE *sphere, double t)
	{
	sphere->c += t;
	}
/*...e*/
/*...sscale_sphere       \45\ scale a sphere equally in all directions:0:*/
/*
The idea is to make the sphere 'factor' times bigger.
Just multiply d by the scale factor.
*/

void	scale_sphere(SPHERE *sphere, double factor)
	{
	sphere->d *= factor;
	}
/*...e*/
/*...srot_x_sphere       \45\ rotate about x axis by given angle:0:*/
void	rot_x_sphere(SPHERE *sphere, double angle)
	{
	double	b  = sphere->b;
	double	c  = sphere->c;
	double	ca = cos(angle);
	double	sa = sin(angle);

	sphere->b = b * ca - c * sa;
	sphere->c = b * sa + c * ca;
	}
/*...e*/
/*...srot_y_sphere       \45\ rotate about y axis by given angle:0:*/
void	rot_y_sphere(SPHERE *sphere, double angle)
	{
	double	c  = sphere->c;
	double	a  = sphere->a;
	double	ca = cos(angle);
	double	sa = sin(angle);

	sphere->c = c * ca - a * sa;
	sphere->a = c * sa + a * ca;
	}
/*...e*/
/*...srot_z_sphere       \45\ rotate about z axis by given angle:0:*/
void	rot_z_sphere(SPHERE *sphere, double angle)
	{
	double	a  = sphere->a;
	double	b  = sphere->b;
	double	ca = cos(angle);
	double	sa = sin(angle);

	sphere->a = a * ca - b * sa;
	sphere->b = a * sa + b * ca;
	}
/*...e*/

/*...sisects_reqd_sphere \45\ max number of isects we will generate:0:*/
int	isects_reqd_sphere(SPHERE *sphere)
	{
	sphere=sphere; /* Suppress 'unref arg' compiler warning */

	return 2;
	}
/*...e*/
/*...sintersect_sphere   \45\ determine intersection range of t:0:*/
/*

Any point along the line we are interested in is of the form p + tq.
							     ~    ~
	     2      2      2  2
Given:	(x-a) +(y-b) +(z-c) -d =0

Subs:	x = x +tx		y = y +ty		z = z +tz
	     p   q		     p   q		     p   q

	  2  2  2  2                      2  2  2  2
Gives:	(x +y +z )t - 2(ix +jy +kz )t + (i +j +k -d ) = 0
	  q  q  q         q   q   q

Where:	i = (a-x ), j = (b-y ), k = (c-z )
	        p           p           p

*/

/*...svalue_of_sphere:0:*/
static double value_of_sphere(void *shapeinfo, VECTOR v)
	{
	SPHERE	*sphere = (SPHERE *) shapeinfo;
	double	i = v.x - sphere->a;
	double	j = v.y - sphere->b;
	double	k = v.z - sphere->c;

	return i * i + j * j + k * k - sphere->d * sphere->d;
	}
/*...e*/

void	intersect_sphere(SPHERE *sphere, VECTOR p, VECTOR q, SIL *sil)
	{
	double	i  = sphere->a - p.x;
	double	j  = sphere->b - p.y;
	double	k  = sphere->c - p.z;
	double	qa = q.x * q.x + q.y * q.y + q.z * q.z;
	double	xx = -(i * q.x + j * q.y + k * q.z);
	double	qb = xx + xx;
	double	qc = i * i + j * j + k * k - sphere->d * sphere->d;

	intersect_quadratic_t(qa, qb, qc, p, q,
			      (void *) sphere, value_of_sphere, sil);
	}
/*...e*/
/*...snormal_to_sphere   \45\ find normal of surface at a given point:0:*/
/*

Use partial derivatives to find normal at a given point.
We will ignore the multiple of 2 on all components, producing a half length,
but correctly pointing vector, which is fair enough, as it will be normalised
later anyhow.

	     2      2      2  2
Given:	(x-a) +(y-b) +(z-c) +d =0

d/dx:	2(x-a)

d/dy:	2(y-b)

d/dz:	2(z-c)

*/

VECTOR	normal_to_sphere(SPHERE *sphere, VECTOR p)
	{
	VECTOR	normal;

	normal.x = p.x - sphere->a;
	normal.y = p.y - sphere->b;
	normal.z = p.z - sphere->c;

	return normal;
	}
/*...e*/
