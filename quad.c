/*

quad.c - Arbitrary quadratic logic

		       2   2   2
Defined solid for :- ax +by +cz +dxy+eyz+fzx+gx+hy+iz+j<=0

Ellisoids, spheres, infinitely long cylinders, double cones and planes are
all special case of this equation.

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
#define	_QUAD_
#include "quad.h"

/*...vrt\46\h:0:*/
/*...vvector\46\h:0:*/
/*...vsil\46\h:0:*/
/*...vquad\46\h:0:*/
/*...e*/

/*...screate_quad        \45\ create any general quadratic:0:*/
QUAD	*create_quad(
	double a, double b, double c,
	double d, double e, double f,
	double g, double h, double i,
	double j
	)
	{
	QUAD	*quad;

	if ( (quad = malloc(sizeof(QUAD))) == NULL )
		return NULL;

	quad->a = a;
	quad->b = b;
	quad->c = c;
	quad->d = d;
	quad->e = e;
	quad->f = f;
	quad->g = g;
	quad->h = h;
	quad->i = i;
	quad->j = j;

	return quad;
	}
/*...e*/
/*...screate_ellipsoid   \45\ create an ellipsoid of a given size at the origin:0:*/
QUAD	*create_ellipsoid(double rx, double ry, double rz)
	{
	return create_quad(1.0 / (rx * rx), 1.0 / (ry * ry), 1.0 / (rz * rz),
			   0.0, 0.0, 0.0,
			   0.0, 0.0, 0.0,
			   -1.0);
	}
/*...e*/
/*...screate_x_ell_cyl   \45\ create elliptical cylinder along x axis of given radii:0:*/
QUAD	*create_x_ell_cyl(double ry, double rz)
	{
	return create_quad(0.0, 1.0 / (ry * ry), 1.0 / (rz * rz),
			   0.0, 0.0, 0.0,
			   0.0, 0.0, 0.0,
			   -1.0);
	}
/*...e*/
/*...screate_y_ell_cyl   \45\ create elliptical cylinder along y axis of given radii:0:*/
QUAD	*create_y_ell_cyl(double rx, double rz)
	{
	return create_quad(1.0 / (rx * rx), 0.0, 1.0 / (rz * rz),
			   0.0, 0.0, 0.0,
			   0.0, 0.0, 0.0,
			   -1.0);
	}
/*...e*/
/*...screate_z_ell_cyl   \45\ create elliptical cylinder along z axis of given radii:0:*/
QUAD	*create_z_ell_cyl(double rx, double ry)
	{
	return create_quad(1.0 / (rx * rx), 1.0 / (ry * ry), 0.0,
			   0.0, 0.0, 0.0,
			   0.0, 0.0, 0.0,
			   -1.0);
	}
/*...e*/
/*...screate_x_cyl       \45\ create cylinder along x axis of given radii:0:*/
QUAD	*create_x_cyl(double r)
	{
	return create_x_ell_cyl(r, r);
	}
/*...e*/
/*...screate_y_cyl       \45\ create cylinder along y axis of given radii:0:*/
QUAD	*create_y_cyl(double r)
	{
	return create_y_ell_cyl(r, r);
	}
/*...e*/
/*...screate_z_cyl       \45\ create cylinder along z axis of given radii:0:*/
QUAD	*create_z_cyl(double r)
	{
	return create_z_ell_cyl(r, r);
	}
/*...e*/
/*...screate_x_ell_cone  \45\ create elliptical cone along x axis of given radii:0:*/
QUAD	*create_x_ell_cone(double ky, double kz)
	{
	return create_quad(-1.0, ky * ky, kz * kz,
			   0.0, 0.0, 0.0,
			   0.0, 0.0, 0.0,
			   0.0);
	}
/*...e*/
/*...screate_y_ell_cone  \45\ create elliptical cone along y axis of given radii:0:*/
QUAD	*create_y_ell_cone(double kx, double kz)
	{
	return create_quad(kx * kx, -1.0, kz * kz,
			   0.0, 0.0, 0.0,
			   0.0, 0.0, 0.0,
			   0.0);
	}
/*...e*/
/*...screate_z_ell_cone  \45\ create elliptical cone along z axis of given radii:0:*/
QUAD	*create_z_ell_cone(double kx, double ky)
	{
	return create_quad(kx * kx, ky * ky, -1.0,
			   0.0, 0.0, 0.0,
			   0.0, 0.0, 0.0,
			   0.0);
	}
/*...e*/
/*...screate_x_cone      \45\ create cone along x axis of given radii:0:*/
QUAD	*create_x_cone(double k)
	{
	return create_x_ell_cone(k, k);
	}
/*...e*/
/*...screate_y_cone      \45\ create cone along y axis of given radii:0:*/
QUAD	*create_y_cone(double k)
	{
	return create_y_ell_cone(k, k);
	}
/*...e*/
/*...screate_z_cone      \45\ create cone along z axis of given radii:0:*/
QUAD	*create_z_cone(double k)
	{
	return create_z_ell_cone(k, k);
	}
/*...e*/
/*...scopy_quad          \45\ make a copy of a quad:0:*/
QUAD	*copy_quad(QUAD *quad)
	{
	QUAD	*copy;

	if ( (copy = malloc(sizeof(QUAD))) == NULL )
		return NULL;

	memcpy(copy, quad, sizeof(QUAD));
	return copy;
	}
/*...e*/
/*...sdestroy_quad       \45\ destroy a quad:0:*/
void	destroy_quad(QUAD *quad)
	{
	free(quad);
	}
/*...e*/

/*...strans_x_quad       \45\ translate by amount in x direction:0:*/
void	trans_x_quad(QUAD *quad, double t)
	{
	quad->j += (quad->a * t * t - quad->g * t);
	quad->i -= quad->f * t;
	quad->h -= quad->d * t;
	quad->g -= 2.0 * quad->a * t;
	}
/*...e*/
/*...strans_y_quad       \45\ translate by amount in y direction:0:*/
void	trans_y_quad(QUAD *quad, double t)
	{
	quad->j += (quad->b * t * t - quad->h * t);
	quad->i -= quad->e * t;
	quad->h -= 2.0 * quad->b * t;
	quad->g -= quad->d * t;
	}
/*...e*/
/*...strans_z_quad       \45\ translate by amount in z direction:0:*/
void	trans_z_quad(QUAD *quad, double t)
	{
	quad->j += (quad->c * t * t - quad->i * t);
	quad->i -= 2.0 * quad->c * t;
	quad->h -= quad->e * t;
	quad->g -= quad->f * t;
	}
/*...e*/
/*...sscale_x_quad       \45\ scale by factor in x direction:0:*/
void	scale_x_quad(QUAD *quad, double factor)
	{
	quad->a /= (factor * factor);
	quad->d /= factor;
	quad->f /= factor;
	quad->g /= factor;
	}
/*...e*/
/*...sscale_y_quad       \45\ scale by factor in y direction:0:*/
void	scale_y_quad(QUAD *quad, double factor)
	{
	quad->b /= (factor * factor);
	quad->d /= factor;
	quad->e /= factor;
	quad->h /= factor;
	}
/*...e*/
/*...sscale_z_quad       \45\ scale by factor in z direction:0:*/
void	scale_z_quad(QUAD *quad, double factor)
	{
	quad->c /= (factor * factor);
	quad->e /= factor;
	quad->f /= factor;
	quad->i /= factor;
	}
/*...e*/
/*...srot_x_quad         \45\ rotate about x axis by given angle:0:*/
void	rot_x_quad(QUAD *quad, double angle)
	{
	double	b    = quad->b;
	double	c    = quad->c;
	double	d    = quad->d;
	double	e    = quad->e;
	double	f    = quad->f;
	double	h    = quad->h;
	double	i    = quad->i;
	double	ca   = cos(angle);
	double	sa   = sin(angle);
	double	caca = ca * ca;
	double	sasa = sa * sa;
	double	saca = sa * ca;

	quad->b = b * caca + c * sasa - e * saca;
	quad->c = b * sasa + c * caca + e * saca;
	quad->d = d * ca - f * sa;
	quad->e = 2.0 * (b - c) * saca + e * (caca - sasa);
	quad->f = d * sa + f * ca;
	quad->h = h * ca - i * sa;
	quad->i = h * sa + i * ca;
	}
/*...e*/
/*...srot_y_quad         \45\ rotate about y axis by given angle:0:*/
void	rot_y_quad(QUAD *quad, double angle)
	{
	double	a    = quad->a;
	double	c    = quad->c;
	double	d    = quad->d;
	double	e    = quad->e;
	double	f    = quad->f;
	double	g    = quad->g;
	double	i    = quad->i;
	double	ca   = cos(angle);
	double	sa   = sin(angle);
	double	caca = ca * ca;
	double	sasa = sa * sa;
	double	saca = sa * ca;

	quad->a = a * caca + c * sasa + f * saca;
	quad->c = a * sasa + c * caca - f * saca;
	quad->d = e * sa + d * ca;
	quad->e = e * ca - d * sa;
	quad->f = 2.0 * (c - a) * saca + f * (caca - sasa);
	quad->g = i * sa + g * ca;
	quad->i = i * ca - g * sa;
	}
/*...e*/
/*...srot_z_quad         \45\ rotate about z axis by given angle:0:*/
void	rot_z_quad(QUAD *quad, double angle)
	{
	double	a    = quad->a;
	double	b    = quad->b;
	double	d    = quad->d;
	double	e    = quad->e;
	double	f    = quad->f;
	double	g    = quad->g;
	double	h    = quad->h;
	double	ca   = cos(angle);
	double	sa   = sin(angle);
	double	caca = ca * ca;
	double	sasa = sa * sa;
	double	saca = sa * ca;

	quad->a = a * caca + b * sasa - d * saca;
	quad->b = a * sasa + b * caca + d * saca;
	quad->d = 2.0 * (a - b) * saca + d * (caca - sasa);
	quad->e = f * sa + e * ca;
	quad->f = f * ca - e * sa;
	quad->g = g * ca - h * sa;
	quad->h = g * sa + h * ca;
	}
/*...e*/

/*...sisects_reqd_quad   \45\ max number of isects we will generate:0:*/
int	isects_reqd_quad(QUAD *quad)
	{
	quad=quad; /* Suppress 'unref arg' compiler warning */

	return 4;
	}
/*...e*/
/*...sintersect_quad     \45\ intersect quadratic:0:*/
/*

Those who don't like maths should page down a few times very rapidly!

Any point along the line we are interested in is of the form p + tq.
							     ~    ~
	  2   2   2
Given:	ax +by +cz +dxy+eyz+fzx+gx+hy+iz+j=0

Subs:	x = x +tx		y = y +ty		z = z +tz
	     p   q		     p   q		     p   q

	   2   2   2                    2
Gives:	(ax +by +cz +dx y +ey z +fz x )t  +
	   q   q   q   q q   q q   q q

	(2[ax x +by y +cz z ]+d[x y +x y ]+e[y z +y z ]+f[z x +z x ]+gx +hy +iz )t +
	     p q   p q   p q     p q  q p     p q  q p     p q  q p    q   q   q

	   2   2   2
	(ax +by +cz +dx y +ey z +fz x +gx +hy +iz +j) = 0
	   p   p   p   p p   p p   p p   p   p   p


	  2   2   2                   
Opt:	ax +by +cz +dx y +ey z +fz x			costs 12x's and 5+'s
	  q   q   q   q q   q q   q q

	(ax +dy )x +(by +ez )y +(cz +fx )z		costs 9x's and 5+'s
	   q   q  q    q   q  q    q   q  q


	  2   2   2
Opt:	ax +by +cz +dx y +ey z +fz x +gx +hy +iz +j	costs 15x's and 9+'s
	  p   p   p   p p   p p   p p   p   p   p

	(ax +dy +g)x +(by +ez +h)y +(cz +fx +i)z +j	costs 9x's and 9+'s
	   p   p    p    p   p    p    p   p    p
*/

/*...svalue_of_quad:0:*/
/*
Use this to determine value of quadratic at given point in space.

		  2   2   2
Opt:	ax +by +cz +dxy+eyz+fzx+gx+hy+iz+j		costs 15x's and 9+'s

	(ax+dy+g)x+(by+ez+h)y+(cz+fx+i)z+j		costs 9x's and 9+'s

*/

static double value_of_quad(void *shapeinfo, VECTOR v)
	{
	QUAD	*quad = (QUAD *) shapeinfo;
	double	x = v.x, y = v.y, z = v.z;

	return (quad->a * x + quad->d * y + quad->g) * x +
	       (quad->b * y + quad->e * z + quad->h) * y +
	       (quad->c * z + quad->f * x + quad->i) * z +
	        quad->j;
	}
/*...e*/

void	intersect_quad(QUAD *quad, VECTOR p, VECTOR q, SIL *sil)
	{
	double	a  = quad->a;
	double	b  = quad->b;
	double	c  = quad->c;
	double	d  = quad->d;
	double	e  = quad->e;
	double	f  = quad->f;
	double	g  = quad->g;
	double	h  = quad->h;
	double	i  = quad->i;
	double	j  = quad->j;
	double	qa = (a * q.x + d * q.y) * q.x +
		     (b * q.y + e * q.z) * q.y +
		     (c * q.z + f * q.x) * q.z;
	double	xx = (a * p.x * q.x + b * p.y * q.y + c * p.z * q.z);
	double	qb = xx + xx +
		     d * (p.x * q.y + q.x * p.y) +
		     e * (p.y * q.z + q.y * p.z) +
		     f * (p.z * q.x + q.z * p.x) +
		     g * q.x + h * q.y + i * q.z;
	double	qc = (a * p.x + d * p.y + g) * p.x +
		     (b * p.y + e * p.z + h) * p.y +
		     (c * p.z + f * p.x + i) * p.z + j;

	intersect_quadratic_t(qa, qb, qc, p, q,
			      (void *) quad, value_of_quad, sil);
	}
/*...e*/
/*...snormal_to_quad     \45\ find normal of surface at a given point:0:*/
/*

Use partial derivatives to find normal at a given point.

Given:	  2   2   2
	ax +by +cz +dxy+eyz+fzx+gx+hy+iz+j=0

d/dx:	2ax+dy+fz+g

d/dy:	2by+dx+ez+h

d/dz:	2cz+ey+fx+i

*/

VECTOR	normal_to_quad(QUAD *quad, VECTOR p)
	{
	double	a = quad->a;
	double	b = quad->b;
	double	c = quad->c;
	double	d = quad->d;
	double	e = quad->e;
	double	f = quad->f;
	double	g = quad->g;
	double	h = quad->h;
	double	i = quad->i;
	double	x = p.x;
	double	y = p.y;
	double	z = p.z;
	double	ax = a * x;
	double	by = b * y;
	double	cz = c * z;
	VECTOR	normal;

	normal.x = ax + ax + d * y + f * z + g;
	normal.y = by + by + d * x + e * z + h;
	normal.z = cz + cz + e * y + f * x + i;
	return normal;
	}
/*...e*/
