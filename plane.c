/*

plane.c - Half-plane logic

Defined solid for :- ax+by+cz+d<=0

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
#define	_PLANE_
#include "plane.h"

/*...vrt\46\h:0:*/
/*...vvector\46\h:0:*/
/*...vsil\46\h:0:*/
/*...vplane\46\h:0:*/
/*...e*/

/*...screate_plane      \45\ create any general plane:0:*/
PLANE	*create_plane(double a, double b, double c, double d)
	{
	PLANE	*plane;

	if ( (plane = malloc(sizeof(PLANE))) == NULL )
		return NULL;

	plane->a = a;
	plane->b = b;
	plane->c = c;
	plane->d = d;
	return plane;
	}
/*...e*/
/*...screate_x_lt_plane \45\ create plane solid for x \60\\61\ some value:0:*/
PLANE	*create_x_lt_plane(double x)
	{
	return create_plane(1.0, 0.0, 0.0, -x);
	}
/*...e*/
/*...screate_x_gt_plane \45\ create plane solid for x \62\\61\ some value:0:*/
PLANE	*create_x_gt_plane(double x)
	{
	return create_plane(-1.0, 0.0, 0.0, x);
	}
/*...e*/
/*...screate_y_lt_plane \45\ create plane solid for y \60\\61\ some value:0:*/
PLANE	*create_y_lt_plane(double y)
	{
	return create_plane(0.0, 1.0, 0.0, -y);
	}
/*...e*/
/*...screate_y_gt_plane \45\ create plane solid for y \62\\61\ some value:0:*/
PLANE	*create_y_gt_plane(double y)
	{
	return create_plane(0.0, -1.0, 0.0, y);
	}
/*...e*/
/*...screate_z_lt_plane \45\ create plane solid for z \60\\61\ some value:0:*/
PLANE	*create_z_lt_plane(double z)
	{
	return create_plane(0.0, 0.0, 1.0, -z);
	}
/*...e*/
/*...screate_z_gt_plane \45\ create plane solid for z \62\\61\ some value:0:*/
PLANE	*create_z_gt_plane(double z)
	{
	return create_plane(0.0, 0.0, -1.0, z);
	}
/*...e*/
/*...scopy_plane        \45\ make a copy of a plane:0:*/
PLANE	*copy_plane(PLANE *plane)
	{
	PLANE	*copy;

	if ( (copy = malloc(sizeof(PLANE))) == NULL )
		return NULL;

	memcpy(copy, plane, sizeof(PLANE));
	return copy;
	}
/*...e*/
/*...sdestroy_plane     \45\ destroy a plane:0:*/
void	destroy_plane(PLANE *plane)
	{
	free(plane);
	}
/*...e*/

/*...strans_x_plane     \45\ translate by amount in x direction:0:*/
void	trans_x_plane(PLANE *plane, double t)
	{
	plane->d -= plane->a * t;
	}
/*...e*/
/*...strans_y_plane     \45\ translate by amount in y direction:0:*/
void	trans_y_plane(PLANE *plane, double t)
	{
	plane->d -= plane->b * t;
	}
/*...e*/
/*...strans_z_plane     \45\ translate by amount in z direction:0:*/
void	trans_z_plane(PLANE *plane, double t)
	{
	plane->d -= plane->c * t;
	}
/*...e*/
/*...sscale_x_plane     \45\ scale by factor in x direction:0:*/
void	scale_x_plane(PLANE *plane, double factor)
	{
	plane->a /= factor;
	}
/*...e*/
/*...sscale_y_plane     \45\ scale by factor in y direction:0:*/
void	scale_y_plane(PLANE *plane, double factor)
	{
	plane->b /= factor;
	}
/*...e*/
/*...sscale_z_plane     \45\ scale by factor in z direction:0:*/
void	scale_z_plane(PLANE *plane, double factor)
	{
	plane->c /= factor;
	}
/*...e*/
/*...srot_x_plane       \45\ rotate about x axis by given angle:0:*/
void	rot_x_plane(PLANE *plane, double angle)
	{
	double	b  = plane->b;
	double	c  = plane->c;
	double	ca = cos(angle);
	double	sa = sin(angle);

	plane->b = b * ca - c * sa;
	plane->c = b * sa + c * ca;
	}
/*...e*/
/*...srot_y_plane       \45\ rotate about y axis by given angle:0:*/
void	rot_y_plane(PLANE *plane, double angle)
	{
	double	c  = plane->c;
	double	a  = plane->a;
	double	ca = cos(angle);
	double	sa = sin(angle);

	plane->c = c * ca - a * sa;
	plane->a = c * sa + a * ca;
	}
/*...e*/
/*...srot_z_plane       \45\ rotate about z axis by given angle:0:*/
void	rot_z_plane(PLANE *plane, double angle)
	{
	double	a  = plane->a;
	double	b  = plane->b;
	double	ca = cos(angle);
	double	sa = sin(angle);

	plane->a = a * ca - b * sa;
	plane->b = a * sa + b * ca;
	}
/*...e*/

/*...sisects_reqd_plane \45\ max number of isects we will generate:0:*/
int	isects_reqd_plane(PLANE *plane)
	{
	plane=plane; /* Suppress 'unref arg' compiler warning */

	return 2;
	}
/*...e*/
/*...sintersect_plane   \45\ determine intersection range of t for plane:0:*/
/*

Any point along the line we are interested in is of the form p + tq.
							     ~    ~
Given:	ax+by+cz+d=0

Subs:	x = x +tx		y = y +ty		z = z +tz
	     p   q		     p   q		     p   q

Gives:	(ax +by +cz )t + ax +by +cz +d = 0
	   q   q   q       p   p   p

*/

/*...svalue_of_plane:0:*/
static double value_of_plane(void *shapeinfo, VECTOR v)
	{
	PLANE	*plane = (PLANE *) shapeinfo;

	return plane->a * v.x +
	       plane->b * v.y +
	       plane->c * v.z +
	       plane->d;
	}
/*...e*/

void	intersect_plane(PLANE *plane, VECTOR p, VECTOR q, SIL *sil)
	{
	double	a = plane->a;
	double	b = plane->b;
	double	c = plane->c;
	double	d = plane->d;
	double	coeff_of_t    = a * q.x + b * q.y + c * q.z;
	double	constant_term = a * p.x + b * p.y + c * p.z + d;

	intersect_linear_t(coeff_of_t, constant_term, p, q,
			   (void *) plane, value_of_plane, sil);
	}
/*...e*/
/*...snormal_to_plane   \45\ find normal of surface at a given point:0:*/
/*

Use partial derivatives to find normal at a given point.

Given:	ax+by+cz+d=0

d/dx:	a

d/dy:	b

d/dz:	c

*/

VECTOR	normal_to_plane(PLANE *plane)
	{
	VECTOR	normal;

	normal.x = plane->a;
	normal.y = plane->b;
	normal.z = plane->c;

	return normal;
	}
/*...e*/
