/*

biplane.c - Half-plane logic

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
#define	_BIPLANE_
#include "biplane.h"

/*...vrt\46\h:0:*/
/*...vvector\46\h:0:*/
/*...vsil\46\h:0:*/
/*...vbiplane\46\h:0:*/
/*...e*/

/*...screate_biplane      \45\ create any general biplane:0:*/
BIPLANE	*create_biplane(double a, double b, double c, double d1, double d2)
	{
	BIPLANE	*biplane;

	if ( (biplane = malloc(sizeof(BIPLANE))) == NULL )
		return NULL;

	biplane->a  = a;
	biplane->b  = b;
	biplane->c  = c;
	biplane->d1 = d1;
	biplane->d2 = d2;
	return biplane;
	}
/*...e*/
/*...screate_x_in_biplane \45\ create biplane solid for x1 \60\\61\ x \60\\61\ x2:0:*/
BIPLANE	*create_x_in_biplane(double x1, double x2)
	{
	return create_biplane(1.0, 0.0, 0.0, -x1, -x2);
	}
/*...e*/
/*...screate_y_in_biplane \45\ create biplane solid for y1 \60\\61\ y \60\\61\ y2:0:*/
BIPLANE	*create_y_in_biplane(double y1, double y2)
	{
	return create_biplane(0.0, 1.0, 0.0, -y1, -y2);
	}
/*...e*/
/*...screate_z_in_biplane \45\ create biplane solid for z1 \60\\61\ z \60\\61\ z2:0:*/
BIPLANE	*create_z_in_biplane(double z1, double z2)
	{
	return create_biplane(0.0, 0.0, 1.0, -z1, -z2);
	}
/*...e*/
/*...scopy_biplane        \45\ make a copy of a biplane:0:*/
BIPLANE	*copy_biplane(BIPLANE *biplane)
	{
	BIPLANE	*copy;

	if ( (copy = malloc(sizeof(BIPLANE))) == NULL )
		return NULL;

	memcpy(copy, biplane, sizeof(BIPLANE));
	return copy;
	}
/*...e*/
/*...sdestroy_biplane     \45\ destroy a biplane:0:*/
void	destroy_biplane(BIPLANE *biplane)
	{
	free(biplane);
	}
/*...e*/

/*...strans_x_biplane     \45\ translate by amount in x direction:0:*/
void	trans_x_biplane(BIPLANE *biplane, double t)
	{
	double at = biplane->a * t;
	biplane->d1 -= at;
	biplane->d2 -= at;
	}
/*...e*/
/*...strans_y_biplane     \45\ translate by amount in y direction:0:*/
void	trans_y_biplane(BIPLANE *biplane, double t)
	{
	double bt = biplane->b * t;
	biplane->d1 -= bt;
	biplane->d2 -= bt;
	}
/*...e*/
/*...strans_z_biplane     \45\ translate by amount in z direction:0:*/
void	trans_z_biplane(BIPLANE *biplane, double t)
	{
	double ct = biplane->c * t;
	biplane->d1 -= ct;
	biplane->d2 -= ct;
	}
/*...e*/
/*...sscale_x_biplane     \45\ scale by factor in x direction:0:*/
void	scale_x_biplane(BIPLANE *biplane, double factor)
	{
	biplane->a /= factor;
	}
/*...e*/
/*...sscale_y_biplane     \45\ scale by factor in y direction:0:*/
void	scale_y_biplane(BIPLANE *biplane, double factor)
	{
	biplane->b /= factor;
	}
/*...e*/
/*...sscale_z_biplane     \45\ scale by factor in z direction:0:*/
void	scale_z_biplane(BIPLANE *biplane, double factor)
	{
	biplane->c /= factor;
	}
/*...e*/
/*...srot_x_biplane       \45\ rotate about x axis by given angle:0:*/
void	rot_x_biplane(BIPLANE *biplane, double angle)
	{
	double	b  = biplane->b;
	double	c  = biplane->c;
	double	ca = cos(angle);
	double	sa = sin(angle);

	biplane->b = b * ca - c * sa;
	biplane->c = b * sa + c * ca;
	}
/*...e*/
/*...srot_y_biplane       \45\ rotate about y axis by given angle:0:*/
void	rot_y_biplane(BIPLANE *biplane, double angle)
	{
	double	c  = biplane->c;
	double	a  = biplane->a;
	double	ca = cos(angle);
	double	sa = sin(angle);

	biplane->c = c * ca - a * sa;
	biplane->a = c * sa + a * ca;
	}
/*...e*/
/*...srot_z_biplane       \45\ rotate about z axis by given angle:0:*/
void	rot_z_biplane(BIPLANE *biplane, double angle)
	{
	double	a  = biplane->a;
	double	b  = biplane->b;
	double	ca = cos(angle);
	double	sa = sin(angle);

	biplane->a = a * ca - b * sa;
	biplane->b = a * sa + b * ca;
	}
/*...e*/

/*...sisects_reqd_biplane \45\ max number of isects we will generate:0:*/
int	isects_reqd_biplane(BIPLANE *biplane)
	{
	biplane=biplane; /* Suppress 'unref arg' compiler warning */

	return 2;
	}
/*...e*/
/*...sintersect_biplane   \45\ determine intersection range of t for biplane:0:*/
/*

Any point along the line we are interested in is of the form p + tq.
							     ~    ~
Given:	ax+by+cz+d1 >0 and
	ax+by+cz+d2<=0

Gives:	-d1<ax+by+cz<=d2

Subs:	x = x +tx		y = y +ty		z = z +tz
	     p   q		     p   q		     p   q

Gives:	(ax +by +cz )t + ax +by +cz +d = 0
	   q   q   q       p   p   p

*/

#define	ZERO(x)		( ((x)>=-1.0e-100) && ((x)<=1.0e-100) )

void	intersect_biplane(BIPLANE *biplane, VECTOR p, VECTOR q, SIL *sil)
	{
	double	a  = biplane->a;
	double	b  = biplane->b;
	double	c  = biplane->c;
	double	d1 = biplane->d1;
	double	d2 = biplane->d2;
	double	coeff_of_t = a * q.x + b * q.y + c * q.z;
	double	axbycz     = a * p.x + b * p.y + c * p.z;

	if ( ZERO(coeff_of_t) )
		/* No intersection with solid surface */
		{
		if ( axbycz > -d1 && axbycz <= -d2 )
			/* Completely in solid */
			{
			sil->n_sis            = 2;
			sil->sis[0].t        = -INFINITE;
			sil->sis[0].entering = TRUE;
			sil->sis[1].t        = INFINITE;
			sil->sis[1].entering = FALSE;
			}
		else
			/* Completely out of solid */
			sil->n_sis            = 0;
		}
	else
		/* Intersects solid exactly twice */
		{
		double	t1 = - (axbycz + d1) / coeff_of_t;
		double	t2 = - (axbycz + d2) / coeff_of_t;

		sil->n_sis            = 2;
		sil->sis[0].entering = TRUE;
		sil->sis[1].entering = FALSE;
		if ( t1 < t2 )
			{
			sil->sis[0].t = t1;
			sil->sis[1].t = t2;
			}
		else
			{
			sil->sis[0].t = t2;
			sil->sis[1].t = t1;
			}
		}
	}
/*...e*/
/*...snormal_to_biplane   \45\ find normal of surface at a given point:0:*/
/*
This is much like the simple plane case except that we know need to know the
point of intersection in order to know which of the 2 half-plane surfaces were
involved. We compute a d value, and if its closest to d2, then great, else
the normal is the other way around.
*/

VECTOR	normal_to_biplane(BIPLANE *biplane, VECTOR p)
	{
	VECTOR	normal;
	double	a = biplane->a;
	double	b = biplane->b;
	double	c = biplane->c;
	double	d = - ( a * p.x + b * p.y + c * p.z );

	if ( fabs(d - biplane->d2) < fabs(d - biplane->d1) )
		{
		normal.x = a;
		normal.y = b;
		normal.z = c;
		}
	else
		{
		normal.x = -a;
		normal.y = -b;
		normal.z = -c;
		}

	return normal;
	}
/*...e*/
