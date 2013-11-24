/*

shape.c - General shapes

*/

/*...sincludes:0:*/
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <malloc.h>
#include <memory.h>
#include <math.h>
#include "rt.h"
#include "fio.h"
#include "tex.h"
#include "vector.h"
#include "rgbvec.h"
#include "col.h"
#include "surf.h"
#include "sil.h"
#include "plane.h"
#include "biplane.h"
#include "sphere.h"
#include "quad.h"
#define	_SHAPE_
#include "shape.h"

/*...vrt\46\h:0:*/
/*...vfio\46\h:0:*/
/*...vtex\46\h:0:*/
/*...vvector\46\h:0:*/
/*...vrgbvec\46\h:0:*/
/*...vcol\46\h:0:*/
/*...vsurf\46\h:0:*/
/*...vsil\46\h:0:*/
/*...vplane\46\h:0:*/
/*...vbiplane\46\h:0:*/
/*...vsphere\46\h:0:*/
/*...vquad\46\h:0:*/
/*...vshape\46\h:0:*/
/*...e*/

/*...screate_plane_shape   \45\ create a shape tree consisting of a single plane:0:*/
SHAPE	*create_plane_shape(PLANE *plane, SURF *surf)
	{
	SHAPE	*shape;

	if ( (shape = malloc(sizeof(SHAPE))) == NULL )
		return NULL;

	shape->stype   = STYPE_PLANE;
	shape->u.plane = plane;
	shape->surf    = surf;
	return shape;
	}
/*...e*/
/*...screate_biplane_shape \45\ create a shape tree consisting of a single plane:0:*/
SHAPE	*create_biplane_shape(BIPLANE *biplane, SURF *surf)
	{
	SHAPE	*shape;

	if ( (shape = malloc(sizeof(SHAPE))) == NULL )
		return NULL;

	shape->stype     = STYPE_BIPLANE;
	shape->u.biplane = biplane;
	shape->surf      = surf;
	return shape;
	}
/*...e*/
/*...screate_sphere_shape  \45\ create a shape tree consisting of a single sphere:0:*/
SHAPE	*create_sphere_shape(SPHERE *sphere, SURF *surf)
	{
	SHAPE	*shape;

	if ( (shape = malloc(sizeof(SHAPE))) == NULL )
		return NULL;

	shape->stype    = STYPE_SPHERE;
	shape->u.sphere = sphere;
	shape->surf     = surf;
	return shape;
	}
/*...e*/
/*...screate_quad_shape    \45\ create a shape tree consisting of a single quadratic:0:*/
SHAPE	*create_quad_shape(QUAD *quad, SURF *surf)
	{
	SHAPE	*shape;

	if ( (shape = malloc(sizeof(SHAPE))) == NULL )
		return NULL;

	shape->stype  = STYPE_QUAD;
	shape->u.quad = quad;
	shape->surf   = surf;
	return shape;
	}
/*...e*/
/*...screate_bin_shape     \45\ create a binary shape tree made of 2 shapes:0:*/
SHAPE	*create_bin_shape(STYPE stype, SHAPE *shape0, SHAPE *shape1)
	{
	SHAPE	*shape;

	if ( (shape = malloc(sizeof(SHAPE))) == NULL )
		return NULL;

	shape->stype       = stype;
	shape->u.shapes[0] = shape0;
	shape->u.shapes[1] = shape1;
	return shape;
	}
/*...e*/
/*...scopy_shape           \45\ make a complete copy of a shape tree:0:*/
SHAPE	*copy_shape(SHAPE *shape)
	{
	SHAPE	*copy;

	if ( (copy = malloc(sizeof(SHAPE))) == NULL )
		return NULL;

	copy->stype = shape->stype;

	if ( shape->stype <= STYPE_QUAD )
		/* Leaf node */
		{
		if ( (copy->surf = copy_surf(shape->surf)) == NULL )
			{
			free(copy);
			return NULL;
			}
		switch ( shape->stype )
			{
/*...sSTYPE_PLANE:24:*/
case STYPE_PLANE:
	if ( (copy->u.plane = copy_plane(shape->u.plane)) == NULL )
		{
		destroy_surf(copy->surf);
		free(copy);
		return NULL;
		}
	break;
/*...e*/
/*...sSTYPE_BIPLANE:24:*/
case STYPE_BIPLANE:
	if ( (copy->u.biplane = copy_biplane(shape->u.biplane)) == NULL )
		{
		destroy_surf(copy->surf);
		free(copy);
		return NULL;
		}
	break;
/*...e*/
/*...sSTYPE_SPHERE:24:*/
case STYPE_SPHERE:
	if ( (copy->u.sphere = copy_sphere(shape->u.sphere)) == NULL )
		{
		destroy_surf(copy->surf);
		free(copy);
		return NULL;
		}
	break;
/*...e*/
/*...sSTYPE_QUAD:24:*/
case STYPE_QUAD:
	if ( (copy->u.quad = copy_quad(shape->u.quad)) == NULL )
		{
		destroy_surf(copy->surf);
		free(copy);
		return NULL;
		}
	break;
/*...e*/
			}
		}
	else
/*...sboolean combinations:16:*/
{
if ( (copy->u.shapes[0] = copy_shape(shape->u.shapes[0])) == NULL )
	{
	free(copy);
	return NULL;
	}
if ( (copy->u.shapes[1] = copy_shape(shape->u.shapes[1])) == NULL )
	{
	free(copy->u.shapes[0]);
	free(copy);
	return NULL;
	}
}
/*...e*/

	return copy;
	}
/*...e*/
/*...sdestroy_shape        \45\ delete a shape tree:0:*/
void	destroy_shape(SHAPE *shape)
	{
	if ( shape->stype <= STYPE_QUAD )
		/* Leaf node */
		{
		destroy_surf(shape->surf);
		switch ( shape->stype )
			{
/*...sSTYPE_PLANE:24:*/
case STYPE_PLANE:
	destroy_plane(shape->u.plane);
	break;
/*...e*/
/*...sSTYPE_BIPLANE:24:*/
case STYPE_BIPLANE:
	destroy_biplane(shape->u.biplane);
	break;
/*...e*/
/*...sSTYPE_SPHERE:24:*/
case STYPE_SPHERE:
	destroy_sphere(shape->u.sphere);
	break;
/*...e*/
/*...sSTYPE_QUAD:24:*/
case STYPE_QUAD:
	destroy_quad(shape->u.quad);
	break;
/*...e*/
			}
		}
	else
		{
		destroy_shape(shape->u.shapes[0]);
		destroy_shape(shape->u.shapes[1]);
		}

	free(shape);
	}
/*...e*/
/*...ssphere_to_quad       \45\ convert sphere to quad:0:*/
/*
Spheres cannot be rescaled in each direction and remain a sphere.
Hence this code will convert a sphere into the equivelent quadratic.

			     2      2      2
	Sphere:		(x-a) +(y-b) +(z-c) -d <=0

			 2      2  2      2  2      2
			x -2ax+a +y -2by+b +z -2cz+c -d <= 0

			  2   2   2
	Quadratic:	ax +by +cz +dxy+eyz+fzx+gx+hy+iz+j<=0

*/

static BOOLEAN sphere_to_quad(SHAPE *shape)
	{
	SPHERE *sphere = shape->u.sphere;
	QUAD *quad;

	if ( (quad = create_quad(1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
	     -2.0 * sphere->a, -2.0 * sphere->b, -2.0 * sphere->c,
	     - sphere->d)) == NULL )
		return FALSE;

	destroy_sphere(shape->u.sphere);
	shape->stype = STYPE_QUAD;
	shape->u.quad = quad;

	return TRUE;
	}
/*...e*/

/*...strans_x              \45\ translate by amount in x direction:0:*/
void	trans_x(SHAPE *shape, double t)
	{
	if ( shape->stype <= STYPE_QUAD )
		{
		switch ( shape->stype )
			{
			case STYPE_PLANE:
				trans_x_plane(shape->u.plane, t);
				break;
			case STYPE_BIPLANE:
				trans_x_biplane(shape->u.biplane, t);
				break;
			case STYPE_SPHERE:
				trans_x_sphere(shape->u.sphere, t);
				break;
			case STYPE_QUAD:
				trans_x_quad(shape->u.quad, t);
				break;
			}
		trans_x_surf(shape->surf, t);
		}
	else
		{
		trans_x(shape->u.shapes[0], t);
		trans_x(shape->u.shapes[1], t);
		}
	}
/*...e*/
/*...strans_y              \45\ translate by amount in y direction:0:*/
void	trans_y(SHAPE *shape, double t)
	{
	if ( shape->stype <= STYPE_QUAD )
		{
		switch ( shape->stype )
			{
			case STYPE_PLANE:
				trans_y_plane(shape->u.plane, t);
				break;
			case STYPE_BIPLANE:
				trans_y_biplane(shape->u.biplane, t);
				break;
			case STYPE_SPHERE:
				trans_y_sphere(shape->u.sphere, t);
				break;
			case STYPE_QUAD:
				trans_y_quad(shape->u.quad, t);
				break;
			}
		trans_y_surf(shape->surf, t);
		}
	else
		{
		trans_y(shape->u.shapes[0], t);
		trans_y(shape->u.shapes[1], t);
		}
	}
/*...e*/
/*...strans_z              \45\ translate by amount in z direction:0:*/
void	trans_z(SHAPE *shape, double t)
	{
	if ( shape->stype <= STYPE_QUAD )
		{
		switch ( shape->stype )
			{
			case STYPE_PLANE:
				trans_z_plane(shape->u.plane, t);
				break;
			case STYPE_BIPLANE:
				trans_z_biplane(shape->u.biplane, t);
				break;
			case STYPE_SPHERE:
				trans_z_sphere(shape->u.sphere, t);
				break;
			case STYPE_QUAD:
				trans_z_quad(shape->u.quad, t);
				break;
			}
		trans_z_surf(shape->surf, t);
		}
	else
		{
		trans_z(shape->u.shapes[0], t);
		trans_z(shape->u.shapes[1], t);
		}
	}
/*...e*/
/*...strans                \45\ translate by vector:0:*/
void	trans(SHAPE *shape, VECTOR v)
	{
	trans_x(shape, v.x);
	trans_y(shape, v.y);
	trans_z(shape, v.z);
	}
/*...e*/
/*...sscale_x              \45\ scale by factor in x direction:0:*/
BOOLEAN	scale_x(SHAPE *shape, double factor)
	{
	if ( factor == 1.0 )
		return TRUE;

	if ( shape->stype <= STYPE_QUAD )
		{
		switch ( shape->stype )
			{
			case STYPE_PLANE:
				scale_x_plane(shape->u.plane, factor);
				break;
			case STYPE_BIPLANE:
				scale_x_biplane(shape->u.biplane, factor);
				break;
			case STYPE_SPHERE:
				if ( !sphere_to_quad(shape) )
					return FALSE;
				scale_x_quad(shape->u.quad, factor);
				break;
			case STYPE_QUAD:
				scale_x_quad(shape->u.quad, factor);
				break;
			}
		scale_x_surf(shape->surf, factor);
		}
	else
		{
		if ( !scale_x(shape->u.shapes[0], factor) ||
		     !scale_x(shape->u.shapes[1], factor) )
			return FALSE;
		}
	return TRUE;
	}
/*...e*/
/*...sscale_y              \45\ scale by factor in y direction:0:*/
BOOLEAN	scale_y(SHAPE *shape, double factor)
	{
	if ( factor == 1.0 )
		return TRUE;

	if ( shape->stype <= STYPE_QUAD )
		{
		switch ( shape->stype )
			{
			case STYPE_PLANE:
				scale_y_plane(shape->u.plane, factor);
				break;
			case STYPE_BIPLANE:
				scale_y_biplane(shape->u.biplane, factor);
				break;
			case STYPE_SPHERE:
				if ( !sphere_to_quad(shape) )
					return FALSE;
				scale_y_quad(shape->u.quad, factor);
				break;
			case STYPE_QUAD:
				scale_y_quad(shape->u.quad, factor);
				break;
			}
		scale_y_surf(shape->surf, factor);
		}
	else
		{
		if ( !scale_y(shape->u.shapes[0], factor) ||
		     !scale_y(shape->u.shapes[1], factor) )
			return FALSE;
		}
	return TRUE;
	}
/*...e*/
/*...sscale_z              \45\ scale by factor in z direction:0:*/
BOOLEAN	scale_z(SHAPE *shape, double factor)
	{
	if ( factor == 1.0 )
		return TRUE;

	if ( shape->stype <= STYPE_QUAD )
		{
		switch ( shape->stype )
			{
			case STYPE_PLANE:
				scale_z_plane(shape->u.plane, factor);
				break;
			case STYPE_BIPLANE:
				scale_z_biplane(shape->u.biplane, factor);
				break;
			case STYPE_SPHERE:
				if ( !sphere_to_quad(shape) )
					return FALSE;
				scale_z_quad(shape->u.quad, factor);
				break;
			case STYPE_QUAD:
				scale_z_quad(shape->u.quad, factor);
				break;
			}
		scale_z_surf(shape->surf, factor);
		}
	else
		{
		if ( !scale_z(shape->u.shapes[0], factor) ||
		     !scale_z(shape->u.shapes[1], factor) )
			return FALSE;
		}
	return TRUE;
	}
/*...e*/
/*...sscale                \45\ scale by vector factor:0:*/
BOOLEAN	scale(SHAPE *shape, VECTOR factor)
	{
	if ( shape->stype == STYPE_SPHERE &&
	     factor.x == factor.y && factor.y == factor.z )
		{
		scale_sphere(shape->u.sphere, factor.x);
		scale_x_surf(shape->surf, factor.x);
		scale_y_surf(shape->surf, factor.y);
		scale_z_surf(shape->surf, factor.z);
		}
	else
		{
		if ( !scale_x(shape, factor.x) ||
		     !scale_y(shape, factor.y) ||
		     !scale_z(shape, factor.z) )
			return FALSE;
		}
	return TRUE;
	}
/*...e*/
/*...srot_x                \45\ rotate about x axis by given angle:0:*/
void	rot_x(SHAPE *shape, double angle)
	{
	if ( shape->stype <= STYPE_QUAD )
		{
		switch ( shape->stype )
			{
			case STYPE_PLANE:
				rot_x_plane(shape->u.plane, angle);
				break;
			case STYPE_BIPLANE:
				rot_x_biplane(shape->u.biplane, angle);
				break;
			case STYPE_SPHERE:
				rot_x_sphere(shape->u.sphere, angle);
				break;
			case STYPE_QUAD:
				rot_x_quad(shape->u.quad, angle);
				break;
			}
		rot_x_surf(shape->surf, angle);
		}
	else
		{
		rot_x(shape->u.shapes[0], angle);
		rot_x(shape->u.shapes[1], angle);
		}
	}
/*...e*/
/*...srot_y                \45\ rotate about y axis by given angle:0:*/
void	rot_y(SHAPE *shape, double angle)
	{
	if ( shape->stype <= STYPE_QUAD )
		{
		switch ( shape->stype )
			{
			case STYPE_PLANE:
				rot_y_plane(shape->u.plane, angle);
				break;
			case STYPE_BIPLANE:
				rot_y_biplane(shape->u.biplane, angle);
				break;
			case STYPE_SPHERE:
				rot_y_sphere(shape->u.sphere, angle);
				break;
			case STYPE_QUAD:
				rot_y_quad(shape->u.quad, angle);
				break;
			}
		rot_y_surf(shape->surf, angle);
		}
	else
		{
		rot_y(shape->u.shapes[0], angle);
		rot_y(shape->u.shapes[1], angle);
		}
	}
/*...e*/
/*...srot_z                \45\ rotate about z axis by given angle:0:*/
void	rot_z(SHAPE *shape, double angle)
	{
	if ( shape->stype <= STYPE_QUAD )
		{
		switch ( shape->stype )
			{
			case STYPE_PLANE:
				rot_z_plane(shape->u.plane, angle);
				break;
			case STYPE_BIPLANE:
				rot_z_biplane(shape->u.biplane, angle);
				break;
			case STYPE_SPHERE:
				rot_z_sphere(shape->u.sphere, angle);
				break;
			case STYPE_QUAD:
				rot_z_quad(shape->u.quad, angle);
				break;
			}
		rot_z_surf(shape->surf, angle);
		}
	else
		{
		rot_z(shape->u.shapes[0], angle);
		rot_z(shape->u.shapes[1], angle);
		}
	}
/*...e*/

/*...sresurf               \45\ change the surface of a shape:0:*/
BOOLEAN resurf(SHAPE *shape, SURF *surf)
	{
	if ( shape->stype <= STYPE_QUAD )
		{
		SURF *surf_copy;

		if ( (surf_copy = copy_surf(surf)) == NULL )
			return FALSE;

		destroy_surf(shape->surf);
		shape->surf = surf_copy;
		}
	else
		{
		if ( !resurf(shape->u.shapes[0], surf) )
			return FALSE;
		if ( !resurf(shape->u.shapes[1], surf) )
			return FALSE;
		}

	return TRUE;
	}
/*...e*/

/*...sis_empty_isectl      \45\ is intersection list empty:0:*/
BOOLEAN	is_empty_isectl(ISECTL *il)
	{
	return il->n_isects == 0;
	}
/*...e*/
/*...sis_solid_isectl      \45\ is intersection list solid:0:*/
/*
We will allow t from -INFINITE onwards or t from -INFINITE to INFINITE.
*/

BOOLEAN	is_solid_isectl(ISECTL *il)
	{
	if ( il->n_isects > 2 )
		return FALSE;

	if ( il->isects[0].t != -INFINITE )
		return FALSE;

	return il->n_isects == 1 || il->isects[1].t == INFINITE;
	}
/*...e*/
/*...st_after_isectl       \45\ eliminate all intersections before a t value:0:*/
void	t_after_isectl(ISECTL *il, double t)
	{
	int	i;

	for ( i = 0; i < il->n_isects; i++ )
		if ( il->isects[i].t >= t )
			break;

	if ( i == il->n_isects )
		/* All behind t */
		il->n_isects = 0;
	else if ( i != 0 )
		/* Some behind and some after t */
		{
		int	j = 0;

		if ( il->isects[i].entering == FALSE )
			/* Had better make an isect case first */
			{
			il->isects[j  ]   = il->isects[i - 1];
			il->isects[j++].t = t;
			}

		while ( i < il->n_isects )
			il->isects[j++] = il->isects[i++];
		il->n_isects = j;
		}
	}
/*...e*/

/*...screate_isectl        \45\ create an isectl:0:*/
ISECTL	*create_isectl(int n_isects)
	{
	return malloc(sizeof(ISECTL) + (n_isects - 1) * sizeof(ISECT));
	}
/*...e*/
/*...sdestroy_isectl       \45\ destroy an isectl:0:*/
void	destroy_isectl(ISECTL *il)
	{
	free(il);
	}
/*...e*/

/*...spreprocess_shape     \45\ preprocess shape tree to save work when tracing:0:*/
/*...sextent_shape:0:*/
/*
If we can find shapes that do not overlap, we can avoid tracing time later.
This type of function is only made possible because the code that implements
spheres etc. export their internal representations of the shapes.
*/

/*...sspheres_overlap:0:*/
/*
2 spheres overlap if the sum of the radii > distance between the centres.
We square both sides of this equation.
*/

static BOOLEAN spheres_overlap(SPHERE *sphere_a, SPHERE *sphere_b)
	{
	double	dx = sphere_a->a - sphere_b->a;
	double	dy = sphere_a->b - sphere_b->b;
	double	dz = sphere_a->c - sphere_b->c;
	double	rr = sphere_a->d + sphere_b->d;

	return rr*rr >= dx*dx + dy*dy + dz*dz;
	}
/*...e*/

typedef struct { VECTOR xmin, xmax; } EXTENT;

/*...sextent_infinite:0:*/
static EXTENT extent_infinite(void)
	{
	EXTENT	extent;

	extent.xmin.x = extent.xmin.y = extent.xmin.z = -INFINITE;
	extent.xmax.x = extent.xmax.y = extent.xmax.z =  INFINITE;
	return extent;
	}
/*...e*/
/*...sextent_of_plane:0:*/
/*

ax+by+cz+d<= 0

Quite a lot of models have their plane faces aligned along x,y y,z or x,z
planes. Hence this code should discover that fact quite often.

*/

static EXTENT extent_of_plane(PLANE *plane)
	{
	double	a = plane->a;
	double	b = plane->b;
	double	c = plane->c;
	double	d = plane->d;
	EXTENT	extent = extent_infinite();

	if ( b == 0.0 && c == 0.0 )
		/* ax<=-d */
		{
		if ( a >= 0.0 )
			extent.xmax.x = -d/a;
		else
			extent.xmin.x = -d/a;
		}
	else if ( a == 0.0 && c == 0.00 )
		/* by<=-d */
		{
		if ( b >= 0.0 )
			extent.xmax.y = -d/b;
		else
			extent.xmin.y = -d/b;
		}
	else if ( a == 0.0 && b == 0.0 )
		/* cz<=-d */
		{
		if ( c >= 0.0 )
			extent.xmax.z = -d/c;
		else
			extent.xmin.z = -d/c;
		}

	return extent;
	}
/*...e*/
/*...sextent_of_biplane:0:*/
/*

ax+by+cz+d1 >  0 and
ax+by+cz+d2 <= 0

Quite a lot of models have their plane faces aligned along x,y y,z or x,z
planes. Hence this code should discover that fact quite often.

*/

static EXTENT extent_of_biplane(BIPLANE *biplane)
	{
	double	a  = biplane->a;
	double	b  = biplane->b;
	double	c  = biplane->c;
	double	d1 = biplane->d1;
	double	d2 = biplane->d2;
	EXTENT	extent = extent_infinite();

	if ( b == 0.0 && c == 0.0 )
		/* ax>-d1 and ax<=-d2 */
		{
		if ( a >= 0.0 )
			{
			extent.xmin.x = -d1/a;
			extent.xmax.x = -d2/a;
			}
		else
			{
			extent.xmin.x = -d2/a;
			extent.xmax.x = -d1/a;
			}
		}
	else if ( a == 0.0 && c == 0.00 )
		/* by>-d1 and by<=-d2 */
		{
		if ( b >= 0.0 )
			{
			extent.xmin.y = -d1/b;
			extent.xmax.y = -d2/b;
			}
		else
			{
			extent.xmin.y = -d2/b;
			extent.xmax.y = -d1/b;
			}
		}
	else if ( a == 0.0 && b == 0.0 )
		/* cz>-d1 and cz<=-d2 */
		{
		if ( c >= 0.0 )
			{
			extent.xmin.z = -d1/c;
			extent.xmax.z = -d2/c;
			}
		else
			{
			extent.xmin.z = -d2/c;
			extent.xmax.z = -d1/c;
			}
		}

	return extent;
	}
/*...e*/
/*...sextent_of_sphere:0:*/
static EXTENT extent_of_sphere(SPHERE *sphere)
	{
	double	a = sphere->a;
	double	b = sphere->b;
	double	c = sphere->c;
	double	d = sphere->d;
	EXTENT	extent;

	extent.xmin.x = a - d; extent.xmax.x = a + d;
	extent.xmin.y = b - d; extent.xmax.y = b + d;
	extent.xmin.z = c - d; extent.xmax.z = c + d;

	return extent;
	}
/*...e*/
/*...sextents_overlap:0:*/
static BOOLEAN extents_overlap(EXTENT *extent_a, EXTENT *extent_b)
	{
	return extent_a->xmax.x > extent_b->xmin.x &&
	       extent_a->xmax.y > extent_b->xmin.y &&
	       extent_a->xmax.z > extent_b->xmin.z &&
	       extent_b->xmax.x > extent_a->xmin.x &&
	       extent_b->xmax.y > extent_a->xmin.y &&
	       extent_b->xmax.z > extent_a->xmin.z ;
	}
/*...e*/
/*...sextents_union:0:*/
static EXTENT extents_union(EXTENT *extent_a, EXTENT *extent_b)
	{
	EXTENT extent;

	extent.xmin.x = min(extent_a->xmin.x, extent_b->xmin.x);
	extent.xmin.y = min(extent_a->xmin.y, extent_b->xmin.y);
	extent.xmin.z = min(extent_a->xmin.z, extent_b->xmin.z);
	extent.xmax.x = max(extent_a->xmax.x, extent_b->xmax.x);
	extent.xmax.y = max(extent_a->xmax.y, extent_b->xmax.y);
	extent.xmax.z = max(extent_a->xmax.z, extent_b->xmax.z);

	return extent;
	}
/*...e*/
/*...sextents_isect:0:*/
static EXTENT extents_isect(EXTENT *extent_a, EXTENT *extent_b)
	{
	EXTENT	extent;

	extent.xmin.x = max(extent_a->xmin.x, extent_b->xmin.x);
	extent.xmin.y = max(extent_a->xmin.y, extent_b->xmin.y);
	extent.xmin.z = max(extent_a->xmin.z, extent_b->xmin.z);
	extent.xmax.x = min(extent_a->xmax.x, extent_b->xmax.x);
	extent.xmax.y = min(extent_a->xmax.y, extent_b->xmax.y);
	extent.xmax.z = min(extent_a->xmax.z, extent_b->xmax.z);

	return extent;
	}
/*...e*/

static EXTENT extent_shape(SHAPE *shape)
	{
	static EXTENT dummy_extent = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };

	if ( shape->stype <= STYPE_QUAD )
		{
		switch ( shape->stype )
			{
/*...sSTYPE_PLANE:24:*/
case STYPE_PLANE:
	return extent_of_plane(shape->u.plane);
/*...e*/
/*...sSTYPE_BIPLANE:24:*/
case STYPE_BIPLANE:
	return extent_of_biplane(shape->u.biplane);
/*...e*/
/*...sSTYPE_QUAD \45\ too difficult:24:*/
case STYPE_QUAD:
	return extent_infinite();
/*...e*/
/*...sSTYPE_SPHERE:24:*/
case STYPE_SPHERE:
	return extent_of_sphere(shape->u.sphere);
/*...e*/
			}
		}
	else if ( shape->stype == STYPE_EXTENT )
		{
		extent_shape(shape->u.shapes[0]);
		/* Performed for side effect of extent_shape() call */

		return extent_shape(shape->u.shapes[1]);
		}
	else
		{
		SHAPE	*shape_a = shape->u.shapes[0];
		SHAPE	*shape_b = shape->u.shapes[1];
		EXTENT	extent_a = extent_shape(shape_a);
		EXTENT	extent_b = extent_shape(shape_b);
		EXTENT	extent;

		switch ( shape->stype )
			{
/*...sSTYPE_UNION\44\ STYPE_SDIFF:24:*/
case STYPE_UNION:
case STYPE_SDIFF:
	extent = extents_union(&extent_a, &extent_b);
	break;
/*...e*/
/*...sSTYPE_ISECT:24:*/
case STYPE_ISECT:
	extent = extents_isect(&extent_a, &extent_b);
	break;
/*...e*/
/*...sSTYPE_DIFF:24:*/
case STYPE_DIFF:
	extent = extent_a;
	break;
/*...e*/
			}

		if ( shape_a->stype == STYPE_SPHERE &&
		     shape_b->stype == STYPE_SPHERE )
			shape->overlap = spheres_overlap(
				shape_a->u.sphere, shape_b->u.sphere);
		else
			shape->overlap = extents_overlap(&extent_a, &extent_b);

		return extent;
		}

	return dummy_extent; /* Keep fussy C compiler happy */
	}
/*...e*/
/*...sn_isectls_reqd_shape:0:*/
/*
We use this function so that we can work out how many intersection lists
we will need to pre-allocate for tracing the shape.
*/

int	n_isectls_reqd_shape(SHAPE *shape)
	{
	int nest0, nest1;

	if ( shape->stype <= STYPE_QUAD )
		return 1;

	nest0 = n_isectls_reqd_shape(shape->u.shapes[0]);
	nest1 = n_isectls_reqd_shape(shape->u.shapes[1]);

	return 2 + max(nest0, nest1);
	}
/*...e*/
/*...sisects_reqd_shape:0:*/
/*
Determine largest needed intersection list.
*/

int	isects_reqd_shape(SHAPE *shape)
	{
	if ( shape->stype <= STYPE_QUAD )
		switch ( shape->stype )
			{
/*...sSTYPE_PLANE:24:*/
case STYPE_PLANE:
	return isects_reqd_plane(shape->u.plane);
/*...e*/
/*...sSTYPE_BIPLANE:24:*/
case STYPE_BIPLANE:
	return isects_reqd_biplane(shape->u.biplane);
/*...e*/
/*...sSTYPE_SPHERE:24:*/
case STYPE_SPHERE:
	return isects_reqd_sphere(shape->u.sphere);
/*...e*/
/*...sSTYPE_QUAD:24:*/
case STYPE_QUAD:
	return isects_reqd_quad(shape->u.quad);
/*...e*/
			}
	else
		{
		int	reqd0 = isects_reqd_shape(shape->u.shapes[0]);
		int	reqd1 = isects_reqd_shape(shape->u.shapes[1]);

		switch ( shape->stype )
			{
/*...sSTYPE_UNION:24:*/
case STYPE_UNION:
	return reqd0 + reqd1;
/*...e*/
/*...sSTYPE_ISECT:24:*/
case STYPE_ISECT:
	return ( shape->overlap ) ? reqd0 + reqd1 : 0;
/*...e*/
/*...sSTYPE_DIFF:24:*/
case STYPE_DIFF:
	return ( shape->overlap ) ? reqd0 + reqd1 : reqd0;
/*...e*/
/*...sSTYPE_SDIFF:24:*/
case STYPE_SDIFF:
	return reqd0 + reqd1;
/*...e*/
/*...sSTYPE_EXTENT:24:*/
case STYPE_EXTENT:
	return max(reqd0, reqd1);
/*...e*/
			}
		}
	return 0; /* Keep fussy C compiler happy */
	}
/*...e*/

void	preprocess_shape(SHAPE *shape, int *n_isectls, int *n_isects)
	{
	extent_shape(shape);
		/* The side effect is to label the shape tree */

	*n_isectls = n_isectls_reqd_shape(shape);
	*n_isects  = isects_reqd_shape(shape);
	}
/*...e*/
/*...sintersect_shape      \45\ intersect with a shape to give a ISECTL:0:*/
/*
Given a shape tree, find the intersection list of a vector with it.
If the shape tree is not a leaf node, then combine the intersection lists of
the subtrees. We can make several optimisations in this area.
*/

/*...smake_empty_isectl:0:*/
static void make_empty_isectl(ISECTL *il)
	{
	il->n_isects = 0;
	}
/*...e*/
/*...scopy_isectl:0:*/
static void copy_isectl(ISECTL *il_dst, ISECTL *il_src)
	{
	int	i;

	il_dst->n_isects = il_src->n_isects;
	for ( i = 0; i < il_dst->n_isects; i++ )
		il_dst->isects[i] = il_src->isects[i];
	}
/*...e*/
/*...scombine_isectl:0:*/
/*
Given 2 intersection lists, produce a new one that is a combination of them.

The in_old flag is used to determine if the point of intersection changes
meaning from in->out to out->in, or vice-versa. If this happens, the
sense of the normal vector must change :-

         ..... .....                .....
        .     .     .              .     .
       .  A  . .  B  .            . A-B .
      .     .   .     .          .     .
    <-.   <-.   .->   .->      <-.     .-> This vector changes direction!
      .     .   .     .          .     .
       .     . .     .            .     .
        .     .     .              .     .
         ..... .....                .....

*/

static void combine_isectl(
	BOOLEAN (*combine)(BOOLEAN in_a, BOOLEAN in_b),
	ISECTL *il_a, ISECTL *il_b, ISECTL *il)
	{
	BOOLEAN	in_a = FALSE;	/* When t = -INFINITE, in_a = FALSE */
	BOOLEAN	in_b = FALSE;	/* When t = -INFINITE, in_b = FALSE */
	BOOLEAN	in = FALSE;	/* Therefore combination starts as FALSE */
	int	ptr_a = 0;
	int	ptr_b = 0;
	BOOLEAN	in_new, in_old;

	il->n_isects = 0;

	/* Work through both a and b, looking at nearest ones first */

	while ( ptr_a < il_a->n_isects && ptr_b < il_b->n_isects )
		{
		ISECT	*isect;
		double	t_a = il_a->isects[ptr_a].t;
		double	t_b = il_b->isects[ptr_b].t;

		if ( t_a < t_b )
			{
			in_old = in_a = !in_a;
			isect = &(il_a->isects[ptr_a++]);
			}
		else if ( t_a > t_b )
			{
			in_old = in_b = !in_b;
			isect = &(il_b->isects[ptr_b++]);
			}
		else
			/* Two surfaces at exactly the same place */
			/* Not a very frequent event, but problematical */
			/* Just label intersection arbitrarily as with B */
			{
			in_a = !in_a; ptr_a++;
			in_old = in_b = !in_b;
			isect = &(il_b->isects[ptr_b++]);
			}

		if ( (in_new = (*combine)(in_a, in_b)) != in )
			/* Need to keep a record of this transition */
			{
			il->isects[il->n_isects] = *isect;
			il->isects[il->n_isects].entering = in = in_new;
			if ( in_new != in_old )
				il->isects[il->n_isects].negate_normal ^= TRUE;
			il->n_isects++;
			}
		}

	/* Either a or b is exhausted, so one of a or b may be left */

	while ( ptr_a < il_a->n_isects )
		{
		in_old = in_a = !in_a;
		if ( (in_new = (*combine)(in_a, in_b)) != in )
			/* Need to keep a record of this transition */
			{
			il->isects[il->n_isects] = il_a->isects[ptr_a];
			il->isects[il->n_isects].entering = in = in_new;
			if ( in_new != in_old )
				il->isects[il->n_isects].negate_normal ^= TRUE;
			il->n_isects++;
			}
		ptr_a++;
		}

	while ( ptr_b < il_b->n_isects )
		{
		in_old = in_b = !in_b;
		if ( (in_new = (*combine)(in_a, in_b)) != in )
			/* Need to keep a record of this transition */
			{
			il->isects[il->n_isects] = il_b->isects[ptr_b];
			il->isects[il->n_isects].entering = in = in_new;
			if ( in_new != in_old )
				il->isects[il->n_isects].negate_normal ^= TRUE;
			il->n_isects++;
			}
		ptr_b++;
		}
	}
/*...e*/
/*...sunion_isectl:0:*/
/*...scombine_union:0:*/
static BOOLEAN combine_union(BOOLEAN in_a, BOOLEAN in_b)
	{
	return in_a || in_b;
	}
/*...e*/

static void union_isectl(ISECTL *il_a, ISECTL *il_b, ISECTL *il)
	{
	combine_isectl(combine_union, il_a, il_b, il);
	}
/*...e*/
/*...sisect_isectl:0:*/
/*...scombine_isect:0:*/
static BOOLEAN combine_isect(BOOLEAN in_a, BOOLEAN in_b)
	{
	return in_a && in_b;
	}
/*...e*/

static void isect_isectl(ISECTL *il_a, ISECTL *il_b, ISECTL *il)
	{
	combine_isectl(combine_isect, il_a, il_b, il);
	}
/*...e*/
/*...sdiff_isectl:0:*/
/*...scombine_diff:0:*/
static BOOLEAN combine_diff(BOOLEAN in_a, BOOLEAN in_b)
	{
	return in_a && !in_b;
	}
/*...e*/

static void diff_isectl(ISECTL *il_a, ISECTL *il_b, ISECTL *il)
	{
	combine_isectl(combine_diff, il_a, il_b, il);
	}
/*...e*/
/*...ssdiff_isectl:0:*/
/*...scombine_sdiff:0:*/
static BOOLEAN combine_sdiff(BOOLEAN in_a, BOOLEAN in_b)
	{
	return in_a ^ in_b;
	}
/*...e*/

static void sdiff_isectl(ISECTL *il_a, ISECTL *il_b, ISECTL *il)
	{
	combine_isectl(combine_sdiff, il_a, il_b, il);
	}
/*...e*/
/*...sconcat_isectl:0:*/
/*
This is a quick case of unioning two intersection lists for use when it is
known that they do not overlap.
*/

static void concat_isectl(ISECTL *il_a, ISECTL *il_b, ISECTL *il)
	{
	ISECTL	*il_rest;
	int	i, j;

	if ( il_a->n_isects == 0 )
		{
		copy_isectl(il, il_b);
		return;
		}

	if ( il_b->n_isects == 0 )
		{
		copy_isectl(il, il_a);
		return;
		}

	if ( il_a->isects[0].t < il_b->isects[0].t )
		{
		copy_isectl(il, il_a);
		il_rest = il_b;
		}
	else
		{
		copy_isectl(il, il_b);
		il_rest = il_a;
		}

	for ( i = 0, j = il->n_isects; i < il_rest->n_isects; i++, j++ )
		il->isects[j] = il_rest->isects[i];

	il->n_isects = j;
	}
/*...e*/

void	intersect_shape(SHAPE *shape, VECTOR p, VECTOR q, ISECTL *ils[])
	{
	if ( shape->stype <= STYPE_QUAD )
		/* Is a leaf node */
		{
		SIL	sil;
		int	i;
	
		switch ( shape->stype )
			{
/*...sSTYPE_PLANE:24:*/
case STYPE_PLANE:
	intersect_plane(shape->u.plane, p, q, &sil);
	break;
/*...e*/
/*...sSTYPE_BIPLANE:24:*/
case STYPE_BIPLANE:
	intersect_biplane(shape->u.biplane, p, q, &sil);
	break;
/*...e*/
/*...sSTYPE_SPHERE:24:*/
case STYPE_SPHERE:
	intersect_sphere(shape->u.sphere, p, q, &sil);
	break;
/*...e*/
/*...sSTYPE_QUAD:24:*/
case STYPE_QUAD:
	intersect_quad(shape->u.quad, p, q, &sil);
	break;
/*...e*/
			}
		ils[0]->n_isects = sil.n_sis;
		for ( i = 0; i < sil.n_sis; i++ )
			{
			ils[0]->isects[i].t             = sil.sis[i].t;
			ils[0]->isects[i].entering      = sil.sis[i].entering;
			ils[0]->isects[i].shape         = shape;
			ils[0]->isects[i].negate_normal = FALSE;
			}
		}
	else
		/* Binary combination of two subtrees */
		{
		SHAPE	*shape_a = shape->u.shapes[0];
		SHAPE	*shape_b = shape->u.shapes[1];

		switch ( shape->stype )
			{
/*...sSTYPE_UNION:24:*/
case STYPE_UNION:
	if ( shape->overlap )
		{
		intersect_shape(shape_b, p, q, ils + 1);
		if ( is_empty_isectl(ils[1]) )
			intersect_shape(shape_a, p, q, ils);
		else if ( is_solid_isectl(ils[1]) )
			copy_isectl(ils[0], ils[1]);
		else
			{
			intersect_shape(shape_a, p, q, ils + 2);
			union_isectl(ils[2], ils[1], ils[0]);
			}
		}
	else
		/* No overlap, treat like concatentation */
		{
		intersect_shape(shape_a, p, q, ils + 1);
		intersect_shape(shape_b, p, q, ils + 2);
		concat_isectl(ils[1], ils[2], ils[0]);
		}
	break;
/*...e*/
/*...sSTYPE_ISECT:24:*/
case STYPE_ISECT:
	if ( shape->overlap )
		{
		intersect_shape(shape_b, p, q, ils + 1);
		if ( is_empty_isectl(ils[1]) )
			make_empty_isectl(ils[0]);
		else if ( is_solid_isectl(ils[1]) )
			intersect_shape(shape_a, p, q, ils);
		else
			{
			intersect_shape(shape_a, p, q, ils + 2);
			isect_isectl(ils[2], ils[1], ils[0]);
			}
		}
	else
		make_empty_isectl(ils[0]);
	break;
/*...e*/
/*...sSTYPE_DIFF:24:*/
case STYPE_DIFF:
	if ( shape->overlap )
		{
		intersect_shape(shape_b, p, q, ils + 1);
		if ( is_empty_isectl(ils[1]) )
			intersect_shape(shape_a, p, q, ils);
		else if ( is_solid_isectl(ils[1]) )
			make_empty_isectl(ils[0]);
		else
			{
			intersect_shape(shape_a, p, q, ils + 2);
			diff_isectl(ils[2], ils[1], ils[0]);
			}
		}
	else
		intersect_shape(shape_a, p, q, ils);
	break;
/*...e*/
/*...sSTYPE_SDIFF:24:*/
case STYPE_SDIFF:
	if ( shape->overlap )
		{
		intersect_shape(shape_b, p, q, ils + 1);
		if ( is_empty_isectl(ils[1]) )
			intersect_shape(shape_a, p, q, ils);
		else
			{
			intersect_shape(shape_a, p, q, ils + 2);
			sdiff_isectl(ils[2], ils[1], ils[0]);
			}
		}
	else
		/* No overlap, treat like concatentation */
		{
		intersect_shape(shape_a, p, q, ils + 1);
		intersect_shape(shape_b, p, q, ils + 2);
		concat_isectl(ils[1], ils[2], ils[0]);
		}
	break;
/*...e*/
/*...sSTYPE_EXTENT:24:*/
case STYPE_EXTENT:
	intersect_shape(shape_b, p, q, ils);
	if ( !is_empty_isectl(ils[0]) )
		intersect_shape(shape_a, p, q, ils);
	break;
/*...e*/
			}
		}
	}
/*...e*/
/*...snormal_to_shape      \45\ find normal at point on edge of shape:0:*/
VECTOR	normal_to_shape(SHAPE *shape, VECTOR p)
	{
	static VECTOR dummy_vector = { 0.0, 0.0, 0.0 };

	switch ( shape->stype )
		{
		case STYPE_PLANE:
			return normal_to_plane(shape->u.plane);
		case STYPE_BIPLANE:
			return normal_to_biplane(shape->u.biplane, p);
		case STYPE_SPHERE:
			return normal_to_sphere(shape->u.sphere, p);
		case STYPE_QUAD:
			return normal_to_quad(shape->u.quad, p);
		}
	return dummy_vector; /* Keep fussy C compiler happy */
	}
/*...e*/
