/*

surf.c - Surface datatype

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
#define	_SURF_
#include "surf.h"

/*...vrt\46\h:0:*/
/*...vfio\46\h:0:*/
/*...vtex\46\h:0:*/
/*...vvector\46\h:0:*/
/*...vrgbvec\46\h:0:*/
/*...vcol\46\h:0:*/
/*...vsurf\46\h:0:*/
/*...e*/

/*...screate_surf:0:*/
SURF	*create_surf(
	double	ka, double kd, double ks, double kt,
	COL	*od, COL *os,
	double	phong,
	double	rinx
	)
	{
	SURF	*surf;

	if ( (surf = malloc(sizeof(SURF))) == NULL )
		return NULL;

	surf->ka    = ka;
	surf->kd    = kd;
	surf->ks    = ks;
	surf->kt    = kt;
	surf->od    = od;
	surf->os    = os;
	surf->phong = phong;
	surf->rinx  = rinx;
	return surf;
	}
/*...e*/
/*...scopy_surf:0:*/
SURF	*copy_surf(SURF *surf)
	{
	SURF	*copy;

	if ( (copy = malloc(sizeof(SURF))) == NULL )
		return NULL;

	if ( (copy->od = copy_col(surf->od)) == NULL )
		{
		free(copy);
		return NULL;
		}

	if ( (copy->os = copy_col(surf->os)) == NULL )
		{
		destroy_col(copy->od);
		free(copy);
		return NULL;
		}

	copy->ka    = surf->ka;
	copy->kd    = surf->kd;
	copy->ks    = surf->ks;
	copy->kt    = surf->kt;
	copy->phong = surf->phong;
	copy->rinx  = surf->rinx;

	return copy;
	}
/*...e*/
/*...sdestroy_surf:0:*/
void	destroy_surf(SURF *surf)
	{
	destroy_col(surf->od);
	destroy_col(surf->os);
	free(surf);
	}
/*...e*/

/*...strans_x_surf:0:*/
void	trans_x_surf(SURF *surf, double t)
	{
	trans_x_col(surf->od, t);
	trans_x_col(surf->os, t);
	}
/*...e*/
/*...strans_y_surf:0:*/
void	trans_y_surf(SURF *surf, double t)
	{
	trans_y_col(surf->od, t);
	trans_y_col(surf->os, t);
	}
/*...e*/
/*...strans_z_surf:0:*/
void	trans_z_surf(SURF *surf, double t)
	{
	trans_z_col(surf->od, t);
	trans_z_col(surf->os, t);
	}
/*...e*/
/*...sscale_x_surf:0:*/
void	scale_x_surf(SURF *surf, double factor)
	{
	scale_x_col(surf->od, factor);
	scale_x_col(surf->os, factor);
	}
/*...e*/
/*...sscale_y_surf:0:*/
void	scale_y_surf(SURF *surf, double factor)
	{
	scale_y_col(surf->od, factor);
	scale_y_col(surf->os, factor);
	}
/*...e*/
/*...sscale_z_surf:0:*/
void	scale_z_surf(SURF *surf, double factor)
	{
	scale_z_col(surf->od, factor);
	scale_z_col(surf->os, factor);
	}
/*...e*/
/*...srot_x_surf:0:*/
void	rot_x_surf(SURF *surf, double angle)
	{
	rot_x_col(surf->od, angle);
	rot_x_col(surf->os, angle);
	}
/*...e*/
/*...srot_y_surf:0:*/
void	rot_y_surf(SURF *surf, double angle)
	{
	rot_y_col(surf->od, angle);
	rot_y_col(surf->os, angle);
	}
/*...e*/
/*...srot_z_surf:0:*/
void	rot_z_surf(SURF *surf, double angle)
	{
	rot_z_col(surf->od, angle);
	rot_z_col(surf->os, angle);
	}
/*...e*/
