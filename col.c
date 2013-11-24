/*

col.c - Colour (within colour fields) code

*/

/*...simprovements:0:*/
/*

	INPUT	PRIMITIVE				PASSED DOWN

	_,_,_	col(a_rgb)
	x,y,_	col_field2d(bx,by,"fn.bmp")
	x,y,z	col_field3d(bx,by,bz,"fn.tex")

	x,y,z	col_interp0(a_col)			x(+1),y,z
	x,y,z	col_interp1(a_col)			x,y(+1),z
	x,y,z	col_interp2(a_col)			x,y,z(+1)

	x,y,z	col_remap(a_xyz_base,
			  a_xyz_v0,
			  a_xyz_v1,a_col)		x',y',z'
	x,y,z	col_cyl(lond,rd,hd,a_col)		lon/lond,r/rd,z/hd
	x,y,z	col_sph(lond,latd,rd,a_col)		lon/lond,lat/latd,r/rd
	x,y,z	col_nomove(a_col)			x,y,z
	a,b,c	col_mat2d(a00,a01,a10,a11)		a',b',c'
	a,b,c	col_mat3d(a00,a01,a02,a10,...,a22)	a',b',c'

+ origin to bitmaps
+ 3d texture files
+ optional interpolation
+ polar stuff isolated
+ more polar possibilities
+ no move with object
- more complicated
- huge matrices
- more primitives

*/
/*...e*/

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

/*...vrt\46\h:0:*/
/*...vfio\46\h:0:*/
/*...vtex\46\h:0:*/
/*...vvector\46\h:0:*/
/*...vrgbvec\46\h:0:*/
/*...e*/

typedef byte CTYPE;
#define	CTYPE_CONST	((CTYPE) 0)
#define	CTYPE_NO_MOVE	((CTYPE) 1)
#define	CTYPE_INTERP0	((CTYPE) 2)
#define	CTYPE_INTERP1	((CTYPE) 3)
#define	CTYPE_INTERP2	((CTYPE) 4)
#define	CTYPE_FIELD2D	((CTYPE) 5)
#define	CTYPE_FIELD3D	((CTYPE) 6)
#define	CTYPE_REMAP	((CTYPE) 7)
#define	CTYPE_CYLPOLAR	((CTYPE) 8)
#define	CTYPE_SPHPOLAR	((CTYPE) 9)
#define	CTYPE_MATRIX2D	((CTYPE) 10)
#define	CTYPE_MATRIX3D	((CTYPE) 11)

typedef struct _COL COL;

typedef struct { int bx,by; BITMAP *bitmap; } FIELD2D;
typedef struct { int bx,by,bz; TEX *tex; } FIELD3D;
typedef struct { VECTOR base,v0,v1,v2; double inv_m[3][3]; COL *col; } REMAP;
typedef struct { double lond,rd,hd; COL *col; } CYLPOLAR;
typedef struct { double lond,latd,rd; COL *col; } SPHPOLAR;
typedef struct { double m[2][2]; COL *col; } MATRIX2D;
typedef struct { double m[3][3]; COL *col; } MATRIX3D;

struct _COL
	{
	CTYPE ctype;
	union
		{
		RGBVEC rgbvec;		/* _CONST */
		COL *col;		/* _INTERP?D or _NO_MOVE */
		FIELD2D field2d;	/* _FIELD2D */
		FIELD3D field3d;	/* _FIELD3D */
		REMAP remap;		/* _REMAP */
		CYLPOLAR cylpolar;	/* _CYLPOLAR */
		SPHPOLAR sphpolar;	/* _SPHPOLAR */
		MATRIX2D matrix2d;	/* _MATRIX2D */
		MATRIX3D matrix3d;	/* _MATRIX3D */
		} u;
	};

/*...sinvert_matrix:0:*/
/*
	 -1     1          T
	M   = ------ cof(M)
	      det(M)
*/

/*...scofactor:0:*/
/*
See Stephenson p277 for defn of a 'minor'.
See Stephenson p307 for defn of a 'cofactor'.
*/

static double cofactor(int row, int col, double m[3][3])
	{
	static int lower[] = { 1, 0, 0 };
	static int upper[] = { 2, 2, 1 };
	int	lower_row = lower[row], upper_row = upper[row];
	int	lower_col = lower[col], upper_col = upper[col];
	double	minor = m[lower_row][lower_col] * m[upper_row][upper_col] -
			m[lower_row][upper_col] * m[upper_row][lower_col];

	return ((row + col) & 1) ? -minor : minor;
	}
/*...e*/

static BOOLEAN invert_matrix(double m[3][3], double inv_m[3][3])
	{
	int	row, col;
	double	det = 0.0;
	double	cof_m[3][3];

	for ( col = 0; col < 3; col++ )
		{
		for ( row = 0; row < 3; row++ )
			cof_m[row][col] = cofactor(row, col, m);
		det += m[0][col] * cof_m[0][col];
		}

	if ( det == 0.0 )
		return FALSE;

	det = 1.0 / det;

	for ( col = 0; col < 3; col++ )
		for ( row = 0; row < 3; row++ )
			inv_m[col][row] = det * cof_m[row][col];

	return TRUE;
	}
/*...e*/
/*...scompile_inv_m:0:*/
static void compile_inv_m(REMAP *remap)
	{
	double m[3][3];
	VECTOR v0, v1, v2;

	v0 = remap->v0;
	v1 = remap->v1;
	v2 = remap->v2;

	m[0][0] = v0.x; m[0][1] = v1.x; m[0][2] = v2.x;
	m[1][0] = v0.y; m[1][1] = v1.y; m[1][2] = v2.y;
	m[2][0] = v0.z; m[2][1] = v1.z; m[2][2] = v2.z;

	invert_matrix(m, remap->inv_m);
	}
/*...e*/
/*...sdetermine_coeffs:0:*/
/*

Point (p) = base (b) + multiples (c1,c2 and c3) of 3 vectors (v0, v1 and v2).

        /                \ /    \     /    \       /                \
        | v0.x v1.x v2.x | | c1 |     | c1 |       | v0.x v1.x v2.x |
p = b + | v0.y v1.y v2.y | | c2 | =>  | c2 | = inv | v0.y v1.y v2.y | ( p - b )
~   ~   | v0.z v1.z v2.z | | c3 |     | c3 |       | v0.z v1.z v2.z |   ~   ~
        \                / \    /     \    /       \                /

The inverse matrix is kept ready computed in the remap structure.

*/

static void determine_coeffs(
	REMAP *remap,
	VECTOR p,
	double *c0, double *c1, double *c2
	)
	{
	VECTOR pb;

	pb = vector_difference(p, remap->base);
	*c0 = remap->inv_m[0][0] * pb.x +
	      remap->inv_m[0][1] * pb.y +
	      remap->inv_m[0][2] * pb.z;
	*c1 = remap->inv_m[1][0] * pb.x +
	      remap->inv_m[1][1] * pb.y +
	      remap->inv_m[1][2] * pb.z;
	*c2 = remap->inv_m[2][0] * pb.x +
	      remap->inv_m[2][1] * pb.y +
	      remap->inv_m[2][2] * pb.z;
	}
/*...e*/

/*...screate_const_col:0:*/
COL *create_const_col(RGBVEC rgbvec)
	{
	COL *col;

	if ( (col = malloc(sizeof(COL))) == NULL )
		return NULL;
	col->ctype    = CTYPE_CONST;
	col->u.rgbvec = rgbvec;
	return col;
	}
/*...e*/
/*...screate_chain_on_col:0:*/
static COL *create_chain_on_col(CTYPE ctype, COL *col_chain_on)
	{
	COL *col;

	if ( (col = malloc(sizeof(COL))) == NULL )
		return NULL;
	col->ctype = ctype;
	col->u.col = col_chain_on;
	return col;
	}
/*...e*/
/*...screate_no_move_col:0:*/
COL *create_no_move_col(COL *col)
	{
	return create_chain_on_col(CTYPE_NO_MOVE, col);
	}
/*...e*/
/*...screate_interp0_col:0:*/
COL *create_interp0_col(COL *col)
	{
	return create_chain_on_col(CTYPE_INTERP0, col);
	}
/*...e*/
/*...screate_interp1_col:0:*/
COL *create_interp1_col(COL *col)
	{
	return create_chain_on_col(CTYPE_INTERP1, col);
	}
/*...e*/
/*...screate_interp2_col:0:*/
COL *create_interp2_col(COL *col)
	{
	return create_chain_on_col(CTYPE_INTERP2, col);
	}
/*...e*/
/*...screate_2d_field_col:0:*/
COL *create_2d_field_col(double bx, double by, BITMAP *bitmap)
	{
	COL *col;

	if ( (col = malloc(sizeof(COL))) == NULL )
		return NULL;
	col->ctype            = CTYPE_FIELD2D;
	col->u.field2d.bx     = (int) bx;
	col->u.field2d.by     = (int) by;
	col->u.field2d.bitmap = bitmap;
	return col;
	}
/*...e*/
/*...screate_3d_field_col:0:*/
COL *create_3d_field_col(double bx, double by, double bz, TEX *tex)
	{
	COL *col;

	if ( (col = malloc(sizeof(COL))) == NULL )
		return NULL;
	col->ctype         = CTYPE_FIELD3D;
	col->u.field3d.bx  = (int) bx;
	col->u.field3d.by  = (int) by;
	col->u.field3d.bz  = (int) bz;
	col->u.field3d.tex = tex;
	return col;
	}
/*...e*/
/*...screate_remap_col:0:*/
COL *create_remap_col(
	VECTOR base,
	VECTOR v0, VECTOR v1, VECTOR v2,
	COL *col_chain_on
	)
	{
	COL *col;

	if ( (col = malloc(sizeof(COL))) == NULL )
		return NULL;
	col->ctype        = CTYPE_REMAP;
	col->u.remap.base = base;
	col->u.remap.v0   = v0;
	col->u.remap.v1   = v1;
	col->u.remap.v2   = v2;
	col->u.remap.col  = col_chain_on;
	compile_inv_m(&(col->u.remap));
	return col;
	}
/*...e*/
/*...screate_cyl_polar_col:0:*/
COL *create_cyl_polar_col(double lond, double rd, double hd, COL *col_chain_on)
	{
	COL *col;

	if ( (col = malloc(sizeof(COL))) == NULL )
		return NULL;
	col->ctype           = CTYPE_CYLPOLAR;
	col->u.cylpolar.lond = lond;
	col->u.cylpolar.rd   = rd;
	col->u.cylpolar.hd   = hd;
	col->u.cylpolar.col  = col_chain_on;
	return col;
	}
/*...e*/
/*...screate_sph_polar_col:0:*/
COL *create_sph_polar_col(double lond, double latd, double rd, COL *col_chain_on)
	{
	COL *col;

	if ( (col = malloc(sizeof(COL))) == NULL )
		return NULL;
	col->ctype           = CTYPE_SPHPOLAR;
	col->u.sphpolar.lond = lond;
	col->u.sphpolar.latd = latd;
	col->u.sphpolar.rd   = rd;
	col->u.sphpolar.col  = col_chain_on;
	return col;
	}
/*...e*/
/*...screate_2d_matrix_col:0:*/
COL *create_2d_matrix_col(double m[2][2], COL *col_chain_on)
	{
	COL *col;

	if ( (col = malloc(sizeof(COL))) == NULL )
		return NULL;
	col->ctype          = CTYPE_MATRIX2D;
	memcpy(col->u.matrix2d.m, m, sizeof(col->u.matrix2d.m));
	col->u.matrix2d.col = col_chain_on;
	return col;
	}
/*...e*/
/*...screate_3d_matrix_col:0:*/
COL *create_3d_matrix_col(double m[3][3], COL *col_chain_on)
	{
	COL *col;

	if ( (col = malloc(sizeof(COL))) == NULL )
		return NULL;
	col->ctype          = CTYPE_MATRIX3D;
	memcpy(col->u.matrix3d.m, m, sizeof(col->u.matrix3d.m));
	col->u.matrix3d.col = col_chain_on;
	return col;
	}
/*...e*/

/*...scopy_col:0:*/
COL *copy_col(COL *col)
	{
	COL *copy;

	if ( (copy = malloc(sizeof(COL))) == NULL )
		return NULL;

	copy->ctype = col->ctype;
	switch ( col->ctype )
		{
		case CTYPE_CONST:
			copy->u.rgbvec = col->u.rgbvec;
			break;
		case CTYPE_NO_MOVE:
		case CTYPE_INTERP0:
		case CTYPE_INTERP1:
		case CTYPE_INTERP2:
			if ( (copy->u.col = copy_col(col->u.col)) == NULL )
				{
				free(col);
				return NULL;
				}
			break;
		case CTYPE_FIELD2D:
			copy->u.field2d.bx = col->u.field2d.bx;
			copy->u.field2d.by = col->u.field2d.by;
			if ( (copy->u.field2d.bitmap = fio_copy_bitmap(col->u.field2d.bitmap)) == NULL )
				{
				free(col);
				return NULL;
				}
			break;
		case CTYPE_FIELD3D:
			copy->u.field3d.bx = col->u.field3d.bx;
			copy->u.field3d.by = col->u.field3d.by;
			copy->u.field3d.bz = col->u.field3d.bz;
			if ( (copy->u.field3d.tex = copy_tex(col->u.field3d.tex)) == NULL )
				{
				free(col);
				return NULL;
				}
			break;
		case CTYPE_REMAP:
			copy->u.remap = col->u.remap;
			if ( (copy->u.remap.col = copy_col(col->u.remap.col)) == NULL )
				{
				free(col);
				return NULL;
				}
			break;
		case CTYPE_SPHPOLAR:
			copy->u.sphpolar = col->u.sphpolar;
			if ( (copy->u.sphpolar.col = copy_col(col->u.sphpolar.col)) == NULL )
				{
				free(col);
				return NULL;
				}
			break;
		case CTYPE_CYLPOLAR:
			copy->u.cylpolar = col->u.cylpolar;
			if ( (copy->u.cylpolar.col = copy_col(col->u.cylpolar.col)) == NULL )
				{
				free(col);
				return NULL;
				}
			break;
		case CTYPE_MATRIX2D:
			copy->u.matrix2d = col->u.matrix2d;
			if ( (copy->u.matrix2d.col = copy_col(col->u.matrix2d.col)) == NULL )
				{
				free(col);
				return NULL;
				}
			break;
		case CTYPE_MATRIX3D:
			copy->u.matrix3d = col->u.matrix3d;
			if ( (copy->u.matrix3d.col = copy_col(col->u.matrix3d.col)) == NULL )
				{
				free(col);
				return NULL;
				}
			break;
		}

	return copy;
	}
/*...e*/
/*...sdestroy_col:0:*/
void destroy_col(COL *col)
	{
	switch ( col->ctype )
		{
		case CTYPE_NO_MOVE:
		case CTYPE_INTERP0:
		case CTYPE_INTERP1:
		case CTYPE_INTERP2:
			destroy_col(col->u.col);
			break;
		case CTYPE_FIELD2D:
			fio_destroy_bitmap(col->u.field2d.bitmap);
			break;
		case CTYPE_FIELD3D:
			destroy_tex(col->u.field3d.tex);
			break;
		case CTYPE_REMAP:
			destroy_col(col->u.remap.col);
			break;
		case CTYPE_SPHPOLAR:
			destroy_col(col->u.sphpolar.col);
			break;
		case CTYPE_CYLPOLAR:
			destroy_col(col->u.cylpolar.col);
			break;
		case CTYPE_MATRIX2D:
			destroy_col(col->u.matrix2d.col);
			break;
		case CTYPE_MATRIX3D:
			destroy_col(col->u.matrix3d.col);
			break;
		}
	free(col);
	}
/*...e*/

/*...strans_x_col:0:*/
void trans_x_col(COL *col, double t)
	{
	switch ( col->ctype )
		{
		case CTYPE_INTERP0:
		case CTYPE_INTERP1:
		case CTYPE_INTERP2:
			trans_x_col(col->u.col, t);
			break;
		case CTYPE_REMAP:
			col->u.remap.base.x += t;
			break;
		case CTYPE_SPHPOLAR:
			trans_x_col(col->u.sphpolar.col, t);
			break;
		case CTYPE_CYLPOLAR:
			trans_x_col(col->u.cylpolar.col, t);
			break;
		case CTYPE_MATRIX2D:
			trans_x_col(col->u.matrix2d.col, t);
			break;
		case CTYPE_MATRIX3D:
			trans_x_col(col->u.matrix3d.col, t);
			break;
		}
	}
/*...e*/
/*...strans_y_col:0:*/
void trans_y_col(COL *col, double t)
	{
	switch ( col->ctype )
		{
		case CTYPE_INTERP0:
		case CTYPE_INTERP1:
		case CTYPE_INTERP2:
			trans_y_col(col->u.col, t);
			break;
		case CTYPE_REMAP:
			col->u.remap.base.y += t;
			break;
		case CTYPE_SPHPOLAR:
			trans_y_col(col->u.sphpolar.col, t);
			break;
		case CTYPE_CYLPOLAR:
			trans_y_col(col->u.cylpolar.col, t);
			break;
		case CTYPE_MATRIX2D:
			trans_y_col(col->u.matrix2d.col, t);
			break;
		case CTYPE_MATRIX3D:
			trans_y_col(col->u.matrix3d.col, t);
			break;
		}
	}
/*...e*/
/*...strans_z_col:0:*/
void trans_z_col(COL *col, double t)
	{
	switch ( col->ctype )
		{
		case CTYPE_INTERP0:
		case CTYPE_INTERP1:
		case CTYPE_INTERP2:
			trans_z_col(col->u.col, t);
			break;
		case CTYPE_REMAP:
			col->u.remap.base.z += t;
			break;
		case CTYPE_SPHPOLAR:
			trans_z_col(col->u.sphpolar.col, t);
			break;
		case CTYPE_CYLPOLAR:
			trans_z_col(col->u.cylpolar.col, t);
			break;
		case CTYPE_MATRIX2D:
			trans_z_col(col->u.matrix2d.col, t);
			break;
		case CTYPE_MATRIX3D:
			trans_z_col(col->u.matrix3d.col, t);
			break;
		}
	}
/*...e*/
/*...sscale_x_col:0:*/
/*
How to scale a colour is quite tricky to define.
If a bitmap is mapped onto a surface of a shape, then if the shape expands,
then we would expect the bitmap to 'expand' to cover the new shape.
To do this we enlarge the basis vectors for a CTYPE_REMAP.
Remapping between coordinate systems and general matrix manipulation of
parameters p0,p1 and p2 remain untouched for now.
*/

void scale_x_col(COL *col, double factor)
	{
	switch ( col->ctype )
		{
		case CTYPE_INTERP0:
		case CTYPE_INTERP1:
		case CTYPE_INTERP2:
			scale_x_col(col->u.col, factor);
			break;
		case CTYPE_REMAP:
			{
			REMAP	*remap = &(col->u.remap);

			remap->base.x *= factor;
			remap->v0.x   *= factor;
			remap->v1.x   *= factor;
			remap->v2.x   *= factor;
			compile_inv_m(remap);
			}
			break;
		case CTYPE_SPHPOLAR:
			scale_x_col(col->u.sphpolar.col, factor);
			break;
		case CTYPE_CYLPOLAR:
			scale_x_col(col->u.cylpolar.col, factor);
			break;
		case CTYPE_MATRIX2D:
			scale_x_col(col->u.matrix2d.col, factor);
			break;
		case CTYPE_MATRIX3D:
			scale_x_col(col->u.matrix3d.col, factor);
			break;
		}
	}
/*...e*/
/*...sscale_y_col:0:*/
void scale_y_col(COL *col, double factor)
	{
	switch ( col->ctype )
		{
		case CTYPE_INTERP0:
		case CTYPE_INTERP1:
		case CTYPE_INTERP2:
			scale_y_col(col->u.col, factor);
			break;
		case CTYPE_REMAP:
			{
			REMAP	*remap = &(col->u.remap);

			remap->base.y *= factor;
			remap->v0.y   *= factor;
			remap->v1.y   *= factor;
			remap->v2.y   *= factor;
			compile_inv_m(remap);
			}
			break;
		case CTYPE_SPHPOLAR:
			scale_y_col(col->u.sphpolar.col, factor);
			break;
		case CTYPE_CYLPOLAR:
			scale_y_col(col->u.cylpolar.col, factor);
			break;
		case CTYPE_MATRIX2D:
			scale_y_col(col->u.matrix2d.col, factor);
			break;
		case CTYPE_MATRIX3D:
			scale_y_col(col->u.matrix3d.col, factor);
			break;
		}
	}
/*...e*/
/*...sscale_z_col:0:*/
void scale_z_col(COL *col, double factor)
	{
	switch ( col->ctype )
		{
		case CTYPE_INTERP0:
		case CTYPE_INTERP1:
		case CTYPE_INTERP2:
			scale_z_col(col->u.col, factor);
			break;
		case CTYPE_REMAP:
			{
			REMAP	*remap = &(col->u.remap);

			remap->base.z *= factor;
			remap->v0.z   *= factor;
			remap->v1.z   *= factor;
			remap->v2.z   *= factor;
			compile_inv_m(remap);
			}
			break;
		case CTYPE_SPHPOLAR:
			scale_z_col(col->u.sphpolar.col, factor);
			break;
		case CTYPE_CYLPOLAR:
			scale_z_col(col->u.cylpolar.col, factor);
			break;
		case CTYPE_MATRIX2D:
			scale_z_col(col->u.matrix2d.col, factor);
			break;
		case CTYPE_MATRIX3D:
			scale_z_col(col->u.matrix3d.col, factor);
			break;
		}
	}
/*...e*/
/*...srot_x_col:0:*/
void rot_x_col(COL *col, double angle)
	{
	switch ( col->ctype )
		{
		case CTYPE_INTERP0:
		case CTYPE_INTERP1:
		case CTYPE_INTERP2:
			rot_x_col(col->u.col, angle);
			break;
		case CTYPE_REMAP:
			{
			REMAP	*remap = &(col->u.remap);

			remap->base = rot_x_vector(remap->base, angle);
			remap->v0   = rot_x_vector(remap->v0  , angle);
			remap->v1   = rot_x_vector(remap->v1  , angle);
			remap->v2   = rot_x_vector(remap->v2  , angle);
			compile_inv_m(remap);
			}
			break;
		case CTYPE_SPHPOLAR:
			rot_x_col(col->u.sphpolar.col, angle);
			break;
		case CTYPE_CYLPOLAR:
			rot_x_col(col->u.cylpolar.col, angle);
			break;
		case CTYPE_MATRIX2D:
			rot_x_col(col->u.matrix2d.col, angle);
			break;
		case CTYPE_MATRIX3D:
			rot_x_col(col->u.matrix3d.col, angle);
			break;
		}
	}
/*...e*/
/*...srot_y_col:0:*/
void rot_y_col(COL *col, double angle)
	{
	switch ( col->ctype )
		{
		case CTYPE_INTERP0:
		case CTYPE_INTERP1:
		case CTYPE_INTERP2:
			rot_y_col(col->u.col, angle);
			break;
		case CTYPE_REMAP:
			{
			REMAP	*remap = &(col->u.remap);

			remap->base = rot_y_vector(remap->base, angle);
			remap->v0   = rot_y_vector(remap->v0  , angle);
			remap->v1   = rot_y_vector(remap->v1  , angle);
			remap->v2   = rot_y_vector(remap->v2  , angle);
			compile_inv_m(remap);
			}
			break;
		case CTYPE_SPHPOLAR:
			rot_y_col(col->u.sphpolar.col, angle);
			break;
		case CTYPE_CYLPOLAR:
			rot_y_col(col->u.cylpolar.col, angle);
			break;
		case CTYPE_MATRIX2D:
			rot_y_col(col->u.matrix2d.col, angle);
			break;
		case CTYPE_MATRIX3D:
			rot_y_col(col->u.matrix3d.col, angle);
			break;
		}
	}
/*...e*/
/*...srot_z_col:0:*/
void rot_z_col(COL *col, double angle)
	{
	switch ( col->ctype )
		{
		case CTYPE_INTERP0:
		case CTYPE_INTERP1:
		case CTYPE_INTERP2:
			rot_z_col(col->u.col, angle);
			break;
		case CTYPE_REMAP:
			{
			REMAP	*remap = &(col->u.remap);

			remap->base = rot_z_vector(remap->base, angle);
			remap->v0   = rot_z_vector(remap->v0  , angle);
			remap->v1   = rot_z_vector(remap->v1  , angle);
			remap->v2   = rot_z_vector(remap->v2  , angle);
			compile_inv_m(remap);
			}
			break;
		case CTYPE_SPHPOLAR:
			rot_z_col(col->u.sphpolar.col, angle);
			break;
		case CTYPE_CYLPOLAR:
			rot_z_col(col->u.cylpolar.col, angle);
			break;
		case CTYPE_MATRIX2D:
			rot_z_col(col->u.matrix2d.col, angle);
			break;
		case CTYPE_MATRIX3D:
			rot_z_col(col->u.matrix3d.col, angle);
			break;
		}
	}
/*...e*/

/*...sevaluate_col:0:*/
/*...spos_fmod:0:*/
static double pos_fmod(double num, double denom)
	{
	if ( num >= 0.0 )
		return fmod(num, denom);
	else
		return denom + fmod(num, denom);
	}
/*...e*/

RGBVEC evaluate_col(COL *col, double p0, double p1, double p2)
	{
	static RGBVEC dummy_rgbvec = { 0.0, 0.0, 0.0 };

	switch ( col->ctype )
		{
/*...sCTYPE_CONST:16:*/
case CTYPE_CONST:
	return col->u.rgbvec;
/*...e*/
/*...sCTYPE_NO_MOVE:16:*/
case CTYPE_NO_MOVE:
	return evaluate_col(col->u.col, p0, p1, p2);
/*...e*/
/*...sCTYPE_INTERP0:16:*/
case CTYPE_INTERP0:
	{
	double fp0 = floor(p0);
	RGBVEC r[2];
	r[0] = evaluate_col(col->u.col, fp0      , p1, p2);
	r[1] = evaluate_col(col->u.col, fp0 + 1.0, p1, p2);
	return rgbvec_interp_1d(r, pos_fmod(p0, 1.0));
	}
/*...e*/
/*...sCTYPE_INTERP1:16:*/
case CTYPE_INTERP1:
	{
	double fp1 = floor(p1);
	RGBVEC r[2];
	r[0] = evaluate_col(col->u.col, p0, fp1      , p2);
	r[1] = evaluate_col(col->u.col, p0, fp1 + 1.0, p2);
	return rgbvec_interp_1d(r, pos_fmod(p1, 1.0));
	}
/*...e*/
/*...sCTYPE_INTERP2:16:*/
case CTYPE_INTERP2:
	{
	double fp2 = floor(p2);
	RGBVEC r[2];
	r[0] = evaluate_col(col->u.col, p0, p1, fp2      );
	r[1] = evaluate_col(col->u.col, p0, p1, fp2 + 1.0);
	return rgbvec_interp_1d(r, pos_fmod(p2, 1.0));
	}
/*...e*/
/*...sCTYPE_FIELD2D:16:*/
case CTYPE_FIELD2D:
	{
	BITMAP *bitmap = col->u.field2d.bitmap;
	int w = fio_width(bitmap);
	int h = fio_height(bitmap);
	int x = (int) pos_fmod(p0 + col->u.field2d.bx, (double) w);
	int y = (int) pos_fmod(p1 + col->u.field2d.by, (double) h);
	byte r, g, b;
	RGBVEC rgbvec;

	fio_get_pixel(bitmap, x, y, &r, &g, &b);
	rgbvec.b = ((double) b) / 255.0;
	rgbvec.g = ((double) g) / 255.0;
	rgbvec.r = ((double) r) / 255.0;
	return rgbvec;
	}
/*...e*/
/*...sCTYPE_FIELD3D:16:*/
case CTYPE_FIELD3D:
	{
	TEX *tex = col->u.field3d.tex;
	int w = width_tex(tex);
	int h = height_tex(tex);
	int d = depth_tex(tex);
	int x = (int) pos_fmod(p0 + col->u.field3d.bx, (double) w);
	int y = (int) pos_fmod(p1 + col->u.field3d.by, (double) h);
	int z = (int) pos_fmod(p2 + col->u.field3d.bz, (double) d);
	byte r, g, b;
	RGBVEC rgbvec;

	get_voxel_tex(tex, x, y, z, &r, &g, &b);
	rgbvec.b = ((double) b) / 255.0;
	rgbvec.g = ((double) g) / 255.0;
	rgbvec.r = ((double) r) / 255.0;
	return rgbvec;
	}
/*...e*/
/*...sCTYPE_REMAP:16:*/
case CTYPE_REMAP:
	{
	VECTOR p;
	double c0, c1, c2;

	p.x = p0;
	p.y = p1;
	p.z = p2;
	determine_coeffs(&col->u.remap, p, &c0, &c1, &c2);

	return evaluate_col(col->u.remap.col, c0, c1, c2);
	}
/*...e*/
/*...sCTYPE_SPHPOLAR:16:*/
case CTYPE_SPHPOLAR:
	{
	double lon = atan2(p1, p0);
	double p01 = p0 * p0 + p1 * p1;
	double lat = atan2(p2, sqrt(p01));
	double p012 = p01 + p2 * p2;

	return evaluate_col(col->u.sphpolar.col,
			    lon        / col->u.sphpolar.lond,
			    lat        / col->u.sphpolar.latd,
			    sqrt(p012) / col->u.sphpolar.rd);
	}
/*...e*/
/*...sCTYPE_CYLPOLAR:16:*/
case CTYPE_CYLPOLAR:
	{
	double lon = atan2(p1, p0);
	double p01 = p0 * p0 + p1 * p1;

	return evaluate_col(col->u.cylpolar.col,
			    lon       / col->u.cylpolar.lond,
			    sqrt(p01) / col->u.cylpolar.rd  ,
			    p2        / col->u.cylpolar.hd  );
	}
/*...e*/
/*...sCTYPE_MATRIX2D:16:*/
case CTYPE_MATRIX2D:
	{
	double (*m)[2] = col->u.matrix2d.m;
	double q0 = m[0][0] * p0 + m[0][1] * p1;
	double q1 = m[1][0] * p0 + m[1][1] * p1;

	return evaluate_col(col->u.matrix2d.col, q0, q1, p2);
	}
/*...e*/
/*...sCTYPE_MATRIX3D:16:*/
case CTYPE_MATRIX3D:
	{
	double (*m)[3] = col->u.matrix3d.m;
	double q0 = m[0][0] * p0 + m[0][1] * p1 + m[0][2] * p2;
	double q1 = m[1][0] * p0 + m[1][1] * p1 + m[1][2] * p2;
	double q2 = m[2][0] * p0 + m[2][1] * p1 + m[2][2] * p2;

	return evaluate_col(col->u.matrix3d.col, q0, q1, q2);
	}
/*...e*/
		}
	return dummy_rgbvec; /* Please fussy C compiler */
	}
/*...e*/
