/*

shape.h - General shape module

*/

typedef byte STYPE;			/* Type (discriminant of union)      */
#define	STYPE_PLANE	((STYPE) 0)	/* Shape is a half-plane             */
#define	STYPE_BIPLANE	((STYPE) 1)	/* Shape is a bi-plane               */
#define	STYPE_SPHERE	((STYPE) 2)	/* Shape is a sphere                 */
#define	STYPE_QUAD	((STYPE) 3)	/* Shape is a quadratic solid        */
#define	STYPE_UNION	((STYPE) 4)	/* Shape is the union of 2 shapes    */
#define	STYPE_ISECT	((STYPE) 5)	/* Shape is the intersection of 2    */
#define	STYPE_DIFF	((STYPE) 6)	/* Shape is the difference of 2      */
#define	STYPE_SDIFF	((STYPE) 7)	/* Shape is the symmetric difference */
#define	STYPE_EXTENT	((STYPE) 8)	/* Shape is 1st if intersects 2nd    */

#define	SHAPE_F		struct shape_struct

typedef struct shape_struct
	{
	STYPE	stype;			/* Shape type                        */
	union
		{
		PLANE	*plane;		/* Used if STYPE_PLANE               */
		BIPLANE	*biplane;	/* Used if STYPE_BIPLANE             */
		SPHERE	*sphere;	/* Used if STYPE_SPHERE              */
		QUAD	*quad;		/* Used if STYPE_QUAD                */
		SHAPE_F	*shapes[2];	/* For boolean combinations          */
		} u;
	BOOLEAN	overlap;		/* If boolean, could they overlap?   */
	SURF	*surf;			/* Surface characteristic of shape   */
	} SHAPE;

typedef struct
	{
	double	t;			/* t value at point of intersection  */
	SHAPE	*shape;			/* Shape intersection is with        */
	BOOLEAN	entering;		/* Are we entering the shape?        */
	BOOLEAN	negate_normal;		/* Is normal to be negated?          */
	} ISECT;

typedef struct
	{
	int	n_isects;		/* How many intersections in list    */
	ISECT	isects[1];		/* Intersections                     */
	} ISECTL;

#ifndef _SHAPE_

extern SHAPE *create_plane_shape(PLANE *plane, SURF *surf);
extern SHAPE *create_biplane_shape(BIPLANE *biplane, SURF *surf);
extern SHAPE *create_sphere_shape(SPHERE *sphere, SURF *surf);
extern SHAPE *create_quad_shape(QUAD *quad, SURF *surf);
extern SHAPE *create_bin_shape(STYPE stype, SHAPE *shape0, SHAPE *shape1);
extern SHAPE *copy_shape(SHAPE *shape);
extern void destroy_shape(SHAPE *shape);

extern void trans_x(SHAPE *shape, double t);
extern void trans_y(SHAPE *shape, double t);
extern void trans_z(SHAPE *shape, double t);
extern void trans(SHAPE *shape, VECTOR t);
extern BOOLEAN scale_x(SHAPE *shape, double factor);
extern BOOLEAN scale_y(SHAPE *shape, double factor);
extern BOOLEAN scale_z(SHAPE *shape, double factor);
extern BOOLEAN scale(SHAPE *shape, VECTOR factor);
extern void rot_x(SHAPE *shape, double angle);
extern void rot_y(SHAPE *shape, double angle);
extern void rot_z(SHAPE *shape, double angle);

extern BOOLEAN is_empty_isectl(ISECTL *il);
extern BOOLEAN is_solid_isectl(ISECTL *il);
extern void t_after_isectl(ISECTL *il, double t);

extern ISECTL *create_isectl(int n_isects);
extern void destroy_isectl(ISECTL *il);

extern void preprocess_shape(SHAPE *shape, int *n_isectls, int *n_isects);
extern void intersect_shape(SHAPE *shape, VECTOR p, VECTOR q, ISECTL *ils[]);
extern VECTOR normal_to_shape(SHAPE *shape, VECTOR p);

extern BOOLEAN resurf(SHAPE *shape, SURF *surf);

#endif
