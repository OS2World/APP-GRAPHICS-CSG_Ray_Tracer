/*

sil.h - Interface to simple intersection lists

*/

typedef struct
	{
	double	t;			/* t value at point of intersection  */
	BOOLEAN	entering;		/* Are we entering the shape?        */
	} SI;

#define	N_SIS		(2+2)		/* Highest order of t is quadratic   */
					/* +2 for possible +/- INFINITE      */

typedef struct
	{
	int	n_sis;			/* How many intersections in list    */
	SI	sis[N_SIS];		/* Simple intersections              */
	} SIL;

#ifndef _SIL_

extern void intersect_linear_t(
	double coeff_of_t, double constant_term,
	VECTOR p, VECTOR q,
	void *shapeinfo, double (*value_of_shape)(void *shapeinfo, VECTOR v),
	SIL *sil
	);

extern void intersect_quadratic_t(
	double qa, double qb, double qc,
	VECTOR p, VECTOR q,
	void *shapeinfo, double (*value_of_shape)(void *shapeinfo, VECTOR v),
	SIL *sil
	);

#endif
