/*

sil.c - Simple intersection lists

*/

/*...sincludes:0:*/
#include <stdlib.h>
#include <stddef.h>
#include <memory.h>
#include <math.h>
#include "rt.h"
#include "vector.h"
#define	_SIL_
#include "sil.h"

/*...vrt\46\h:0:*/
/*...vvector\46\h:0:*/
/*...vsil\46\h:0:*/
/*...e*/

#define	ZERO(x) ( ((x)>=-1.0e-10) && ((x)<=1.0e-10) )
#define	SIDESTEP (1.0e6)

/*...sintersect_linear_t:0:*/
void	intersect_linear_t(
	double coeff_of_t, double constant_term,
	VECTOR p, VECTOR q,
	void *shapeinfo, double (*value_of_shape)(void *shapeinfo, VECTOR v),
	SIL *sil
	)
	{
	if ( ZERO(coeff_of_t) )
		/* No intersection with solid surface */
		{
		if ( (*value_of_shape)(shapeinfo, p) <= 0.0 )
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
		/* Intersects solid exactly once */
		{
		double t = - constant_term / coeff_of_t;

		sil->n_sis            = 2;
		sil->sis[0].entering = TRUE;
		sil->sis[1].entering = FALSE;

		if ( (*value_of_shape)(shapeinfo, t_along_pq(p, q, t+SIDESTEP)) <= 0.0 )
			/* Entering solid */
			{
			sil->sis[0].t = t;
			sil->sis[1].t = INFINITE;
			}
		else
			/* Leaving solid */
			{
			sil->sis[0].t = -INFINITE;
			sil->sis[1].t = t;
			}
		}
	}
/*...e*/
/*...sintersect_quadratic_t:0:*/
void	intersect_quadratic_t(
	double qa, double qb, double qc,
	VECTOR p, VECTOR q,
	void *shapeinfo, double (*value_of_shape)(void *shapeinfo, VECTOR v),
	SIL *sil
	)
	{
	double	qs;

	if ( ZERO(qa) )
/*...sdo linear case:16:*/
intersect_linear_t(qb, qc, p, q, shapeinfo, value_of_shape, sil);
/*...e*/
	else if ( (qs = qb*qb - 4.0*qa*qc) < 0.0 )
/*...sno real roots:16:*/
/* No real roots => no intersections, all solid or all empty */
{
if ( (*value_of_shape)(shapeinfo, p) <= 0.0 )
	{
	sil->n_sis            = 2;
	sil->sis[0].t        = -INFINITE;
	sil->sis[0].entering = TRUE;
	sil->sis[1].t        = INFINITE;
	sil->sis[1].entering = FALSE;
	}
else
	sil->n_sis = 0;
}
/*...e*/
	else if ( ZERO(qs) )
/*...sone double root:16:*/
{
double	t = -qb / (qa + qa);
BOOLEAN from_inside = ( (*value_of_shape)(shapeinfo, t_along_pq(p, q, t-SIDESTEP)) <= 0.0 );
BOOLEAN to_inside   = ( (*value_of_shape)(shapeinfo, t_along_pq(p, q, t+SIDESTEP)) <= 0.0 );

#define	DUAL_SWITCH(f,t)	(((f)<<1)+(t))

switch ( DUAL_SWITCH(from_inside, to_inside) )
	{
	case DUAL_SWITCH(FALSE, FALSE):		/* Grazing surface */
		sil->n_sis            = 0;
		break;
	case DUAL_SWITCH(FALSE,  TRUE):		/* Entering solid area */
		sil->n_sis            = 2;
		sil->sis[0].t        = t;
		sil->sis[0].entering = TRUE;
		sil->sis[1].t        = INFINITE;
		sil->sis[1].entering = FALSE;
		break;
	case DUAL_SWITCH( TRUE, FALSE):		/* Leaving solid area */
		sil->n_sis            = 2;
		sil->sis[0].t        = -INFINITE;
		sil->sis[0].entering = TRUE;
		sil->sis[1].t        = t;
		sil->sis[1].entering = FALSE;
		break;
	case DUAL_SWITCH( TRUE,  TRUE):		/* Never left solid area */
		sil->n_sis            = 2;
		sil->sis[0].t        = -INFINITE;
		sil->sis[0].entering = TRUE;
		sil->sis[1].t        = INFINITE;
		sil->sis[1].entering = FALSE;
		break;
	}
}
/*...e*/
	else
/*...stwo roots:16:*/
/*
Where I say t1-SIDESTEP, I used to say the more obvious (t1+t2)*0.5.
I changed it because if t1 and t2 are very close, inaccuracies in the
arithmetic can cause the value of shape evaluated at (t1+t2)*0.5 to be just the
wrong side of 0.0. The newer test is much more likely to give the right result
as we are now evaluating the shape function further away from any roots.
*/

{
double	rooted = sqrt(qs);
double	t1 = (-qb - rooted) / (qa + qa);
double	t2 = (-qb + rooted) / (qa + qa);

if ( t1 > t2 )
	/* Ensure t1 is lower than t2 */
	{ double t = t1; t1 = t2; t2 = t; }

if ( (*value_of_shape)(shapeinfo, t_along_pq(p, q, t1-SIDESTEP)) > 0.0 )
	/* Middle part is the solid part */
	{
	sil->n_sis            = 2;
	sil->sis[0].t        = t1;
	sil->sis[0].entering = TRUE;
	sil->sis[1].t        = t2;
	sil->sis[1].entering = FALSE;
	}
else
	/* Regions before and after middle part are the solid parts */
	{
	sil->n_sis            = 4;
	sil->sis[0].t        = -INFINITE;
	sil->sis[0].entering = TRUE;
	sil->sis[1].t        = t1;
	sil->sis[1].entering = FALSE;
	sil->sis[2].t        = t2;
	sil->sis[2].entering = TRUE;
	sil->sis[3].t        = INFINITE;
	sil->sis[3].entering = FALSE;
	}
}
/*...e*/
	}
/*...e*/
