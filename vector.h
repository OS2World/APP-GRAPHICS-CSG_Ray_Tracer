/*

vector.h - Interface to vector datatype

*/

typedef struct { double	x, y, z; } VECTOR;

#ifndef _VECTOR_

extern double magnitude(VECTOR p);
extern VECTOR unit_vector(VECTOR p);
extern VECTOR vector_product(VECTOR p, VECTOR q);
extern double scalar_product(VECTOR p, VECTOR q);
extern VECTOR vector_sum(VECTOR p, VECTOR q);
extern VECTOR negate_vector(VECTOR p);
extern VECTOR vector_difference(VECTOR p, VECTOR q);
extern VECTOR scale_vector(VECTOR p, double scalar);
extern VECTOR inv_scale_vector(VECTOR p, double inv_scalar);
extern VECTOR t_along_pq(VECTOR p, VECTOR q, double t);
extern VECTOR rot_x_vector(VECTOR v, double angle);
extern VECTOR rot_y_vector(VECTOR v, double angle);
extern VECTOR rot_z_vector(VECTOR v, double angle);

#endif
