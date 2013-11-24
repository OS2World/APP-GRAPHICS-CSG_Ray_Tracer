/*

sphere.h - Interface to sphere logic

*/

typedef struct { double a, b, c, d; } SPHERE;

#ifndef _SPHERE_

extern SPHERE *create_sphere(double r);
extern SPHERE *copy_sphere(SPHERE *sphere);
extern void destroy_sphere(SPHERE *sphere);

extern void trans_x_sphere(SPHERE *sphere, double t);
extern void trans_y_sphere(SPHERE *sphere, double t);
extern void trans_z_sphere(SPHERE *sphere, double t);
extern void scale_sphere(SPHERE *sphere, double factor);
extern void rot_x_sphere(SPHERE *sphere, double angle);
extern void rot_y_sphere(SPHERE *sphere, double angle);
extern void rot_z_sphere(SPHERE *sphere, double angle);

extern int isects_reqd_sphere(SPHERE *sphere);
extern void intersect_sphere(SPHERE *sphere, VECTOR p, VECTOR q, SIL *sil);
extern VECTOR normal_to_sphere(SPHERE *sphere, VECTOR p);

#endif
