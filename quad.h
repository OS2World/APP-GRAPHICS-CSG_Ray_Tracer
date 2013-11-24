/*

quad.h - Interface to arbitrary quadratic logic

*/

typedef struct { double a, b, c, d, e, f, g, h, i, j; } QUAD;

#ifndef _QUAD_

extern QUAD *create_quad(
	double a, double b, double c,
	double d, double e, double f,
	double g, double h, double i,
	double j
	);

extern QUAD *create_ellipsoid(double rx, double ry, double rz);

extern QUAD *create_x_ell_cyl(double ry, double rz);
extern QUAD *create_y_ell_cyl(double rx, double rz);
extern QUAD *create_z_ell_cyl(double rx, double ry);
extern QUAD *create_x_cyl(double r);
extern QUAD *create_y_cyl(double r);
extern QUAD *create_z_cyl(double r);

extern QUAD *create_x_ell_cone(double ky, double kz);
extern QUAD *create_y_ell_cone(double kx, double kz);
extern QUAD *create_z_ell_cone(double kx, double ky);
extern QUAD *create_x_cone(double k);
extern QUAD *create_y_cone(double k);
extern QUAD *create_z_cone(double k);

extern QUAD *copy_quad(QUAD *quad);
extern void destroy_quad(QUAD *quad);

extern void trans_x_quad(QUAD *quad, double t);
extern void trans_y_quad(QUAD *quad, double t);
extern void trans_z_quad(QUAD *quad, double t);
extern void scale_x_quad(QUAD *quad, double factor);
extern void scale_y_quad(QUAD *quad, double factor);
extern void scale_z_quad(QUAD *quad, double factor);
extern void rot_x_quad(QUAD *quad, double angle);
extern void rot_y_quad(QUAD *quad, double angle);
extern void rot_z_quad(QUAD *quad, double angle);

extern void extent_xyz_quad(QUAD *quad, VECTOR *min, VECTOR *max);
extern int isects_reqd_quad(QUAD *quad);
extern void intersect_quad(QUAD *quad, VECTOR p, VECTOR q, SIL *sil);
extern VECTOR normal_to_quad(QUAD *quad, VECTOR p);

#endif
