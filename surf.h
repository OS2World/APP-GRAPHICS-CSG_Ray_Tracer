/*

surf.h - Surface datatype

*/

typedef struct
	{
	double	ka, kd, ks, kt;		/* Ambient, diffuse, specular, trans.*/
	COL	*od, *os;		/* Diffuse and specular colours      */
	double	phong;			/* Phong power number                */
	double	rinx;			/* Refractive index of inside        */
	} SURF;

#ifndef _SURF_

extern SURF *create_surf(
	double ka, double kd, double ks, double kt,
	COL *od, COL *os,
	double phong,
	double rinx
	);

extern SURF *copy_surf(SURF *surf);

extern void destroy_surf(SURF *surf);

extern void trans_x_surf(SURF *surf, double t);
extern void trans_y_surf(SURF *surf, double t);
extern void trans_z_surf(SURF *surf, double t);
extern void scale_x_surf(SURF *surf, double factor);
extern void scale_y_surf(SURF *surf, double factor);
extern void scale_z_surf(SURF *surf, double factor);
extern void rot_x_surf(SURF *surf, double angle);
extern void rot_y_surf(SURF *surf, double angle);
extern void rot_z_surf(SURF *surf, double angle);

#endif
