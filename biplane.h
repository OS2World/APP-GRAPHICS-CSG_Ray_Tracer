/*

biplane.h - Interface to intersection of two half-planes logic

	Biplane:	ax+by+cz+d1> 0 and
			ax+by+cz+d2<=0

*/

typedef struct { double a, b, c, d1, d2; } BIPLANE;

#ifndef _BIPLANE_

extern BIPLANE *create_biplane(double a, double b, double c, double d1, double d2);
extern BIPLANE *create_x_in_biplane(double x1, double x2);
extern BIPLANE *create_y_in_biplane(double y1, double y2);
extern BIPLANE *create_z_in_biplane(double z1, double z2);

extern BIPLANE *copy_biplane(BIPLANE *biplane);
extern void destroy_biplane(BIPLANE *biplane);

extern void trans_x_biplane(BIPLANE *biplane, double t);
extern void trans_y_biplane(BIPLANE *biplane, double t);
extern void trans_z_biplane(BIPLANE *biplane, double t);
extern void scale_x_biplane(BIPLANE *biplane, double factor);
extern void scale_y_biplane(BIPLANE *biplane, double factor);
extern void scale_z_biplane(BIPLANE *biplane, double factor);
extern void rot_x_biplane(BIPLANE *biplane, double angle);
extern void rot_y_biplane(BIPLANE *biplane, double angle);
extern void rot_z_biplane(BIPLANE *biplane, double angle);

extern int isects_reqd_biplane(BIPLANE *biplane);
extern void intersect_biplane(BIPLANE *biplane, VECTOR p, VECTOR q, SIL *sil);
extern VECTOR normal_to_biplane(BIPLANE *biplane, VECTOR p);

#endif
