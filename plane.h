/*

plane.h - Interface to half-plane logic

*/

typedef struct { double a, b, c, d; } PLANE;

#ifndef _PLANE_

extern PLANE *create_plane(double a, double b, double c, double d);
extern PLANE *create_x_lt_plane(double x);
extern PLANE *create_x_gt_plane(double x);
extern PLANE *create_y_lt_plane(double y);
extern PLANE *create_y_gt_plane(double y);
extern PLANE *create_z_lt_plane(double z);
extern PLANE *create_z_gt_plane(double z);

extern PLANE *copy_plane(PLANE *plane);
extern void destroy_plane(PLANE *plane);

extern void trans_x_plane(PLANE *plane, double t);
extern void trans_y_plane(PLANE *plane, double t);
extern void trans_z_plane(PLANE *plane, double t);
extern void scale_x_plane(PLANE *plane, double factor);
extern void scale_y_plane(PLANE *plane, double factor);
extern void scale_z_plane(PLANE *plane, double factor);
extern void rot_x_plane(PLANE *plane, double angle);
extern void rot_y_plane(PLANE *plane, double angle);
extern void rot_z_plane(PLANE *plane, double angle);

extern int isects_reqd_plane(PLANE *plane);
extern void intersect_plane(PLANE *plane, VECTOR p, VECTOR q, SIL *sil);
extern VECTOR normal_to_plane(PLANE *plane);

#endif
