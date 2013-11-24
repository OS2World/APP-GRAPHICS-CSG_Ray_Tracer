/*

col.h - Interface to colour code

When evaluating a colour, 3 values are passed as p0,p1,p2.
Some types of colour modify these and call on recursively.

*/

typedef void COL;

extern COL *create_const_col(RGBVEC rgbvec);
extern COL *create_no_move_col(COL *col);
extern COL *create_interp0_col(COL *col);
extern COL *create_interp1_col(COL *col);
extern COL *create_interp2_col(COL *col);
extern COL *create_2d_field_col(double bx, double by, BITMAP *bitmap);
extern COL *create_3d_field_col(double bx, double by, double bz, TEX *tex);
extern COL *create_remap_col(VECTOR base, VECTOR v0, VECTOR v1, VECTOR v2, COL *col);
extern COL *create_cyl_polar_col(double lond, double rd, double hd, COL *col);
extern COL *create_sph_polar_col(double lond, double latd, double rd, COL *col);
extern COL *create_2d_matrix_col(double m[2][2], COL *col);
extern COL *create_3d_matrix_col(double m[3][3], COL *col);

extern COL *copy_col(COL *col);

extern void destroy_col(COL *col);

extern void trans_x_col(COL *col, double t);
extern void trans_y_col(COL *col, double t);
extern void trans_z_col(COL *col, double t);
extern void scale_x_col(COL *col, double factor);
extern void scale_y_col(COL *col, double factor);
extern void scale_z_col(COL *col, double factor);
extern void rot_x_col(COL *col, double angle);
extern void rot_y_col(COL *col, double angle);
extern void rot_z_col(COL *col, double angle);

extern RGBVEC evaluate_col(COL *col, double p0, double p1, double p2);
