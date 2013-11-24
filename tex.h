/*

tex.h - Interface to 3D texture map handling code

*/

typedef void TEX;

extern TEX *read_tex(char *fn);
extern TEX *copy_tex(TEX *tex);
extern void destroy_tex(TEX *tex);

extern void get_voxel_tex(
	TEX *tex,
	int x, int y, int z,
	byte *r, byte *g, byte *b
	);

extern int width_tex(TEX *tex);
extern int height_tex(TEX *tex);
extern int depth_tex(TEX *tex);
