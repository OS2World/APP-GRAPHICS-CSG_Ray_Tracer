/*

fio.h - Interface to bitmap file IO

*/

typedef void BITMAP;

extern void fio_init(void);
extern void fio_deinit(void);
extern BITMAP *fio_create_bitmap(int w, int h);
extern BITMAP *fio_copy_bitmap(BITMAP *bitmap);
extern void fio_destroy_bitmap(BITMAP *bitmap);
extern BITMAP *fio_read_bitmap(char *fn);
extern BOOLEAN fio_write_bitmap(BITMAP *bitmap, char *fn);

extern void fio_get_pixel(
	BITMAP *bitmap,
	int x, int y,
	byte *r, byte *g, byte *b
	);

extern void fio_set_pixel(
	BITMAP *bitmap,
	int x, int y,
	byte r, byte g, byte b
	);

extern int fio_width(BITMAP *bitmap);
extern int fio_height(BITMAP *bitmap);
