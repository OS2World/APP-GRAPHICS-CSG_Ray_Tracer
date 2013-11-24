/*

fio.c - Bitmap file IO

This module uses the Generalised Bitmap Module to do the I/O.
This code is also portable, therefore the whole raytracer still is.

*/

/*...sincludes:0:*/
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <string.h>
#include <memory.h>
#include <malloc.h>
#if defined(AIX) || defined(LINUX)
#include <unistd.h>
#else
#include <io.h>
#endif
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "rt.h"
#include "gbm.h"

#ifndef O_BINARY
#define	O_BINARY 0
#endif
/*...e*/

typedef struct
	{
	GBM gbm;
	GBMRGB gbmrgb[0x100];
	byte *data;
	int stride;
	int ref_count;
	} BITMAP;

/*...sfio_init:0:*/
void fio_init(void)
	{
	gbm_init();
	}
/*...e*/
/*...sfio_deinit:0:*/
void fio_deinit(void)
	{
	gbm_deinit();
	}
/*...e*/
/*...sfio_create_bitmap:0:*/
BITMAP *fio_create_bitmap(int w, int h)
	{
	BITMAP *bitmap;

	if ( (bitmap = malloc(sizeof(BITMAP))) == NULL )
		return NULL;

	bitmap->gbm.w   = w;
	bitmap->gbm.h   = h;
	bitmap->gbm.bpp = 24;
	bitmap->stride  = ((w * 24 + 31)/32)*4;
	if ( (bitmap->data = malloc(bitmap->stride * h)) == NULL )
		{ free(bitmap); return NULL; }

	/* Initialse surface to medium grey */
	memset(bitmap->data, 0x80, bitmap->stride * h);

	bitmap->ref_count = 1;

	return bitmap;
	}
/*...e*/
/*...sfio_copy_bitmap:0:*/
BITMAP *fio_copy_bitmap(BITMAP *bitmap)
	{
	(bitmap->ref_count)++;
	return bitmap;
	}
/*...e*/
/*...sfio_destroy_bitmap:0:*/
void fio_destroy_bitmap(BITMAP *bitmap)
	{
	if ( --(bitmap->ref_count) == 0 )
		{
		free(bitmap->data);
		free(bitmap);
		}
	}
/*...e*/
/*...sfio_read_bitmap:0:*/
BITMAP *fio_read_bitmap(char *fn)
	{
	BITMAP *bitmap;
	char *opt;
	int fd, ft;

	if ( (opt = strchr(fn, ',')) != NULL )
		*opt++ = '\0';
	else
		opt = "";

	if ( gbm_guess_filetype(fn, &ft) != GBM_ERR_OK )
		return NULL;

	if ( (fd = gbm_io_open(fn, O_RDONLY | O_BINARY)) == -1 )
		return NULL;

	if ( (bitmap = malloc(sizeof(BITMAP))) == NULL )
		{ gbm_io_close(fd); return NULL; }

	if ( gbm_read_header(fn, fd, ft, &(bitmap->gbm), opt) != GBM_ERR_OK )
		{ free(bitmap); gbm_io_close(fd); return NULL; }

	if ( gbm_read_palette(fd, ft, &(bitmap->gbm), bitmap->gbmrgb) != GBM_ERR_OK )
		{ free(bitmap); gbm_io_close(fd); return NULL; }

	bitmap->stride = ((bitmap->gbm.w * bitmap->gbm.bpp + 31)/32)*4;
	if ( (bitmap->data = malloc(bitmap->stride * bitmap->gbm.h)) == NULL )
		{ free(bitmap); gbm_io_close(fd); return NULL; }

	if ( gbm_read_data(fd, ft, &(bitmap->gbm), bitmap->data) )
		{ free(bitmap); gbm_io_close(fd); return NULL; }

	gbm_io_close(fd);

	bitmap->ref_count = 1;

	return bitmap;
	}
/*...e*/
/*...sfio_write_bitmap:0:*/
BOOLEAN fio_write_bitmap(BITMAP *bitmap, char *fn)
	{
	char *opt;
	int fd, ft;
	GBMFT gbmft;

	if ( (opt = strchr(fn, ',')) != NULL )
		*opt++ = '\0';
	else
		opt = "";

	if ( gbm_guess_filetype(fn, &ft) != GBM_ERR_OK )
		return FALSE;

	gbm_query_filetype(ft, &gbmft);

	if ( (gbmft.flags & GBM_FT_W24) == 0 )
		return FALSE;

	if ( (fd = gbm_io_create(fn, O_WRONLY|O_BINARY)) == -1 )
		return FALSE;

	if ( gbm_write(fn, fd, ft, &(bitmap->gbm), bitmap->gbmrgb,
			bitmap->data, opt) != GBM_ERR_OK )
		{ gbm_io_close(fd); unlink(fn); return FALSE; }

	gbm_io_close(fd);

	return TRUE;
	}
/*...e*/
/*...sfio_get_pixel:0:*/
void fio_get_pixel(
	BITMAP *bitmap,
	int x, int y,
	byte *r, byte *g, byte *b
	)
	{
	byte inx, *data = bitmap->data + y * bitmap->stride;

	switch ( bitmap->gbm.bpp )
		{
/*...s1:16:*/
case 1:
	inx = data[x >> 3];
	inx >>= ( 7 - (x & 7) );
	inx &= 1;	
	*b = bitmap->gbmrgb[inx].b;
	*g = bitmap->gbmrgb[inx].g;
	*r = bitmap->gbmrgb[inx].r;
	break;
/*...e*/
/*...s4:16:*/
case 4:
	inx = data[x >> 1];
	if ( x & 1 )
		inx &= 0x0f;
	else
		inx >>= 4;
	*b = bitmap->gbmrgb[inx].b;
	*g = bitmap->gbmrgb[inx].g;
	*r = bitmap->gbmrgb[inx].r;
	break;
/*...e*/
/*...s8:16:*/
case 8:
	inx = data[x];
	*b = bitmap->gbmrgb[inx].b;
	*g = bitmap->gbmrgb[inx].g;
	*r = bitmap->gbmrgb[inx].r;
	break;
/*...e*/
/*...s24:16:*/
case 24:
	data += (x * 3);
	*b = *data++;
	*g = *data++;
	*r = *data;
	break;
/*...e*/
		}
	}
/*...e*/
/*...sfio_set_pixel:0:*/
void fio_set_pixel(
	BITMAP *bitmap,
	int x, int y,
	byte r, byte g, byte b
	)
	{
	int stride = bitmap->stride;
	byte *data = bitmap->data + y * stride + x * 3;

	*data++ = b;
	*data++ = g;
	*data   = r;
	}
/*...e*/
/*...sfio_width:0:*/
int fio_width(BITMAP *bitmap)
	{
	return bitmap->gbm.w;
	}
/*...e*/
/*...sfio_height:0:*/
int fio_height(BITMAP *bitmap)
	{
	return bitmap->gbm.h;
	}
/*...e*/
