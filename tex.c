/*

tex.c - Handle 3D texture map data

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

#ifndef O_BINARY
#define	O_BINARY 0
#endif
/*...e*/

typedef struct { byte b, g, r, dummy; } PAL;

typedef struct
	{
	dword w, h, d, bpv, stride;
	PAL pals[0x100];
	byte *data;
	int ref_count;
	} TEX;

#define	TEX_MAGIC	((dword) 0x1a584554)

/*...sread_tex:0:*/
/*...sread_dword:0:*/
static dword read_dword(int fd)
	{
	byte b[4];

	read(fd, b, 4);
	return   (dword) b[0]        +
		((dword) b[1] <<  8) +
		((dword) b[2] << 16) +
		((dword) b[3] << 24) ;
	}
/*...e*/

TEX *read_tex(char *fn)
	{
	TEX *tex;
	int fd, i;
	dword bytes;

	if ( (fd = open(fn, O_RDONLY | O_BINARY)) == -1 )
		return NULL;

	if ( (tex = malloc(sizeof(TEX))) == NULL )
		{ close(fd); return NULL; }

	if ( read_dword(fd) != TEX_MAGIC )
		{ free(tex); close(fd); return NULL; }

	tex->w   = read_dword(fd);
	tex->h   = read_dword(fd);
	tex->d   = read_dword(fd);
	tex->bpv = read_dword(fd);

	if ( tex->bpv != 1 && tex->bpv != 4 && tex->bpv != 8 && tex->bpv != 24 )
		{ free(tex); close(fd); return NULL; }

	for ( i = 0; i < ((1 << tex->bpv) & 0x1ff); i++ )
		if ( read(fd, (char *) &(tex->pals[i]), sizeof(PAL)) != sizeof(PAL) )
			{ free(tex); close(fd); return NULL; }

	tex->stride = ((tex->w * tex->bpv + 31)/32)*4;
	bytes = tex->stride * tex->h * tex->d;

	if ( (tex->data = malloc((size_t) bytes)) == NULL )
		{ free(tex); close(fd); return NULL; }

	if ( read(fd, tex->data, (unsigned int) bytes) != (int) bytes )
		{ free(tex->data); free(tex); close(fd); return NULL; }

	close(fd);

	tex->ref_count = 1;

	return tex;
	}
/*...e*/
/*...scopy_tex:0:*/
TEX *copy_tex(TEX *tex)
	{
	(tex->ref_count)++;
	return tex;
	}
/*...e*/
/*...sdestroy_tex:0:*/
void destroy_tex(TEX *tex)
	{
	if ( --(tex->ref_count) == 0 )
		{
		free(tex->data);
		free(tex);
		}
	}
/*...e*/

/*...sget_voxel_tex:0:*/
void get_voxel_tex(
	TEX *tex,
	int x, int y, int z,
	byte *r, byte *g, byte *b
	)
	{
	byte inx, *data = tex->data + (z * tex->h + y) * tex->stride;

	switch ( tex->bpv )
		{
/*...s1:16:*/
case 1:
	inx = data[x >> 3];
	inx >>= ( 7 - (x & 7) );
	inx &= 1;	
	*b = tex->pals[inx].b;
	*g = tex->pals[inx].g;
	*r = tex->pals[inx].r;
	break;
/*...e*/
/*...s4:16:*/
case 4:
	inx = data[x >> 1];
	if ( x & 1 )
		inx &= 0x0f;
	else
		inx >>= 4;
	*b = tex->pals[inx].b;
	*g = tex->pals[inx].g;
	*r = tex->pals[inx].r;
	break;
/*...e*/
/*...s8:16:*/
case 8:
	inx = data[x];
	*b = tex->pals[inx].b;
	*g = tex->pals[inx].g;
	*r = tex->pals[inx].r;
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
/*...swidth_tex:0:*/
int width_tex(TEX *tex)
	{
	return (int) tex->w;
	}
/*...e*/
/*...sheight_tex:0:*/
int height_tex(TEX *tex)
	{
	return (int) tex->h;
	}
/*...e*/
/*...sdepth_tex:0:*/
int depth_tex(TEX *tex)
	{
	return (int) tex->d;
	}
/*...e*/
