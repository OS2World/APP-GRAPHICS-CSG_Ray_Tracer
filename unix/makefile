#
# CSG Raytracer
# UNIX
#

GBM =		../../gbm

ifdef AIXxlc
CC =		cc
CFLAGS =	-c -DUNIX -DAIX
else
ifdef AIX
CC =		gcc
CFLAGS =	-c -DUNIX -DAIX
else
ifdef HP
CC =		cc
CFLAGS =	-c -DUNIX -DHP    -Aa
else
ifdef SUN
CC =		cc
CFLAGS =	-c -DUNIX -DSUN   -Xa
else
ifdef XNX
CC =		cc
CFLAGS =	-c -DUNIX -DXNX
else
CC =		cc
CFLAGS =	-c -DUNIX -DLINUX
endif
endif
endif
endif
endif

.SUFFIXES:	.c .o

.c.o:
		$(CC) $(CFLAGS) -I$(GBM) $*.c

#

rt:		fio.o tex.o vector.o rgbvec.o col.o surf.o sil.o plane.o biplane.o sphere.o quad.o shape.o rt.o $(GBM)/gbm.a
		$(CC) -o rt fio.o tex.o vector.o rgbvec.o col.o surf.o sil.o plane.o biplane.o sphere.o quad.o shape.o rt.o $(GBM)/gbm.a -lm

fio.o:		fio.c $(GBM)/gbm.h

tex.o:		tex.c

vector.o:	vector.c vector.h

rgbvec.o:	rgbvec.c rgbvec.h

col.o:		col.c rt.h fio.h tex.h vector.h rgbvec.h

surf.o:		surf.c rt.h fio.h tex.h vector.h rgbvec.h col.h surf.h

sil.o:		sil.c rt.h vector.h sil.h

plane.o:	plane.c rt.h vector.h sil.h plane.h

biplane.o:	biplane.c rt.h vector.h sil.h biplane.h

sphere.o:	sphere.c rt.h vector.h sil.h sphere.h

quad.o:		quad.c rt.h vector.h sil.h quad.h

shape.o:	shape.c rt.h fio.h tex.h vector.h rgbvec.h col.h surf.h sil.h plane.h biplane.h sphere.h quad.h shape.h

rt.o:		rt.c rt.h fio.h tex.h vector.h rgbvec.h col.h surf.h sil.h plane.h biplane.h sphere.h quad.h shape.h

#

clean:
		@rm -f *.o

#

package:
		zip -q -r rt.zip *
