#
# CSG Raytracer
#
# Win32
# MSVC
#

GBM =		..\..\gbm\gbm

CFLAGS =	/c /W3 /DWIN32 /G5 /Gs /Oit /nologo /WX /MT /I$(GBM)
LFLAGS =	/NOLOGO /INCREMENTAL:NO /STACK:0x30000

.SUFFIXES:	.c .obj

.c.obj:
		cl $(CFLAGS) $*.c

#

rt.exe:		fio.obj tex.obj vector.obj rgbvec.obj col.obj surf.obj \
		sil.obj plane.obj biplane.obj sphere.obj quad.obj shape.obj \
		rt.obj $(GBM)\gbm.lib
		link $(LFLAGS) /OUT:$@ $**

fio.obj:	fio.c $(GBM)\gbm.h

tex.obj:	tex.c

vector.obj:	vector.c vector.h

rgbvec.obj:	rgbvec.c rgbvec.h

col.obj:	col.c rt.h fio.h tex.h vector.h rgbvec.h

surf.obj:	surf.c rt.h fio.h tex.h vector.h rgbvec.h col.h surf.h

sil.obj:	sil.c rt.h vector.h sil.h

plane.obj:	plane.c rt.h vector.h sil.h plane.h

biplane.obj:	biplane.c rt.h vector.h sil.h biplane.h

sphere.obj:	sphere.c rt.h vector.h sil.h sphere.h

quad.obj:	quad.c rt.h vector.h sil.h quad.h

shape.obj:	shape.c rt.h fio.h tex.h vector.h rgbvec.h col.h surf.h sil.h plane.h biplane.h sphere.h quad.h shape.h

rt.obj:		rt.c rt.h fio.h tex.h vector.h rgbvec.h col.h surf.h sil.h plane.h biplane.h sphere.h quad.h shape.h

#

clean:
		-del *.obj 2> nul
