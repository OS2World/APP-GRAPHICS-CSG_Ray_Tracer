#
# CSG Raytracer
#
# 32 bit OS/2
# IBM VisualAge C++, with CTC306, CTD302 and CTU304
#
# Incorrect output seen with /O+ /Os-.
# Indeed with earlier levels of compiler, floating point exceptions were seen.
#

GBM =		..\..\gbm

CWARNS =	/W3 /Wcmp+cnd+dcl+ord+par+use+
CFLAGS =	/C+ /I$(GBM) /G4 /Gd-e+m+ /Gf+ $(CWARNS) /Q+ /DOS2 /O-
LFLAGS =	/NOI /ALIGN:4 /BASE:0x10000 /STACK:0x30000 /PMTYPE:VIO /NOLOGO

.SUFFIXES:	.c .obj

.c.obj:
		icc $(CFLAGS) $*.c

#

rt.exe:		fio.obj tex.obj vector.obj rgbvec.obj col.obj surf.obj \
		sil.obj plane.obj biplane.obj sphere.obj quad.obj shape.obj \
		rt.obj $(GBM)\gbm.lib
		ilink $(LFLAGS) /OUT:$@ $**

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
