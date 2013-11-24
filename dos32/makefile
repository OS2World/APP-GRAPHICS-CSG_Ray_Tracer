#
# CSG Raytracer
# 32 bit DOS
# Watcom C/C++ 10.6 + CauseWay 1.3
#

GBM =		..\..\gbm

CFLAGS =	-bt=DOS -dDOS32 -oit -4r -s -w3 -we -zp4 -mf -zq -fr -i=$(GBM)

.SUFFIXES:	.c .obj

.c.obj:
		wcc386 $(CFLAGS) $*.c

#

rt.exe:		fio.obj tex.obj vector.obj rgbvec.obj col.obj surf.obj \
		sil.obj plane.obj biplane.obj sphere.obj quad.obj shape.obj \
		rt.obj $(GBM)\gbm.lib
		wlink @<<
			System CauseWay
			Name $@
			File fio.obj,tex.obj,vector.obj,rgbvec.obj,col.obj
			File surf.obj,sil.obj,plane.obj,biplane.obj,sphere.obj
			File quad.obj,shape.obj,rt.obj
			Library $(GBM)\gbm.lib
			Option Quiet
			Option Stack = 0x6000
<<

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
