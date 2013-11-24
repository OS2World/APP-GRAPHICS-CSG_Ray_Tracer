/*

rt.h - Stuff global to raytracer

*/

#ifndef RT_H
#define	RT_H

#ifndef BOOLEAN_DEFINED
#define	BOOLEAN_DEFINED
typedef	int BOOLEAN;
#define	TRUE  1
#define	FALSE 0
#endif

#ifndef BASICTYPES_DEFINED
#define	BASICTYPES_DEFINED
typedef unsigned  char  byte;
typedef unsigned short  word;
typedef unsigned  long dword;
#endif

#define	INFINITE	(1.0e+20)	/* A hoooge number                   */
#ifndef PI
#define	PI		3.14159265359	/* Pi                                */
#endif

#ifndef min
#define	min(a,b)	(((a)<(b))?(a):(b))
#endif

#ifndef max
#define	max(a,b)	(((a)>(b))?(a):(b))
#endif

#endif
