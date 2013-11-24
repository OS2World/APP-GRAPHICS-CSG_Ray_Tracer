/*

rt.c - CSG Ray Tracer

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
#include <math.h>
/*...e*/
/*...sisalnum fix:0:*/
#ifdef LINUX
/* On Slackware 3.4 /lib/libc.so.5 -> libc.5.4.33
   isalnum(c) is implemented as (__ctype_b(c)&_ISalnum)
   Programs compiled on Slackware therefore look for this bit
   On RedHat 6.0, /usr/i486-linux-libc5/lib/libc.so.5 -> libc.so.5.3.12
   isalnum(c) is probably done as (__ctype_b(c)&(_ISalnum|_ISdigit))
   So the __ctype_b array doesn't have this bit.
   So programs compiled, using isalnum on Slackware 3.4, don't work
   when run on RedHat 6.0. Best to avoid it. */
#undef isalnum
#define isalnum(c) (isalpha(c)||isdigit(c))
#endif
/*...e*/

static char progname[] = "rt";

/*...suseful:0:*/
/*...sfatal:0:*/
static void fatal(const char *fmt, ...)
	{
	va_list	vars;
	char	s[256+1];

	va_start(vars, fmt);
	vsprintf(s, fmt, vars);
	va_end(vars);
	fprintf(stderr, "%s: %s\n", progname, s);
	exit(1);
	}
/*...e*/
/*...smemcheck:0:*/
static void *memcheck(void *p)
	{
	if ( p == NULL )
		fatal("out of memory");
	return p;
	}
/*...e*/
/*...sstrsave:0:*/
static char *strsave(char *s)
	{
	int	len = strlen(s);
	char	*t = memcheck(malloc(len + 1));

	return strcpy(t, s);
	}
/*...e*/
/*...e*/
/*...smain:0:*/
/*...sincludes:0:*/
#include "rt.h"
#include "fio.h"
#include "tex.h"
#include "vector.h"
#include "rgbvec.h"
#include "col.h"
#include "surf.h"
#include "sil.h"
#include "plane.h"
#include "biplane.h"
#include "sphere.h"
#include "quad.h"
#include "shape.h"

/*...vrt\46\h:0:*/
/*...vfio\46\h:0:*/
/*...vtex\46\h:0:*/
/*...vvector\46\h:0:*/
/*...vrgbvec\46\h:0:*/
/*...vcol\46\h:0:*/
/*...vsurf\46\h:0:*/
/*...vsil\46\h:0:*/
/*...vplane\46\h:0:*/
/*...vbiplane\46\h:0:*/
/*...vsphere\46\h:0:*/
/*...vquad\46\h:0:*/
/*...vshape\46\h:0:*/
/*...e*/

#ifdef OS2
#include <float.h>
#endif

/*...slighting:0:*/
/*
Surfaces do not emit light of their own. They only reflect and refract light
from the light sources in the scene. The lighting/shading model used is the
one defined in F,vD,F,H p734.

Global parameters :-
	Ia = ambient intensity (rgb vector)

Parameters defined per surface :-
	ka               = ambient-reflection-coefficient
	kd               = diffuse coefficient, range 0.0 to 1.0
	ks               = specular coefficient, range 0.0 to 1.0
	kt               = transmission coefficent, range 0.0 to 1.0
	od               = diffuse colour (rgb vector)
	os               = specular colour (rgb vector)
	phong_number     = Phong power
	rinx = inside relative to outside of surface

Ambient contribution is :-
	Ia * ka * od

Diffuse for light of intensity Ip, attenuated by distance by fatt, using
Lamberts law gives :-
	fatt * Ip * kd * od * cos_nl

Specular for light of intensity Ip, attenuated by distance by fatt, using
the Phong shading model gives :-
	fatt * Ip * ks * os * pow(cos_rv, phong_number)

Reflective component (specular only), using recursion (if depth not too deep) :-
	intensity_from_looking_along_reflected_ray * ks * os

Refractive component (transmissive), given no total internal reflection,
using recursion (if depth not too deep) :-
	intensity_from_refractive_ray * kt

Before returning intensity, scale down by distance to origin of ray (eye).
*/

typedef struct { VECTOR posn; RGBVEC i; } LIGHT;

#define	N_LIGHTS	10
static LIGHT lights[N_LIGHTS];
static int n_lights = 0;

static RGBVEC i_background = { 0.0, 0.0, 0.0 };
static RGBVEC i_ambient = { 0.0, 0.0, 0.0 };
static double af1 = 1.0, af2 = 0.9;

/*...sscale_by_distance:0:*/
/*
The further light travels (from light source to surface, or from surface to
eye etc.) the fainter it gets. Therefore as t increases, the amount intensity
is reduced by should increase.
*/

static RGBVEC scale_by_distance(RGBVEC i, double t)
	{
	double	scalar = af1 * pow(af2, t);

	if ( scalar > 1.0 )
		scalar = 1.0;

	i.r *= scalar;
	i.g *= scalar;
	i.b *= scalar;

	return i;
	}
/*...e*/
/*...e*/
/*...strace:0:*/
#define	EPSILON		(1.0e-8)		/* A tiny number             */

/*...sshadow_calc:0:*/
/*

In a CSG system without refraction, no rays can pass through any solid.
Hence this code would return 1.0 or 0.0 depending on whether il is empty.
Actually we are only interested in the range of t from 0 to dist_l being empty.
This is because the light is dist_l away, things further don't cause shadows.

With refraction, things can be lit through glass-like solids.
Let me admit, from the start, that refraction in a CSG system is a kludge!
Glass-like solids transmit kt of the light coming from their other side.
So we look along the intersections for transmissive surfaces and work out the
combined transmissiveness.

We assume that pathalogical shapes (transmissive on one side, non-transmissive
on the the other) do not exist. We also assume the kt coefficients match for
a given shape on its entry and exit. We only use kt once (on entry) per shape.

Actually, to correctly handle illumination through refractive shapes is a
nightmare, and requires that we calculate a full multi-part path back to the
light. This is horrendous, and we will assume that refraction does not alter
the path (much) back to the light.

In doing so, we will prevent glass causing shadows, which is the main goal.

*/

static double shadow_calc(ISECTL *il, double dist_l)
	{
	int	j = 0;
	double	kt_accum = 1.0;

	while ( j < il->n_isects && il->isects[j].t < dist_l )
		{
		double	kt = il->isects[j].shape->surf->kt;

		if ( kt == 0.0 )
			return 0.0;
		kt_accum *= kt;
		j += 2;
		}

	return kt_accum;
	}
/*...e*/

/*...sreflect:0:*/
static VECTOR reflect(VECTOR unit_v, VECTOR unit_n, double cos_vn)
	{
	return vector_difference(scale_vector(vector_sum(unit_n, unit_n), cos_vn), unit_v);
	}
/*...e*/
/*...srefract:0:*/
/*
          ^ N             Given an incident unit vector I, approaching a
          | -             surface with unit normal vector N, compute the
  \       |               transmitted ray T. This ray will not necessarily
   \      |               be a unit vector.
    \     |               
     \  0 |               If the term under the square root is negative
      \  i|               then this indicates total internal reflection.
       \  |               
     I  \||               n_it is the relative refractive index of the
     -  -\|               incident medium relative the the transmissive
----------+----------> S  medium. Thus for air (1.0) entering crown glass
          |\           -  (1.5) this number would be 1.0/1.5 = 0.66 approx.
          | \             
          |  \            We use the equation given in Byte Magazine Dec 90.
          | 0 \           
          |  t \          
          |     \         
          |   T  \|       
          |   -  -\       
          |               
*/

static BOOLEAN refract(
	VECTOR unit_i,
	VECTOR unit_n,
	double n_it,
	VECTOR *t
	)
	{
	double	cos_ni = -scalar_product(unit_n, unit_i);
	double	under_root = 1.0 + n_it*n_it * (cos_ni*cos_ni - 1.0);
	double	n_comp;

	if ( under_root < 0.0 )
		return FALSE; /* Total internal reflection */

	n_comp = n_it * cos_ni - sqrt(under_root);

	*t = unit_vector(vector_sum(scale_vector(unit_i, n_it),
			            scale_vector(unit_n, n_comp)));

	return TRUE;
	}
/*...e*/

static RGBVEC trace(
	SHAPE *root_shape,
	VECTOR start, VECTOR direction,
	int depth,
	ISECTL *ils[]
	)
	{
	ISECTL	*il = *ils;
	double	t;
	SHAPE	*shape;
	SURF	*surf;
	RGBVEC	i, od, os;
	VECTOR	unit_direction, unit_v, isect_posn, unit_n;
	int	j;

	unit_direction = unit_vector(direction);
	unit_v = negate_vector(unit_direction);

	intersect_shape(root_shape, start, unit_direction, ils);
	t_after_isectl(il, EPSILON);

	if ( is_empty_isectl(il) )
		return i_background;

	/* Hit something */

	t          = il->isects[0].t;
	shape      = il->isects[0].shape;
	surf       = shape->surf;
	isect_posn = t_along_pq(start, unit_direction, t);
	unit_n     = unit_vector(normal_to_shape(shape, isect_posn));

	if ( il->isects[0].negate_normal )
		unit_n = negate_vector(unit_n);

	/* Calculate colours at intersection position */

	od = evaluate_col(surf->od, isect_posn.x, isect_posn.y, isect_posn.z);
	os = evaluate_col(surf->os, isect_posn.x, isect_posn.y, isect_posn.z);

	/* Ambient light */

	i.r = i_ambient.r * surf->ka * od.r;
	i.g = i_ambient.g * surf->ka * od.g;
	i.b = i_ambient.b * surf->ka * od.b;

	/* For each light source */

	for ( j = 0; j < n_lights; j++ )
/*...shandle contribution of this light source:16:*/
/*
l is the vector from the intersection point to the light.
dist_l is the distance to l.
unit_l is a unit vector pointing to the light.
We can reuse the intersection list used to hit the object for the shadow calc.
*/

{
VECTOR	l, unit_l;
double	dist_l, kt_accum;

l      = vector_difference(lights[j].posn, isect_posn);
dist_l = magnitude(l);
unit_l = inv_scale_vector(l, dist_l);

/* Can we see the light from the point of intersection */

intersect_shape(root_shape, isect_posn, unit_l, ils);
t_after_isectl(il, EPSILON);

if ( (kt_accum = shadow_calc(il, dist_l)) > 0.0 )
	{
	RGBVEC	i_light;
	VECTOR	unit_r;
	double	cos_ln, cos_rv;

	i_light = scale_rgbvec(scale_by_distance(lights[j].i, dist_l), kt_accum);

	/* Diffuse lighting, using Lambert's law */

	if ( (cos_ln = scalar_product(unit_l, unit_n)) > 0.0 )
		{
		double	kd_cos_ln = surf->kd * cos_ln;

		i.r += i_light.r * od.r * kd_cos_ln;
		i.g += i_light.g * od.g * kd_cos_ln;
		i.b += i_light.b * od.b * kd_cos_ln;
		}

	/* Specular lighting by light source, using Phong model */

	unit_r = reflect(unit_l, unit_n, cos_ln);

	if ( (cos_rv = scalar_product(unit_r, unit_v)) > 0.0 )
		{
		double	ks_cos_rv_n = surf->ks * pow(cos_rv, surf->phong);

		i.r += i_light.r * os.r * ks_cos_rv_n;
		i.g += i_light.g * os.g * ks_cos_rv_n;
		i.b += i_light.b * os.b * ks_cos_rv_n;
		}
	}
}
/*...e*/

	if ( depth > 0 )
		{
		if ( surf->ks > 0.0 )
/*...sreflection:24:*/
{
double	cos_nv;

if ( (cos_nv = scalar_product(unit_n, unit_v)) > 0.0 )
	{
	VECTOR unit_r;
	RGBVEC i_r;

	unit_r = reflect(unit_v, unit_n, cos_nv);
	i_r    = trace(root_shape, isect_posn, unit_r, depth - 1, ils);

	i.r += surf->ks * os.r * i_r.r;
	i.g += surf->ks * os.g * i_r.g;
	i.b += surf->ks * os.b * i_r.b;
	}
}
/*...e*/
		if ( surf->kt > 0.0 )
/*...srefraction:24:*/
/*
Refractive index of outside medium is assumed to be 1.0.
Bend the ray into the medium, work out where it comes out, and bend it again.
Then trace the emerging ray and accumulate its effect.
If total internal reflection occurs, entering the shape, treat as reflection.
If it occurs leaving, bounce the ray inside the shape and try to exit again.
*/

#define	MAX_BOUNCE	10

{
VECTOR	unit_r;				/* Unit refracted ray at entry       */

if ( refract(unit_direction, unit_n, 1.0 / surf->rinx, &unit_r) )
	/* Refraction has occurred */
	{
	double	t_out;			/* Value of t where leave solid      */
	SHAPE	*shape_out;		/* Shape where leave solid           */
	SURF	*surf_out;		/* Surface of shape where leave      */
	VECTOR	isect_posn_out;		/* Place where leave solid           */
	VECTOR	unit_n_out;		/* Unit inward pointing normal       */
	VECTOR	unit_r_out;		/* Unit leaving refracted ray        */
	double	kt_comp = surf->kt;	/* Composite scale down factor       */
	double	t_comp = 0.0;		/* Composite/total distance inside   */
	int	n_bounce = 0;		/* # of internal bounces             */

	isect_posn_out = isect_posn;
	unit_r_out     = unit_r;
	for ( ;; )
		{
		double	cos_rn;

		intersect_shape(root_shape, isect_posn_out, unit_r_out, ils);
		t_after_isectl(il, EPSILON);

		t_out          = il->isects[1].t;
		shape_out      = il->isects[1].shape;
		surf_out       = shape_out->surf;
		isect_posn_out = t_along_pq(isect_posn_out, unit_r_out, t_out);
		unit_n_out     = unit_vector(normal_to_shape(shape_out, isect_posn_out));

		if ( !il->isects[1].negate_normal )
			unit_n_out = negate_vector(unit_n_out);

		t_comp += t_out;

		if ( refract(unit_r_out, unit_n_out, surf_out->rinx, &unit_r_out) )
			break; /* Refracted out of solid */

		/* Total internal reflection trying to leave solid */
		/* This implies the solid has an index > 1.0 */

		if ( ++n_bounce == MAX_BOUNCE )
			break; /* Reached our bounce limit, give up! */

		cos_rn = scalar_product(unit_r_out, unit_n_out);
		unit_r_out = reflect(negate_vector(unit_r_out), unit_n_out, -cos_rn);
		kt_comp *= surf_out->kt; /* Accumulate this as we effectively re-enter solid */
		}

	if ( n_bounce < MAX_BOUNCE )
		{
		RGBVEC i_r;

		i_r = trace(root_shape, isect_posn_out, unit_r_out, depth - 1, ils);
		i_r = scale_by_distance(i_r, t_comp);

		i.r += kt_comp * i_r.r;
		i.g += kt_comp * i_r.g;
		i.b += kt_comp * i_r.b;
		}
	}
else
	/* Total internal reflection trying to enter solid */
	/* This implies the solid has an index < 1.0 */
	/* This is not actually very likely (glass is 1.5, water 1.3 etc) */
	{
	double	cos_nv;

	if ( (cos_nv = scalar_product(unit_n, unit_v)) > 0.0 )
		{
		VECTOR u_r;
		RGBVEC i_r;

		u_r = reflect(unit_v, unit_n, cos_nv);
		i_r = trace(root_shape, isect_posn, u_r, depth - 1, ils);

		i.r += surf->kt * i_r.r;
		i.g += surf->kt * i_r.g;
		i.b += surf->kt * i_r.b;
		}
	}
}
/*...e*/
		}

	return scale_by_distance(i, t);
	}
/*...e*/
/*...slogging:0:*/
#include <time.h>

static time_t t_start;
static double m = 0.0; /* Maximum brightness so far */

static void bs(int n)
	{
	while ( n-- )
		fputc('\b', stdout);
	}

static char *str_of_time(time_t t)
	{
	static char buf[30+1];

	strcpy(buf, ctime(&t));
	buf[19] = '\0';
	return buf + 11;
	}

static void log_init(char *fn)
	{
	t_start = time(NULL);

	printf("%s, started %8s,   0.0%% done, etc ________, max int 0.00 ", fn, str_of_time(t_start));
	fflush(stdout);
	m = 0.0;
	}

static void log_bright(RGBVEC rgb)
	{
	if ( rgb.r > m ) m = rgb.r;
	if ( rgb.g > m ) m = rgb.g;
	if ( rgb.b > m ) m = rgb.b;
	}

static void log_done_so_far(int done, int total)
	{
	double	percent = (100.0 * (double) done) / (double) total;
	time_t	t_fin = t_start + (((time(NULL) - t_start) * total) / done);

	bs(40);
	printf("%5.1lf%% done, etc %8s, max int %4.2lf ",
		percent, str_of_time(t_fin), m);
	fflush(stdout);
	}

static void log_done(void)
	{
	time_t	t_now = time(NULL);

	bs(40);
	printf("done %8s, took %ds, max int %4.2lf         \n",
		str_of_time(t_now), t_now - t_start, m);
	fflush(stdout);
	}
/*...e*/
/*...smake_r_g_b:0:*/
/*
Given a RGBVEC value computed for a ray,
work out a suitable set of 3 R,G and B bytes to write into the file.
Mostly the RGBVEC should be supplied to us in a good range.
Occasionally we must clamp the value.
We ought to do something more clever here.
*/

static void make_r_g_b(RGBVEC rgb, byte *r, byte *g, byte *b)
	{
	if ( rgb.r <= 1.0 ) *r = (byte) (rgb.r * 255.0); else *r = 255;
	if ( rgb.g <= 1.0 ) *g = (byte) (rgb.g * 255.0); else *g = 255;
	if ( rgb.b <= 1.0 ) *b = (byte) (rgb.b * 255.0); else *b = 255;
	}
/*...e*/
static double least_vis_diff = 1.0 / 255.0;
/*...srender:0:*/
typedef	VECTOR (*RAY_FUNC)(
		int h, int v, int hpixels, int vpixels,
		double hangle, double vangle,
		VECTOR unit_forward, VECTOR unit_right, VECTOR unit_up
		);

/*...scalc_normal_ray:0:*/
/*
Calculates a ray given usual perspective rules.
Need to know which pixel of how many we are tracing.
Need to know right and up vectors.
*/

static VECTOR calc_normal_ray(
	int h, int v, int hpixels, int vpixels,
	double hangle, double vangle,
	VECTOR unit_forward, VECTOR unit_right, VECTOR unit_up
	)
	{
	double hpixels2   = (double) (hpixels >> 1);
	double vpixels2   = (double) (vpixels >> 1);
	double vfactor    = (v - vpixels2) / vpixels2;
	double hfactor    = (h - hpixels2) / hpixels2;
	VECTOR vis_right  = scale_vector(unit_right, tan(hangle));
	VECTOR vis_up     = scale_vector(unit_up   , tan(vangle));
	VECTOR comp_up    = scale_vector(vis_up, vfactor);
	VECTOR ray_plane  = vector_sum(unit_forward, comp_up);
	VECTOR comp_right = scale_vector(vis_right, hfactor);
	return vector_sum(ray_plane, comp_right);
	}
/*...e*/
/*...scalc_escher_ray:0:*/
/*
Calculates a ray given using the twisted perspective Escher used in his
"Above and Below" and "House of Stairs" works.
*/

static VECTOR calc_escher_ray(
	int h, int v, int hpixels, int vpixels,
	double hangle, double vangle,
	VECTOR unit_forward, VECTOR unit_right, VECTOR unit_up
	)
	{
	double hpixels2     = (double) (hpixels >> 1);
	double vpixels2     = (double) (vpixels >> 1);
	double vfactor      = (v - vpixels2) / vpixels2;
	double hfactor      = (h - hpixels2) / hpixels2;
	VECTOR vis_right    = scale_vector(unit_right, tan(hangle));
	VECTOR comp_right   = scale_vector(vis_right, hfactor);
	double va           = vangle * vfactor;
	VECTOR comp_forward = scale_vector(unit_forward, cos(va));
	VECTOR comp_up      = scale_vector(unit_up, sin(va));
	return vector_sum(comp_forward, vector_sum(comp_right, comp_up));
	}
/*...e*/
/*...sbasic_render:0:*/
/*
Traditional basic render operation.
Send one ray for each pixel (excluding extras for reflection etc.).
ie: Sample once for each pixel, do NOT supersample.
Also, do not perform any 'jittering' of rays.
*/

static void basic_render(
	char *fn,
	BITMAP *bitmap,
	ISECTL **ils,
	SHAPE *shape,
	VECTOR eye, VECTOR unit_forward, VECTOR unit_right, VECTOR unit_up,
	int hpixels, int vpixels, double hangle, double vangle,
	RAY_FUNC calc_ray,
	int depth
	)
	{
	int h, v;
	log_init(fn);
	for ( v = 0; v < vpixels; v++ )
		{
		for ( h = 0; h < hpixels; h++ )
			{
			VECTOR ray = (*calc_ray)(h, v, hpixels, vpixels, hangle, vangle, unit_forward, unit_right, unit_up);
			RGBVEC rgb = trace(shape, eye, ray, depth, ils);
			byte r, g, b;
			log_bright(rgb);
			make_r_g_b(rgb, &r, &g, &b);
			fio_set_pixel(bitmap, h, v, r, g, b);
			}
		log_done_so_far(v + 1, vpixels);
		}
	log_done();
	}
/*...e*/
/*...swhitted_render:0:*/
/*
Render using Whitted adaptive supersampling.
Send rays at the 4 corners of a pixel.
Ask function whitted for a combined pixel value.
*/

typedef struct { VECTOR ray; RGBVEC rgb; BOOLEAN traced; } SAMPLE;

#define	MAX_SS	4		/* Supersample by upto 4 each way */

static long normal = 0, extra = 0;
/*...swhitted:0:*/
/*
Return a combined weighting of the 4 corner samples.
If they differ 'notably', break the area into 4 sub-cells.
Generate the 5 extra samples needed for 4 sub-cells.
Then call ourself recursively to get colours for each cell, and return average.
*/

/*...saverage_of_4:0:*/
static RGBVEC average_of_4(RGBVEC a, RGBVEC b, RGBVEC c, RGBVEC d)
	{
	RGBVEC av;
	av.r = (a.r + b.r + c.r + d.r) * 0.25;
	av.g = (a.g + b.g + c.g + d.g) * 0.25;
	av.b = (a.b + b.b + c.b + d.b) * 0.25;
	return av;
	}
/*...e*/
/*...smidpoint_of_2:0:*/
static VECTOR midpoint_of_2(VECTOR a, VECTOR b)
	{
	VECTOR mid;
	mid.x = (a.x + b.x) * 0.5;
	mid.y = (a.y + b.y) * 0.5;
	mid.z = (a.z + b.z) * 0.5;
	return mid;
	}
/*...e*/
/*...smidpoint_of_4:0:*/
static VECTOR midpoint_of_4(VECTOR a, VECTOR b, VECTOR c, VECTOR d)
	{
	VECTOR mid;
	mid.x = (a.x + b.x + c.x + d.x) * 0.25;
	mid.y = (a.y + b.y + c.y + d.y) * 0.25;
	mid.z = (a.z + b.z + c.z + d.z) * 0.25;
	return mid;
	}
/*...e*/
/*...sis_close:0:*/
/*
2 colours are not close when their red green and blue components differ
more than the least visible difference.
*/

static BOOLEAN is_close(RGBVEC a, RGBVEC b)
	{
	return fabs(a.r - b.r) < least_vis_diff &&
	       fabs(a.g - b.g) < least_vis_diff &&
	       fabs(a.b - b.b) < least_vis_diff ;
	}
/*...e*/

static RGBVEC whitted(
	SAMPLE *samples[MAX_SS+1],
	int h,				/* Scaled by ss */
	int ss,				/* Sample spacing */
	SHAPE *shape,
	VECTOR eye,
	int depth,
	ISECTL **ils
	)
	{
	SAMPLE *sa = &(samples[ 0][h     ]);
	SAMPLE *sb = &(samples[ 0][h + ss]);
	SAMPLE *sc = &(samples[ss][h     ]);
	SAMPLE *sd = &(samples[ss][h + ss]);
	RGBVEC m_rgb = average_of_4(sa->rgb, sb->rgb, sc->rgb, sd->rgb);

	if ( ss == 1 )
		return m_rgb;

	/* Are all 4 corners close to the middle */

	if ( is_close(m_rgb, sa->rgb) &&
	     is_close(m_rgb, sb->rgb) &&
	     is_close(m_rgb, sc->rgb) &&
	     is_close(m_rgb, sd->rgb) )
		return m_rgb;
	else
		{
		VECTOR a = sa->ray, b = sb->ray, c = sc->ray, d = sd->ray;
		int ss2 = (ss >> 1);
		SAMPLE *sab   = &(samples[0  ][h + ss2]);
		SAMPLE *scd   = &(samples[ss ][h + ss2]);
		SAMPLE *sac   = &(samples[ss2][h      ]);
		SAMPLE *sbd   = &(samples[ss2][h + ss ]);
		SAMPLE *sabcd = &(samples[ss2][h + ss2]);

		/* Trace any that are not already traced */

		if ( ! sab->traced )
			{
			sab->ray = midpoint_of_2(a, b);
			sab->rgb = trace(shape, eye, sab->ray, depth, ils);
			sab->traced = TRUE;
			extra++;
			}

		if ( ! scd->traced )
			{
			scd->ray = midpoint_of_2(c, d);
			scd->rgb = trace(shape, eye, scd->ray, depth, ils);
			scd->traced = TRUE;
			extra++;
			}

		if ( ! sac->traced )
			{
			sac->ray = midpoint_of_2(a, c);
			sac->rgb = trace(shape, eye, sac->ray, depth, ils);
			sac->traced = TRUE;
			extra++;
			}

		if ( ! sbd->traced )
			{
			sbd->ray = midpoint_of_2(b, d);
			sbd->rgb = trace(shape, eye, sbd->ray, depth, ils);
			sbd->traced = TRUE;
			extra++;
			}

		if ( ! sabcd->traced )
			{
			sabcd->ray = midpoint_of_4(a, b, c, d);
			sabcd->rgb = trace(shape, eye, sabcd->ray, depth, ils);
			sabcd->traced = TRUE;
			extra++;
			}

		/* Now obtain average of 4 nested whitted values */

		return average_of_4(
			whitted(samples    , h    , ss2, shape, eye, depth, ils),
			whitted(samples+ss2, h    , ss2, shape, eye, depth, ils),
			whitted(samples    , h+ss2, ss2, shape, eye, depth, ils),
			whitted(samples+ss2, h+ss2, ss2, shape, eye, depth, ils)
			);
		}
	}
/*...e*/

static void whitted_render(
	char *fn,
	BITMAP *bitmap,
	ISECTL **ils,
	SHAPE *shape,
	VECTOR eye, VECTOR unit_forward, VECTOR unit_right, VECTOR unit_up,
	int hpixels, int vpixels, double hangle, double vangle,
	RAY_FUNC calc_ray,
	int depth
	)
	{
	int h, v, ssh, ssv;
	int hsamples = (MAX_SS * hpixels) + 1;
	SAMPLE *samples[MAX_SS+1], *t;

	for ( ssv = 0; ssv <= MAX_SS; ssv++ )
		{
		if ( (samples[ssv] = malloc(hsamples * sizeof(SAMPLE))) == NULL)
			fatal("out of memory");
		for ( ssh = 0; ssh < hsamples; ssh++ )
			samples[ssv][ssh].traced = FALSE;
		}

	log_init(fn);

	/* Work out the top row */
	for ( h = 0; h <= hpixels; h++ )
		{
		VECTOR ray = (*calc_ray)(h, 0, hpixels, vpixels, hangle, vangle, unit_forward, unit_right, unit_up);
		samples[0][h * MAX_SS].ray = ray;
		samples[0][h * MAX_SS].rgb = trace(shape, eye, ray, depth, ils);
		samples[0][h * MAX_SS].traced = TRUE;
		normal++;
		}

	log_done_so_far(1, vpixels+1);

	for ( v = 0; v < vpixels; v++ )
		{
		/* Work out bottom row (for this scan line of pixels) */
		for ( h = 0; h <= hpixels; h++ )
			{
			VECTOR ray = (*calc_ray)(h, v+1, hpixels, vpixels, hangle, vangle, unit_forward, unit_right, unit_up);
			samples[MAX_SS][h * MAX_SS].ray = ray;
			samples[MAX_SS][h * MAX_SS].rgb = trace(shape, eye, ray, depth, ils);
			samples[MAX_SS][h * MAX_SS].traced = TRUE;
			normal++;
			}
		for ( h = 0; h < hpixels; h++ )
			{
			RGBVEC rgb = whitted(samples, h * MAX_SS, MAX_SS, shape, eye, depth, ils);
			byte r, g, b;
			log_bright(rgb);
			make_r_g_b(rgb, &r, &g, &b);
			fio_set_pixel(bitmap, h, v, r, g, b);
			}
		t = samples[0]; samples[0] = samples[MAX_SS]; samples[MAX_SS] = t;
		for ( ssv = 1; ssv <= MAX_SS; ssv++ )
			for ( ssh = 0; ssh < hsamples; ssh++ )
				samples[ssv][ssh].traced = FALSE;
		log_done_so_far(v + 1, vpixels+1);
		}
	log_done();
	printf("%ld normal rays, %ld extra rays\n", normal, extra);

	for ( ssv = 0; ssv <= MAX_SS; ssv++ )
		free(samples[ssv]);
	}
/*...e*/

static void render(
	SHAPE *shape,				/* Root of shape tree        */
	VECTOR eye,				/* Eye position vector       */
	VECTOR forward,				/* Forward vector            */
	VECTOR up,				/* Up vector                 */
	double hangle, double vangle,		/* Veiwing angles            */
	int hpixels, int vpixels,		/* No of pixels to render    */
	int depth,				/* Depth to recurse to       */
	int render_type,			/* Rendering type            */
	char *fn				/* Output filename           */
	)
	{
	VECTOR unit_forward = unit_vector(forward);
	VECTOR unit_up      = unit_vector(up);
	VECTOR unit_right   = vector_product(unit_forward, unit_up);
	BITMAP	*bitmap;
	int	n_isectls, n_isects, i;
	ISECTL	**ils;
	RAY_FUNC calc_ray;

	if ( (bitmap = fio_create_bitmap(hpixels, vpixels)) == NULL )
		fatal("out of memory for bitmap %dx%d", hpixels, vpixels);

	preprocess_shape(shape, &n_isectls, &n_isects);

	/* Now allocate intersection lists, for use in tracing */

	printf("Require %d intersection lists, each of %d intersections\n",
		n_isectls, n_isects);
	if ( (ils = malloc(n_isectls * sizeof(ISECTL *))) == NULL )
		fatal("out of memory");
	for ( i = 0; i < n_isectls; i++ )
		if ( (ils[i] = create_isectl(n_isects)) == NULL )
			fatal("out of memory");

	switch ( render_type / 10 )
		{
		case 0:
			calc_ray = calc_normal_ray;
			break;
		case 1:
			calc_ray = calc_escher_ray;
			break;
		default:
			fatal("unknown projection");
		}

	switch ( render_type % 10 )
		{
		case 0:
			basic_render(fn, bitmap, ils, shape, eye, unit_forward, unit_right, unit_up, hpixels, vpixels, hangle, vangle, calc_ray, depth);
			break;
		case 1:
			whitted_render(fn, bitmap, ils, shape, eye, unit_forward, unit_right, unit_up, hpixels, vpixels, hangle, vangle, calc_ray, depth);
			break;
		default:
			fatal("unknown render type");
		}

	for ( i = 0; i < n_isectls; i++ )
		destroy_isectl(ils[i]);
	free(ils);

	if ( !fio_write_bitmap(bitmap, fn) )
		fatal("unable to write %s", fn);

	fio_destroy_bitmap(bitmap);
	}
/*...e*/
/*...sread_data_file:0:*/
/*...slexical stuff:0:*/
#define	S_EOF		0
#define	S_ID		1
#define	S_VALUE		2
#define	S_STRING	3
#define	S_COMMA		4
#define	S_LPAR		5
#define	S_RPAR		6
#define	S_RAD		7
#define	S_VECTOR	8
#define	S_RGBVEC	9
#define	S_COL_CONST	10
#define	S_COL_NO_MOVE	11
#define	S_COL_INTERP0	12
#define	S_COL_INTERP1	13
#define	S_COL_INTERP2	14
#define	S_COL_FIELD2D	15
#define	S_COL_FIELD3D	16
#define	S_COL_REMAP	17
#define	S_COL_CYLPOLAR	18
#define	S_COL_SPHPOLAR	19
#define	S_COL_MATRIX2D	20
#define	S_COL_MATRIX3D	21
#define	S_SURF		22
#define	S_RESURF	23
#define	S_PLANE		24
#define	S_X_LT		25
#define	S_X_GT		26
#define	S_Y_LT		27
#define	S_Y_GT		28
#define	S_Z_LT		29
#define	S_Z_GT		30
#define	S_BIPLANE	31
#define	S_X_IN		32
#define	S_Y_IN		33
#define	S_Z_IN		34
#define	S_QUAD		35
#define	S_ELLIPSOID	36
#define	S_SPHERE	37
#define	S_X_ELL_CYL	38
#define	S_Y_ELL_CYL	39
#define	S_Z_ELL_CYL	40
#define	S_X_CYL		41
#define	S_Y_CYL		42
#define	S_Z_CYL		43
#define	S_X_ELL_CONE	44
#define	S_Y_ELL_CONE	45
#define	S_Z_ELL_CONE	46
#define	S_X_CONE	47
#define	S_Y_CONE	48
#define	S_Z_CONE	49
#define	S_TRANS		50
#define	S_TRANS_X	51
#define	S_TRANS_Y	52
#define	S_TRANS_Z	53
#define	S_SCALE		54
#define	S_SCALE_X	55
#define	S_SCALE_Y	56
#define	S_SCALE_Z	57
#define	S_ROTATE_X	58
#define	S_ROTATE_Y	59
#define	S_ROTATE_Z	60
#define	S_UNION		61
#define	S_ISECT		62
#define	S_DIFF		63
#define	S_SDIFF		64
#define	S_EXTENT	65
#define	S_SET_VALUE	66
#define	S_SET_VECTOR	67
#define	S_SET_RGBVEC	68
#define	S_SET_COL	69
#define	S_SET_SURF	70
#define	S_SET_SHAPE	71
#define	S_SET_BKGND	72
#define	S_SET_AMBIENT	73
#define	S_SET_ATTEN	74
#define	S_ADD_LIGHT	75
#define	S_RENDER	76
#define	S_VISDIFF	77
#define	S_INCLUDE	78

typedef struct { char *id_name; int id; } RESERVED_WORD;

static RESERVED_WORD my_reserved_words[] =
	{
	"rad",		S_RAD,
	"xyz",		S_VECTOR,
	"rgb",		S_RGBVEC,
	"col",		S_COL_CONST,
	"col_nomove",	S_COL_NO_MOVE,
	"col_interp0",	S_COL_INTERP0,
	"col_interp1",	S_COL_INTERP1,
	"col_interp2",	S_COL_INTERP2,
	"col_field2d",	S_COL_FIELD2D,
	"col_field3d",	S_COL_FIELD3D,
	"col_remap",	S_COL_REMAP,
	"col_cyl",	S_COL_CYLPOLAR,
	"col_sph",	S_COL_SPHPOLAR,
	"col_mat2d",	S_COL_MATRIX2D,
	"col_mat3d",	S_COL_MATRIX3D,
	"surf",		S_SURF,
	"resurf",	S_RESURF,
	"plane",	S_PLANE,
	"x_lt",		S_X_LT,
	"x_gt",		S_X_GT,
	"y_lt",		S_Y_LT,
	"y_gt",		S_Y_GT,
	"z_lt",		S_Z_LT,
	"z_gt",		S_Z_GT,
	"biplane",	S_BIPLANE,
	"x_in",		S_X_IN,
	"y_in",		S_Y_IN,
	"z_in",		S_Z_IN,
	"quad",		S_QUAD,
	"ellipsoid",	S_ELLIPSOID,
	"sphere",	S_SPHERE,
	"x_ell_cyl", 	S_X_ELL_CYL,
	"y_ell_cyl", 	S_Y_ELL_CYL,
	"z_ell_cyl", 	S_Z_ELL_CYL,
	"x_cyl", 	S_X_CYL,
	"y_cyl", 	S_Y_CYL,
	"z_cyl", 	S_Z_CYL,
	"x_ell_cone", 	S_X_ELL_CONE,
	"y_ell_cone", 	S_Y_ELL_CONE,
	"z_ell_cone", 	S_Z_ELL_CONE,
	"x_cone", 	S_X_CONE,
	"y_cone", 	S_Y_CONE,
	"z_cone", 	S_Z_CONE,
	"trans",	S_TRANS,
	"trans_x",	S_TRANS_X,
	"trans_y",	S_TRANS_Y,
	"trans_z",	S_TRANS_Z,
	"scale",	S_SCALE,
	"scale_x",	S_SCALE_X,
	"scale_y",	S_SCALE_Y,
	"scale_z",	S_SCALE_Z,
	"rot_x",	S_ROTATE_X,
	"rot_y",	S_ROTATE_Y,
	"rot_z",	S_ROTATE_Z,
	"union",	S_UNION,
	"isect",	S_ISECT,
	"diff",		S_DIFF,
	"sdiff",	S_SDIFF,
	"extent",	S_EXTENT,
	"set_value",	S_SET_VALUE,
	"set_xyz",	S_SET_VECTOR,
	"set_rgb",	S_SET_RGBVEC,
	"set_col",	S_SET_COL,
	"set_surf",	S_SET_SURF,
	"set_shape",	S_SET_SHAPE,
	"set_background",S_SET_BKGND,
	"set_ambient",	S_SET_AMBIENT,
	"set_attenuation",S_SET_ATTEN,
	"add_light",	S_ADD_LIGHT,
	"render",	S_RENDER,
	"visdiff",	S_VISDIFF,
	"include",	S_INCLUDE,
	};

#define	N_RESERVED_WORDS (sizeof(my_reserved_words)/sizeof(my_reserved_words[0]))

typedef struct
	{
	FILE *fp;
	char fn[500+1];
	unsigned long line_num;
	int chr;
	char str[100+1];
	char id_name[100+1];
	double id_value;
	} F;

/*...sopen_stream:0:*/
static F *open_stream(FILE *fp, char *fn)
	{
	F *f = (F *) memcheck(malloc(sizeof(F)));

	f->fp = fp;
	f->line_num = 1UL;
	strcpy(f->fn, fn);
	f->chr = getc(f->fp);
	return f;
	}
/*...e*/
/*...sclose_stream:0:*/
static FILE *close_stream(F *f)
	{
	FILE	*fp = f->fp;

	free(f);
	return fp;
	}
/*...e*/
/*...sreaderr:0:*/
static void readerr(F *f, const char *fmt, ...)
	{
	va_list	vars;
	char	s[256+1];

	va_start(vars, fmt);
	vsprintf(s, fmt, vars);
	va_end(vars);
	fprintf(stderr, "%s(%lu): %s\n", f->fn, f->line_num, s);
	exit(1);
	}
/*...e*/
/*...sgetsym:0:*/
static int getsym(F *f)
	{
	int	i = 0;

	for ( ;; )
		{
		while ( f->chr != EOF && isspace(f->chr) )
			{
			if ( f->chr == '\n' )
				f->line_num++;
			f->chr = getc(f->fp);
			}

		if ( f->chr == EOF )
			return S_EOF;

		if ( f->chr != ';' )
			break;

		while ( (f->chr = getc(f->fp)) != EOF && f->chr != '\n' )
			if ( f->chr == '\n' )
				f->line_num++;
		}

	if ( f->chr == ',' )
		{
		f->chr = getc(f->fp);
		return S_COMMA;
		}

	if ( f->chr == '(' )
		{
		f->chr = getc(f->fp);
		return S_LPAR;
		}

	if ( f->chr == ')' )
		{
		f->chr = getc(f->fp);
		return S_RPAR;
		}

	if ( isdigit(f->chr) || f->chr == '+' || f->chr == '-' || f->chr == '.' )
		{
		char	num[50+1];

		do
			{
			num[i++] = (char) f->chr;
			f->chr = getc(f->fp);
			}
		while ( isdigit(f->chr) || f->chr == '+' || f->chr == '-' ||
			f->chr == '.' || f->chr == 'e' || f->chr == 'E' );
		num[i] = (char) '\0';

		sscanf(num, "%lf", &(f->id_value));

		return S_VALUE;
		}

	if ( f->chr == '"' )
		{
		int	j = 0;

		while ( (f->chr = getc(f->fp)) != EOF && f->chr != '"' )
			f->str[j++] = (char) f->chr;
		f->str[j] = '\0';
		if ( f->chr == '"' )
			f->chr = getc(f->fp);
		return S_STRING;
		}

	if ( !(isalnum(f->chr) || f->chr != '_') )
		readerr(f, "character 0x%02x not expected", f->chr);

	while ( f->chr != EOF && (isalnum(f->chr) || f->chr == '_') )
		{
		f->id_name[i++] = (char) f->chr;
		f->chr = getc(f->fp);
		}
	f->id_name[i] = '\0';

	for ( i = 0; i < N_RESERVED_WORDS; i++ )
		if ( !strcmp(f->id_name, my_reserved_words[i].id_name) )
			return my_reserved_words[i].id;

	return S_ID;
	}
/*...e*/
/*...sskip:0:*/
static void skip(F *f, int symbol, char *symbol_name)
	{
	if ( getsym(f) != symbol )
		readerr(f, "expected %s", symbol_name);
	}
/*...e*/
/*...e*/
/*...suser variables:0:*/
typedef byte VTYPE;
#define	VTYPE_VALUE	((VTYPE) 0)
#define	VTYPE_VECTOR	((VTYPE) 1)
#define	VTYPE_RGBVEC	((VTYPE) 2)
#define	VTYPE_COL	((VTYPE) 3)
#define	VTYPE_SURF	((VTYPE) 4)
#define	VTYPE_SHAPE	((VTYPE) 5)

static char *vtype_names[] =
	{
	"value",
	"xyz_vector",
	"rgb_vector",
	"colour",
	"surface",
	"shape",
	};

typedef struct
	{
	char	*name;
	VTYPE	vtype;
	union
		{
		double	value;
		VECTOR	vector;
		RGBVEC	rgb;
		COL	*col;
		SURF	*surf;
		SHAPE	*shape;
		} u;
	} VAR;

#define	N_VARS		500
static VAR vars[N_VARS];
static int n_vars = 0;

/*...slookup_var:0:*/
static int lookup_var(char *name)
	{
	int	i;

	for ( i = 0; i < n_vars; i++ )
		if ( !strcmp(name, vars[i].name) )
			return i;
	return -1;
	}
/*...e*/
/*...slookup_defined_var:0:*/
static int lookup_defined_var(F *f, char *name)
	{
	int	i;

	if ( (i = lookup_var(name)) == -1 )
		readerr(f, "undefined variable %s", name);
	return i;
	}
/*...e*/
/*...slookup_defined_var_vtype:0:*/
static int lookup_defined_var_vtype(F *f, char *name, VTYPE vtype)
	{
	int	i = lookup_defined_var(f, name);

	if ( vars[i].vtype != vtype )
		readerr(f, "expected %s variable", vtype_names[vtype]);
	return i;
	}
/*...e*/

/*...sadd_value_var:0:*/
static void add_value_var(F *f, char *name, double value)
	{
	if ( n_vars == N_VARS )
		readerr(f, "too many variables");

	vars[n_vars  ].name    = strsave(name);
	vars[n_vars  ].vtype   = VTYPE_VALUE;
	vars[n_vars++].u.value = value;
	}
/*...e*/
/*...sadd_vector_var:0:*/
static void add_vector_var(F *f, char *name, VECTOR vector)
	{
	if ( n_vars == N_VARS )
		readerr(f, "too many variables");

	vars[n_vars  ].name     = strsave(name);
	vars[n_vars  ].vtype    = VTYPE_VECTOR;
	vars[n_vars++].u.vector = vector;
	}
/*...e*/
/*...sadd_rgb_var:0:*/
static void add_rgb_var(F *f, char *name, RGBVEC rgb)
	{
	if ( n_vars == N_VARS )
		readerr(f, "too many variables");

	vars[n_vars  ].name  = strsave(name);
	vars[n_vars  ].vtype = VTYPE_RGBVEC;
	vars[n_vars++].u.rgb = rgb;
	}
/*...e*/
/*...sadd_col_var:0:*/
static void add_col_var(F *f, char *name, COL *col)
	{
	if ( n_vars == N_VARS )
		readerr(f, "too many variables");

	vars[n_vars  ].name  = strsave(name);
	vars[n_vars  ].vtype = VTYPE_COL;
	vars[n_vars++].u.col = col;
	}
/*...e*/
/*...sadd_surf_var:0:*/
static void add_surf_var(F *f, char *name, SURF *surf)
	{
	if ( n_vars == N_VARS )
		readerr(f, "too many variables");

	vars[n_vars  ].name   = strsave(name);
	vars[n_vars  ].vtype  = VTYPE_SURF;
	vars[n_vars++].u.surf = surf;
	}
/*...e*/
/*...sadd_shape_var:0:*/
static void add_shape_var(F *f, char *name, SHAPE *shape)
	{
	if ( n_vars == N_VARS )
		readerr(f, "too many variables");

	vars[n_vars  ].name    = strsave(name);
	vars[n_vars  ].vtype   = VTYPE_SHAPE;
	vars[n_vars++].u.shape = shape;
	}
/*...e*/
/*...e*/

static void read_data_file(char *fn);

/*...sget_file:0:*/
/*...sget_new_id:0:*/
static void get_new_id(F *f, char *name)
	{
	if ( getsym(f) != S_ID )
		readerr(f, "expected new identifier");

	strncpy(name, f->id_name, 30);

	if ( lookup_var(name) != -1 )
		readerr(f, "redefinition of identifier %s", name);
	}
/*...e*/
/*...sget_value:0:*/
static double get_value(F *f)
	{
	switch ( getsym(f) )
		{
		case S_VALUE:
			return f->id_value;
		case S_ID:
			return vars[lookup_defined_var_vtype(f, f->id_name, VTYPE_VALUE)].u.value;
		case S_RAD:
			{
			double	value;

			skip(f, S_LPAR, "(");
			value = get_value(f);
			skip(f, S_RPAR, ")");
			return value * PI / 180.0;
			}
		default:
			readerr(f, "number or numeric variable expected");
		}
	return 0.0; /* Keep fussy C compiler happy */
	}
/*...e*/
/*...sget_vector:0:*/
static VECTOR get_vector(F *f)
	{
	static VECTOR dummy_vector = { 0.0, 0.0, 0.0 };

	switch ( getsym(f) )
		{
/*...sS_VECTOR:16:*/
case S_VECTOR:
	{
	VECTOR	vector;

	skip(f, S_LPAR , "("); vector.x = get_value(f);
	skip(f, S_COMMA, ","); vector.y = get_value(f);
	skip(f, S_COMMA, ","); vector.z = get_value(f);
	skip(f, S_RPAR , ")");
	return vector;
	}
/*...e*/
/*...sS_TRANS:16:*/
case S_TRANS:
	{
	VECTOR v, t;

	skip(f, S_LPAR , "("); v = get_vector(f);
	skip(f, S_COMMA, ","); t = get_vector(f);
	skip(f, S_RPAR , ")");
	return vector_sum(v, t);
	}
/*...e*/
/*...sS_TRANS_X:16:*/
case S_TRANS_X:
	{
	VECTOR v;
	double t;

	skip(f, S_LPAR , "("); v = get_vector(f);
	skip(f, S_COMMA, ","); t = get_value(f);
	skip(f, S_RPAR , ")");
	v.x += t;
	return v;
	}
/*...e*/
/*...sS_TRANS_Y:16:*/
case S_TRANS_Y:
	{
	VECTOR v;
	double t;

	skip(f, S_LPAR , "("); v = get_vector(f);
	skip(f, S_COMMA, ","); t = get_value(f);
	skip(f, S_RPAR , ")");
	v.y += t;
	return v;
	}
/*...e*/
/*...sS_TRANS_Z:16:*/
case S_TRANS_Z:
	{
	VECTOR v;
	double t;

	skip(f, S_LPAR , "("); v = get_vector(f);
	skip(f, S_COMMA, ","); t = get_value(f);
	skip(f, S_RPAR , ")");
	v.z += t;
	return v;
	}
/*...e*/
/*...sS_SCALE:16:*/
case S_SCALE:
	{
	VECTOR v, factor;

	skip(f, S_LPAR , "("); v      = get_vector(f);
	skip(f, S_COMMA, ","); factor = get_vector(f);
	skip(f, S_RPAR , ")");
	v.x *= factor.x;
	v.y *= factor.y;
	v.z *= factor.z;
	return v;
	}
/*...e*/
/*...sS_SCALE_X:16:*/
case S_SCALE_X:
	{
	VECTOR v;
	double factor;

	skip(f, S_LPAR , "("); v      = get_vector(f);
	skip(f, S_COMMA, ","); factor = get_value(f);
	skip(f, S_RPAR , ")");
	v.x *= factor;
	return v;
	}
/*...e*/
/*...sS_SCALE_Y:16:*/
case S_SCALE_Y:
	{
	VECTOR v;
	double factor;

	skip(f, S_LPAR , "("); v      = get_vector(f);
	skip(f, S_COMMA, ","); factor = get_value(f);
	skip(f, S_RPAR , ")");
	v.y *= factor;
	return v;
	}
/*...e*/
/*...sS_SCALE_Z:16:*/
case S_SCALE_Z:
	{
	VECTOR v;
	double factor;

	skip(f, S_LPAR , "("); v      = get_vector(f);
	skip(f, S_COMMA, ","); factor = get_value(f);
	skip(f, S_RPAR , ")");
	v.z *= factor;
	return v;
	}
/*...e*/
/*...sS_ROTATE_X:16:*/
case S_ROTATE_X:
	{
	VECTOR v;
	double angle;

	skip(f, S_LPAR , "("); v = get_vector(f);
	skip(f, S_COMMA, ","); angle = get_value(f);
	skip(f, S_RPAR , ")");
	return rot_x_vector(v, angle);
	}
/*...e*/
/*...sS_ROTATE_Y:16:*/
case S_ROTATE_Y:
	{
	VECTOR v;
	double angle;

	skip(f, S_LPAR , "("); v = get_vector(f);
	skip(f, S_COMMA, ","); angle = get_value(f);
	skip(f, S_RPAR , ")");
	return rot_y_vector(v, angle);
	}
/*...e*/
/*...sS_ROTATE_Z:16:*/
case S_ROTATE_Z:
	{
	VECTOR v;
	double angle;

	skip(f, S_LPAR , "("); v = get_vector(f);
	skip(f, S_COMMA, ","); angle = get_value(f);
	skip(f, S_RPAR , ")");
	return rot_z_vector(v, angle);
	}
/*...e*/
/*...sS_ID:16:*/
case S_ID:
	return vars[lookup_defined_var_vtype(f, f->id_name, VTYPE_VECTOR)].u.vector;
/*...e*/
/*...sdefault:16:*/
default:
	readerr(f, "vector expected");
/*...e*/
		}
	return dummy_vector; /* Keep fussy C compiler happy */
	}
/*...e*/
/*...sget_rgb:0:*/
static RGBVEC get_rgb(F *f)
	{
	static RGBVEC dummy_rgbvec = { 0.0, 0.0, 0.0 };

	switch ( getsym(f) )
		{
		case S_RGBVEC:
			{
			RGBVEC	rgb;

			skip(f, S_LPAR , "("); rgb.r = get_value(f);
			skip(f, S_COMMA, ","); rgb.g = get_value(f);
			skip(f, S_COMMA, ","); rgb.b = get_value(f);
			skip(f, S_RPAR , ")");
			return rgb;
			}
		case S_ID:
			return vars[lookup_defined_var_vtype(f, f->id_name, VTYPE_RGBVEC)].u.rgb;
		default:
			readerr(f, "rgb colour vector expected");
		}

	return dummy_rgbvec; /* Keep fussy C compiler happy */
	}
/*...e*/
/*...sget_col:0:*/
/*
This returns a colour.
If it comes from a user variable, then a copy of the variable is returned.
If not, then the new datastructure is returned.
*/

static COL *get_col(F *f)
	{
	switch ( getsym(f) )
		{
/*...sS_COL_CONST:16:*/
case S_COL_CONST:
	{
	RGBVEC rgbvec;

	skip(f, S_LPAR , "("); rgbvec = get_rgb(f);
	skip(f, S_RPAR , ")");
	return memcheck(create_const_col(rgbvec));
	}
/*...e*/
/*...sS_COL_NO_MOVE:16:*/
case S_COL_NO_MOVE:
	{
	COL *col;

	skip(f, S_LPAR , "("); col = get_col(f);
	skip(f, S_RPAR , ")");

	return memcheck(create_no_move_col(col));
	}
/*...e*/
/*...sS_COL_INTERP0:16:*/
case S_COL_INTERP0:
	{
	COL *col;

	skip(f, S_LPAR , "("); col = get_col(f);
	skip(f, S_RPAR , ")");

	return memcheck(create_interp0_col(col));
	}
/*...e*/
/*...sS_COL_INTERP1:16:*/
case S_COL_INTERP1:
	{
	COL *col;

	skip(f, S_LPAR , "("); col = get_col(f);
	skip(f, S_RPAR , ")");

	return memcheck(create_interp1_col(col));
	}
/*...e*/
/*...sS_COL_INTERP2:16:*/
case S_COL_INTERP2:
	{
	COL *col;

	skip(f, S_LPAR , "("); col = get_col(f);
	skip(f, S_RPAR , ")");

	return memcheck(create_interp2_col(col));
	}
/*...e*/
/*...sS_COL_FIELD2D:16:*/
case S_COL_FIELD2D:
	{
	double bx, by;
	BITMAP *bitmap;

	skip(f, S_LPAR , "("); bx = get_value(f);
	skip(f, S_COMMA, ","); by = get_value(f);
	skip(f, S_COMMA, ",");

	if ( getsym(f) != S_STRING )
		readerr(f, "expected input bitmap filename");

	if ( (bitmap = fio_read_bitmap(f->str)) == NULL )
		readerr(f, "unable to read bitmap %s", f->str);

	skip(f, S_RPAR , ")");
	return memcheck(create_2d_field_col(bx, by, bitmap));
	}
/*...e*/
/*...sS_COL_FIELD3D:16:*/
case S_COL_FIELD3D:
	{
	double bx, by, bz;
	TEX *tex;

	skip(f, S_LPAR , "("); bx = get_value(f);
	skip(f, S_COMMA, ","); by = get_value(f);
	skip(f, S_COMMA, ","); bz = get_value(f);
	skip(f, S_COMMA, ",");

	if ( getsym(f) != S_STRING )
		readerr(f, "expected input 3d texture-map filename");

	if ( (tex = read_tex(f->str)) == NULL )
		readerr(f, "unable to read 3d texture-map %s", f->str);

	skip(f, S_RPAR , ")");
	return memcheck(create_3d_field_col(bx, by, bz, tex));
	}
/*...e*/
/*...sS_COL_REMAP:16:*/
case S_COL_REMAP:
	{
	VECTOR base, v0, v1, v2;
	COL *col;

	skip(f, S_LPAR , "("); base = get_vector(f);
	skip(f, S_COMMA, ","); v0   = get_vector(f);
	skip(f, S_COMMA, ","); v1   = get_vector(f);
	skip(f, S_COMMA, ","); v2   = get_vector(f);
	skip(f, S_COMMA, ","); col  = get_col(f);
	skip(f, S_RPAR , ")");
	return memcheck(create_remap_col(base, v0, v1, v2, col));
	}
/*...e*/
/*...sS_COL_CYLPOLAR:16:*/
case S_COL_CYLPOLAR:
	{
	double lond, rd, hd;
	COL *col;

	skip(f, S_LPAR , "("); lond = get_value(f);
	skip(f, S_COMMA, ","); rd   = get_value(f);
	skip(f, S_COMMA, ","); hd   = get_value(f);
	skip(f, S_COMMA, ","); col  = get_col(f);
	skip(f, S_RPAR , ")");

	return memcheck(create_cyl_polar_col(lond, rd, hd, col));
	}
/*...e*/
/*...sS_COL_SPHPOLAR:16:*/
case S_COL_SPHPOLAR:
	{
	double lond, latd, rd;
	COL *col;

	skip(f, S_LPAR , "("); lond = get_value(f);
	skip(f, S_COMMA, ","); latd = get_value(f);
	skip(f, S_COMMA, ","); rd   = get_value(f);
	skip(f, S_COMMA, ","); col  = get_col(f);
	skip(f, S_RPAR , ")");

	return memcheck(create_sph_polar_col(lond, latd, rd, col));
	}
/*...e*/
/*...sS_COL_MATRIX2D:16:*/
case S_COL_MATRIX2D:
	{
	double mat[2][2];
	int i, j;
	COL *col;

	skip(f, S_LPAR , "(");
	for ( j = 0; j < 2; j++ )
		for ( i = 0; i < 2; i++ )
			{
			mat[j][i] = get_value(f);
			skip(f, S_COMMA, ",");
			}
	col = get_col(f);
	skip(f, S_RPAR , ")");

	return memcheck(create_2d_matrix_col(mat, col));
	}
/*...e*/
/*...sS_COL_MATRIX3D:16:*/
case S_COL_MATRIX3D:
	{
	double mat[3][3];
	int i, j;
	COL *col;

	skip(f, S_LPAR , "(");
	for ( j = 0; j < 3; j++ )
		for ( i = 0; i < 3; i++ )
			{
			mat[j][i] = get_value(f);
			skip(f, S_COMMA, ",");
			}
	col = get_col(f);
	skip(f, S_RPAR , ")");

	return memcheck(create_3d_matrix_col(mat, col));
	}
/*...e*/
/*...sS_ID:16:*/
case S_ID:
	{
	COL *col = vars[lookup_defined_var_vtype(f, f->id_name, VTYPE_COL)].u.col;

	return memcheck(copy_col(col));
	}
/*...e*/
/*...sdefault:16:*/
default:
	readerr(f, "colour definition or colour variable expected");
/*...e*/
		}
	return NULL; /* Keep fussy C compiler happy */
	}
/*...e*/
/*...sget_surf:0:*/
/*
This returns a surface.
If it comes from a user variable, then a copy of the variable is returned.
If not, then the new datastructure is returned.
*/

static SURF *get_surf(F *f)
	{
	switch ( getsym(f) )
		{
		case S_SURF:
			{
			double	ka, kd, ks, kt;
			COL	*od, *os;
			double	phong, rinx;

			skip(f, S_LPAR , "("); ka    = get_value(f);
			skip(f, S_COMMA, ","); kd    = get_value(f);
			skip(f, S_COMMA, ","); ks    = get_value(f);
			skip(f, S_COMMA, ","); kt    = get_value(f);
			skip(f, S_COMMA, ","); od    = get_col(f);
			skip(f, S_COMMA, ","); os    = get_col(f);
			skip(f, S_COMMA, ","); phong = get_value(f);
			skip(f, S_COMMA, ","); rinx  = get_value(f);
			skip(f, S_RPAR , ")");

			return memcheck(create_surf(ka, kd, ks, kt, od, os, phong, rinx));
			}
		case S_ID:
			{
			SURF *surf = vars[lookup_defined_var_vtype(f, f->id_name, VTYPE_SURF)].u.surf;

			return memcheck(copy_surf(surf));
			}
		default:
			readerr(f, "surface definition or surface variable expected");
		}
	return NULL; /* Keep fussy C compiler happy */
	}
/*...e*/
/*...sget_opt_surf:0:*/
/*
Handle both the case that ", surface)" is supplied and that surface is
to be returned, and also the case ")" is supplied, and a default is to be used.
*/

static SURF *get_opt_surf(F *f)
	{
	switch ( getsym(f) )
		{
		case S_COMMA:
			{
			SURF *surf = get_surf(f);
			skip(f, S_RPAR, ")");
			return surf;
			}
		case S_RPAR:
			{
			static RGBVEC rgbvec_g = { 0.5, 0.5, 0.5 };
			COL *col_g = memcheck(create_const_col(rgbvec_g));
			return memcheck(create_surf(0.8,0.8,0.8,0.0,col_g,col_g,100.0,1.0));
			}
		default:
			readerr(f, "surface definition, surface variable or ) expected");
		}
	return NULL; /* Keep fussy C compiler happy */
	}
/*...e*/
/*...sget_shape:0:*/
/*
This returns a shape.
If it comes from a user variable, then a copy of the variable is returned.
If not, then the new datastructure is returned.
*/

static SHAPE *get_shape(F *f)
	{
	switch ( getsym(f) )
		{
/*...sS_PLANE:16:*/
case S_PLANE:
	{
	double	a, b, c, d;
	SURF	*surf;

	skip(f, S_LPAR , "("); a = get_value(f);
	skip(f, S_COMMA, ","); b = get_value(f);
	skip(f, S_COMMA, ","); c = get_value(f);
	skip(f, S_COMMA, ","); d = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_plane_shape(memcheck(create_plane(a, b, c, d)), surf));
	}
/*...e*/
/*...sS_X_LT:16:*/
case S_X_LT:
	{
	double	x;
	SURF	*surf;

	skip(f, S_LPAR , "("); x = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_plane_shape(memcheck(create_x_lt_plane(x)), surf));
	}
/*...e*/
/*...sS_X_GT:16:*/
case S_X_GT:
	{
	double	x;
	SURF	*surf;

	skip(f, S_LPAR , "("); x = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_plane_shape(memcheck(create_x_gt_plane(x)), surf));
	}
/*...e*/
/*...sS_Y_LT:16:*/
case S_Y_LT:
	{
	double	y;
	SURF	*surf;

	skip(f, S_LPAR , "("); y = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_plane_shape(memcheck(create_y_lt_plane(y)), surf));
	}
/*...e*/
/*...sS_Y_GT:16:*/
case S_Y_GT:
	{
	double	y;
	SURF	*surf;

	skip(f, S_LPAR , "("); y = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_plane_shape(memcheck(create_y_gt_plane(y)), surf));
	}
/*...e*/
/*...sS_Z_LT:16:*/
case S_Z_LT:
	{
	double	z;
	SURF	*surf;

	skip(f, S_LPAR , "("); z = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_plane_shape(memcheck(create_z_lt_plane(z)), surf));
	}
/*...e*/
/*...sS_Z_GT:16:*/
case S_Z_GT:
	{
	double	z;
	SURF	*surf;

	skip(f, S_LPAR , "("); z = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_plane_shape(memcheck(create_z_gt_plane(z)), surf));
	}
/*...e*/
/*...sS_BIPLANE:16:*/
case S_BIPLANE:
	{
	double	a, b, c, d1, d2;
	SURF	*surf;

	skip(f, S_LPAR , "("); a  = get_value(f);
	skip(f, S_COMMA, ","); b  = get_value(f);
	skip(f, S_COMMA, ","); c  = get_value(f);
	skip(f, S_COMMA, ","); d1 = get_value(f);
	skip(f, S_COMMA, ","); d2 = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_biplane_shape(memcheck(create_biplane(a, b, c, d1, d2)), surf));
	}
/*...e*/
/*...sS_X_IN:16:*/
case S_X_IN:
	{
	double	x1, x2;
	SURF	*surf;

	skip(f, S_LPAR , "("); x1 = get_value(f);
	skip(f, S_COMMA, ","); x2 = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_biplane_shape(memcheck(create_x_in_biplane(x1, x2)), surf));
	}
/*...e*/
/*...sS_Y_IN:16:*/
case S_Y_IN:
	{
	double	yy1, yy2;
	SURF	*surf;

	skip(f, S_LPAR , "("); yy1 = get_value(f);
	skip(f, S_COMMA, ","); yy2 = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_biplane_shape(memcheck(create_y_in_biplane(yy1, yy2)), surf));
	}
/*...e*/
/*...sS_Z_IN:16:*/
case S_Z_IN:
	{
	double	z1, z2;
	SURF	*surf;

	skip(f, S_LPAR , "("); z1 = get_value(f);
	skip(f, S_COMMA, ","); z2 = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_biplane_shape(memcheck(create_z_in_biplane(z1, z2)), surf));
	}
/*...e*/
/*...sS_SPHERE:16:*/
case S_SPHERE:
	{
	double	r;
	SURF	*surf;

	skip(f, S_LPAR , "("); r = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_sphere_shape(memcheck(create_sphere(r)), surf));
	}
/*...e*/
/*...sS_QUAD:16:*/
case S_QUAD:
	{
	double	a, b, c, d, e, ff, g, h, i, j;
	SURF	*surf;

	skip(f, S_LPAR , "("); a  = get_value(f);
	skip(f, S_COMMA, ","); b  = get_value(f);
	skip(f, S_COMMA, ","); c  = get_value(f);
	skip(f, S_COMMA, ","); d  = get_value(f);
	skip(f, S_COMMA, ","); e  = get_value(f);
	skip(f, S_COMMA, ","); ff = get_value(f);
	skip(f, S_COMMA, ","); g  = get_value(f);
	skip(f, S_COMMA, ","); h  = get_value(f);
	skip(f, S_COMMA, ","); i  = get_value(f);
	skip(f, S_COMMA, ","); j  = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_quad_shape(memcheck(create_quad(a, b, c, d, e, ff, g, h, i, j)), surf));
	}
/*...e*/
/*...sS_ELLIPSOID:16:*/
case S_ELLIPSOID:
	{
	double	rx, ry, rz;
	SURF	*surf;

	skip(f, S_LPAR , "("); rx = get_value(f);
	skip(f, S_COMMA, ","); ry = get_value(f);
	skip(f, S_COMMA, ","); rz = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_quad_shape(memcheck(create_ellipsoid(rx, ry, rz)), surf));
	}
/*...e*/
/*...sS_X_ELL_CYL:16:*/
case S_X_ELL_CYL:
	{
	double	ry, rz;
	SURF	*surf;

	skip(f, S_LPAR , "("); ry = get_value(f);
	skip(f, S_COMMA, ","); rz = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_quad_shape(memcheck(create_x_ell_cyl(ry, rz)), surf));
	}
/*...e*/
/*...sS_Y_ELL_CYL:16:*/
case S_Y_ELL_CYL:
	{
	double	rx, rz;
	SURF	*surf;

	skip(f, S_LPAR , "("); rx = get_value(f);
	skip(f, S_COMMA, ","); rz = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_quad_shape(memcheck(create_y_ell_cyl(rx, rz)), surf));
	}
/*...e*/
/*...sS_Z_ELL_CYL:16:*/
case S_Z_ELL_CYL:
	{
	double	rx, ry;
	SURF	*surf;

	skip(f, S_LPAR , "("); rx = get_value(f);
	skip(f, S_COMMA, ","); ry = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_quad_shape(memcheck(create_z_ell_cyl(rx, ry)), surf));
	}
/*...e*/
/*...sS_X_CYL:16:*/
case S_X_CYL:
	{
	double	r;
	SURF	*surf;

	skip(f, S_LPAR , "("); r = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_quad_shape(memcheck(create_x_cyl(r)), surf));
	}
/*...e*/
/*...sS_Y_CYL:16:*/
case S_Y_CYL:
	{
	double	r;
	SURF	*surf;

	skip(f, S_LPAR , "("); r = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_quad_shape(memcheck(create_y_cyl(r)), surf));
	}
/*...e*/
/*...sS_Z_CYL:16:*/
case S_Z_CYL:
	{
	double	r;
	SURF	*surf;

	skip(f, S_LPAR , "("); r = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_quad_shape(memcheck(create_z_cyl(r)), surf));
	}
/*...e*/
/*...sS_X_ELL_CONE:16:*/
case S_X_ELL_CONE:
	{
	double	ky, kz;
	SURF	*surf;

	skip(f, S_LPAR , "("); ky = get_value(f);
	skip(f, S_COMMA, ","); kz = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_quad_shape(memcheck(create_x_ell_cone(ky, kz)), surf));
	}
/*...e*/
/*...sS_Y_ELL_CONE:16:*/
case S_Y_ELL_CONE:
	{
	double	kx, kz;
	SURF	*surf;

	skip(f, S_LPAR , "("); kx = get_value(f);
	skip(f, S_COMMA, ","); kz = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_quad_shape(memcheck(create_y_ell_cone(kx, kz)), surf));
	}
/*...e*/
/*...sS_Z_ELL_CONE:16:*/
case S_Z_ELL_CONE:
	{
	double	kx, ky;
	SURF	*surf;

	skip(f, S_LPAR , "("); kx = get_value(f);
	skip(f, S_COMMA, ","); ky = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_quad_shape(memcheck(create_z_ell_cone(kx, ky)), surf));
	}
/*...e*/
/*...sS_X_CONE:16:*/
case S_X_CONE:
	{
	double	k;
	SURF	*surf;

	skip(f, S_LPAR , "("); k = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_quad_shape(memcheck(create_x_cone(k)), surf));
	}
/*...e*/
/*...sS_Y_CONE:16:*/
case S_Y_CONE:
	{
	double	k;
	SURF	*surf;

	skip(f, S_LPAR , "("); k = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_quad_shape(memcheck(create_y_cone(k)), surf));
	}
/*...e*/
/*...sS_Z_CONE:16:*/
case S_Z_CONE:
	{
	double	k;
	SURF	*surf;

	skip(f, S_LPAR , "("); k = get_value(f);
	surf = get_opt_surf(f);
	return memcheck(create_quad_shape(memcheck(create_z_cone(k)), surf));
	}
/*...e*/
/*...sS_TRANS:16:*/
case S_TRANS:
	{
	SHAPE	*shape;
	VECTOR	t;

	skip(f, S_LPAR , "("); shape = get_shape(f);
	skip(f, S_COMMA, ","); t = get_vector(f);
	skip(f, S_RPAR , ")");
	trans(shape, t);
	return shape;
	}
/*...e*/
/*...sS_TRANS_X:16:*/
case S_TRANS_X:
	{
	SHAPE	*shape;
	double	t;

	skip(f, S_LPAR , "("); shape = get_shape(f);
	skip(f, S_COMMA, ","); t = get_value(f);
	skip(f, S_RPAR , ")");
	trans_x(shape, t);
	return shape;
	}
/*...e*/
/*...sS_TRANS_Y:16:*/
case S_TRANS_Y:
	{
	SHAPE	*shape;
	double	t;

	skip(f, S_LPAR , "("); shape = get_shape(f);
	skip(f, S_COMMA, ","); t = get_value(f);
	skip(f, S_RPAR , ")");
	trans_y(shape, t);
	return shape;
	}
/*...e*/
/*...sS_TRANS_Z:16:*/
case S_TRANS_Z:
	{
	SHAPE	*shape;
	double	t;

	skip(f, S_LPAR , "("); shape = get_shape(f);
	skip(f, S_COMMA, ","); t = get_value(f);
	skip(f, S_RPAR , ")");
	trans_z(shape, t);
	return shape;
	}
/*...e*/
/*...sS_SCALE:16:*/
case S_SCALE:
	{
	SHAPE	*shape;
	VECTOR	factor;

	skip(f, S_LPAR , "("); shape = get_shape(f);
	skip(f, S_COMMA, ","); factor = get_vector(f);
	if ( factor.x == 0.0 || factor.y == 0.0 || factor.z == 0 )
		readerr(f, "at least one component of factor is zero");
	skip(f, S_RPAR , ")");
	if ( !scale(shape, factor) )
		readerr(f, "can't scale shape");
	return shape;
	}
/*...e*/
/*...sS_SCALE_X:16:*/
case S_SCALE_X:
	{
	SHAPE	*shape;
	double	factor;

	skip(f, S_LPAR , "("); shape = get_shape(f);
	skip(f, S_COMMA, ","); factor = get_value(f);
	if ( factor == 0.0 )
		readerr(f, "factor is zero");
	skip(f, S_RPAR , ")");
	if ( !scale_x(shape, factor) )
		readerr(f, "can't scale_x shape");
	return shape;
	}
/*...e*/
/*...sS_SCALE_Y:16:*/
case S_SCALE_Y:
	{
	SHAPE	*shape;
	double	factor;

	skip(f, S_LPAR , "("); shape = get_shape(f);
	skip(f, S_COMMA, ","); factor = get_value(f);
	if ( factor == 0.0 )
		readerr(f, "factor is zero");
	skip(f, S_RPAR , ")");
	if ( !scale_y(shape, factor) )
		readerr(f, "can't scale_y shape");
	return shape;
	}
/*...e*/
/*...sS_SCALE_Z:16:*/
case S_SCALE_Z:
	{
	SHAPE	*shape;
	double	factor;

	skip(f, S_LPAR , "("); shape = get_shape(f);
	skip(f, S_COMMA, ","); factor = get_value(f);
	if ( factor == 0.0 )
		readerr(f, "factor is zero");
	skip(f, S_RPAR , ")");
	if ( !scale_z(shape, factor) )
		readerr(f, "can't scale_z shape");
	return shape;
	}
/*...e*/
/*...sS_ROTATE_X:16:*/
case S_ROTATE_X:
	{
	SHAPE	*shape;
	double	angle;

	skip(f, S_LPAR , "("); shape = get_shape(f);
	skip(f, S_COMMA, ","); angle = get_value(f);
	skip(f, S_RPAR , ")");
	rot_x(shape, angle);
	return shape;
	}
/*...e*/
/*...sS_ROTATE_Y:16:*/
case S_ROTATE_Y:
	{
	SHAPE	*shape;
	double	angle;

	skip(f, S_LPAR , "("); shape = get_shape(f);
	skip(f, S_COMMA, ","); angle = get_value(f);
	skip(f, S_RPAR , ")");
	rot_y(shape, angle);
	return shape;
	}
/*...e*/
/*...sS_ROTATE_Z:16:*/
case S_ROTATE_Z:
	{
	SHAPE	*shape;
	double	angle;

	skip(f, S_LPAR , "("); shape = get_shape(f);
	skip(f, S_COMMA, ","); angle = get_value(f);
	skip(f, S_RPAR , ")");
	rot_z(shape, angle);
	return shape;
	}
/*...e*/
/*...sS_UNION:16:*/
case S_UNION:
	{
	SHAPE	*shape;
	int	sym;

	skip(f, S_LPAR, "(");
	shape = get_shape(f);
	while ( (sym = getsym(f)) == S_COMMA )
		shape = memcheck(create_bin_shape(STYPE_UNION, shape, get_shape(f)));
	if ( sym != S_RPAR )
		readerr(f, "expected )");

	return shape;
	}
/*...e*/
/*...sS_ISECT:16:*/
case S_ISECT:
	{
	SHAPE	*shape;
	int	sym;

	skip(f, S_LPAR, "(");
	shape = get_shape(f);
	while ( (sym = getsym(f)) == S_COMMA )
		shape = memcheck(create_bin_shape(STYPE_ISECT, shape, get_shape(f)));
	if ( sym != S_RPAR )
		readerr(f, "expected )");

	return shape;
	}
/*...e*/
/*...sS_DIFF:16:*/
case S_DIFF:
	{
	SHAPE	*shape;
	int	sym;

	skip(f, S_LPAR, "(");
	shape = get_shape(f);
	while ( (sym = getsym(f)) == S_COMMA )
		shape = memcheck(create_bin_shape(STYPE_DIFF, shape, get_shape(f)));
	if ( sym != S_RPAR )
		readerr(f, "expected )");

	return shape;
	}
/*...e*/
/*...sS_SDIFF:16:*/
case S_SDIFF:
	{
	SHAPE	*shape;
	int	sym;

	skip(f, S_LPAR, "(");
	shape = get_shape(f);
	while ( (sym = getsym(f)) == S_COMMA )
		shape = memcheck(create_bin_shape(STYPE_SDIFF, shape, get_shape(f)));
	if ( sym != S_RPAR )
		readerr(f, "expected )");

	return shape;
	}
/*...e*/
/*...sS_EXTENT:16:*/
case S_EXTENT:
	{
	SHAPE	*shape;
	int	sym;

	skip(f, S_LPAR, "(");
	shape = get_shape(f);
	while ( (sym = getsym(f)) == S_COMMA )
		shape = memcheck(create_bin_shape(STYPE_EXTENT, shape, get_shape(f)));
	if ( sym != S_RPAR )
		readerr(f, "expected )");

	return shape;
	}
/*...e*/
/*...sS_RESURF:16:*/
case S_RESURF:
	{
	SHAPE	*shape;
	SURF	*surf;

	skip(f, S_LPAR , "("); shape = get_shape(f);
	surf = get_opt_surf(f);

	if ( !resurf(shape, surf) )
		readerr(f, "unable to resurface shape");
	destroy_surf(surf);

	return shape;
	}
/*...e*/
/*...sS_ID:16:*/
case S_ID:
	{
	SHAPE *shape = vars[lookup_defined_var_vtype(f, f->id_name, VTYPE_SHAPE)].u.shape;

	return memcheck(copy_shape(shape));
	}
/*...e*/
/*...sdefault:16:*/
default:
	readerr(f, "expected a shape specification");
/*...e*/
		}
	return NULL; /* Keep fussy C compiler happy */
	}
/*...e*/
/*...sget_set_value:0:*/
static void get_set_value(F *f)
	{
	char	name[30+1];

	get_new_id(f, name);
	add_value_var(f, name, get_value(f));
	}

/*...e*/
/*...sget_set_vector:0:*/
static void get_set_vector(F *f)
	{
	char	name[30+1];

	get_new_id(f, name);
	add_vector_var(f, name, get_vector(f));
	}
/*...e*/
/*...sget_set_rgb:0:*/
static void get_set_rgb(F *f)
	{
	char	name[30+1];

	get_new_id(f, name);
	add_rgb_var(f, name, get_rgb(f));
	}
/*...e*/
/*...sget_set_col:0:*/
static void get_set_col(F *f)
	{
	char	name[30+1];

	get_new_id(f, name);
	add_col_var(f, name, get_col(f));
	}
/*...e*/
/*...sget_set_surf:0:*/
static void get_set_surf(F *f)
	{
	char	name[30+1];

	get_new_id(f, name);
	add_surf_var(f, name, get_surf(f));
	}
/*...e*/
/*...sget_set_shape:0:*/
static void get_set_shape(F *f)
	{
	char	name[30+1];

	get_new_id(f, name);
	add_shape_var(f, name, get_shape(f));
	}
/*...e*/
/*...sget_set_background:0:*/
static void get_set_background(F *f)
	{
	i_background = get_rgb(f);
	}
/*...e*/
/*...sget_set_ambient:0:*/
static void get_set_ambient(F *f)
	{
	i_ambient = get_rgb(f);
	}
/*...e*/
/*...sget_set_atten:0:*/
static void get_set_atten(F *f)
	{
	af1 = get_value(f);
	af2 = get_value(f);
	}
/*...e*/
/*...sget_add_light:0:*/
static void get_add_light(F *f)
	{
	if ( n_lights == N_LIGHTS )
		readerr(f, "maximum of %d lights exceeded", N_LIGHTS);

	lights[n_lights  ].posn = get_vector(f);
	lights[n_lights++].i    = get_rgb(f);
	}
/*...e*/
/*...sget_render:0:*/
static void get_render(F *f)
	{
	SHAPE	*shape;
	VECTOR	eye, forward, up;
	double	hangle, vangle;
	int	hpixels, vpixels, depth, type;

	shape   = get_shape(f);
	eye     = get_vector(f);
	forward = get_vector(f);
	up      = get_vector(f);
	hangle  = get_value(f);
	vangle  = get_value(f);
	hpixels = (int) get_value(f);
	vpixels = (int) get_value(f);
	depth   = (int) get_value(f);
	type    = (int) get_value(f);

	if ( getsym(f) != S_STRING )
		readerr(f, "expected output filename");

	render(shape, eye, forward, up, hangle, vangle, hpixels, vpixels, depth, type, f->str);

	destroy_shape(shape);
	}
/*...e*/
/*...sget_visdiff:0:*/
static void get_visdiff(F *f)
	{
	least_vis_diff = get_value(f);
	}
/*...e*/
/*...sget_include:0:*/
static void get_include(F *f)
	{
	if ( getsym(f) != S_STRING )
		readerr(f, "expected filename");

	read_data_file(f->str);
	}
/*...e*/

static void get_file(F *f)
	{
	int	sym;

	while ( (sym = getsym(f)) != S_EOF )
		switch ( sym )
			{
			case S_SET_VALUE:	get_set_value(f);	break;
			case S_SET_VECTOR:	get_set_vector(f);	break;
			case S_SET_RGBVEC:	get_set_rgb(f);		break;
			case S_SET_COL:		get_set_col(f);		break;
			case S_SET_SURF:	get_set_surf(f);	break;
			case S_SET_SHAPE:	get_set_shape(f);	break;
			case S_SET_BKGND:	get_set_background(f);	break;
			case S_SET_AMBIENT:	get_set_ambient(f);	break;
			case S_SET_ATTEN:	get_set_atten(f);	break;
			case S_ADD_LIGHT:	get_add_light(f);	break;
			case S_RENDER:		get_render(f);		break;
			case S_VISDIFF:		get_visdiff(f);		break;
			case S_INCLUDE:		get_include(f);		break;
			default:		readerr(f, "expected assignment, set_ambient, set_attenuation, add_light, render or include statement");
			}
	}
/*...e*/

static void read_data_file(char *fn)
	{
	FILE	*fp;
	F	*f;

	if ( (fp = fopen(fn, "r")) == NULL )
		fatal("can't open %s", fn);

	f = open_stream(fp, fn);
	get_file(f);
	fp = close_stream(f);
	fclose(fp);
	}
/*...e*/

int main(int argc, char *argv[])
	{
	if ( argc != 2 )
		exit(1);

#ifdef OS2
	/* Prevent numeric exceptions from terminating this program */
	_control87(EM_UNDERFLOW|EM_DENORMAL, EM_UNDERFLOW|EM_DENORMAL);
#endif

	read_data_file(argv[1]);

	return 0;
	}
/*...e*/
