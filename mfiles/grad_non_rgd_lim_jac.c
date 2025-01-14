
/*********************************************************************************************
 *    B-spline based nonrigid transformation with transformed image and image gradients
 *
 *    Date: June 13, 2003
 *    Programmer: Jeongtae Kim
 *                (EECS:Systems, U of Michigan)
 *    Note: 1. Parts of the program for image interpolation are originated from M. Unser's group.
 *          2. Penalty function and boundary conditions are not finalized yet.
 *   
 *    Bug report: jeongtae@umich.edu, +1-734-647-8390
 *********************************************************************************************/

/*****************************************************************************
 *	System includes
 ****************************************************************************/
#include	<math.h>
#include	<stddef.h>
#include	<stdio.h>
#include	<stdlib.h>

#include        "mex.h"
#include        "matrix.h"

/*****************************************************************************
 *	Other includes
 ****************************************************************************/
#include	"interpol3d.h"
#include        "nonrigid.h"
 
#define out_val 45
#define out_grad_val 0
#define rpi 3.141592
#define kpi 0.0175



/*****************************************************************************
 *	Definition of extern procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
extern double	InterpolatedValue
				(
				        double	*Bcoeff,   /* input B-spline array of coefficients */

					double   *xgradient,   /* gradient x-value */
					double   *ygradient,   /* gradient y-value */
					double   *zgradient,   /* gradient z-value */

					long	Width,	   /* width of the image */
					long	Height,	   /* height of the image */
					long    Slice,
					double	x,	   /* x coordinate where to interpolate */
					double	y,	   /* y coordinate where to interpolate */
					double  z,
					int     ind3d,
					long	SplineDegree  /* degree of the spline model */
			
				)
{ /* begin InterpolatedValue */
	double	*p;
	double	xWeight[6], yWeight[6], zWeight[6];
	double	gxWeight[6], gyWeight[6], gzWeight[6];
	double	interpolated, grad_x, grad_y, grad_z;
	double	w, wx, w2, w2x, w2y, w4, t, t0, t1;
	long	xIndex[6], yIndex[6], zIndex[6];
	long	Width2 = 2L * Width - 2L, Height2 = 2L * Height - 2L, Slice2 = 2L * Slice - 2L;
	long	i, j, k, l;
	/* compute the interpolation indexes */
	if (SplineDegree & 1L) {
		i = (long)floor(x) - SplineDegree / 2L;
		j = (long)floor(y) - SplineDegree / 2L;
		k = (long)floor(z) - SplineDegree / 2L;
		for (l = 0L; l <= SplineDegree; l++) {
			xIndex[l] = i++;
			yIndex[l] = j++;
			zIndex[l] = k++;
		}
	}
	else {
		i = (long)floor(x + 0.5) - SplineDegree / 2L;
		j = (long)floor(y + 0.5) - SplineDegree / 2L;
		k = (long)floor(z + 0.5) - SplineDegree / 2L;
		for (l = 0L; l <= SplineDegree; l++) {
			xIndex[l] = i++;
			yIndex[l] = j++;
			zIndex[l] = k++;
		}
	}
	/* compute the interpolation weights */
	switch (SplineDegree) {
		case 2L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[1] = 3.0 / 4.0 - w * w;
			xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
			xWeight[0] = 1.0 - xWeight[1] - xWeight[2];
			/* y */
			w = y - (double)yIndex[1];
			yWeight[1] = 3.0 / 4.0 - w * w;
			yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
			yWeight[0] = 1.0 - yWeight[1] - yWeight[2];
			/* z */
			if (ind3d ==1){
			  w = z - (double)zIndex[1];
			  zWeight[1] = 3.0 / 4.0 - w * w;
			  zWeight[2] = (1.0 / 2.0) * (w - zWeight[1] + 1.0);
			  zWeight[0] = 1.0 - zWeight[1] - zWeight[2];
			}
			break;
		case 3L:
			/* x-gradient */
			w = x - (double)xIndex[1];
			gxWeight[3] = (1.0/ 2.0) * w * w ;
			gxWeight[0] = (-1.0 / 2.0) + w - gxWeight[3];
			gxWeight[2] = 3.0 * gxWeight[0] + 2.0 * (1 - w);
			gxWeight[1] = 3.0 * gxWeight[3] - 2.0 * w;
			/* printf("Came here\n");*/

			/* x */
			xWeight[3] = (1.0 / 6.0) * w * w * w;
			xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - xWeight[3];
			xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
			xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];


			/* y-gradient */
			w = y - (double)yIndex[1];
			gyWeight[3] = (1.0/ 2.0) * w * w ;
			gyWeight[0] = (-1.0 / 2.0) + w - gyWeight[3];
			gyWeight[2] = 3.0 * gyWeight[0] + 2.0 * (1 - w);
			gyWeight[1] = 3.0 * gyWeight[3] - 2.0 * w;

                        /* y */
			yWeight[3] = (1.0 / 6.0) * w * w * w;
			yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - yWeight[3];
			yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
			yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];
			/* z */
			if (ind3d == 1){
			  /* z- gradient */
			  w = z - (double)zIndex[1];
			  gzWeight[3] = (1.0/ 2.0) * w * w ;
			  gzWeight[0] = (-1.0 / 2.0) + w - gzWeight[3];
			  gzWeight[2] = 3.0 * gzWeight[0] + 2.0 * (1 - w);
			  gzWeight[1] = 3.0 * gzWeight[3] - 2.0 * w; 

			  /* z */
			  zWeight[3] = (1.0 / 6.0) * w * w * w;
			  zWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - zWeight[3];
			  zWeight[2] = w + zWeight[0] - 2.0 * zWeight[3];
			  zWeight[1] = 1.0 - zWeight[0] - zWeight[2] - zWeight[3];
			}
			break;
		case 4L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= (1.0 / 24.0) * xWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			xWeight[1] = t1 + t0;
			xWeight[3] = t1 - t0;
			xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
			xWeight[2] = 1.0 - xWeight[0] - xWeight[1] - xWeight[3] - xWeight[4];
			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= (1.0 / 24.0) * yWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			yWeight[1] = t1 + t0;
			yWeight[3] = t1 - t0;
			yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
			yWeight[2] = 1.0 - yWeight[0] - yWeight[1] - yWeight[3] - yWeight[4];
			/* z */
			if (ind3d == 1){
			  w = z - (double)zIndex[2];
			  w2 = w * w;
			  t = (1.0 / 6.0) * w2;
			  zWeight[0] = 1.0 / 2.0 - w;
			  zWeight[0] *= zWeight[0];
			  zWeight[0] *= (1.0 / 24.0) * zWeight[0];
			  t0 = w * (t - 11.0 / 24.0);
			  t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			  zWeight[1] = t1 + t0;
			  zWeight[3] = t1 - t0;
			  zWeight[4] = zWeight[0] + t0 + (1.0 / 2.0) * w;
			  zWeight[2] = 1.0 - zWeight[0] - zWeight[1] - zWeight[3] - zWeight[4];
			}
			break;
		case 5L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - xWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			xWeight[2] = t0 + t1;
			xWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			xWeight[1] = t0 + t1;
			xWeight[4] = t0 - t1;
			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - yWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			yWeight[2] = t0 + t1;
			yWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			yWeight[1] = t0 + t1;
			yWeight[4] = t0 - t1;
			/* z */
			if (ind3d ==1){
			  w = z - (double)zIndex[2];
			  w2 = w * w;
			  zWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			  w2 -= w;
			  w4 = w2 * w2;
			  w -= 1.0 / 2.0;
			  t = w2 * (w2 - 3.0);
			  zWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - zWeight[5];
			  t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			  t1 = (-1.0 / 12.0) * w * (t + 4.0);
			  zWeight[2] = t0 + t1;
			  zWeight[3] = t0 - t1;
			  t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			  t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			  zWeight[1] = t0 + t1;
			  zWeight[4] = t0 - t1;
			}
			break;
		default:
			printf("Invalid spline degree\n");
			return(0.0);
	}
	/* apply the mirror boundary conditions  */
	for (k = 0L; k <= SplineDegree; k++) {
		xIndex[k] = (Width == 1L) ? (0L) : ((xIndex[k] < 0L) ?
			(-xIndex[k] - Width2 * ((-xIndex[k]) / Width2))
			: (xIndex[k] - Width2 * (xIndex[k] / Width2)));
		if (Width <= xIndex[k]) {
			xIndex[k] = Width2 - xIndex[k];
		}
		yIndex[k] = (Height == 1L) ? (0L) : ((yIndex[k] < 0L) ?
			(-yIndex[k] - Height2 * ((-yIndex[k]) / Height2))
			: (yIndex[k] - Height2 * (yIndex[k] / Height2)));
		if (Height <= yIndex[k]) {
			yIndex[k] = Height2 - yIndex[k];
		}
		if (ind3d ==1){
		  zIndex[k] = (Slice == 1L) ? (0L) : ((zIndex[k] < 0L) ?
			(-zIndex[k] - Slice2 * ((-zIndex[k]) / Slice2))
			: (zIndex[k] - Slice2 * (zIndex[k] / Slice2)));
		  if (Slice <= zIndex[k]) {
			zIndex[k] = Slice2 - zIndex[k];
		  }
		}		
	     }
	/* perform interpolation */
	interpolated = 0.0;
	grad_x = 0.0;
	grad_y = 0.0;
	grad_z = 0.0;

	if (ind3d != 1){
	  	for (j = 0L; j <= SplineDegree; j++) {
		  p = Bcoeff + (ptrdiff_t)(yIndex[j] * Width);
		  w = 0.0;
		  wx = 0.0;
		  for (i = 0L; i <= SplineDegree; i++) {
		    w += xWeight[i] * p[xIndex[i]];
		    wx += gxWeight[i] * p[xIndex[i]];
		  }
		  interpolated += yWeight[j] * w;
		  grad_x += yWeight[j] * wx;
		  grad_y += gyWeight[j] * w;
		}
	}
	else {
	  for (k = 0L; k <= SplineDegree; k++) {
	    w2 = 0.0;
	    w2x= 0.0;
	    w2y =0.0;
	    for (j = 0L; j <= SplineDegree; j++) {
	            w = 0.0;
                    wx = 0.0;
		for (i = 0L; i <= SplineDegree; i++) {
		        p = Bcoeff + (ptrdiff_t)(yIndex[j] * Width + zIndex[k] * Width * Height);
			w += xWeight[i] * p[xIndex[i]];
			wx += gxWeight[i] * p[xIndex[i]];
		}
		w2 += yWeight[j] * w;
		w2x += yWeight[j] * wx;
		w2y += gyWeight[j] * w;
	    }
	    interpolated += zWeight[k] * w2;
	    grad_x += zWeight[k] * w2x;
	    grad_y += zWeight[k] * w2y;
	    grad_z += gzWeight[k] * w2;
	  }
	}
	
	/*   printf("%f\n", grad_y); */
       	
	*xgradient = grad_x;
	*ygradient = grad_y;
	*zgradient = grad_z;  
	return(interpolated);
} /* end InterpolatedValue */

/*****************************************************************************
 *	Definition of extern procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
extern double	InterpolatedValue_nogr
				(
				        double	*Bcoeff,   /* input B-spline array of coefficients */

					long	Width,	   /* width of the image */
					long	Height,	   /* height of the image */
					long    Slice,
					double	x,	   /* x coordinate where to interpolate */
					double	y,	   /* y coordinate where to interpolate */
					double  z,
					int     ind3d,
					long	SplineDegree  /* degree of the spline model */
			
				)
{ /* begin InterpolatedValue */
	double	*p;
	double	xWeight[6], yWeight[6], zWeight[6];
	double	interpolated;
	double	w, wx, w2, w2x, w2y, w4, t, t0, t1;
	long	xIndex[6], yIndex[6], zIndex[6];
	long	Width2 = 2L * Width - 2L, Height2 = 2L * Height - 2L, Slice2 = 2L * Slice - 2L;
	long	i, j, k, l;
	/* compute the interpolation indexes */
	if (SplineDegree & 1L) {
		i = (long)floor(x) - SplineDegree / 2L;
		j = (long)floor(y) - SplineDegree / 2L;
		k = (long)floor(z) - SplineDegree / 2L;
		for (l = 0L; l <= SplineDegree; l++) {
			xIndex[l] = i++;
			yIndex[l] = j++;
			zIndex[l] = k++;
		}
	}
	else {
		i = (long)floor(x + 0.5) - SplineDegree / 2L;
		j = (long)floor(y + 0.5) - SplineDegree / 2L;
		k = (long)floor(z + 0.5) - SplineDegree / 2L;
		for (l = 0L; l <= SplineDegree; l++) {
			xIndex[l] = i++;
			yIndex[l] = j++;
			zIndex[l] = k++;
		}
	}
	/* compute the interpolation weights */
	switch (SplineDegree) {
		case 2L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[1] = 3.0 / 4.0 - w * w;
			xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
			xWeight[0] = 1.0 - xWeight[1] - xWeight[2];
			/* y */
			w = y - (double)yIndex[1];
			yWeight[1] = 3.0 / 4.0 - w * w;
			yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
			yWeight[0] = 1.0 - yWeight[1] - yWeight[2];
			/* z */
			if (ind3d ==1){
			  w = z - (double)zIndex[1];
			  zWeight[1] = 3.0 / 4.0 - w * w;
			  zWeight[2] = (1.0 / 2.0) * (w - zWeight[1] + 1.0);
			  zWeight[0] = 1.0 - zWeight[1] - zWeight[2];
			}
			break;
		case 3L:
			/* x-gradient */
			w = x - (double)xIndex[1];
						/* x */
			xWeight[3] = (1.0 / 6.0) * w * w * w;
			xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - xWeight[3];
			xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
			xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];


			/* y-gradient */
			w = y - (double)yIndex[1];
		
                        /* y */
			yWeight[3] = (1.0 / 6.0) * w * w * w;
			yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - yWeight[3];
			yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
			yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];
			/* z */
			if (ind3d == 1){
			  /* z- gradient */
			  w = z - (double)zIndex[1];
			 
			  /* z */
			  zWeight[3] = (1.0 / 6.0) * w * w * w;
			  zWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - zWeight[3];
			  zWeight[2] = w + zWeight[0] - 2.0 * zWeight[3];
			  zWeight[1] = 1.0 - zWeight[0] - zWeight[2] - zWeight[3];
			}
			break;
		case 4L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= (1.0 / 24.0) * xWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			xWeight[1] = t1 + t0;
			xWeight[3] = t1 - t0;
			xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
			xWeight[2] = 1.0 - xWeight[0] - xWeight[1] - xWeight[3] - xWeight[4];
			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= (1.0 / 24.0) * yWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			yWeight[1] = t1 + t0;
			yWeight[3] = t1 - t0;
			yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
			yWeight[2] = 1.0 - yWeight[0] - yWeight[1] - yWeight[3] - yWeight[4];
			/* z */
			if (ind3d == 1){
			  w = z - (double)zIndex[2];
			  w2 = w * w;
			  t = (1.0 / 6.0) * w2;
			  zWeight[0] = 1.0 / 2.0 - w;
			  zWeight[0] *= zWeight[0];
			  zWeight[0] *= (1.0 / 24.0) * zWeight[0];
			  t0 = w * (t - 11.0 / 24.0);
			  t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			  zWeight[1] = t1 + t0;
			  zWeight[3] = t1 - t0;
			  zWeight[4] = zWeight[0] + t0 + (1.0 / 2.0) * w;
			  zWeight[2] = 1.0 - zWeight[0] - zWeight[1] - zWeight[3] - zWeight[4];
			}
			break;
		case 5L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - xWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			xWeight[2] = t0 + t1;
			xWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			xWeight[1] = t0 + t1;
			xWeight[4] = t0 - t1;
			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - yWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			yWeight[2] = t0 + t1;
			yWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			yWeight[1] = t0 + t1;
			yWeight[4] = t0 - t1;
			/* z */
			if (ind3d ==1){
			  w = z - (double)zIndex[2];
			  w2 = w * w;
			  zWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			  w2 -= w;
			  w4 = w2 * w2;
			  w -= 1.0 / 2.0;
			  t = w2 * (w2 - 3.0);
			  zWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - zWeight[5];
			  t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			  t1 = (-1.0 / 12.0) * w * (t + 4.0);
			  zWeight[2] = t0 + t1;
			  zWeight[3] = t0 - t1;
			  t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			  t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			  zWeight[1] = t0 + t1;
			  zWeight[4] = t0 - t1;
			}
			break;
		default:
			printf("Invalid spline degree\n");
			return(0.0);
	}
	/* apply the mirror boundary conditions  */
	for (k = 0L; k <= SplineDegree; k++) {
		xIndex[k] = (Width == 1L) ? (0L) : ((xIndex[k] < 0L) ?
			(-xIndex[k] - Width2 * ((-xIndex[k]) / Width2))
			: (xIndex[k] - Width2 * (xIndex[k] / Width2)));
		if (Width <= xIndex[k]) {
			xIndex[k] = Width2 - xIndex[k];
		}
		yIndex[k] = (Height == 1L) ? (0L) : ((yIndex[k] < 0L) ?
			(-yIndex[k] - Height2 * ((-yIndex[k]) / Height2))
			: (yIndex[k] - Height2 * (yIndex[k] / Height2)));
		if (Height <= yIndex[k]) {
			yIndex[k] = Height2 - yIndex[k];
		}
		if (ind3d ==1){
		  zIndex[k] = (Slice == 1L) ? (0L) : ((zIndex[k] < 0L) ?
			(-zIndex[k] - Slice2 * ((-zIndex[k]) / Slice2))
			: (zIndex[k] - Slice2 * (zIndex[k] / Slice2)));
		  if (Slice <= zIndex[k]) {
			zIndex[k] = Slice2 - zIndex[k];
		  }
		}		
	     }
	/* perform interpolation */
	interpolated = 0.0;
	if (ind3d != 1){
	  	for (j = 0L; j <= SplineDegree; j++) {
		  p = Bcoeff + (ptrdiff_t)(yIndex[j] * Width);
		  w = 0.0;
		  wx = 0.0;
		  for (i = 0L; i <= SplineDegree; i++) {
		    w += xWeight[i] * p[xIndex[i]];
		  }
		  interpolated += yWeight[j] * w;
		}
	}
	else {
	  for (k = 0L; k <= SplineDegree; k++) {
	    w2 = 0.0;
	    w2x= 0.0;
	    w2y =0.0;
	    for (j = 0L; j <= SplineDegree; j++) {
	            w = 0.0;
                    wx = 0.0;
		for (i = 0L; i <= SplineDegree; i++) {
		        p = Bcoeff + (ptrdiff_t)(yIndex[j] * Width + zIndex[k] * Width * Height);
			w += xWeight[i] * p[xIndex[i]];
		}
		w2 += yWeight[j] * w;
	    }
	    interpolated += zWeight[k] * w2;
	  }
	}
	
	/*   printf("%f\n", grad_y); */
     
	return(interpolated);
} /* end InterpolatedValue */


     
/************************************************************/
/*              Mex function interface                      */ 
/************************************************************/

void mexFunction(
      int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])   
{
      double           *img_coeff;
      double           *img_ref;

      double           *xcoeff;
      double           *ycoeff;
      double           *zcoeff;

      double           *xloc_x;
      double           *xloc_y;
      double           *xloc_z;

      double           *yloc_x;
      double           *yloc_y;
      double           *yloc_z;

      double           *zloc_x;
      double           *zloc_y;
      double           *zloc_z;

      float           *img_out;
      float           *xgrad;
      float           *ygrad;
      float           *zgrad;

      double           *grx_pnt;
      double           *gry_pnt;
      double           *grz_pnt;

      float           *x_tran;
      float           *y_tran;
      float           *z_tran;

      float           *x_penx;
      float           *x_peny;
      float           *x_penz;
/*
      double           *x_penpen_xx;
      double           *x_penpen_yy;
      double           *x_penpen_zz;
      double           *x_penpen_xy;
      double           *x_penpen_yz;
      double           *x_penpen_zx;
*/
      float           *y_penx;
      float           *y_peny;
      float           *y_penz;
/*
      double           *y_penpen_xx;
      double           *y_penpen_yy;
      double           *y_penpen_zz;
      double           *y_penpen_xy;
      double           *y_penpen_yz;
      double           *y_penpen_zx;
*/

      float           *z_penx;
      float           *z_peny;
      float           *z_penz;
/*
      double           *z_penpen_xx;
      double           *z_penpen_yy;
      double           *z_penpen_zz;
      double           *z_penpen_xy;
      double           *z_penpen_yz;
      double           *z_penpen_zx;
*/

      float           *jacobian;
      double            jac_thr;
      float           *grad_xcoeff, *grad_ycoeff, *grad_zcoeff;

      double           *obj_val, *pen_val;
      double           *obj_grad_x, *obj_grad_y, *obj_grad_z, *hess_obj_x, *hess_obj_y, *hess_obj_z, *pen_grad_x, *pen_grad_y, *pen_grad_z, *pen_hess_x, *pen_hess_y, *pen_hess_z;
      double           *hess_diag;

      double           *xi, *yi, *zi;
      double           xi_ax, yi_ax, zi_ax, xt, yt, zt;
    
      double           mse, img_pix,  wx, wy, wz, gwx, gwy, gwz, ggwx, ggwy, ggwz;
      double           fdouble_subi, fdouble_subj, fdouble_subk, double_subi, double_subj, double_subk;
      double           isgn, jsgn, ksgn;

      double           penalty, grad_penx, grad_peny, grad_penz, reg_param, jac_pen, slop;
      double           x_pen_grad_theta, y_pen_grad_theta, z_pen_grad_theta;
  
      long             width, height, slice,t_size, loc, locx, locy, locz, num_spl_y, num_pix, out_reg, splinedegree, sub_smp_xax, sub_smp_yax, sub_smp_zax, x_limit, y_limit, z_limit;   
      int              ndim, ii, i, j, k, test, ind3d;            
      const    int     *dim, *x_size, *y_size, *z_size, *x_con_size, *y_con_size, *z_con_size;
      
      if(nrhs!=27){
            mexErrMsgTxt("Number of Input argument is not correct"); 
       }
      if(nlhs!=22){
	    mexErrMsgTxt("Number of output argument is not correct");
      }

      img_coeff = mxGetPr(prhs[0]);
      img_ref = mxGetPr(prhs[1]);

      xi = mxGetPr(prhs[2]);
      yi = mxGetPr(prhs[3]);
      zi = mxGetPr(prhs[4]);

      xcoeff = mxGetPr(prhs[5]);
      xloc_x = mxGetPr(prhs[6]);
      xloc_y = mxGetPr(prhs[7]);
      xloc_z = mxGetPr(prhs[8]);

      ycoeff = mxGetPr(prhs[9]);
      yloc_x = mxGetPr(prhs[10]);
      yloc_y = mxGetPr(prhs[11]);
      yloc_z = mxGetPr(prhs[12]);

      zcoeff = mxGetPr(prhs[13]);
      zloc_x = mxGetPr(prhs[14]);
      zloc_y = mxGetPr(prhs[15]);
      zloc_z = mxGetPr(prhs[16]);
      
      splinedegree = mxGetScalar(prhs[17]);
      
      sub_smp_xax = mxGetScalar(prhs[18]);
      sub_smp_yax = mxGetScalar(prhs[19]);
      sub_smp_zax = mxGetScalar(prhs[20]);


      x_limit = mxGetScalar(prhs[21]);
      y_limit = mxGetScalar(prhs[22]);
      z_limit = mxGetScalar(prhs[23]);
 

      reg_param = mxGetScalar(prhs[24]);
      slop = mxGetScalar(prhs[25]);
      jac_thr=mxGetScalar(prhs[26]);
      
      ndim = mxGetNumberOfDimensions(prhs[1]);
      dim  = mxGetDimensions(prhs[1]);

     
      
      y_con_size = mxGetDimensions(prhs[9]);     
      num_spl_y = y_con_size[1];
      
      x_size = mxGetDimensions(prhs[2]);
      y_size = mxGetDimensions(prhs[3]);
      z_size = mxGetDimensions(prhs[4]);
      
      width  = x_size[1];
      height = y_size[1];
      slice  = z_size[1];

  
      t_size = width * height * slice;
      ind3d=1;

      if (slice==1) {
	ind3d=0;
	printf(" Input image is 2-Dimensional!!\n");
      }

      printf("%d, %d, %d\n",width, height, slice); 
      printf("Number of control point is %d\n", num_spl_y);      
      printf("Sub sampling is %d, %d, %d\n", sub_smp_xax,sub_smp_yax,sub_smp_zax);
      printf("Bounds are %d, %d, %d\n", x_limit, y_limit, z_limit);

      plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);  
      obj_val = mxGetData(plhs[0]);
    
      plhs[1] = mxCreateNumericArray(ndim, dim, mxSINGLE_CLASS, mxREAL);         
      img_out = mxGetData(plhs[1]);

      plhs[2] = mxCreateNumericArray(ndim, dim, mxSINGLE_CLASS, mxREAL);         
      xgrad   = mxGetData(plhs[2]);

      plhs[3] = mxCreateNumericArray(ndim, dim, mxSINGLE_CLASS, mxREAL);              
      ygrad   = mxGetData(plhs[3]);
     
      plhs[4] = mxCreateNumericArray(ndim, dim, mxSINGLE_CLASS, mxREAL);              
      zgrad   = mxGetData(plhs[4]);

      plhs[5] = mxCreateNumericArray(ndim, dim, mxSINGLE_CLASS, mxREAL);         
      x_tran = mxGetData(plhs[5]);

      plhs[6] = mxCreateNumericArray(ndim, dim, mxSINGLE_CLASS, mxREAL);              
      y_tran = mxGetData(plhs[6]);
     
      plhs[7] = mxCreateNumericArray(ndim, dim, mxSINGLE_CLASS, mxREAL);              
      z_tran = mxGetData(plhs[7]);

      plhs[8] = mxCreateDoubleMatrix(num_spl_y,1, mxREAL);           
      obj_grad_x = mxGetData(plhs[8]);

      plhs[9] = mxCreateDoubleMatrix(num_spl_y,1, mxREAL);           
      obj_grad_y = mxGetData(plhs[9]);

      plhs[10] = mxCreateDoubleMatrix(num_spl_y,1, mxREAL);           
      obj_grad_z = mxGetData(plhs[10]);

      plhs[11] = mxCreateDoubleMatrix(num_spl_y,1, mxREAL);           
      hess_obj_x = mxGetData(plhs[11]);

      plhs[12] = mxCreateDoubleMatrix(num_spl_y,1, mxREAL);           
      hess_obj_y = mxGetData(plhs[12]);

      plhs[13] = mxCreateDoubleMatrix(num_spl_y,1, mxREAL);           
      hess_obj_z = mxGetData(plhs[13]);

      plhs[14] = mxCreateDoubleMatrix(1, 1, mxREAL);         
      pen_val = mxGetData(plhs[14]);

      plhs[15] = mxCreateDoubleMatrix(num_spl_y,1, mxREAL);           
      pen_grad_x = mxGetData(plhs[15]);

      plhs[16] = mxCreateDoubleMatrix(num_spl_y,1, mxREAL);           
      pen_grad_y = mxGetData(plhs[16]);

      plhs[17] = mxCreateDoubleMatrix(num_spl_y,1, mxREAL);           
      pen_grad_z = mxGetData(plhs[17]);

      plhs[18] = mxCreateDoubleMatrix(num_spl_y,1, mxREAL);           
      pen_hess_x = mxGetData(plhs[18]);

      plhs[19] = mxCreateDoubleMatrix(num_spl_y,1, mxREAL);           
      pen_hess_y = mxGetData(plhs[19]);

      plhs[20] = mxCreateDoubleMatrix(num_spl_y,1, mxREAL);           
      pen_hess_z = mxGetData(plhs[20]);

      plhs[21] = mxCreateNumericArray(ndim, dim, mxSINGLE_CLASS, mxREAL);              
      jacobian = mxGetData(plhs[21]);


    
      grx_pnt = mxCalloc(1,sizeof(double));
      gry_pnt = mxCalloc(1,sizeof(double));
      grz_pnt = mxCalloc(1,sizeof(double));

      grad_xcoeff = mxCalloc(t_size, sizeof(float));
      grad_ycoeff = mxCalloc(t_size, sizeof(float));
      grad_zcoeff = mxCalloc(t_size, sizeof(float));

      x_penx = mxCalloc(t_size, sizeof(float));
      x_peny = mxCalloc(t_size, sizeof(float));
      x_penz = mxCalloc(t_size, sizeof(float));
/*
      x_penpen_xx = mxCalloc(t_size, sizeof(double));
      x_penpen_yy = mxCalloc(t_size, sizeof(double));
      x_penpen_zz = mxCalloc(t_size, sizeof(double));
      x_penpen_xy = mxCalloc(t_size, sizeof(double));
      x_penpen_yz = mxCalloc(t_size, sizeof(double));
      x_penpen_zx = mxCalloc(t_size, sizeof(double));

*/
      y_penx = mxCalloc(t_size, sizeof(float));
      y_peny = mxCalloc(t_size, sizeof(float));
      y_penz = mxCalloc(t_size, sizeof(float));
/*
      y_penpen_xx = mxCalloc(t_size, sizeof(double));
      y_penpen_yy = mxCalloc(t_size, sizeof(double));
      y_penpen_zz = mxCalloc(t_size, sizeof(double));
      y_penpen_xy = mxCalloc(t_size, sizeof(double));
      y_penpen_yz = mxCalloc(t_size, sizeof(double));
      y_penpen_zx = mxCalloc(t_size, sizeof(double));

*/
      z_penx = mxCalloc(t_size, sizeof(float));
      z_peny = mxCalloc(t_size, sizeof(float));
      z_penz = mxCalloc(t_size, sizeof(float));
/*
      z_penpen_xx = mxCalloc(t_size, sizeof(double));
      z_penpen_yy = mxCalloc(t_size, sizeof(double));
      z_penpen_zz = mxCalloc(t_size, sizeof(double));
      z_penpen_xy = mxCalloc(t_size, sizeof(double));
      z_penpen_yz = mxCalloc(t_size, sizeof(double));
      z_penpen_zx = mxCalloc(t_size, sizeof(double));
  */    
      printf("allocation, came here?\n");

      /*******************************************************************/
      /*                                                                 */
      /*   Computing transformation and its gradient w.r.t x, y, z       */
      /*                                                                 */
      /*******************************************************************/

      for (ii=0; ii<num_spl_y; ++ii){
	
		  
	for (i = 0; i < (4 * sub_smp_xax -1); ++i){
	  for (j = 0; j < (4 * sub_smp_yax -1); ++j){
	    for (k = 0; k < (4 * sub_smp_zax -1); ++k){ 
	     
	     locx = (i - 2 * sub_smp_xax + 1 + (long)yloc_x[ii]);
             locy = (j - 2 * sub_smp_yax + 1 + (long)yloc_y[ii]);
             locz = (k - 2 * sub_smp_zax + 1 + (long)yloc_z[ii]);

             loc = locx + locy * width + locz * width *height *ind3d;
	
	     if ((locx >= 0) && (locx < width) && (locy >=0) && (locy < height) && (locz >=0) && (locz < slice)){
	             
               double_subi = (double)(i - 2 * sub_smp_xax + 1)/(double)(sub_smp_xax);
	       double_subj = (double)(j - 2 * sub_smp_yax + 1)/(double)(sub_smp_yax);
	       double_subk = (double)(k - 2 * sub_smp_zax + 1)/(double)(sub_smp_zax);

	       fdouble_subi = fabs(double_subi);
	       fdouble_subj = fabs(double_subj);
	       fdouble_subk = fabs(double_subk);

	       if ( fdouble_subi > 0) {
		 isgn = double_subi / fdouble_subi;
	       }
	       else{
		 isgn = 0.0;
	       }

	       if ( fdouble_subj > 0) {
		 jsgn = double_subj / fdouble_subj;
	       }
	       else{
		 jsgn = 0.0;
	       }

	       if ( fdouble_subk > 0) {
		 ksgn = double_subk / fdouble_subk;
	       }
	       else{
		 ksgn = 0.0;
	       }


	       if ( fdouble_subi < 1.0 ) {
		 wx  = (2.0 /3.0) - fdouble_subi * fdouble_subi + (1.0/2.0) * fdouble_subi * fdouble_subi * fdouble_subi;
		 gwx = isgn * ((3.0 /2.0) * fdouble_subi * fdouble_subi -  2.0 * fdouble_subi)/sub_smp_xax;
		 ggwx = ((3.0) * fdouble_subi - 2.0)/sub_smp_xax/sub_smp_xax;
	       }
	       else{
		 wx  = (2.0 - fdouble_subi) * (2.0 - fdouble_subi) * (2.0 - fdouble_subi) / 6.0;
		 gwx = isgn * ((-1.0 /2.0) * ( 2.0 - fdouble_subi) *(2.0 - fdouble_subi))/sub_smp_xax;
		 ggwx = (2.0 - fdouble_subi)/sub_smp_xax/sub_smp_xax; 
	       }

	       if ( fdouble_subj <  1.0) {
		 wy =  (2.0 /3.0) - fdouble_subj * fdouble_subj + (1.0/2.0) * fdouble_subj * fdouble_subj * fdouble_subj;
		 gwy = jsgn * ((3.0 /2.0) * fdouble_subj * fdouble_subj -  2.0 * fdouble_subj)/sub_smp_yax;
		 ggwy = ((3.0) * fdouble_subj - 2.0)/sub_smp_yax/sub_smp_yax;
	       }
	       else{
		 wy = (2.0 - fdouble_subj) * (2.0 - fdouble_subj) * (2.0 - fdouble_subj)/ 6.0;
		 gwy = jsgn * ((-1.0 /2.0) * ( 2.0 - fdouble_subj) *(2.0 - fdouble_subj))/sub_smp_yax;
		 ggwy = (2.0 - fdouble_subj)/sub_smp_yax/sub_smp_yax;
	       }

	       if ( fdouble_subk <  1.0) {
		 wz =  (2.0 /3.0) - fdouble_subk * fdouble_subk + (1.0/2.0) * fdouble_subk * fdouble_subk * fdouble_subk;
		 gwz = ksgn * ((3.0 /2.0) * fdouble_subk * fdouble_subk - 2.0 * fdouble_subk)/sub_smp_zax;
		 ggwz = ((3.0) * fdouble_subk - 2.0)/sub_smp_zax/sub_smp_zax;
	       }
	       else{
		 wz = (2.0 - fdouble_subk) * (2.0 - fdouble_subk) * (2.0 - fdouble_subk)/ 6.0;
		 gwz = ksgn * ((-1.0 /2.0) * ( 2.0 - fdouble_subk) *(2.0 - fdouble_subk))/sub_smp_zax;
		 ggwz = (2.0 - fdouble_subk)/sub_smp_zax/sub_smp_zax;
	       }
		  
	 
	       x_tran[loc] += xcoeff[ii] * wx * wy * wz;
	       y_tran[loc] += ycoeff[ii] * wx * wy * wz;
	       z_tran[loc] += zcoeff[ii] * wx * wy * wz;

	       x_penx[loc] += xcoeff[ii] * gwx * wy * wz; 
	       x_peny[loc] += xcoeff[ii] * wx * gwy * wz;
	       x_penz[loc] += xcoeff[ii] * wx * wy * gwz;
/*
	       x_penpen_xx[loc] += xcoeff[ii] * ggwx * wy * wz;
	       x_penpen_yy[loc] += xcoeff[ii] * wx * ggwy * wz;
	       x_penpen_zz[loc] += xcoeff[ii] * wx * wy * ggwz;
	       x_penpen_xy[loc] += xcoeff[ii] * gwx * gwy * wz;
	       x_penpen_zx[loc] += xcoeff[ii] * gwx * wy * gwz;
	       x_penpen_yz[loc] += xcoeff[ii] * wx * gwy * gwz;
*/
	       y_penx[loc] += ycoeff[ii] * gwx * wy * wz;
	       y_peny[loc] += ycoeff[ii] * wx * gwy * wz;
	       y_penz[loc] += ycoeff[ii] * wx * wy * gwz;
  /*        
	       y_penpen_xx[loc] += ycoeff[ii] * ggwx * wy * wz;
	       y_penpen_yy[loc] += ycoeff[ii] * wx * ggwy * wz;
	       y_penpen_zz[loc] += ycoeff[ii] * wx * wy * ggwz;
	       y_penpen_xy[loc] += ycoeff[ii] * gwx * gwy * wz;
	       y_penpen_zx[loc] += ycoeff[ii] * gwx * wy * gwz;
	       y_penpen_yz[loc] += ycoeff[ii] * wx * gwy * gwz;

*/
	       z_penx[loc] += zcoeff[ii] * gwx * wy * wz;
	       z_peny[loc] += zcoeff[ii] * wx * gwy * wz;
	       z_penz[loc] += zcoeff[ii] * wx * wy * gwz;
/*
	       z_penpen_xx[loc] += zcoeff[ii] * ggwx * wy * wz;
	       z_penpen_yy[loc] += zcoeff[ii] * wx * ggwy * wz;
	       z_penpen_zz[loc] += zcoeff[ii] * wx * wy * ggwz;
	       z_penpen_xy[loc] += zcoeff[ii] * gwx * gwy * wz;
	       z_penpen_zx[loc] += zcoeff[ii] * gwx * wy * gwz;
	       z_penpen_yz[loc] += zcoeff[ii] * wx * gwy * gwz; 
*/
	     }
	    }
	  }
	}
	 
      }
    
      printf("came here?\n");

      /*****************************************/
      /*                                       */
      /*  Computing Mean Square Error          */
      /*                                       */
      /*****************************************/


     mse = 0.0;
     penalty = 0.0;
     jac_pen = 0.0;

     for (k= 0 ; k < slice  ; ++k){
	for (j= 0 ; j < height ; ++j){
	  for (i= 0 ; i < width ; ++i){

	      loc = i + j * width  + ind3d * k * height * width;
	
              xi_ax= xi[i] + x_tran[loc];  
	      yi_ax= yi[j] + y_tran[loc]; 
	      zi_ax= zi[k] + z_tran[loc];      
       
	      /*      printf("xi is %f, yi is %f, zi is %f\n", xi_ax, yi_ax, zi_ax); */  
              img_pix =  InterpolatedValue(img_coeff, grx_pnt, gry_pnt, grz_pnt, width + 2 * x_limit, height + 2 * y_limit, slice + 2 * z_limit, xi_ax, yi_ax, zi_ax, ind3d, splinedegree);

	      *img_out++ = img_pix; 

	      mse += (img_pix - img_ref[loc]) * (img_pix - img_ref[loc]);
	      
	      /*	      penalty += x_penx[loc] * x_penx[loc] + x_peny[loc] * x_peny[loc] + x_penz[loc] * x_penz[loc] + y_penx[loc] * y_penx[loc] + y_peny[loc] * y_peny[loc] + y_penz[loc] * y_penz[loc] + z_penx[loc] * z_penx[loc] + z_peny[loc] * z_peny[loc] + z_penz[loc] * z_penz[loc] ; 

	      penalty += x_penpen_xx[loc] * x_penpen_xx[loc] + x_penpen_yy[loc] * x_penpen_yy[loc] + x_penpen_zz[loc] * x_penpen_zz[loc] + y_penpen_xx[loc] * y_penpen_xx[loc] + y_penpen_yy[loc] * y_penpen_yy[loc] + y_penpen_zz[loc] * y_penpen_zz[loc] + z_penpen_xx[loc] * z_penpen_xx[loc] + z_penpen_yy[loc] * z_penpen_yy[loc] + z_penpen_zz[loc] * z_penpen_zz[loc];
	      penalty += (2.0) * x_penpen_xy[loc] * x_penpen_xy[loc] + (2.0) * x_penpen_zx[loc] * x_penpen_zx[loc] + (2.0) * x_penpen_yz[loc] * x_penpen_yz[loc] + (2.0) * y_penpen_xy[loc] * y_penpen_xy[loc] + (2.0) * y_penpen_zx[loc] * y_penpen_zx[loc] + (2.0) * y_penpen_yz[loc] * y_penpen_yz[loc] + (2.0) * z_penpen_xy[loc] * z_penpen_xy[loc] + (2.0) * z_penpen_zx[loc] * z_penpen_zx[loc] + (2.0) * z_penpen_yz[loc] * z_penpen_yz[loc];  */

	      *xgrad++ = grx_pnt[0];
	      *ygrad++ = gry_pnt[0];
	      *zgrad++ = grz_pnt[0];
	      
	      jacobian[loc] = (1.0 + x_penx[loc]) * ( (1.0 + y_peny[loc]) * (1.0 + z_penz[loc]) - y_penz[loc] * z_peny[loc]) - y_penx[loc] * (x_peny[loc] * (1.0 + z_penz[loc]) - x_penz[loc] * z_peny[loc]) + z_penx[loc] * (x_peny[loc] * y_penz[loc] - x_penz[loc] * (1.0 + y_peny[loc]));
          
	      if ( jacobian[loc] < jac_thr) {
		        jac_pen +=  0.5 * (jacobian[loc]  - jac_thr) * (jacobian[loc] - jac_thr);
        		/*printf("Low jacobian found!\n");*/
              }
	      else{
    		jac_pen += 0;	
	      }
	   /*   jac_pen += exp( -1.0 * slop * jacobian[loc]);*/
	      }
	    }
	}

     obj_val[0] = (1.0/2.0) * mse / (double)(width * height * slice);
     /*  pen_val[0] = (1.0/2.0) * reg_param * penalty / (double)(width * height * slice); */
     pen_val[0] = reg_param * jac_pen/(double)(width * height * slice);
     
      printf("came here?\n");

     /*******************************************/
     /*                                         */
     /*    Computing Gradient  & Hessian        */
     /*                                         */
     /*******************************************/


     xgrad -= (width * height * slice);
     ygrad -= (width * height * slice);
     zgrad -= (width * height * slice); 
     img_out -= (width * height * slice);

     printf("successful?\n");

     out_reg=0;

    
     for (ii=0; ii < num_spl_y; ++ii){	 
 

       for (i = 0; i < (4 * sub_smp_xax -1); ++i){
	 for (j = 0; j < (4 * sub_smp_yax -1); ++j){
	   for (k = 0; k < (4 * sub_smp_zax-1); ++k){ 

	     locx = (i - 2 * sub_smp_xax + 1 + (long)yloc_x[ii]);
             locy = (j - 2 * sub_smp_yax + 1 + (long)yloc_y[ii]);
             locz = (k - 2 * sub_smp_zax + 1 + (long)yloc_z[ii]);

             loc = locx + locy * width + locz * width *height *ind3d;

	     /*   printf("locx is %d, locy is %d, locz is %d, loc is %d\n",locx,locy,locz,loc);*/
	     
	     if ((locx >= 0 ) && (locx < width ) && (locy >= 0 ) && (locy < height ) && (locz >= 0 ) && (locz < slice)){

	       double_subi = (double)(i - 2 * sub_smp_xax + 1)/(double)(sub_smp_xax);
	       double_subj = (double)(j - 2 * sub_smp_yax + 1)/(double)(sub_smp_yax);
	       double_subk = (double)(k - 2 * sub_smp_zax + 1)/(double)(sub_smp_zax);

	       fdouble_subi = fabs(double_subi);
	       fdouble_subj = fabs(double_subj);
	       fdouble_subk = fabs(double_subk);

	       if ( fdouble_subi > 0) {
		 isgn = double_subi / fdouble_subi;
	       }
	       else{
		 isgn = 0.0;
	       }

	       if ( fdouble_subj > 0) {
		 jsgn = double_subj / fdouble_subj;
	       }
	       else{
		 jsgn = 0.0;
	       }

	       if ( fdouble_subk > 0) {
		 ksgn = double_subk / fdouble_subk;
	       }
	       else{
		 ksgn = 0.0;
	       }
	       
	       
	       if ( fdouble_subi < 1.0 ) {
		 wx  = (2.0 /3.0) - fdouble_subi * fdouble_subi + (1.0/2.0) * fdouble_subi * fdouble_subi * fdouble_subi;
		 gwx = isgn * ((3.0 /2.0) * fdouble_subi * fdouble_subi -  2.0 * fdouble_subi)/sub_smp_xax;
		 ggwx = ((3.0) * fdouble_subi - 2.0)/sub_smp_xax/sub_smp_xax;
	       }
	       else{
		 wx  = (2.0 - fdouble_subi) * (2.0 - fdouble_subi) * (2.0 - fdouble_subi) / 6.0;
		 gwx = isgn * ((-1.0 /2.0) * ( 2.0 - fdouble_subi) *(2.0 - fdouble_subi))/sub_smp_xax;
		 ggwx = (2.0 - fdouble_subi)/sub_smp_xax/sub_smp_xax; 
	       }

	       if ( fdouble_subj <  1.0) {
		 wy =  (2.0 /3.0) - fdouble_subj * fdouble_subj + (1.0/2.0) * fdouble_subj * fdouble_subj * fdouble_subj;
		 gwy = jsgn * ((3.0 /2.0) * fdouble_subj * fdouble_subj -  2.0 * fdouble_subj)/sub_smp_yax;
		 ggwy = ((3.0) * fdouble_subj - 2.0)/sub_smp_yax/sub_smp_yax;
	       }
	       else{
		 wy = (2.0 - fdouble_subj) * (2.0 - fdouble_subj) * (2.0 - fdouble_subj)/ 6.0;
		 gwy = jsgn * ((-1.0 /2.0) * ( 2.0 - fdouble_subj) *(2.0 - fdouble_subj))/sub_smp_yax;
		 ggwy = (2.0 - fdouble_subj)/sub_smp_yax/sub_smp_yax;
	       }

	       if ( fdouble_subk <  1.0) {
		 wz =  (2.0 /3.0) - fdouble_subk * fdouble_subk + (1.0/2.0) * fdouble_subk * fdouble_subk * fdouble_subk;
		 gwz = ksgn * ((3.0 /2.0) * fdouble_subk * fdouble_subk - 2.0 * fdouble_subk)/sub_smp_zax;
		 ggwz = ((3.0) * fdouble_subk - 2.0)/sub_smp_zax/sub_smp_zax;
	       }
	       else{
		 wz = (2.0 - fdouble_subk) * (2.0 - fdouble_subk) * (2.0 - fdouble_subk)/ 6.0;
		 gwz = ksgn * ((-1.0 /2.0) * ( 2.0 - fdouble_subk) *(2.0 - fdouble_subk))/sub_smp_zax;
		 ggwz = (2.0 - fdouble_subk)/sub_smp_zax/sub_smp_zax;
	       }
		
	        
	       grad_xcoeff[loc] = xgrad[loc] * wx * wy * wz;
	       grad_ycoeff[loc] = ygrad[loc] * wx * wy * wz;
	       grad_zcoeff[loc] = zgrad[loc] * wx * wy * wz;
	
	       /* This penalty is for bending energy */
	       /*  x_pen_grad_theta = x_penpen_xx[loc] * ggwx * wy * wz + x_penpen_yy[loc] * wx * ggwy *wz + x_penpen_zz[loc] * wx *wy * ggwz + (2.0) * x_penpen_xy[loc] * gwx *gwy *wz + (2.0) * x_penpen_yz[loc] * wx *gwy *gwz + (2.0) * x_penpen_zx[loc] *gwx * wy * gwz ;
	       y_pen_grad_theta = y_penpen_xx[loc] * ggwx * wy * wz + y_penpen_yy[loc] * wx * ggwy *wz + y_penpen_zz[loc] * wx *wy * ggwz + (2.0) * y_penpen_xy[loc] * gwx *gwy *wz + (2.0) * y_penpen_yz[loc] * wx *gwy *gwz + (2.0) * y_penpen_zx[loc] *gwx * wy * gwz ;
	       z_pen_grad_theta = z_penpen_xx[loc] * ggwx * wy * wz + z_penpen_yy[loc] * wx * ggwy *wz + z_penpen_zz[loc] * wx *wy * ggwz + (2.0) * z_penpen_xy[loc] * gwx *gwy *wz + (2.0) * z_penpen_yz[loc] * wx *gwy *gwz + (2.0) * z_penpen_zx[loc] *gwx * wy * gwz ; */
	      
	       /* This penalty is for negative Jacobian */
	       x_pen_grad_theta =  (gwx * wy * wz) * ( (1.0 + y_peny[loc]) * (1.0 + z_penz[loc]) - y_penz[loc] * z_peny[loc]) - y_penx[loc] * (wx * gwy * wz * (1.0 + z_penz[loc]) - wx * wy * gwz * z_peny[loc]) + z_penx[loc] * (wx * gwy * wz * y_penz[loc] - wx * wy * gwz * (1.0 + y_peny[loc]));

	       y_pen_grad_theta = (1.0 + x_penx[loc]) * ( wx * gwy * wz * (1.0 + z_penz[loc]) - wx * wy * gwz * z_peny[loc]) - gwx * wy * wz * (x_peny[loc] * (1.0 + z_penz[loc]) - x_penz[loc] * z_peny[loc]) + z_penx[loc] * (x_peny[loc] * wx * wy * gwz - x_penz[loc] * wx * gwy * wz);

	       z_pen_grad_theta = (1.0 + x_penx[loc]) * ( (1.0 + y_peny[loc]) *  wx * wy * gwz - y_penz[loc] * wx * gwy * wz) - y_penx[loc] * (x_peny[loc] * wx * wy * gwz - x_penz[loc] * wx * gwy * wz) + gwx * wy * wz * (x_peny[loc] * y_penz[loc] - x_penz[loc] * (1.0 + y_peny[loc]));

	       obj_grad_x[ii] += (img_out[loc] - img_ref[loc]) * grad_xcoeff[loc]; 
	       obj_grad_y[ii] += (img_out[loc] - img_ref[loc]) * grad_ycoeff[loc];
	       obj_grad_z[ii] += (img_out[loc] - img_ref[loc]) * grad_zcoeff[loc];

	       /*
	       pen_grad_x[ii] += (-1.0) * slop * reg_param * (x_pen_grad_theta) * exp( -1.0 * slop * jacobian[loc]);
	       pen_grad_y[ii] += (-1.0) * slop * reg_param * (y_pen_grad_theta) * exp( -1.0 * slop * jacobian[loc]);
	       pen_grad_z[ii] += (-1.0) * slop * reg_param * (z_pen_grad_theta) * exp( -1.0 * slop * jacobian[loc]);
	       */
	        
	     
	       hess_obj_x[ii] += grad_xcoeff[loc] * grad_xcoeff[loc];
	       hess_obj_y[ii] += grad_ycoeff[loc] * grad_ycoeff[loc];
	       hess_obj_z[ii] += grad_zcoeff[loc] * grad_zcoeff[loc];
	   
           /*
	       pen_hess_x[ii] += reg_param * slop * slop * x_pen_grad_theta * x_pen_grad_theta * exp( -1.0 * slop * jacobian[loc]);
	       pen_hess_y[ii] += reg_param * slop * slop * y_pen_grad_theta * y_pen_grad_theta * exp( -1.0 * slop * jacobian[loc]); 
	       pen_hess_z[ii] += reg_param * slop * slop * z_pen_grad_theta * z_pen_grad_theta * exp( -1.0 * slop * jacobian[loc]);
	       */
	       
	       if ( jacobian[loc] < jac_thr){
        		 pen_grad_x[ii] +=  reg_param * x_pen_grad_theta *  (jacobian[loc] - jac_thr);
        		 pen_grad_y[ii] +=  reg_param * y_pen_grad_theta *  (jacobian[loc] - jac_thr);
        		 pen_grad_z[ii] +=  reg_param * z_pen_grad_theta *  (jacobian[loc] - jac_thr);

        		 pen_hess_x[ii] += reg_param * ( x_pen_grad_theta * x_pen_grad_theta * jacobian[loc]);
	             pen_hess_y[ii] += reg_param * ( y_pen_grad_theta * y_pen_grad_theta * jacobian[loc]);
    	         pen_hess_z[ii] += reg_param * ( y_pen_grad_theta * y_pen_grad_theta * jacobian[loc]);
	       }
	       else{
		        pen_grad_x[ii] += 0;
        	    pen_grad_y[ii] += 0;
        		pen_grad_z[ii] += 0;

		        pen_hess_x[ii] += reg_param * slop *  x_pen_grad_theta * x_pen_grad_theta * (jacobian[loc] *jacobian[loc] * jacobian[loc] - 3.0 * jacobian[loc] + 3.0) / (jacobian[loc] * jacobian[loc] * jacobian[loc] * jacobian[loc]);
	            pen_hess_y[ii] += reg_param * slop * y_pen_grad_theta * y_pen_grad_theta * (jacobian[loc] *jacobian[loc] * jacobian[loc] - 3.0 * jacobian[loc] + 3.0) / (jacobian[loc] * jacobian[loc] * jacobian[loc] * jacobian[loc]);
	            pen_hess_z[ii] += reg_param * slop * z_pen_grad_theta * z_pen_grad_theta * (jacobian[loc] *jacobian[loc] * jacobian[loc] - 3.0 * jacobian[loc] + 3.0) / (jacobian[loc] * jacobian[loc] * jacobian[loc] * jacobian[loc]);

	       }


	     }

	     else{
	       out_reg+=1; 
	       /*	       printf("locx is %d, locy is %d, locz is %d\n",locx,locy,locz); */
	     }
	   }
	 }
       }
       

       obj_grad_x[ii] = obj_grad_x[ii]/(double)(width * height * slice);
       obj_grad_y[ii] = obj_grad_y[ii]/(double)(width * height * slice);
       obj_grad_z[ii] = obj_grad_z[ii]/(double)(width * height * slice);

       pen_grad_x[ii] = pen_grad_x[ii]/(double)(width * height * slice) ;
       pen_grad_y[ii] = pen_grad_y[ii]/(double)(width * height * slice);
       pen_grad_z[ii] = pen_grad_z[ii]/(double)(width * height * slice);

       hess_obj_x[ii] = hess_obj_x[ii]/(double)(width * height *slice); 
       hess_obj_y[ii] = hess_obj_y[ii]/(double)(width * height *slice); 
       hess_obj_z[ii] = hess_obj_z[ii]/(double)(width * height *slice); 

       pen_hess_x[ii] = pen_hess_x[ii]/(double)(width * height * slice); 
       pen_hess_y[ii] = pen_hess_y[ii]/(double)(width * height * slice); 
       pen_hess_z[ii] = pen_hess_z[ii]/(double)(width * height * slice); 
       
     }
     printf("Number of out region is %d\n",out_reg);

     mxFree(grx_pnt);mxFree(gry_pnt);mxFree(grz_pnt);mxFree(grad_xcoeff);mxFree(grad_ycoeff);mxFree(grad_zcoeff);
     mxFree(x_penx);mxFree(x_peny);mxFree(x_penz);mxFree(y_penx);mxFree(y_peny);mxFree(y_penz);mxFree(z_penx);mxFree(z_peny);mxFree(z_penz);
/*     mxFree(x_penpen_xx);mxFree(x_penpen_yy);mxFree(x_penpen_zz);mxFree(x_penpen_xy);mxFree(x_penpen_yz);mxFree(x_penpen_zx);
     mxFree(y_penpen_xx);mxFree(y_penpen_yy);mxFree(y_penpen_zz);mxFree(y_penpen_xy);mxFree(y_penpen_yz);mxFree(y_penpen_zx);
     mxFree(z_penpen_xx);mxFree(z_penpen_yy);mxFree(z_penpen_zz);mxFree(z_penpen_xy);mxFree(z_penpen_yz);mxFree(z_penpen_zx);
*/
 }


   














