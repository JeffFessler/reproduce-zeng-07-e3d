
/*********************************************************************************************
 *    3D Bspline interpolation.
 *    Programmer: Rongping Zeng
 *                (EECS:Systems, U of Michigan)
 *    05/2005
 *    Note: 1. Parts of the program for image interpolation are originated from M. Unser's group.
 *         
 *   
 *   
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
      double           *coeff;
      double           *xi, *yi, *zi;
      long             splinedegree;
      int              ind3d;

      float           *out;
      float           *xgrad, *ygrad,*zgrad;
      double	       img,*grx,*gry,*grz;
      double           hx, hy, hz;
      int               ndim;
      const int           *size, *pntsize;
      long           width, height, slice, vol;
      long            lpnt;
      int               i;
      
      if(nrhs!=8){
            mexErrMsgTxt("Number of Input argument is not correct"); 
       }
      if(nlhs!=4){
	    mexErrMsgTxt("Number of output argument is not correct");
      }

      coeff = mxGetPr(prhs[0]);
      xi = mxGetPr(prhs[1]);
      yi = mxGetPr(prhs[2]);
      zi = mxGetPr(prhs[3]);
      splinedegree = mxGetScalar(prhs[4]);      
      hx = mxGetScalar(prhs[5]);
      hy = mxGetScalar(prhs[6]);
      hz = mxGetScalar(prhs[7]);

      ndim = mxGetNumberOfDimensions(prhs[0]);
      size = mxGetDimensions(prhs[0]);
      width = size[0];
      height = size[1];
      slice = 1;
  
      pntsize = mxGetDimensions(prhs[1]);
      lpnt = pntsize[0];

      ind3d=0;
      if (ndim==3) {
	        ind3d=1;
	        slice = size[2];
/*	    	printf(" Input image is 3-Dimensional!!\n");*/
      }
      

/*      printf("%d, %d, %d\n",width, height, slice); 
      printf("Sub sampling is %d, %d, %d\n", hx,hy,hz);
*/

      plhs[0] = mxCreateNumericMatrix(lpnt,1,mxSINGLE_CLASS, mxREAL);           
      out = mxGetData(plhs[0]);
      plhs[1] = mxCreateNumericMatrix(lpnt,1,mxSINGLE_CLASS, mxREAL);           
      xgrad = mxGetData(plhs[1]);
      plhs[2] = mxCreateNumericMatrix(lpnt,1,mxSINGLE_CLASS, mxREAL);           
      ygrad = mxGetData(plhs[2]);
      plhs[3] = mxCreateNumericMatrix(lpnt,1,mxSINGLE_CLASS, mxREAL);           
      zgrad = mxGetData(plhs[3]);
      
      grx=mxCalloc(1, sizeof(double));
      gry=mxCalloc(1, sizeof(double));
      grz=mxCalloc(1, sizeof(double));


      /****---  interpolating -------------****/
      for(i=0;i<lpnt;i++){
            img  =  InterpolatedValue(coeff, grx, gry, grz, width, height, slice, (*xi)/hx, (*yi)/hy, (*zi)/hz, ind3d, splinedegree);
            *out=(float)img;
            out++; xi++; yi++; zi++; 
            *xgrad=(float)(*grx/hx); *ygrad=(float)(*gry/hy); *zgrad=(float)(*grz/hz);
            xgrad++; ygrad++; zgrad++;
            }
      out -= lpnt;
      xi -= lpnt;
      yi -= lpnt;
      zi -= lpnt;
      xgrad -= lpnt;
      ygrad -= lpnt;
      zgrad -= lpnt;

      mxFree(grx); mxFree(gry); mxFree(grz);
      }  
