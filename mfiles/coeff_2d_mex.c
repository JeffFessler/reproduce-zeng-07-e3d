/******************************************************************************
/*	System includes
 ****************************************************************************/
#include	<float.h>
#include	<math.h>
#include	<stddef.h>
#include	<stdio.h>
#include	<stdlib.h>
#include        "mex.h"
#include        "matrix.h"

/*****************************************************************************
*	Other includes
****************************************************************************/
#include	"coeff2d.h"
 
/*****************************************************************************
 *	Declaration of static procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
static void		ConvertToInterpolationCoefficients
				(
					double	c[],		/* input samples --> output coefficients */
					long	DataLength,	/* number of samples or coefficients */
					double	z[],		/* poles */
					long	NbPoles,	/* number of poles */
					double	Tolerance	/* admissible relative error */
				);
/*--------------------------------------------------------------------------*/
static void		GetColumn
				(
					float	*Image,		/* input image array */
					long	Width,		/* width of the image */
					long	x,			/* x coordinate of the selected line */
					double	Line[],		/* output linear array */
					long	Height		/* length of the line and height of the image */
				);
/*--------------------------------------------------------------------------*/
static void		GetRow
				(
					float	*Image,		/* input image array */
					long	y,			/* y coordinate of the selected line */
					double	Line[],		/* output linear array */
					long	Width		/* length of the line and width of the image */
				);
/*--------------------------------------------------------------------------*/
static double	InitialCausalCoefficient
				(
					double	c[],		/* coefficients */
					long	DataLength,	/* number of coefficients */
					double	z,			/* actual pole */
					double	Tolerance	/* admissible relative error */
				);
/*--------------------------------------------------------------------------*/
static double	InitialAntiCausalCoefficient
				(
					double	c[],		/* coefficients */
					long	DataLength,	/* number of samples or coefficients */
					double	z			/* actual pole */
				);
/*--------------------------------------------------------------------------*/
static void		PutColumn
				(
					float	*Image,		/* output image array */
					long	Width,		/* width of the image */
					long	x,			/* x coordinate of the selected line */
					double	Line[],		/* input linear array */
					long	Height		/* length of the line and height of the image */
				);
/*--------------------------------------------------------------------------*/
static void		PutRow
				(
					float	*Image,		/* output image array */
					long	y,			/* y coordinate of the selected line */
					double	Line[],		/* input linear array */
					long	Width		/* length of the line and width of the image */
				);
/*****************************************************************************
 *	Definition of static procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
static void		ConvertToInterpolationCoefficients
				(
					double	c[],		/* input samples --> output coefficients */
					long	DataLength,	/* number of samples or coefficients */
					double	z[],		/* poles */
					long	NbPoles,	/* number of poles */
					double	Tolerance	/* admissible relative error */
				)
{ /* begin ConvertToInterpolationCoefficients */
	double	Lambda = 1.0;
	long	n, k;
	/* special case required by mirror boundaries */
	if (DataLength == 1L) {
		return;
	}
	/* compute the overall gain */
	for (k = 0L; k < NbPoles; k++) {
		Lambda = Lambda * (1.0 - z[k]) * (1.0 - 1.0 / z[k]);
	}
	/* apply the gain */
	for (n = 0L; n < DataLength; n++) {
		c[n] *= Lambda;
	}
	/* loop over all poles */
	for (k = 0L; k < NbPoles; k++) {
		/* causal initialization */
		c[0] = InitialCausalCoefficient(c, DataLength, z[k], Tolerance);
		/* causal recursion */
		for (n = 1L; n < DataLength; n++) {
			c[n] += z[k] * c[n - 1L];
		}
		/* anticausal initialization */
		c[DataLength - 1L] = InitialAntiCausalCoefficient(c, DataLength, z[k]);
		/* anticausal recursion */
		for (n = DataLength - 2L; 0 <= n; n--) {
			c[n] = z[k] * (c[n + 1L] - c[n]);
		}
	}
} /* end ConvertToInterpolationCoefficients */
/*--------------------------------------------------------------------------*/
static double	InitialCausalCoefficient
				(
					double	c[],		/* coefficients */
					long	DataLength,	/* number of coefficients */
					double	z,			/* actual pole */
					double	Tolerance	/* admissible relative error */
				)
{ /* begin InitialCausalCoefficient */
	double	Sum, zn, z2n, iz;
	long	n, Horizon;
	/* this initialization corresponds to mirror boundaries */
	Horizon = DataLength;
	if (Tolerance > 0.0) {
		Horizon = (long)ceil(log(Tolerance) / log(fabs(z)));
	}
	if (Horizon < DataLength) {
		/* accelerated loop */
		zn = z;
		Sum = c[0];
		for (n = 1L; n < Horizon; n++) {
			Sum += zn * c[n];
			zn *= z;
		}
		return(Sum);
	}
	else {
		/* full loop */
		zn = z;
		iz = 1.0 / z;
		z2n = pow(z, (double)(DataLength - 1L));
		Sum = c[0] + z2n * c[DataLength - 1L];
		z2n *= z2n * iz;
		for (n = 1L; n <= DataLength - 2L; n++) {
			Sum += (zn + z2n) * c[n];
			zn *= z;
			z2n *= iz;
		}
		return(Sum / (1.0 - zn * zn));
	}
}
 /* end InitialCausalCoefficient */

/*--------------------------------------------------------------------------*/
static void		GetColumn
				(
					float	*Image,		/* input image array */
					long	Width,		/* width of the image */
					long	x,			/* x coordinate of the selected line */
					double	Line[],		/* output linear array */
					long	Height		/* length of the line */
				)
{ /* begin GetColumn */
	long	y;
	Image = Image + (ptrdiff_t)x;
	for (y = 0L; y < Height; y++) {
		Line[y] = (double)*Image;
		Image += (ptrdiff_t)Width;
	}
} /* end GetColumn */
/*--------------------------------------------------------------------------*/
static void		GetRow
				(
					float	*Image,		/* input image array */
					long	y,			/* y coordinate of the selected line */
					double	Line[],		/* output linear array */
					long	Width		/* length of the line */
				)
{ /* begin GetRow */
	long	x;
	Image = Image + (ptrdiff_t)(y * Width);
	for (x = 0L; x < Width; x++) {
		Line[x] = (double)*Image++;
	}
} /* end GetRow */
/*--------------------------------------------------------------------------*/
static double	InitialAntiCausalCoefficient
				(
					double	c[],		/* coefficients */
					long	DataLength,	/* number of samples or coefficients */
					double	z			/* actual pole */
				)
{ /* begin InitialAntiCausalCoefficient */
	/* this initialization corresponds to mirror boundaries */
	return((z / (z * z - 1.0)) * (z * c[DataLength - 2L] + c[DataLength - 1L]));
} /* end InitialAntiCausalCoefficient */
/*--------------------------------------------------------------------------*/
static void		PutColumn
				(
					float	*Image,		/* output image array */
					long	Width,		/* width of the image */
					long	x,			/* x coordinate of the selected line */
					double	Line[],		/* input linear array */
					long	Height		/* length of the line and height of the image */
				)
{ /* begin PutColumn */
	long	y;
	Image = Image + (ptrdiff_t)x;
	for (y = 0L; y < Height; y++) {
		*Image = (float)Line[y];
		Image += (ptrdiff_t)Width;
	}
} /* end PutColumn */
/*--------------------------------------------------------------------------*/
static void		PutRow
				(
					float	*Image,		/* output image array */
					long	y,			/* y coordinate of the selected line */
					double	Line[],		/* input linear array */
					long	Width		/* length of the line and width of the image */
				)
{ /* begin PutRow */
	long	x;
	Image = Image + (ptrdiff_t)(y * Width);
	for (x = 0L; x < Width; x++) {
		*Image++ = (float)Line[x];
	}
} /* end PutRow */

/*****************************************************************************
 *	Definition of extern procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
extern int		SamplesToCoefficients
				(
					float	*Image,		/* in-place processing */
					long	Width,		/* width of the image */
					long	Height,		/* height of the image */
					long	SplineDegree/* degree of the spline model */
				)
{ /* begin SamplesToCoefficients */
	double	*Line;
	double	Pole[2];
	long	NbPoles;
	long	x, y;
	/* recover the poles from a lookup table */
	switch (SplineDegree) {
		case 2L:
			NbPoles = 1L;
			Pole[0] = sqrt(8.0) - 3.0;
			break;
		case 3L:
			NbPoles = 1L;
			Pole[0] = sqrt(3.0) - 2.0;
			break;
		case 4L:
			NbPoles = 2L;
			Pole[0] = sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0;
			Pole[1] = sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0;
			break;
		case 5L:
			NbPoles = 2L;
			Pole[0] = sqrt(135.0 / 2.0 - sqrt(17745.0 / 4.0)) + sqrt(105.0 / 4.0)
				- 13.0 / 2.0;
			Pole[1] = sqrt(135.0 / 2.0 + sqrt(17745.0 / 4.0)) - sqrt(105.0 / 4.0)
				- 13.0 / 2.0;
			break;
		default:
			printf("Invalid spline degree\n");
			return(1);
	}
	/* convert the image samples into interpolation coefficients */
	/* in-place separable process, along x */
	Line = (double *)malloc((size_t)(Width * (long)sizeof(double)));
	if (Line == (double *)NULL) {
		printf("Row allocation failed\n");
		return(1);
	}
	for (y = 0L; y < Height; y++) {
		GetRow(Image, y, Line, Width);
		ConvertToInterpolationCoefficients(Line, Width, Pole, NbPoles, DBL_EPSILON);
		PutRow(Image, y, Line, Width);
	}
	free(Line);
	/* in-place separable process, along y */
	Line = (double *)malloc((size_t)(Height * (long)sizeof(double)));
	if (Line == (double *)NULL) {
		printf("Column allocation failed\n");
		return(1);
	}
	for (x = 0L; x < Width; x++) {
		GetColumn(Image, Width, x, Line, Height);
		ConvertToInterpolationCoefficients(Line, Height, Pole, NbPoles, DBL_EPSILON);
		PutColumn(Image, Width, x, Line, Height);
	}
	free(Line);
	return(0);
} /* end SamplesToCoefficients */



void mexFunction(
      int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])   
{
      const    double          *im_in;
      float                    *im_out;

      const    int     *im_size;
      long     width, height, num_pix,  splinedegree;
      int      dim_image, ii;


      if(nrhs!=2){
            mexErrMsgTxt("Number of Input argument is not correct"); 
      }
      if(nlhs!=1){
	    mexErrMsgTxt("Number of output argument is not correct");
      }

      im_in=mxGetPr(prhs[0]);
      splinedegree=mxGetScalar(prhs[1]);
   
      
      dim_image=mxGetNumberOfDimensions(prhs[0]);
      if (dim_image!=2)mexErrMsgTxt("Image should be 2-dimensional");
      im_size=mxGetDimensions(prhs[0]);


      width=im_size[0];
      height=im_size[1];
      num_pix=width*height;
 

      plhs[0]=mxCreateNumericArray(dim_image, im_size, mxSINGLE_CLASS, mxREAL);     
      im_out=mxGetData(plhs[0]);
      
      for(ii=0; ii<num_pix; ii++){
	im_out[ii]=(float)(*im_in);
	im_in++;
      }

      SamplesToCoefficients(im_out, width, height, splinedegree);
     
}
