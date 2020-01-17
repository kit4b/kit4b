// BKPlots.cpp : Defines the entry point for the console application.
//

// On windows the standard resources are used to hold the PLplot palletes, fonts and shapes
// On linux the PLplot palletes, font and shapes are loaded linked in binary images

// Linux example
// a file named foo.bar with some example text is 'imported' into 
// an object file using the following command:
//
//      ld -r -b binary -o foo.bar.o foo.bar
//
// That creates an bject file named "foo.bar.o" with the following 
// symbols:
//
//      _binary_foo_bar_start
//      _binary_foo_bar_end
//      _binary_foo_bar_size
//
// Note that the symbols are addresses (so for example, to get the 
// size value, you have to get the address of the _binary_foo_bar_size
// symbol).
//
// In my example, foo.bar is a simple text file, and this program will
// dump the contents of that file which has been linked in by specifying
// foo.bar.o as an object file input to the linker when the progrma is built
//extern "C" {
//   extern char binary_foo_bar_start[];
//   extern char binary_foo_bar_end[];
//   extern uintptr_t binary_foo_bar_size;
// }


//int main(void)
//{
//    printf( "address of start: %p\n", &binary_foo_bar_start);
//    printf( "address of end: %p\n", &binary_foo_bar_end);

//    for (char* p = binary_foo_bar_start; p != binary_foo_bar_end; ++p) {
//        putchar( *p);
//    }

//    return 0;
//}

#include "stdafx.h"
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "plplot.h"




//int _tmain(int argc, _TCHAR* argv[])
#define XPTS             35     // Data points in x
#define YPTS             46     // Data points in y

#define XSPA             2. / ( XPTS - 1 )
#define YSPA             2. / ( YPTS - 1 )

// polar plot data
#define PERIMETERPTS     100
#define RPTS             40
#define THETAPTS         40

// potential plot data
#define PPERIMETERPTS    100
#define PRPTS            40
#define PTHETAPTS        64
#define PNLEVEL          20

static PLFLT clevel[11] =
{ -1., -.8, -.6, -.4, -.2, 0, .2, .4, .6, .8, 1. };

// Transformation function

PLFLT tr[6] =
{ XSPA, 0.0, -1.0, 0.0, YSPA, -1.0 };

// pltr_data is unused so mark it with the PL_UNUSED macro
static void
mypltr( PLFLT x, PLFLT y, PLFLT *tx, PLFLT *ty, void * PL_UNUSED( pltr_data ) )
{
    *tx = tr[0] * x + tr[1] * y + tr[2];
    *ty = tr[3] * x + tr[4] * y + tr[5];
}

static void polar( void )
//polar contour plot example.
{
    int      i, j;
    PLcGrid2 cgrid2;
    PLFLT    **z;
    PLFLT    px[PERIMETERPTS], py[PERIMETERPTS];
    PLFLT    t, r, theta;
    PLFLT    lev[10];

    plenv( -1., 1., -1., 1., 0, -2 );
    plcol0( 1 );

//Perimeter
    for ( i = 0; i < PERIMETERPTS; i++ )
    {
        t     = ( 2. * M_PI / ( PERIMETERPTS - 1 ) ) * (double) i;
        px[i] = cos( t );
        py[i] = sin( t );
    }
    plline( PERIMETERPTS, px, py );

//create data to be contoured.
    plAlloc2dGrid( &cgrid2.xg, RPTS, THETAPTS );
    plAlloc2dGrid( &cgrid2.yg, RPTS, THETAPTS );
    plAlloc2dGrid( &z, RPTS, THETAPTS );
    cgrid2.nx = RPTS;
    cgrid2.ny = THETAPTS;

    for ( i = 0; i < RPTS; i++ )
    {
        r = i / (double) ( RPTS - 1 );
        for ( j = 0; j < THETAPTS; j++ )
        {
            theta           = ( 2. * M_PI / (double) ( THETAPTS - 1 ) ) * (double) j;
            cgrid2.xg[i][j] = r * cos( theta );
            cgrid2.yg[i][j] = r * sin( theta );
            z[i][j]         = r;
        }
    }

    for ( i = 0; i < 10; i++ )
    {
        lev[i] = 0.05 + 0.10 * (double) i;
    }

    plcol0( 2 );
    plcont( (const PLFLT * const *) z, RPTS, THETAPTS, 1, RPTS, 1, THETAPTS, lev, 10,
        pltr2, (void *) &cgrid2 );
    plcol0( 1 );
    pllab( "", "", "Polar Contour Plot" );
    plFree2dGrid( z, RPTS, THETAPTS );
    plFree2dGrid( cgrid2.xg, RPTS, THETAPTS );
    plFree2dGrid( cgrid2.yg, RPTS, THETAPTS );
}

//--------------------------------------------------------------------------
// f2mnmx
//
// Returns min & max of input 2d array.
//--------------------------------------------------------------------------

static void
f2mnmx( PLFLT **f, PLINT nx, PLINT ny, PLFLT *fnmin, PLFLT *fnmax )
{
    int i, j;

    *fnmax = f[0][0];
    *fnmin = *fnmax;

    for ( i = 0; i < nx; i++ )
    {
        for ( j = 0; j < ny; j++ )
        {
            *fnmax = MAX( *fnmax, f[i][j] );
            *fnmin = MIN( *fnmin, f[i][j] );
        }
    }
}

static void potential( void )
//shielded potential contour plot example.
{
    int      i, j;
    PLcGrid2 cgrid2;
    PLFLT    rmax, xmin, xmax, x0, ymin, ymax, y0, zmin, zmax;
    PLFLT    peps, xpmin, xpmax, ypmin, ypmax;
    PLFLT    eps, q1, d1, q1i, d1i, q2, d2, q2i, d2i;
    PLFLT    div1, div1i, div2, div2i;
    PLFLT    **z;
    PLINT    nlevelneg, nlevelpos;
    PLFLT    dz, clevel2, clevelneg[PNLEVEL], clevelpos[PNLEVEL];
    PLINT    ncollin, ncolbox, ncollab;
    PLFLT    px[PPERIMETERPTS], py[PPERIMETERPTS];
    PLFLT    t, r, theta;

//create data to be contoured.
    plAlloc2dGrid( &cgrid2.xg, PRPTS, PTHETAPTS );
    plAlloc2dGrid( &cgrid2.yg, PRPTS, PTHETAPTS );
    plAlloc2dGrid( &z, PRPTS, PTHETAPTS );
    cgrid2.nx = PRPTS;
    cgrid2.ny = PTHETAPTS;

    for ( i = 0; i < PRPTS; i++ )
    {
        r = 0.5 + (double) i;
        for ( j = 0; j < PTHETAPTS; j++ )
        {
            theta           = ( 2. * M_PI / (double) ( PTHETAPTS - 1 ) ) * ( 0.5 + (double) j );
            cgrid2.xg[i][j] = r * cos( theta );
            cgrid2.yg[i][j] = r * sin( theta );
        }
    }

    rmax = r;
    f2mnmx( cgrid2.xg, PRPTS, PTHETAPTS, &xmin, &xmax );
    f2mnmx( cgrid2.yg, PRPTS, PTHETAPTS, &ymin, &ymax );
    x0 = ( xmin + xmax ) / 2.;
    y0 = ( ymin + ymax ) / 2.;

    // Expanded limits
    peps  = 0.05;
    xpmin = xmin - fabs( xmin ) * peps;
    xpmax = xmax + fabs( xmax ) * peps;
    ypmin = ymin - fabs( ymin ) * peps;
    ypmax = ymax + fabs( ymax ) * peps;

    // Potential inside a conducting cylinder (or sphere) by method of images.
    // Charge 1 is placed at (d1, d1), with image charge at (d2, d2).
    // Charge 2 is placed at (d1, -d1), with image charge at (d2, -d2).
    // Also put in smoothing term at small distances.
    //

    eps = 2.;

    q1 = 1.;
    d1 = rmax / 4.;

    q1i = -q1 * rmax / d1;
    d1i = pow( rmax, 2. ) / d1;

    q2 = -1.;
    d2 = rmax / 4.;

    q2i = -q2 * rmax / d2;
    d2i = pow( rmax, 2. ) / d2;

    for ( i = 0; i < PRPTS; i++ )
    {
        for ( j = 0; j < PTHETAPTS; j++ )
        {
            div1    = sqrt( pow( cgrid2.xg[i][j] - d1, 2. ) + pow( cgrid2.yg[i][j] - d1, 2. ) + pow( eps, 2. ) );
            div1i   = sqrt( pow( cgrid2.xg[i][j] - d1i, 2. ) + pow( cgrid2.yg[i][j] - d1i, 2. ) + pow( eps, 2. ) );
            div2    = sqrt( pow( cgrid2.xg[i][j] - d2, 2. ) + pow( cgrid2.yg[i][j] + d2, 2. ) + pow( eps, 2. ) );
            div2i   = sqrt( pow( cgrid2.xg[i][j] - d2i, 2. ) + pow( cgrid2.yg[i][j] + d2i, 2. ) + pow( eps, 2. ) );
            z[i][j] = q1 / div1 + q1i / div1i + q2 / div2 + q2i / div2i;
        }
    }
    f2mnmx( z, PRPTS, PTHETAPTS, &zmin, &zmax );
//   printf("%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g \n",
//        q1, d1, q1i, d1i, q2, d2, q2i, d2i);
// printf("%.15g %.15g %.15g %.15g %.15g %.15g \n",
//        xmin, xmax, ymin, ymax, zmin, zmax);

    // Positive and negative contour levels.
    dz        = ( zmax - zmin ) / (double) PNLEVEL;
    nlevelneg = 0;
    nlevelpos = 0;
    for ( i = 0; i < PNLEVEL; i++ )
    {
        clevel2 = zmin + ( (double) i + 0.5 ) * dz;
        if ( clevel2 <= 0. )
            clevelneg[nlevelneg++] = clevel2;
        else
            clevelpos[nlevelpos++] = clevel2;
    }
    // Colours!
    ncollin = 11;
    ncolbox = 1;
    ncollab = 2;

    // Finally start plotting this page!
    pladv( 0 );
    plcol0( ncolbox );

    plvpas( 0.1, 0.9, 0.1, 0.9, 1.0 );
    plwind( xpmin, xpmax, ypmin, ypmax );
    plbox( "", 0., 0, "", 0., 0 );

    plcol0( ncollin );
    if ( nlevelneg > 0 )
    {
        // Negative contours
        pllsty( 2 );
        plcont( (const PLFLT * const *) z, PRPTS, PTHETAPTS, 1, PRPTS, 1, PTHETAPTS,
            clevelneg, nlevelneg, pltr2, (void *) &cgrid2 );
    }

    if ( nlevelpos > 0 )
    {
        // Positive contours
        pllsty( 1 );
        plcont( (const PLFLT * const *) z, PRPTS, PTHETAPTS, 1, PRPTS, 1, PTHETAPTS,
            clevelpos, nlevelpos, pltr2, (void *) &cgrid2 );
    }

    // Draw outer boundary
    for ( i = 0; i < PPERIMETERPTS; i++ )
    {
        t     = ( 2. * M_PI / ( PPERIMETERPTS - 1 ) ) * (double) i;
        px[i] = x0 + rmax*cos( t );
        py[i] = y0 + rmax*sin( t );
    }

    plcol0( ncolbox );
    plline( PPERIMETERPTS, px, py );

    plcol0( ncollab );
    pllab( "", "", "Shielded potential of charges in a conducting sphere" );

    plFree2dGrid( z, PRPTS, PTHETAPTS );
    plFree2dGrid( cgrid2.xg, PRPTS, PTHETAPTS );
    plFree2dGrid( cgrid2.yg, PRPTS, PTHETAPTS );
}


//--------------------------------------------------------------------------
// main
//
// Does several contour plots using different coordinate mappings.
//--------------------------------------------------------------------------

int
main( int argc, const char *argv[] )
{
    int          i, j;
    PLFLT        xx, yy, argx, argy, distort;
    static PLINT mark = 1500, space = 1500;

    PLFLT        **z, **w;
    PLFLT        xg1[XPTS], yg1[YPTS];
    PLcGrid      cgrid1;
    PLcGrid2     cgrid2;

// Parse and process command line arguments

    (void) plparseopts( &argc, argv, PL_PARSE_FULL );

// Initialize plplot

    plinit();

// Set up function arrays

    plAlloc2dGrid( &z, XPTS, YPTS );
    plAlloc2dGrid( &w, XPTS, YPTS );

    for ( i = 0; i < XPTS; i++ )
    {
        xx = (double) ( i - ( XPTS / 2 ) ) / (double) ( XPTS / 2 );
        for ( j = 0; j < YPTS; j++ )
        {
            yy      = (double) ( j - ( YPTS / 2 ) ) / (double) ( YPTS / 2 ) - 1.0;
            z[i][j] = xx * xx - yy * yy;
            w[i][j] = 2 * xx * yy;
        }
    }

// Set up grids

    cgrid1.xg = xg1;
    cgrid1.yg = yg1;
    cgrid1.nx = XPTS;
    cgrid1.ny = YPTS;

    plAlloc2dGrid( &cgrid2.xg, XPTS, YPTS );
    plAlloc2dGrid( &cgrid2.yg, XPTS, YPTS );
    cgrid2.nx = XPTS;
    cgrid2.ny = YPTS;

    for ( i = 0; i < XPTS; i++ )
    {
        for ( j = 0; j < YPTS; j++ )
        {
            mypltr( (PLFLT) i, (PLFLT) j, &xx, &yy, NULL );

            argx    = xx * M_PI / 2;
            argy    = yy * M_PI / 2;
            distort = 0.4;

            cgrid1.xg[i] = xx + distort * cos( argx );
            cgrid1.yg[j] = yy - distort * cos( argy );

            cgrid2.xg[i][j] = xx + distort * cos( argx ) * cos( argy );
            cgrid2.yg[i][j] = yy - distort * cos( argx ) * cos( argy );
        }
    }

// Plot using identity transform
//
//  plenv(-1.0, 1.0, -1.0, 1.0, 0, 0);
//  plcol0(2);
//  plcont(z, XPTS, YPTS, 1, XPTS, 1, YPTS, clevel, 11, mypltr, NULL);
//  plstyl(1, &mark, &space);
//  plcol0(3);
//  plcont(w, XPTS, YPTS, 1, XPTS, 1, YPTS, clevel, 11, mypltr, NULL);
//  plstyl(0, &mark, &space);
//  plcol0(1);
//  pllab("X Coordinate", "Y Coordinate", "Streamlines of flow");
//
    pl_setcontlabelformat( 4, 3 );
    pl_setcontlabelparam( 0.006, 0.3, 0.1, 1 );
    plenv( -1.0, 1.0, -1.0, 1.0, 0, 0 );
    plcol0( 2 );
    plcont( (const PLFLT * const *) z, XPTS, YPTS, 1, XPTS, 1, YPTS, clevel, 11, mypltr, NULL );
    plstyl( 1, &mark, &space );
    plcol0( 3 );
    plcont( (const PLFLT * const *) w, XPTS, YPTS, 1, XPTS, 1, YPTS, clevel, 11, mypltr, NULL );
    plstyl( 0, &mark, &space );
    plcol0( 1 );
    pllab( "X Coordinate", "Y Coordinate", "Streamlines of flow" );
    pl_setcontlabelparam( 0.006, 0.3, 0.1, 0 );

// Plot using 1d coordinate transform

    plenv( -1.0, 1.0, -1.0, 1.0, 0, 0 );
    plcol0( 2 );
    plcont( (const PLFLT * const *) z, XPTS, YPTS, 1, XPTS, 1, YPTS, clevel, 11,
        pltr1, (void *) &cgrid1 );

    plstyl( 1, &mark, &space );
    plcol0( 3 );
    plcont( (const PLFLT * const *) w, XPTS, YPTS, 1, XPTS, 1, YPTS, clevel, 11,
        pltr1, (void *) &cgrid1 );
    plstyl( 0, &mark, &space );
    plcol0( 1 );
    pllab( "X Coordinate", "Y Coordinate", "Streamlines of flow" );
    //
    // pl_setcontlabelparam(0.006, 0.3, 0.1, 1);
    // plenv(-1.0, 1.0, -1.0, 1.0, 0, 0);
    // plcol0(2);
    // plcont((const PLFLT **) z, XPTS, YPTS, 1, XPTS, 1, YPTS, clevel, 11,
    //     pltr1, (void *) &cgrid1);
    //
    // plstyl(1, &mark, &space);
    // plcol0(3);
    // plcont((const PLFLT **) w, XPTS, YPTS, 1, XPTS, 1, YPTS, clevel, 11,
    //     pltr1, (void *) &cgrid1);
    // plstyl(0, &mark, &space);
    // plcol0(1);
    // pllab("X Coordinate", "Y Coordinate", "Streamlines of flow");
    // pl_setcontlabelparam(0.006, 0.3, 0.1, 0);
    //
// Plot using 2d coordinate transform

    plenv( -1.0, 1.0, -1.0, 1.0, 0, 0 );
    plcol0( 2 );
    plcont( (const PLFLT * const *) z, XPTS, YPTS, 1, XPTS, 1, YPTS, clevel, 11,
        pltr2, (void *) &cgrid2 );

    plstyl( 1, &mark, &space );
    plcol0( 3 );
    plcont( (const PLFLT * const *) w, XPTS, YPTS, 1, XPTS, 1, YPTS, clevel, 11,
        pltr2, (void *) &cgrid2 );
    plstyl( 0, &mark, &space );
    plcol0( 1 );
    pllab( "X Coordinate", "Y Coordinate", "Streamlines of flow" );
    //
    // pl_setcontlabelparam(0.006, 0.3, 0.1, 1);
    // plenv(-1.0, 1.0, -1.0, 1.0, 0, 0);
    // plcol0(2);
    // plcont((const PLFLT **) z, XPTS, YPTS, 1, XPTS, 1, YPTS, clevel, 11,
    //     pltr2, (void *) &cgrid2);
    //
    // plstyl(1, &mark, &space);
    // plcol0(3);
    // plcont((const PLFLT **) w, XPTS, YPTS, 1, XPTS, 1, YPTS, clevel, 11,
    //     pltr2, (void *) &cgrid2);
    // plstyl(0, &mark, &space);
    // plcol0(1);
    // pllab("X Coordinate", "Y Coordinate", "Streamlines of flow");
    //
    pl_setcontlabelparam( 0.006, 0.3, 0.1, 0 );
    polar();
    //
    // pl_setcontlabelparam(0.006, 0.3, 0.1, 1);
    // polar();
    //
    pl_setcontlabelparam( 0.006, 0.3, 0.1, 0 );
    potential();
    //
    // pl_setcontlabelparam(0.006, 0.3, 0.1, 1);
    // potential();
    //

// Clean up

    plFree2dGrid( z, XPTS, YPTS );
    plFree2dGrid( w, XPTS, YPTS );
    plFree2dGrid( cgrid2.xg, XPTS, YPTS );
    plFree2dGrid( cgrid2.yg, XPTS, YPTS );

    plend();

    exit( 0 );
}
