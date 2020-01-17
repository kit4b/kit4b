// $Id: plctrl.c 12966 2014-01-29 00:01:50Z airwin $
//
//      Misc. control routines, like begin, end, exit, change graphics/text
//      mode, change color.  Includes some spillage from plcore.c.  If you
//      don't know where it should go, put it here.
//
// Copyright (C) 2004  Joao Cardoso
// Copyright (C) 2004  Rafael Laboissiere
// Copyright (C) 2008  Hazen Babcock
// Copyright (C) 2009-2013  Alan W. Irwin
// Copyright (C) 2011  Hezekiah M. Carty
//
// This file is part of PLplot.
//
// PLplot is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published
// by the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// PLplot is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with PLplot; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
//
// ----------------------------------------------------------------------------------------------------------------------------
// Modified by Stuart Stephen, March 2014, enabling loading of embedded palettes
// Contains the colour information for the standard palettes distributed with PLPlot as embedded initialised structures in "EmbeddedPalettes.h" 
// Saves having to install the original palette files and configuring the plotting application for location of these palette files
// Note that processing firstly looks for the original palette file (allows users to modify these to suit application requirements) and
// only if the specific palette file can't be loaded will then attempt to load the embedded palette.
// Please contact stuart.stephen@csiro.au for details

//! @file
//!
//! Part 1: Color map routines.
//! Part 2: "A grab-bag of various control routines".
//!
#include "stdafx.h"
#include "plplot_config.h"

#define DEBUG

#define NEED_PLDEBUG
#include "plplotP.h"
#ifdef macintosh
#include "mac.h"
// for plMacLibOpen prototype; used in plLibOpen
#endif

#ifdef DJGPP                    // dos386/djgpp
#ifdef __unix
#undef __unix
#endif
#endif

#ifdef __unix
#include <sys/types.h>
#include <sys/stat.h>
#ifdef PL_HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <errno.h>
#endif

#include "pdf.h"

// Random number generator (Mersenne Twister)
#include "mt19937ar.h"

// Embedded palette suport
#include "EmbeddedPalettes.h"

#define BUFFER_SIZE    256
#define COLLEN         30
#define PALLEN         160
#define MSGLEN         1024

// small epsilon for fuzzy range checks that is still large enough to
// work even in the single precision floating point case.
#define FUZZ_EPSILON    1.e-4

// Static functions

// Used by any external init code to suggest a path
char *plplotLibDir = 0;

static void
color_set( PLINT i, U_CHAR r, U_CHAR g, U_CHAR b, PLFLT a, const char *name );

static void
strcat_delim( char *dirspec );

static int
( *exit_handler )( const char *errormsg );

static void
( *abort_handler )( const char *errormsg );

static void
plcmap0_def( int imin, int imax );

static void
plcmap1_def( void );

static PLFLT
value( double n1, double n2, double hue );

static char *
read_line( char *buffer, int length, FILE *fp );

static void
cmap0_palette_read( const char *filename,
                    int *number_colors, unsigned int **r, unsigned int **g,
                    unsigned int **b, double **a );

// An additional hardwired location for lib files.
// I have no plans to change these again, ever.

#if defined ( DJGPP )
#ifndef PLLIBDEV
#define PLLIBDEV    "c:/plplot/lib"
#endif

#elif defined ( MSDOS )
#ifndef PLLIBDEV
#define PLLIBDEV    "c:\\plplot\\lib"
#endif

#else

// Anything else is assumed to be Unix

#ifndef PLLIBDEV
#define PLLIBDEV    "/usr/local/plplot/lib"
#endif

#endif

//--------------------------------------------------------------------------
//  Routines that deal with colors & color maps.
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// plcol0()
//
//! Set color, map 0.  Argument is a integer between 0 and m_plsc->ncol0.
//!
//! @param icol0 The index of the color map 0 color to use as the current
//! color. (0 - m_plsc->ncol0).

void
c_plcol0( PLINT icol0 )
{
    if ( m_plsc->level < 1 )
    {
        plabort( "plcol0: Please call plinit first" );
        return;
    }
    if ( icol0 < 0 || icol0 >= m_plsc->ncol0 )
    {
        char buffer[BUFFER_SIZE];
        snprintf( buffer, BUFFER_SIZE, "plcol0: Invalid color map entry: %d", (int) icol0 );
        plabort( buffer );
        return;
    }

    m_plsc->icol0      = icol0;
    m_plsc->curcolor.r = m_plsc->cmap0[icol0].r;
    m_plsc->curcolor.g = m_plsc->cmap0[icol0].g;
    m_plsc->curcolor.b = m_plsc->cmap0[icol0].b;
    m_plsc->curcolor.a = m_plsc->cmap0[icol0].a;

    m_plsc->curcmap = 0;
    plP_state( PLSTATE_COLOR0 );
}

//--------------------------------------------------------------------------
// plcol1()
//
//! Set color, map 1.  Argument is a float between 0. and 1.
//!
//! @param col1 The index of the color map 1 color to use as the current
//! color. (0.0 - 1.0)

void
c_plcol1( PLFLT col1 )
{
    PLINT icol1;

    if ( m_plsc->level < 1 )
    {
        plabort( "plcol1: Please call plinit first" );
        return;
    }
    if ( col1 < 0 || col1 > 1 || isnan( col1 ) )
    {
        char buffer[BUFFER_SIZE];
        snprintf( buffer, BUFFER_SIZE, "plcol1: Invalid color map position: %f", (PLFLT) col1 );
        plabort( buffer );
        return;
    }

    icol1 = (PLINT) ( col1 * m_plsc->ncol1 );
    icol1 = MIN( icol1, m_plsc->ncol1 - 1 );

    m_plsc->icol1      = icol1;
    m_plsc->curcolor.r = m_plsc->cmap1[m_plsc->icol1].r;
    m_plsc->curcolor.g = m_plsc->cmap1[m_plsc->icol1].g;
    m_plsc->curcolor.b = m_plsc->cmap1[m_plsc->icol1].b;
    m_plsc->curcolor.a = m_plsc->cmap1[m_plsc->icol1].a;

    m_plsc->curcmap = 1;
    plP_state( PLSTATE_COLOR1 );
}

//--------------------------------------------------------------------------
// plscolbg()
//
//! Set the background color (cmap0[0]) by 8 bit RGB value
//!
//! @param r Red value of the background color (0 - 255).
//! @param g Green value of the background color (0 - 255).
//! @param b Blue value of the background color (0 - 255).

void
c_plscolbg( PLINT r, PLINT g, PLINT b )
{
    plscol0( 0, r, g, b );
}

//--------------------------------------------------------------------------
// plscolbga()
//
//! Set the background color (cmap0[0]) by 8 bit RGB value and alpha value
//!
//! @param r Red value of the background color (0 - 255).
//! @param g Green value of the background color (0 - 255).
//! @param b Blue value of the background color (0 - 255).
//! @param alpha Alpha (transparency) value of the background color
//! (0.0 - 1.0).

//--------------------------------------------------------------------------

void
c_plscolbga( PLINT r, PLINT g, PLINT b, PLFLT alpha )
{
    plscol0a( 0, r, g, b, alpha );
}

//--------------------------------------------------------------------------
// plgcolbg()
//
//! Returns the background color (cmap0[0]) by 8 bit RGB value
//!
//! @param r Current red value of the background color.
//! @param g Current green value of the background color.
//! @param b Current blue value of the background color.

void
c_plgcolbg( PLINT *r, PLINT *g, PLINT *b )
{
    plgcol0( 0, r, g, b );
}

//--------------------------------------------------------------------------
// plgcolbga()
//
//! Returns the background color (cmap0[0]) by 8 bit RGB value and alpha value
//!
//! @param r Current red value of the background color.
//! @param g Current green value of the background color.
//! @param b Current blue value of the background color.
//! @param alpha Current alpha value of the background color.

void
c_plgcolbga( PLINT *r, PLINT *g, PLINT *b, PLFLT *alpha )
{
    plgcol0a( 0, r, g, b, alpha );
}

//--------------------------------------------------------------------------
// plscol0()
//
//! Set a given color from color map 0 by 8 bit RGB value
//! Does not result in any additional cells to be allocated.
//!
//! @param icol0 index of the color to set (0 - m_plsc->ncol0)
//! @param r Red value of the color (0 - 255).
//! @param g Green value of the color (0 - 255).
//! @param b Blue value of the color (0 - 255).

void
c_plscol0( PLINT icol0, PLINT r, PLINT g, PLINT b )
{
    if ( m_plsc->cmap0 == NULL )
        plscmap0n( 0 );
    if ( icol0 < 0 || icol0 >= m_plsc->ncol0 )
    {
        char buffer[BUFFER_SIZE];
        snprintf( buffer, BUFFER_SIZE, "plscol0: Illegal color table value: %d", (int) icol0 );
        plabort( buffer );
        return;
    }
    if ( ( r < 0 || r > 255 ) || ( g < 0 || g > 255 ) || ( b < 0 || b > 255 ) )
    {
        char buffer[BUFFER_SIZE];
        snprintf( buffer, BUFFER_SIZE, "plscol0: Invalid RGB color: %d, %d, %d",
            (int) r, (int) g, (int) b );
        plabort( buffer );
        return;
    }

    plscol0a( icol0, r, g, b, 1.0 );
}

//--------------------------------------------------------------------------
// plscol0a()
//
//! Set a given color from color map 0 by 8 bit RGB value and alpha value.
//! Does not result in any additional cells to be allocated.
//!
//! @param icol0 index of the color to set (0 - m_plsc->ncol0)
//! @param r Red value of the color (0 - 255).
//! @param g Green value of the color (0 - 255).
//! @param b Blue value of the color (0 - 255).
//! @param alpha Alpha value of the color (0.0 - 1.0).

void
c_plscol0a( PLINT icol0, PLINT r, PLINT g, PLINT b, PLFLT alpha )
{
    if ( m_plsc->cmap0 == NULL )
        plscmap0n( 0 );
    if ( icol0 < 0 || icol0 >= m_plsc->ncol0 )
    {
        char buffer[BUFFER_SIZE];
        snprintf( buffer, BUFFER_SIZE, "plscol0a: Illegal color table value: %d", (int) icol0 );
        plabort( buffer );
        return;
    }
    if ( ( r < 0 || r > 255 ) || ( g < 0 || g > 255 ) || ( b < 0 || b > 255 ) || ( alpha < 0. || alpha > 1.0 ) )
    {
        char buffer[BUFFER_SIZE];
        snprintf( buffer, BUFFER_SIZE, "plscol0a: Invalid RGB color: %d, %d, %d, %f",
            (int) r, (int) g, (int) b, (double) alpha );
        plabort( buffer );
        return;
    }

    m_plsc->cmap0[icol0].r = (unsigned char) r;
    m_plsc->cmap0[icol0].g = (unsigned char) g;
    m_plsc->cmap0[icol0].b = (unsigned char) b;
    m_plsc->cmap0[icol0].a = alpha;

    if ( m_plsc->level > 0 )
        plP_state( PLSTATE_CMAP0 );
}

//--------------------------------------------------------------------------
// plgcol0()
//
//! Returns 8 bit RGB values for given color from color map 0
//! Values are negative if an invalid color id is given
//!
//! @param icol0 Index of the color to be return (0 - m_plsc->ncol0).
//! @param r Current red value of the color.
//! @param g Current green value of the color.
//! @param b Current blue value of the color.

void
c_plgcol0( PLINT icol0, PLINT *r, PLINT *g, PLINT *b )
{
    if ( m_plsc->cmap0 == NULL )
        plscmap0n( 0 );

    *r = -1;
    *g = -1;
    *b = -1;

    if ( icol0 < 0 || icol0 > m_plsc->ncol0 )
    {
        char buffer[BUFFER_SIZE];
        snprintf( buffer, BUFFER_SIZE, "plgcol0: Invalid color index: %d", (int) icol0 );
        plabort( buffer );
        return;
    }

    *r = m_plsc->cmap0[icol0].r;
    *g = m_plsc->cmap0[icol0].g;
    *b = m_plsc->cmap0[icol0].b;

    return;
}

//--------------------------------------------------------------------------
// plgcol0a()
//
//! Returns 8 bit RGB values for given color from color map 0 and alpha value
//! Values are negative if an invalid color id is given
//!
//! @param icol0 Index of the color to be return (0 - m_plsc->ncol0).
//! @param r Current red value of the color.
//! @param g Current green value of the color.
//! @param b Current blue value of the color.
//! @param alpha Current alpha value of the color.

void
c_plgcol0a( PLINT icol0, PLINT *r, PLINT *g, PLINT *b, PLFLT *alpha )
{
    if ( m_plsc->cmap0 == NULL )
        plscmap0n( 0 );

    *r     = -1;
    *g     = -1;
    *b     = -1;
    *alpha = -1.0;

    if ( icol0 < 0 || icol0 > m_plsc->ncol0 )
    {
        char buffer[BUFFER_SIZE];
        snprintf( buffer, BUFFER_SIZE, "plgcol0: Invalid color index: %d", (int) icol0 );
        plabort( buffer );
        return;
    }

    *r     = m_plsc->cmap0[icol0].r;
    *g     = m_plsc->cmap0[icol0].g;
    *b     = m_plsc->cmap0[icol0].b;
    *alpha = m_plsc->cmap0[icol0].a;

    return;
}

//--------------------------------------------------------------------------
// plscmap0()
//
//! Set color map 0 colors by 8 bit RGB values.  This sets the entire color
//! map -- only as many colors as specified will be allocated.
//!
//! @param r Array of red values.
//! @param g Array of green values.
//! @param b Array of blue values.
//! @param ncol0 Total number of RGB values.

void
c_plscmap0( const PLINT *r, const PLINT *g, const PLINT *b, PLINT ncol0 )
{
    int i;

    plscmap0n( ncol0 );

    for ( i = 0; i < m_plsc->ncol0; i++ )
    {
        if ( ( r[i] < 0 || r[i] > 255 ) ||
             ( g[i] < 0 || g[i] > 255 ) ||
             ( b[i] < 0 || b[i] > 255 ) )
        {
            char buffer[BUFFER_SIZE];
            snprintf( buffer, BUFFER_SIZE, "plscmap0: Invalid RGB color: %d, %d, %d",
                (int) r[i], (int) g[i], (int) b[i] );
            plabort( buffer );
            return;
        }

        m_plsc->cmap0[i].r = (unsigned char) r[i];
        m_plsc->cmap0[i].g = (unsigned char) g[i];
        m_plsc->cmap0[i].b = (unsigned char) b[i];
        m_plsc->cmap0[i].a = 1.0;
    }

    if ( m_plsc->level > 0 )
        plP_state( PLSTATE_CMAP0 );
}

//--------------------------------------------------------------------------
// plscmap0a()
//
//! Set color map 0 colors by 8 bit RGB and alpha value.  This sets the
//! entire color map -- only as many colors as specified will be allocated.
//!
//! @param r Array of red values.
//! @param g Array of green values.
//! @param b Array of blue values.
//! @param alpha Array of alpha values.
//! @param ncol0 Total number of RGBA values.

void
c_plscmap0a( const PLINT *r, const PLINT *g, const PLINT *b, const PLFLT *alpha, PLINT ncol0 )
{
    int i;

    plscmap0n( ncol0 );

    for ( i = 0; i < m_plsc->ncol0; i++ )
    {
        if ( ( r[i] < 0 || r[i] > 255 ) ||
             ( g[i] < 0 || g[i] > 255 ) ||
             ( b[i] < 0 || b[i] > 255 ) ||
             ( alpha[i] < 0.0 || alpha[i] > 1.0 ) )
        {
            char buffer[BUFFER_SIZE];
            snprintf( buffer, BUFFER_SIZE, "plscmap0a: Invalid RGB color: %d, %d, %d, %f",
                (int) r[i], (int) g[i], (int) b[i], (double) alpha[i] );
            plabort( buffer );
            return;
        }

        m_plsc->cmap0[i].r = (unsigned char) r[i];
        m_plsc->cmap0[i].g = (unsigned char) g[i];
        m_plsc->cmap0[i].b = (unsigned char) b[i];
        m_plsc->cmap0[i].a = alpha[i];
    }

    if ( m_plsc->level > 0 )
        plP_state( PLSTATE_CMAP0 );
}

//--------------------------------------------------------------------------
// plscmap1()
//
//! Set color map 1 colors by 8 bit RGB values
//! This also sets the number of colors.
//!
//! @param r Array of red values.
//! @param g Array of green values.
//! @param b Array of blue values.
//! @param ncol1 Total number of RGB values.

void
c_plscmap1( const PLINT *r, const PLINT *g, const PLINT *b, PLINT ncol1 )
{
    int i;

    plscmap1n( ncol1 );

    for ( i = 0; i < m_plsc->ncol1; i++ )
    {
        if ( ( r[i] < 0 || r[i] > 255 ) ||
             ( g[i] < 0 || g[i] > 255 ) ||
             ( b[i] < 0 || b[i] > 255 ) )
        {
            char buffer[BUFFER_SIZE];
            snprintf( buffer, BUFFER_SIZE, "plscmap1: Invalid RGB color: %d, %d, %d",
                (int) r[i], (int) g[i], (int) b[i] );
            plabort( buffer );
            return;
        }
        m_plsc->cmap1[i].r = (unsigned char) r[i];
        m_plsc->cmap1[i].g = (unsigned char) g[i];
        m_plsc->cmap1[i].b = (unsigned char) b[i];
        m_plsc->cmap1[i].a = 1.0;
    }

    if ( m_plsc->level > 0 )
        plP_state( PLSTATE_CMAP1 );
}

//--------------------------------------------------------------------------
// plscmap1a()
//
//! Set color map 1 colors by 8 bit RGB and alpha values
//! This also sets the number of colors.
//!
//! @param r Array of red values.
//! @param g Array of green values.
//! @param b Array of blue values.
//! @param alpha Array of alpha values.
//! @param ncol1 Total number of RGBA values.

void
c_plscmap1a( const PLINT *r, const PLINT *g, const PLINT *b, const PLFLT *alpha, PLINT ncol1 )
{
    int i;

    plscmap1n( ncol1 );

    for ( i = 0; i < m_plsc->ncol1; i++ )
    {
        if ( ( r[i] < 0 || r[i] > 255 ) ||
             ( g[i] < 0 || g[i] > 255 ) ||
             ( b[i] < 0 || b[i] > 255 ) ||
             ( alpha[i] < 0.0 || alpha[i] > 1.0 ) )
        {
            char buffer[BUFFER_SIZE];
            snprintf( buffer, BUFFER_SIZE, "plscmap1a: Invalid RGB color: %d, %d, %d, %f",
                (int) r[i], (int) g[i], (int) b[i], (double) alpha[i] );
            plabort( buffer );
            return;
        }
        m_plsc->cmap1[i].r = (unsigned char) r[i];
        m_plsc->cmap1[i].g = (unsigned char) g[i];
        m_plsc->cmap1[i].b = (unsigned char) b[i];
        m_plsc->cmap1[i].a = alpha[i];
    }

    if ( m_plsc->level > 0 )
        plP_state( PLSTATE_CMAP1 );
}

//--------------------------------------------------------------------------
// plscmap1l()
//
//! Set color map 1 colors using a piece-wise linear relationship between
//! position in the color map (from 0 to 1) and position in HLS or RGB color
//! space.  May be called at any time.
//!
//! The idea here is to specify a number of control points that specify the
//! mapping between HLS (or RGB or CMY) and palette 1 value.  Between these
//! points, linear interpolation is used.  By mapping position in the color
//! map to function value, this gives a smooth variation of color with
//! intensity.  Any number of control points may be specified, located at
//! arbitrary positions (intensities), although typically 2 - 4 are enough.
//! Another way of stating this is that we are traversing a given number of
//! lines through HLS (or RGB) space as we move through cmap 1 entries.  The
//! control points at the minimum and maximum intensity (0 and 1) must
//! always be specified.  By adding more control points you can get more
//! variation.  One good technique for plotting functions that vary about
//! some expected average is to use an additional 2 control points in the
//! center (intensity ~= 0.5) that are the same color as the background
//! (typically white for paper output, black for crt), and same hue as the
//! boundary control points.  This allows the highs and lows to be very
//! easily distinguished.
//!
//! Each control point must specify the position in cmap 1 as well as three
//! coordinates in HLS or RGB space.  The first point MUST correspond to
//! position = 0, and the last to position = 1.
//!
//! Every change in hue from one control point to the next can be linearly
//! interpolated in two ways.  The usual (alt_hue_path[i] false) method for the ith interval
//! uses the dh = h[i+1] - h[i] interval for interpolation.  The alternate (alt_hue_path true) method for the ith interval uses the dh = (h[i+1] - h[i]) - 360 if (h[i+1] - h[i]) is positive or dh = 360 - (h[i+1] - h[i]) if (h[i+1] - h[i]) is negative interval for the interpolation.  Thus, alt_hue_path true interpolation intervals always include hue = 0.
//! Specifying
//! alt_hue_path=NULL is equivalent to setting alt_hue_path[]=false for every control point.
//!
//! Bounds on RGB coordinates:
//!	R,G,B		[0, 1]		magnitude
//!
//! Bounds on HLS coordinates:
//!	hue		[0, 360]	degrees
//!	lightness	[0, 1]		magnitude
//!	saturation	[0, 1]		magnitude
//!
//! The inputs are:
//! @param itype 0: HLS, 1: RGB
//! @param npts	number of control points
//! @param intensity[] intensity index for each control point
//! @param coord1[] first coordinate for each control point
//! @param coord2[] second coordinate for each control point
//! @param coord3[] third coordinate for each control point
//! @param alt_hue_path[] if true, use alternative hue interpolation path
//! for the associated interval.

void
c_plscmap1l( PLINT itype, PLINT npts, const PLFLT *intensity,
             const PLFLT *coord1, const PLFLT *coord2, const PLFLT *coord3, const PLINT *alt_hue_path )
{
    int   n;
    PLFLT h, l, s, r, g, b;

    if ( npts < 2 )
    {
        plabort( "plscmap1l: Must specify at least two control points" );
        return;
    }

    if ( ( intensity[0] != 0 ) || ( intensity[npts - 1] != 1 ) )
    {
        plabort( "plscmap1l: First, last control points must lie on boundary" );
        return;
    }

    if ( npts > PL_MAX_CMAP1CP )
    {
        plabort( "plscmap1l: exceeded maximum number of control points" );
        return;
    }

// Allocate if not done yet

    if ( m_plsc->cmap1 == NULL )
        plscmap1n( 0 );

// Save control points

    m_plsc->ncp1 = npts;

    for ( n = 0; n < npts; n++ )
    {
        if ( itype == 0 )
        {
            h = coord1[n];
            l = coord2[n];
            s = coord3[n];
        }
        else
        {
            r = coord1[n];
            g = coord2[n];
            b = coord3[n];
            c_plrgbhls( r, g, b, &h, &l, &s );
        }

        m_plsc->cmap1cp[n].h = h;
        m_plsc->cmap1cp[n].l = l;
        m_plsc->cmap1cp[n].s = s;
        m_plsc->cmap1cp[n].p = intensity[n];
        m_plsc->cmap1cp[n].a = 1.0;

        if ( alt_hue_path == NULL )
            m_plsc->cmap1cp[n].alt_hue_path = 0;
        else
        if ( n != npts - 1 )
            m_plsc->cmap1cp[n].alt_hue_path = alt_hue_path[n];
        else
            // Note final element is unused, so we set to zero for completeness.
            m_plsc->cmap1cp[n].alt_hue_path = 0;
    }

// Calculate and set color map

    plcmap1_calc();
}

//--------------------------------------------------------------------------
// plscmap1la()
//
//! This is the same as plscmap1l, but also allows alpha value interpolation.
//!
//! @param itype 0: HLS, 1: RGB
//! @param npts	number of control points
//! @param intensity[] intensity index for each control point
//! @param coord1[] first coordinate for each control point
//! @param coord2[] second coordinate for each control point
//! @param coord3[] third coordinate for each control point
//! @param alpha[] alpha value for each control point
//! @param alt_hue_path[] if true, use alternative hue interpolation path
//! for the associated interval.

void
c_plscmap1la( PLINT itype, PLINT npts, const PLFLT *intensity,
              const PLFLT *coord1, const PLFLT *coord2, const PLFLT *coord3, const PLFLT *alpha, const PLINT *alt_hue_path )
{
    int   n;
    PLFLT h, l, s, r, g, b;

    if ( npts < 2 )
    {
        plabort( "plscmap1la: Must specify at least two control points" );
        return;
    }

    if ( ( intensity[0] != 0 ) || ( intensity[npts - 1] != 1 ) )
    {
        plabort( "plscmap1la: First, last control points must lie on boundary" );
        return;
    }

    if ( npts > PL_MAX_CMAP1CP )
    {
        plabort( "plscmap1la: exceeded maximum number of control points" );
        return;
    }

// Allocate if not done yet

    if ( m_plsc->cmap1 == NULL )
        plscmap1n( 0 );

// Save control points

    m_plsc->ncp1 = npts;

    for ( n = 0; n < npts; n++ )
    {
        if ( itype == 0 )
        {
            h = coord1[n];
            l = coord2[n];
            s = coord3[n];
        }
        else
        {
            r = coord1[n];
            g = coord2[n];
            b = coord3[n];
            c_plrgbhls( r, g, b, &h, &l, &s );
        }

        m_plsc->cmap1cp[n].h = h;
        m_plsc->cmap1cp[n].l = l;
        m_plsc->cmap1cp[n].s = s;
        m_plsc->cmap1cp[n].p = intensity[n];
        m_plsc->cmap1cp[n].a = alpha[n];

        if ( alt_hue_path == NULL )
            m_plsc->cmap1cp[n].alt_hue_path = 0;
        else
        if ( n != npts - 1 )
            m_plsc->cmap1cp[n].alt_hue_path = alt_hue_path[n];
        else
            // Note final element is unused, so we set to zero for completeness.
            m_plsc->cmap1cp[n].alt_hue_path = 0;
    }

// Calculate and set color map

    plcmap1_calc();
}

//--------------------------------------------------------------------------
// plcmap1_calc()
//
//! Bin up cmap 1 space and assign colors to make inverse mapping easy.
//! Always do interpolation in HLS space.

void
plcmap1_calc( void )
{
    int   i, n;
    PLFLT delta, dp, dh, dl, ds, da;
    PLFLT h, l, s, p, r, g, b, a;

// Loop over all control point pairs

    for ( n = 0; n < m_plsc->ncp1 - 1; n++ )
    {
        if ( m_plsc->cmap1cp[n].p == m_plsc->cmap1cp[n + 1].p )
            continue;

        // Differences in p, h, l, s between ctrl pts

        dp = m_plsc->cmap1cp[n + 1].p - m_plsc->cmap1cp[n].p;
        dh = m_plsc->cmap1cp[n + 1].h - m_plsc->cmap1cp[n].h;
        dl = m_plsc->cmap1cp[n + 1].l - m_plsc->cmap1cp[n].l;
        ds = m_plsc->cmap1cp[n + 1].s - m_plsc->cmap1cp[n].s;
        da = m_plsc->cmap1cp[n + 1].a - m_plsc->cmap1cp[n].a;

        // Adjust dh if we are to go around "the back side"

        if ( m_plsc->cmap1cp[n].alt_hue_path )
            dh = ( dh > 0 ) ? dh - 360 : dh + 360;

        // Loop over all color cells.  Only interested in cells located (in
        // cmap1 space)  between n_th and n+1_th control points

        for ( i = 0; i < m_plsc->ncol1; i++ )
        {
            p = (double) i / ( m_plsc->ncol1 - 1.0 );
            if ( ( p < m_plsc->cmap1cp[n].p ) ||
                 ( p > m_plsc->cmap1cp[n + 1].p ) )
                continue;

            // Interpolate based on position of color cell in cmap1 space

            delta = ( p - m_plsc->cmap1cp[n].p ) / dp;

            // Linearly interpolate to get color cell h, l, s values

            h = m_plsc->cmap1cp[n].h + dh * delta;
            l = m_plsc->cmap1cp[n].l + dl * delta;
            s = m_plsc->cmap1cp[n].s + ds * delta;
            a = m_plsc->cmap1cp[n].a + da * delta;

            while ( h >= 360. )
                h -= 360.;

            while ( h < 0. )
                h += 360.;

            c_plhlsrgb( h, l, s, &r, &g, &b );

            m_plsc->cmap1[i].r = (unsigned char) MAX( 0, MIN( 255, (int) ( 256. * r ) ) );
            m_plsc->cmap1[i].g = (unsigned char) MAX( 0, MIN( 255, (int) ( 256. * g ) ) );
            m_plsc->cmap1[i].b = (unsigned char) MAX( 0, MIN( 255, (int) ( 256. * b ) ) );
            m_plsc->cmap1[i].a = a;
        }
    }

    if ( m_plsc->level > 0 )
        plP_state( PLSTATE_CMAP1 );
}

//--------------------------------------------------------------------------
//! Set the color map 1 value range to use in continuous color plots.
//!
//! @param min_color Specifies the minimum color to use.  A value of 0.0 or
//! less indicates that the range should start at the lowest color map 1
//! value available.
//! @param max_color Specifies the maximum color to use.  A value of 1.0 or
//! greater indicates that the range should exten to the highest color map 1
//! value available.
//!
//! If min_color > max_color or min_color is greater than 1.0 or max_color is
//! less than 0.0 then no change is made.
//--------------------------------------------------------------------------

void
c_plscmap1_range( PLFLT min_color, PLFLT max_color )
{
    if ( min_color > max_color || max_color < 0.0 || min_color > 1.0 )
    {
        plwarn( "plscmap1_range called with invalid color range" );
        return;
    }
    if ( min_color < 0.0 )
    {
        plwarn( "plscmap1_range called with a negative minimum color value" );
        min_color = 0.0;
    }
    if ( max_color > 1.0 )
    {
        plwarn( "plscmap1_range called with an out of range maximum color value" );
        max_color = 1.0;
    }
    m_plsc->cmap1_min = min_color;
    m_plsc->cmap1_max = max_color;
}

//--------------------------------------------------------------------------
//! Get the color map 1 value range used in continuous color plots.
//!
//! @param min_color Specifies the minimum color used.
//! @param max_color Specifies the maximum color used.
//--------------------------------------------------------------------------

void
c_plgcmap1_range( PLFLT *min_color, PLFLT *max_color )
{
    *min_color = m_plsc->cmap1_min;
    *max_color = m_plsc->cmap1_max;
}

//--------------------------------------------------------------------------
// plscmap0n()
//
//! Set number of colors in cmap 0, (re-)allocate cmap 0, and fill with
//! default values for those colors not previously allocated (and less
//! than index 15, after that you just get grey).
//!
//! The driver is not guaranteed to support all of these.
//!
//! @param ncol0 Total number of colors.

void
c_plscmap0n( PLINT ncol0 )
{
    int ncol, size, imin, imax;

// No change

    if ( ncol0 > 0 && m_plsc->ncol0 == ncol0 )
        return;

// Handle all possible startup conditions

    if ( m_plsc->ncol0 <= 0 && ncol0 <= 0 )
        ncol = 16;
    else if ( ncol0 <= 0 )
        ncol = m_plsc->ncol0;
    else
        ncol = ncol0;

    imax = ncol - 1;
    size = ncol * (int) sizeof ( PLColor );

// Allocate the space

    if ( m_plsc->cmap0 == NULL )
    {
        if ( ( m_plsc->cmap0 = (PLColor *) calloc( 1, (size_t) size ) ) == NULL )
        {
            plexit( "c_plscmap0n: Insufficient memory" );
        }
        imin = 0;
    }
    else
    {
        if ( ( m_plsc->cmap0 = (PLColor *) realloc( m_plsc->cmap0, (size_t) size ) ) == NULL )
        {
            plexit( "c_plscmap0n: Insufficient memory" );
        }
        imin = m_plsc->ncol0;
    }

// Fill in default entries

    m_plsc->ncol0 = ncol;
    plcmap0_def( imin, imax );

    if ( m_plsc->level > 0 )
        plP_state( PLSTATE_CMAP0 );
}

//--------------------------------------------------------------------------
// color_set()
//
//! Initializes color table 0 entry by RGB values.
//!
//! @param i Index of the color.
//! @param r Red value of the color.
//! @param g Green value of the color.
//! @param b Blue value of the color.
//! @param a Alpha value of the color.
//! @param name The name of the color.

void
color_set( PLINT i, U_CHAR r, U_CHAR g, U_CHAR b, PLFLT a, const char *name )
{
    m_plsc->cmap0[i].r    = r;
    m_plsc->cmap0[i].g    = g;
    m_plsc->cmap0[i].b    = b;
    m_plsc->cmap0[i].a    = a;
    m_plsc->cmap0[i].name = name;
}

#define color_def( i, r, g, b, a, n ) \
    if ( i >= imin && i <= imax ) color_set( i, r, g, b, a, n );

//--------------------------------------------------------------------------
// plcmap0_def()
//
//! Initializes specified color map 0 color entry to its default for
//! index range from imin to imax.
//!
//! @param imin Index of the first color to set to its default.
//! @param imax Index of the last color to set to its default.

void
plcmap0_def( int imin, int imax )
{
    int          i;
    unsigned int *r, *g, *b;
    double       *a;
    int          number_colors;
    if ( imin <= imax )
    {
        cmap0_palette_read( "", &number_colors, &r, &g, &b, &a );
        for ( i = imin; i <= MIN( ( number_colors - 1 ), imax ); i++ )
            color_def( i, (U_CHAR) r[i], (U_CHAR) g[i], (U_CHAR) b[i], a[i],
                "colors defined by default cmap0 palette file" );
        free( r );
        free( g );
        free( b );
        free( a );
    }
    else
    {
        number_colors = 0;
    }

    // Initialize all colours undefined by the default colour palette file
    // to opaque red as a warning.
    for ( i = MAX( number_colors, imin ); i <= imax; i++ )
        color_def( i, 255, 0, 0, 1.0,
            "opaque red colour to mark not defined by palette file" );
}

//--------------------------------------------------------------------------
// plscmap1n()
//
//! Set number of colors in cmap 1, (re-)allocate cmap 1, and set default
//! values if this is the first allocation.
//!
//! Note that the driver is allowed to disregard this number.
//! In particular, most use fewer than we use internally.
//!
//! @param ncol1 The number of colors in cmap1.

void
c_plscmap1n( PLINT ncol1 )
{
    PLINT  ncol;
    size_t size;

// No change

    if ( ncol1 > 0 && m_plsc->ncol1 == ncol1 )
        return;

// Handle all possible startup conditions

    if ( m_plsc->ncol1 <= 0 && ncol1 <= 0 )
        ncol = 128;
    else if ( ncol1 <= 0 )
        ncol = m_plsc->ncol1;
    else
        ncol = ncol1;

    size = (size_t) ncol * sizeof ( PLColor );

// Allocate the space

    if ( m_plsc->ncol1 > 0 )
    {
        if ( ( m_plsc->cmap1 = (PLColor *) realloc( m_plsc->cmap1, size ) ) == NULL )
        {
            plexit( "c_plscmap1n: Insufficient memory" );
        }
    }
    else
    {
        if ( ( m_plsc->cmap1 = (PLColor *) calloc( (size_t) ncol, sizeof ( PLColor ) ) ) == NULL )
        {
            plexit( "c_plscmap1n: Insufficient memory" );
        }
    }

// Fill in default entries

    m_plsc->ncol1 = ncol;
    if ( m_plsc->ncp1 == 0 )
        plcmap1_def();
    else
        plcmap1_calc();
}

//--------------------------------------------------------------------------
// plcmap1_def()
//
//! Initializes color map 1.
//!
//! The default initialization uses 6 control points in HLS space, the inner
//! ones being very close to one of the vertices of the HLS double cone.  The
//! vertex used (black or white) is chosen to be the closer to the background
//! color.  The 6 points were chosen over the older 4 points in order to make
//! weaker structures more easily visible, and give more control through the
//! palette editor.  If you don't like these settings.. change them!
//!

void
plcmap1_def( void )
{
    PLFLT i[6], h[6], l[6], s[6], midpt = 0., vertex = 0.;

// Positions of control points

    i[0] = 0;           // left boundary
    i[1] = 0.44;        // a little left of center
    i[2] = 0.50;        // at center
    i[3] = 0.50;        // at center
    i[4] = 0.56;        // a little right of center
    i[5] = 1;           // right boundary

// For center control points, pick black or white, whichever is closer to bg
// Be careful to pick just short of top or bottom else hue info is lost

    if ( m_plsc->cmap0 != NULL )
        vertex = ( (PLFLT) m_plsc->cmap0[0].r +
                   (PLFLT) m_plsc->cmap0[0].g +
                   (PLFLT) m_plsc->cmap0[0].b ) / 3. / 255.;

    if ( vertex < 0.5 )
    {
        vertex = 0.01;
        midpt  = 0.10;
    }
    else
    {
        vertex = 0.99;
        midpt  = 0.90;
    }

// Set hue

    h[0] = 260;         // low: blue-violet
    h[1] = 260;         // only change as we go over vertex
    h[2] = 260;         // only change as we go over vertex
    h[3] = 0;           // high: red
    h[4] = 0;           // high: red
    h[5] = 0;           // keep fixed

// Set lightness

    l[0] = 0.5;         // low
    l[1] = midpt;       // midpoint value
    l[2] = vertex;      // bg
    l[3] = vertex;      // bg
    l[4] = midpt;       // midpoint value
    l[5] = 0.5;         // high

// Set saturation -- keep at maximum

    s[0] = 1;
    s[1] = 1;
    s[2] = 1;
    s[3] = 1;
    s[4] = 1;
    s[5] = 1;

    c_plscmap1l( 0, 6, i, h, l, s, NULL );

    if ( m_plsc->level > 0 )
        plP_state( PLSTATE_CMAP1 );
}

//--------------------------------------------------------------------------
// plscolor()
//
//! Used to globally turn color output on/off
//!
//! @param color 0 = no color, Not zero = color.
//--------------------------------------------------------------------------

void
c_plscolor( PLINT color )
{
    m_plsc->colorset = 1;
    m_plsc->color    = color;
}

//--------------------------------------------------------------------------
// void value()
//
//! Auxiliary function used by c_plhlsrgb().
//!
//! @param n1 Lightness/saturation value 1.
//! @param n2 Lightness/saturation value 2.
//! @param hue hue (0.0 - 360.0).
//--------------------------------------------------------------------------

PLFLT
value( double n1, double n2, double hue )
{
    PLFLT val;

    while ( hue >= 360. )
        hue -= 360.;
    while ( hue < 0. )
        hue += 360.;

    if ( hue < 60. )
        val = n1 + ( n2 - n1 ) * hue / 60.;
    else if ( hue < 180. )
        val = n2;
    else if ( hue < 240. )
        val = n1 + ( n2 - n1 ) * ( 240. - hue ) / 60.;
    else
        val = n1;

    return ( val );
}

//--------------------------------------------------------------------------
// void c_plhlsrgb()
//
//! Convert HLS color to RGB color.
//! Bounds on HLS (input):
//!	hue		[0., 360.]	degrees
//!	lightness	[0., 1.]	magnitude
//!	saturation	[0., 1.]	magnitude
//!
//! Hue is always mapped onto the interval [0., 360.] regardless of input.
//! Bounds on RGB (output) is always [0., 1.].  Convert to RGB color values
//! by multiplying by 2**nbits (nbits typically 8).
//!
//! @param h hue in HLS color scheme (0.0 - 360.0)
//! @param l lightness in HLS color scheme (0.0 - 1.0)
//! @param s saturation in HLS color scheme (0.0 - 1.0)
//! @param p_r red value of the HLS color
//! @param p_g green value of the HLS color
//! @param p_b blue value of the HLS color

void
c_plhlsrgb( PLFLT h, PLFLT l, PLFLT s, PLFLT *p_r, PLFLT *p_g, PLFLT *p_b )
{
    PLFLT m1, m2;

    if ( l <= .5 )
        m2 = l * ( s + 1. );
    else
        m2 = l + s - l * s;

    m1 = 2 * l - m2;

    *p_r = value( m1, m2, h + 120. );
    *p_g = value( m1, m2, h );
    *p_b = value( m1, m2, h - 120. );
}

//--------------------------------------------------------------------------
// void c_plrgbhls()
//
//! Convert RGB color to HLS color.
//! Bounds on RGB (input) is always [0., 1.].
//! Bounds on HLS (output):
//!	hue		[0., 360.]	degrees
//!	lightness	[0., 1.]	magnitude
//!	saturation	[0., 1.]	magnitude
//! @param r red in RGB color scheme (0.0 - 1.0)
//! @param g green in RGB color scheme (0.0 - 1.0)
//! @param b blue in RGB color scheme (0.0 - 1.0)
//! @param p_h hue value of the RGB color.
//! @param p_l lightness value of the RGB color.
//! @param p_s saturation value of the RGB color.

void
c_plrgbhls( PLFLT r, PLFLT g, PLFLT b, PLFLT *p_h, PLFLT *p_l, PLFLT *p_s )
{
    PLFLT h, l, s, d, rc, gc, bc, rgb_min, rgb_max;

    rgb_min = MIN( r, MIN( g, b ) );
    rgb_max = MAX( r, MAX( g, b ) );

    l = ( rgb_min + rgb_max ) / 2.0;

    if ( rgb_min == rgb_max )
    {
        s = 0;
        h = 0;
    }
    else
    {
        d = rgb_max - rgb_min;
        if ( l < 0.5 )
            s = 0.5 * d / l;
        else
            s = 0.5 * d / ( 1. - l );

        rc = ( rgb_max - r ) / d;
        gc = ( rgb_max - g ) / d;
        bc = ( rgb_max - b ) / d;

        if ( r == rgb_max )
            h = bc - gc;
        else if ( g == rgb_max )
            h = rc - bc + 2;
        else
            h = gc - rc - 2;

        h = h * 60;
        if ( h < 0 )
            h = h + 360;
        else if ( h >= 360 )
            h = h - 360;
    }
    *p_h = h;
    *p_l = l;
    *p_s = s;
}

//--------------------------------------------------------------------------
// read_line()
//
//! Read a complete line and fill the buffer with its contents up to
//! capacity. Then sanitize the string - no control characters, no
//! trailing blanks
//!
//! @param buffer storage for the line of text.
//! @param length size of the buffer.
//! @param fp open file pointer from which the line of text should be read.
//!
//! @returns The sanitized line from the file.

static char *
read_line( char *buffer, int length, FILE *fp )
{
    char *pchr;

    // Read the string
    if ( fgets( buffer, length, fp ) == NULL )
    {
        return NULL;
    }

    // Sanitize the string we read - it may contain EOL characters
    // Make sure file reading starts at the next line
    pchr = strchr( buffer, '\n' );
    if ( pchr != NULL )
    {
        *pchr = '\0';
    }
    else
    {
        if ( fscanf( fp, "%*[^\n]\n" ) == EOF && ferror( fp ) )
        {
            return NULL;
        }
    }


    pchr = strchr( buffer, '\r' );
    if ( pchr != NULL )
    {
        *pchr = '\0';
    }

    // Remove trailing blanks
    pchr = buffer + strlen( buffer ) - 1;
    while ( pchr != buffer && *pchr == ' ' )
    {
        *pchr = '\0';
        pchr--;
    }

    return buffer;
}

//--------------------------------------------------------------------------
// cmap0_palette_read()
//
//! Read and check r, g, b, a data from a cmap0*.pal format file.
//! The caller must free the returned malloc'ed space for r, g, b, and a.
//!
//! @param filename name of the cmap0 palette file.
//! @param number_colors number of color found in the palette file.
//! @param r red value of each color in the palette file.
//! @param g green value of each color in the palette file.
//! @param b blue value of each color in the palette file.
//! @param a alpha value of each color in the palette file.


tsEmbeddedCol0Palette *			// returned ptr to embedded palette, NULL if unable to locate
GetEmbeddedCol0Palette(const char *pszPaletteName)	// which embedded palette to load - name is same as original file name
{
int Idx;
tsEmbeddedCol0Palette *pSelPalette;
pSelPalette = &EmbeddedCol0Palettes[0];
for(Idx = 0; Idx < cNumEmbeddedCol0Pallets; Idx++,pSelPalette+=1)
	if(!stricmp(pSelPalette->pszPaletteName,pszPaletteName))
		return(pSelPalette);
return(NULL);
}

void
cmap0_palette_read( const char *filename,
                    int *number_colors, unsigned int **r, unsigned int **g, unsigned int **b, double **a )
{
tsEmbeddedCol0Palette *pEmbeddedPalette = NULL;
    int  i, err = 0;
    char color_info[COLLEN];
    char msgbuf[MSGLEN];
    FILE *fp;
    char * save_locale = plsave_set_locale();

    if ( strlen( filename ) == 0 )
		{
        fp = plLibOpen( PL_DEFAULT_CMAP0_FILE );
        if ( fp == NULL )
			{
			if((pEmbeddedPalette = GetEmbeddedCol0Palette(PL_DEFAULT_CMAP0_FILE)) == NULL)
				{
				snprintf( msgbuf, MSGLEN, "Unable to open cmap0 file %s\n", PL_DEFAULT_CMAP0_FILE );
				plwarn( msgbuf );
				err = 1;
				}
			}
		}
    else
		{
        fp = plLibOpen( filename );
        if ( fp == NULL )
			{
			if((pEmbeddedPalette = GetEmbeddedCol0Palette(filename)) == NULL)
				{
				snprintf( msgbuf, MSGLEN, "Unable to open cmap0 file %s\n", filename );
				plwarn( msgbuf );
				err = 1;
				}
			}
		}

	if(pEmbeddedPalette != NULL)
		*number_colors = pEmbeddedPalette->NumColors;
	else
		{
	    if ( !err && ( fscanf( fp, "%d\n", number_colors ) != 1 || *number_colors < 1 ) )
			{
			fclose( fp );
			snprintf( msgbuf, MSGLEN, "Unrecognized cmap0 header\n" );
			plwarn( msgbuf );
			err = 1;
			}
		}
	

    if ( !err )
		{
        // Allocate arrays to hold r, g, b, and a data for calling routine.
        // The caller must free these after it is finished with them.
        if ( ( ( *r = (unsigned int *) malloc( (size_t) ( *number_colors ) * sizeof ( unsigned int ) ) ) == NULL ) ||
             ( ( *g = (unsigned int *) malloc( (size_t) ( *number_colors ) * sizeof ( unsigned int ) ) ) == NULL ) ||
             ( ( *b = (unsigned int *) malloc( (size_t) ( *number_colors ) * sizeof ( unsigned int ) ) ) == NULL ) ||
             ( ( *a = (double *) malloc( (size_t) ( *number_colors ) * sizeof ( double ) ) ) == NULL ) )
			{
			if(fp != NULL)
				fclose( fp );
            plexit( "cmap0_palette_read: insufficient memory" );
			}

		if(pEmbeddedPalette != NULL)
			{
			tsCol0ColorInfo *pColorInfo;

			pColorInfo = &pEmbeddedPalette->ColorInfo[0];
			for(i = 0; i < *number_colors; i+=1, pColorInfo += 1)
				{
				*(*r+i) = pColorInfo->Red;
				*(*b+i) = pColorInfo->Blue;
				*(*g+i) = pColorInfo->Green;
				*(*a+i) = pColorInfo->Alpha;
				}
			plrestore_locale( save_locale );
			return;
			}

        for ( i = 0; i < *number_colors; i++ )
			{
            if ( read_line( color_info, COLLEN, fp ) == NULL )
				{
                err = 1;
                break;
				}

            // Get the color data
            if ( strlen( color_info ) == 7 )
            {
                if ( sscanf( color_info, "#%2x%2x%2x",
                         (unsigned int *) ( *r + i ), (unsigned int *) ( *g + i ),
                         (unsigned int *) ( *b + i ) ) != 3 )
                {
                    err = 1;
                    break;
                }
                *( *a + i ) = 1.0;
            }
            else if ( strlen( color_info ) > 9 )
            {
                if ( sscanf( color_info, "#%2x%2x%2x %lf",
                         (unsigned int *) ( *r + i ), (unsigned int *) ( *g + i ),
                         (unsigned int *) ( *b + i ), (double *) ( *a + i ) ) != 4 )
                {
                    err = 1;
                    break;
                }
                // fuzzy range check.
                if ( *( *a + i ) < -FUZZ_EPSILON || *( *a + i ) > ( 1. + FUZZ_EPSILON ) )
                {
                    err = 1;
                    break;
                }
                else if ( *( *a + i ) < 0. )
                {
                    *( *a + i ) = 0.;
                }
                else if ( *( *a + i ) > 1. )
                {
                    *( *a + i ) = 1.;
                }
            }
            else
            {
                err = 1;
                break;
            }
        }
        fclose( fp );
        if ( err )
        {
            snprintf( msgbuf, MSGLEN, "Unrecognized cmap0 format data line.  Line is %s\n",
                color_info );
            plwarn( msgbuf );
            free( *r );
            free( *g );
            free( *b );
            free( *a );
        }
    }
    // Fall back to opaque red on opaque white as visual warning of any
    // error above.
    if ( err )
    {
        *number_colors = 16;
        if ( ( ( *r = (unsigned int *) malloc( (size_t) ( *number_colors ) * sizeof ( int ) ) ) == NULL ) ||
             ( ( *g = (unsigned int *) malloc( (size_t) ( *number_colors ) * sizeof ( unsigned int ) ) ) == NULL ) ||
             ( ( *b = (unsigned int *) malloc( (size_t) ( *number_colors ) * sizeof ( unsigned int ) ) ) == NULL ) ||
             ( ( *a = (double *) malloc( (size_t) ( *number_colors ) * sizeof ( double ) ) ) == NULL ) )
        {
            plexit( "cmap0_palette_read: insufficient memory" );
        }
        **r = 255;
        **g = 255;
        **b = 255;
        **a = 1.;
        for ( i = 1; i < *number_colors; i++ )
        {
            *( *r + i ) = 255;
            *( *g + i ) = 0;
            *( *b + i ) = 0;
            *( *a + i ) = 1.0;
        }
    }

    plrestore_locale( save_locale );
}

//--------------------------------------------------------------------------
// void c_plspal0(filename)
//
//! Set the palette for color map 0 using a cmap0*.pal format file.
//! filename: the name of the cmap0*.pal file to use.
//!
//! @param filename name of the cmap0 palette file.

void
c_plspal0( const char *filename )
{
    int          i;
    unsigned int *r, *g, *b;
    double       *a;
    int          number_colors;
    cmap0_palette_read( filename, &number_colors, &r, &g, &b, &a );
    // Allocate default number of cmap0 colours if cmap0 allocation not
    // done already.
    plscmap0n( 0 );
    // Allocate sufficient cmap0 colours to contain present data.
    if ( number_colors > m_plsc->ncol0 )
    {
        plscmap0n( number_colors );
    }
    for ( i = 0; i < number_colors; i++ )
    {
        c_plscol0a( i, (PLINT) r[i], (PLINT) g[i], (PLINT) b[i], a[i] );
    }
    free( r );
    free( g );
    free( b );
    free( a );
}

//! This code fragment used a lot in plspal1 to deal with
//! floating-point range checking of a value and the adjustment of that
//! value when close to the range when there is floating-point errors.
//!
//! @param value The value to range check.
//! @param min The minimum allowable value.
//! @param max The maximum allowable value.
//! @param fuzz Amount of slop to allow beyond the range defined by min & max.
//! @param err_number The error number.
#define fuzzy_range_check( value, min, max, fuzz, err_number )                                                                        \
    if ( value < ( min - fuzz ) || value > ( max + fuzz ) ) {                                                                         \
        snprintf( msgbuf, MSGLEN, "Unrecognized cmap1 format data line.  Error number is %d. Line is %s\n", err_number, color_info ); \
        plwarn( msgbuf );                                                                                                             \
        err = 1;                                                                                                                      \
        break;                                                                                                                        \
    } else if ( value < min ) {                                                                                                       \
        value = min;                                                                                                                  \
    } else if ( value > max ) {                                                                                                       \
        value = max;                                                                                                                  \
    }


tsEmbeddedCol1Palette *			// returned ptr to embedded palette, NULL if unable to locate
GetEmbeddedCol1Palette(const char *pszPaletteName)	// which embedded palette to load - name is same as original file name
{
int Idx;
tsEmbeddedCol1Palette *pSelPalette;
pSelPalette = &EmbeddedCol1Palettes[0];
for(Idx = 0; Idx < cNumEmbeddedCol1Pallets; Idx++,pSelPalette+=1)
	if(!stricmp(pSelPalette->pszPaletteName,pszPaletteName))
		return(pSelPalette);
return(NULL);
}

//--------------------------------------------------------------------------
// void c_plspal1(filename)
//
//! Set the palette for color map 1 using a cmap1*.pal format file.
//! filename: the name of the cmap1*.pal file to use.
//!
//! @param filename name of the cmap1 palette file.
//! @param interpolate interpolate between control points.

void
c_plspal1( const char *filename, PLBOOL interpolate )
{
tsEmbeddedCol1Palette *pSelPalette = NULL;

    int          i;
    int          number_colors;
    int          format_version, err;
    PLBOOL       rgb;
    char         color_info[PALLEN];
    unsigned int r_i, g_i, b_i;
    int          pos_i, alt_hue_path_i;
    double       r_d, g_d, b_d, a_d, pos_d;
    PLFLT        *r, *g, *b, *a, *pos;
    PLINT        *ri, *gi, *bi;
    PLBOOL       *alt_hue_path;
    FILE         *fp;
    char         msgbuf[MSGLEN];
    char         * save_locale = plsave_set_locale();

    rgb            = TRUE;
    err            = 0;
    format_version = 0;
    if ( strlen( filename ) == 0 )
		{
        fp = plLibOpen( PL_DEFAULT_CMAP1_FILE );
        if ( fp == NULL )
			{
			if((pSelPalette = GetEmbeddedCol1Palette(PL_DEFAULT_CMAP1_FILE))==NULL)
				{
				snprintf( msgbuf, MSGLEN, "Unable to open cmap1 .pal file %s\n", PL_DEFAULT_CMAP1_FILE );
				plwarn( msgbuf );
				goto finish;
				}
			}
		}
    else
		{
        fp = plLibOpen( filename );
        if ( fp == NULL )
			{
			if((pSelPalette = GetEmbeddedCol1Palette(filename))==NULL)
				{
				snprintf( msgbuf, MSGLEN, "Unable to open cmap1 .pal file %s\n", filename );
				plwarn( msgbuf );
				goto finish;
				}
			}
		}

	if(pSelPalette != NULL)
		{
		rgb = pSelPalette->Format == 1 ? FALSE : TRUE;
		number_colors = pSelPalette->NumColors;
		}
	else
		{
			// Check for new file format
		if ( read_line( color_info, PALLEN, fp ) == NULL )
			{
			snprintf( msgbuf, MSGLEN, "Error reading cmap1 .pal file %s\n", filename );
			plwarn( msgbuf );
			fclose( fp );
			goto finish;
			}
		if ( strncmp( color_info, "v2 ", 2 ) == 0 )
			{
			format_version = 1;
			if ( strncmp( &color_info[3], "hls", 3 ) == 0 )
				rgb = FALSE;
			else if ( strncmp( &color_info[3], "rgb", 3 ) == 0 )
				rgb = TRUE;
			else
				{
				snprintf( msgbuf, MSGLEN, "Invalid color space %s - assuming RGB\n", &color_info[3] );
				plwarn( msgbuf );
				rgb = TRUE;
				}
			if ( read_line( color_info, PALLEN, fp ) == NULL )
				{
				snprintf( msgbuf, MSGLEN, "Error reading cmap1 .pal file %s\n", filename );
				plwarn( msgbuf );
				fclose( fp );
				goto finish;
				}
			}

		if ( sscanf( color_info, "%d\n", &number_colors ) != 1 || number_colors < 2 )
			{
			snprintf( msgbuf, MSGLEN, "Unrecognized cmap1 format (wrong number of colors) %s\n", color_info );
			plwarn( msgbuf );
			fclose( fp );
			goto finish;
			}
		}

    r            = (PLFLT *) malloc( (size_t) number_colors * sizeof ( PLFLT ) );
    g            = (PLFLT *) malloc( (size_t) number_colors * sizeof ( PLFLT ) );
    b            = (PLFLT *) malloc( (size_t) number_colors * sizeof ( PLFLT ) );
    ri           = (PLINT *) malloc( (size_t) number_colors * sizeof ( PLINT ) );
    gi           = (PLINT *) malloc( (size_t) number_colors * sizeof ( PLINT ) );
    bi           = (PLINT *) malloc( (size_t) number_colors * sizeof ( PLINT ) );
    a            = (PLFLT *) malloc( (size_t) number_colors * sizeof ( PLFLT ) );
    pos          = (PLFLT *) malloc( (size_t) number_colors * sizeof ( PLFLT ) );
    alt_hue_path = (PLBOOL *) malloc( (size_t) ( number_colors - 1 ) * sizeof ( PLBOOL ) );

	if(pSelPalette != NULL)
		{
		tsCol1ColorInfo *pColorInfo;
		pColorInfo = pSelPalette->ColorInfo;
		for(i=0;i<number_colors;i++,pColorInfo++)
			{
			r[i]   = pColorInfo->r;
			g[i]   = pColorInfo->g;
			b[i]   = pColorInfo->b;
			a[i]   = pColorInfo->a;
			pos[i] = pColorInfo->pos;
			if(i < (number_colors-1))
				alt_hue_path[i] = pColorInfo->alt_hue_path;
			}
		}
	else
		{
		if (format_version == 0)
			{
			int return_sscanf = -1, return_sscanf_old = 0;
			// Old tk file format
			for ( i = 0; i < number_colors; i++ )
				{
				if ( read_line( color_info, PALLEN, fp ) == NULL )
					{
					snprintf( msgbuf, MSGLEN, "Error reading cmap1 .pal file %s\n", filename );
					plwarn( msgbuf );
					fclose( fp );
					goto finish;
					}
				// Ensure string is null terminated if > 160 characters
				color_info[PALLEN - 1] = '\0';
				return_sscanf          = sscanf( color_info, "#%2x%2x%2x %d %d", &r_i, &g_i, &b_i, &pos_i, &alt_hue_path_i );
				if ( return_sscanf < 4 || ( return_sscanf_old != 0 && return_sscanf != return_sscanf_old ) )
					{
					snprintf( msgbuf, MSGLEN, "Unrecognized cmap1 format (wrong number of items for version 1 of format) %s\n", color_info );
					plwarn( msgbuf );
					err = 1;
					break;
					}
				return_sscanf_old = return_sscanf;
				// For old format, input colours range from 0 to 255 and
				// need to be renormalized to the range from 0. to 1..
				r[i]   = (PLFLT) r_i / 255.;
				g[i]   = (PLFLT) g_i / 255.;
				b[i]   = (PLFLT) b_i / 255.;
				a[i]   = 1.0;
				pos[i] = 0.01 * (PLFLT) pos_i;
				fuzzy_range_check( r[i], 0., 1., FUZZ_EPSILON, 1 );
				fuzzy_range_check( g[i], 0., 1., FUZZ_EPSILON, 2 );
				fuzzy_range_check( b[i], 0., 1., FUZZ_EPSILON, 3 );
				fuzzy_range_check( pos[i], 0., 1., FUZZ_EPSILON, 4 );
				if ( ( return_sscanf == 5 ) && ( i != number_colors - 1 ) )
					{
					// Next to oldest tk format with alt_hue_path specified.
					alt_hue_path[i] = (PLBOOL) alt_hue_path_i;
					}
				}
			if ( return_sscanf == 4 )
				{
				// Oldest tk format.  No alt_hue_path specified.
				free( alt_hue_path );
				alt_hue_path = NULL;
				}
			}
		else
			{
				// New floating point file version with support for alpha and alt_hue_path values
			for ( i = 0; i < number_colors; i++ )
				{
				if ( read_line( color_info, PALLEN, fp ) == NULL )
					{
					snprintf( msgbuf, MSGLEN, "Error reading cmap1 .pal file %s\n", filename );
					plwarn( msgbuf );
					fclose( fp );
					goto finish;
					}
				if ( sscanf( color_info, "%lf %lf %lf %lf %lf %d", &pos_d, &r_d, &g_d, &b_d, &a_d, &alt_hue_path_i ) != 6 )
					{
					snprintf( msgbuf, MSGLEN, "Unrecognized cmap1 format (wrong number of items for version 2 of format) %s\n", color_info );
					plwarn( msgbuf );
					err = 1;
					break;
					}

				r[i]   = (PLFLT) r_d;
				g[i]   = (PLFLT) g_d;
				b[i]   = (PLFLT) b_d;
				a[i]   = (PLFLT) a_d;
				pos[i] = (PLFLT) pos_d;
					// Check that all rgba and pos data within range from 0. to
					// 1. except for the hls colour space case where the first
					// coordinate is checked within range from 0. to 360.
				if ( rgb )
					{
					fuzzy_range_check( r[i], 0., 1., FUZZ_EPSILON, 5 );
					}
				else
					{
					fuzzy_range_check( r[i], 0., 360., ( 360. * FUZZ_EPSILON ), 6 );
					}
				fuzzy_range_check( g[i], 0., 1., FUZZ_EPSILON, 7 );
				fuzzy_range_check( b[i], 0., 1., FUZZ_EPSILON, 8 );
				fuzzy_range_check( a[i], 0., 1., FUZZ_EPSILON, 9 );
				fuzzy_range_check( pos[i], 0., 1., FUZZ_EPSILON, 10 );

				if ( i != number_colors - 1 )
					alt_hue_path[i] = (PLBOOL) alt_hue_path_i;
				}
			}
		fclose( fp );
		}

    if ( !err )
		{
        if ( interpolate )
			{
            c_plscmap1la( rgb, number_colors, pos, r, g, b, a, alt_hue_path );
			}
        else
			{
            for ( i = 0; i < number_colors; i++ )
				{
                ri[i] = (PLINT) ( r[i] * 255.0 );
                gi[i] = (PLINT) ( g[i] * 255.0 );
                bi[i] = (PLINT) ( b[i] * 255.0 );
				}
            c_plscmap1a( ri, gi, bi, a, number_colors );
			}
		}
    else
		{
        // Fall back to red scale as visual warning if some problem occurred
        // above.
        free( r );
        free( g );
        free( b );
        free( pos );
        number_colors = 2;
        r             = (PLFLT *) malloc( (size_t) number_colors * sizeof ( PLFLT ) );
        g             = (PLFLT *) malloc( (size_t) number_colors * sizeof ( PLFLT ) );
        b             = (PLFLT *) malloc( (size_t) number_colors * sizeof ( PLFLT ) );
        pos           = (PLFLT *) malloc( (size_t) number_colors * sizeof ( PLFLT ) );
        r[0]          = 0.;
        r[1]          = 1.;
        g[0]          = 0.;
        g[1]          = 0.;
        b[0]          = 0.;
        b[1]          = 0.;
        pos[0]        = 0.;
        pos[1]        = 1.;
        c_plscmap1l( TRUE, number_colors, pos, r, g, b, NULL );
		}

    free( r );
    free( g );
    free( b );
    free( ri );
    free( gi );
    free( bi );
    free( a );
    free( pos );
    free( alt_hue_path );

finish: plrestore_locale( save_locale );
}

//--------------------------------------------------------------------------
// A grab-bag of various control routines.
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// void plwarn()
//
//! A handy way to issue warnings, if need be.
//!
//! @param errormsg The error message.

void
plwarn( const char *errormsg )
{
    int was_gfx = 0;

    if ( m_plsc->graphx == 1 )
    {
        was_gfx = 1;
        pltext();
    }

    fprintf( stderr, "\n*** PLPLOT WARNING ***\n" );
    if ( *errormsg != '\0' )
        fprintf( stderr, "%s\n", errormsg );

    if ( was_gfx == 1 )
        plgra();
}

//--------------------------------------------------------------------------
// void plabort()
//
//! Much the same as plwarn(), but appends ", aborting operation" to the
//! error message.  Helps to keep source code uncluttered and provides a
//! convention for error aborts.
//!
//! If cleanup needs to be done in the main program, the user should write
//! his/her own exit handler and pass it in via plsabort().
//!
//! @param errormsg The error message.

void
plabort( const char *errormsg )
{
    if ( abort_handler != NULL )
        ( *abort_handler )( errormsg );

    if ( m_plsc->errcode != NULL )
        *( m_plsc->errcode ) = 1;

    if ( m_plsc->errmsg != NULL )
    {
        sprintf( m_plsc->errmsg, "\n*** PLPLOT ERROR, ABORTING OPERATION ***\n" );
        if ( *errormsg != '\0' )
            sprintf( m_plsc->errmsg, "%s, aborting operation\n", errormsg );
    }
    else
    {
        int was_gfx = 0;

        if ( m_plsc->graphx == 1 )
        {
            was_gfx = 1;
            pltext();
        }

        fprintf( stderr, "\n*** PLPLOT ERROR, ABORTING OPERATION ***\n" );
        if ( *errormsg != '\0' )
            fprintf( stderr, "%s, aborting operation\n", errormsg );

        if ( was_gfx == 1 )
            plgra();
    }
}


//--------------------------------------------------------------------------
// void plsabort()
//
//! Sets an optional user abort handler.
//!
//! @param handler A function that takes a const char * argument that will
//! be called in the event of a abort.
//--------------------------------------------------------------------------

void
plsabort( void ( *handler )( const char * ) )
{
    abort_handler = handler;
}

//--------------------------------------------------------------------------
// void plexit()
//
//! In case of an abort this routine is called.  It just prints out an error
//! message and tries to clean up as much as possible.  It's best to turn
//! off pause and then restore previous setting before returning.
//!
//! If cleanup needs to be done in the main program, the user should write
//! his/her own exit handler and pass it in via plsexit().  This function
//! should should either call plend() before exiting, or simply return.
//!
//! @param errormsg The error message.
//--------------------------------------------------------------------------

void
plexit( const char *errormsg )
{
    int status = 1;

    if ( exit_handler != NULL )
        status = ( *exit_handler )( errormsg );

    m_plsc->nopause = 1;
    if ( *errormsg != '\0' )
    {
        fprintf( stderr, "\n*** PLPLOT ERROR, IMMEDIATE EXIT ***\n" );
        fprintf( stderr, "%s\n", errormsg );
    }
    plend();

    fprintf( stderr, "Program aborted\n" );
    exit( status );
}

//--------------------------------------------------------------------------
// void plsexit()
//
//! Sets an optional user exit handler.
//!
//! @param handler A function that takes a const char * argument that will
//! will be called in the event of a exit.
//--------------------------------------------------------------------------

void
plsexit( int ( *handler )( const char * ) )
{
    exit_handler = handler;
}

//--------------------------------------------------------------------------
// void plgra()
//
//! Switches to graphics screen.
//!
//! Here and in pltext() it's a good idea to return silently if plinit()
//! hasn't yet been called, since plwarn() calls pltext() and plgra(), and
//! plwarn() may be called at any time.
//--------------------------------------------------------------------------

void
c_plgra( void )
{
    if ( m_plsc->level > 0 )
        plP_esc( PLESC_GRAPH, NULL );
}

//--------------------------------------------------------------------------
// void plxormod()
//
//! Set xor mode? FIXME: Not really sure what this function does.
//!
//! @param mode Boolean.
//! @param status 1 if successful, 0 otherwise.

void
c_plxormod( PLINT mode, PLINT *status )   // xor mode
{
    static int ostate = 0;

    if ( !m_plsc->dev_xor )
    {
        *status = 0;
        return;
    }

    if ( m_plsc->level > 0 )
    {
        plP_esc( PLESC_XORMOD, &mode );
        if ( mode )
        {
            ostate            = m_plsc->plbuf_write;
            m_plsc->plbuf_write = 0;
        }
        else
            m_plsc->plbuf_write = ostate;
    }
    *status = 1;
}

//--------------------------------------------------------------------------
//! Set drawing mode (depends on device support!)
//!
//! @param mode This determines which drawing mode to use.
//!
void
c_plsdrawmode( PLINT mode )
{
    if ( !m_plsc->dev_modeset )
    {
        plwarn( "plsdrawmode: Mode setting is not supported by this device" );
    }
    else if ( m_plsc->level > 0 )
    {
        plP_esc( PLESC_MODESET, &mode );
    }
    else
    {
        plwarn( "plsdrawmode: Initialize PLplot first" );
    }
    return;
}

//--------------------------------------------------------------------------
//! Get drawing mode (depends on device support!)
//!
//! @returns Current drawing mode
//!
PLINT
c_plgdrawmode( void )
{
    PLINT mode;

    if ( !m_plsc->dev_modeset )
    {
        plwarn( "plgdrawmode: Mode getting is not supported by this device" );
        mode = PL_DRAWMODE_UNKNOWN;
    }
    else if ( m_plsc->level > 0 )
    {
        plP_esc( PLESC_MODEGET, &mode );
    }
    else
    {
        plwarn( "plsdrawmode: Initialize PLplot first" );
        mode = PL_DRAWMODE_UNKNOWN;
    }

    return ( mode );
}

//--------------------------------------------------------------------------
// void pltext()
//
//! Switches to text screen.
//--------------------------------------------------------------------------

void
c_pltext( void )
{
    if ( m_plsc->level > 0 )
        plP_esc( PLESC_TEXT, NULL );
}

//--------------------------------------------------------------------------
// void pl_cmd()
//
//! Front-end to driver escape function.
//! In principle this can be used to pass just about anything directly
//! to the driver.
//!
//! @param op A PLESC command to pass to the driver.
//! @param ptr Data associated with the op command.
//--------------------------------------------------------------------------

void
pl_cmd( PLINT op, void *ptr )
{
    plP_esc( op, ptr );
}

//--------------------------------------------------------------------------
// char *plFindCommand
//
//! Looks for the specified executable file.  Search path:
//!      if command invoked in the build tree:
//!         build_tree/tk (plserver lies there - needed for the tk driver)
//!         source_tree/scripts (plpr lies there - needed for the tk driver)
//!      else
//!	PLPLOT_BIN_ENV = $(PLPLOT_BIN)
//!	current directory
//!	PLPLOT_HOME_ENV/bin = $(PLPLOT_HOME)/bin
//!	BIN_DIR
//!
//! The caller must free the returned pointer (points to malloc'ed memory)
//! when finished with it.
//!
//! @param fn Name of the executable(?).
//!
//! @returns The location of the executable file.
//--------------------------------------------------------------------------

char *
plFindCommand( const char *fn )
{
    char *fs = NULL, *dn;

    //*** see if in build tree **
    if ( plInBuildTree() == 1 )
    {
        plGetName( BUILD_DIR, "bindings/tk", fn, &fs );
        if ( !plFindName( fs ) )
            return fs;
        else
        {
            plGetName( SOURCE_DIR, "scripts", fn, &fs );
            if ( !plFindName( fs ) )
                return fs;
        }
    }

// PLPLOT_BIN_ENV = $(PLPLOT_BIN)

#if defined ( PLPLOT_BIN_ENV )
    if ( ( dn = getenv( PLPLOT_BIN_ENV ) ) != NULL )
    {
        plGetName( dn, "", fn, &fs );
        if ( !plFindName( fs ) )
            return fs;
        fprintf( stderr, PLPLOT_BIN_ENV "=\"%s\"\n", dn ); // what IS set?
    }
#endif  // PLPLOT_BIN_ENV

// Current directory

    plGetName( ".", "", fn, &fs );
    if ( !plFindName( fs ) )
        return fs;

// PLPLOT_HOME_ENV/bin = $(PLPLOT_HOME)/bin

#if defined ( PLPLOT_HOME_ENV )
    if ( ( dn = getenv( PLPLOT_HOME_ENV ) ) != NULL )
    {
        plGetName( dn, "bin", fn, &fs );
        if ( !plFindName( fs ) )
            return fs;
        fprintf( stderr, PLPLOT_HOME_ENV "=\"%s\"\n", dn ); // what IS set?
    }
#endif  // PLPLOT_HOME_ENV

// BIN_DIR

#if defined ( BIN_DIR )
    plGetName( BIN_DIR, "", fn, &fs );
    if ( !plFindName( fs ) )
        return fs;
#endif

// Crapped out

    free_mem( fs );
    fprintf( stderr, "plFindCommand: cannot locate command: %s\n", fn );
#if defined ( BIN_DIR )
    fprintf( stderr, "bin dir=\"" BIN_DIR "\"\n" );      // what WAS set?
#endif  // BIN_DIR
    return NULL;
}

//--------------------------------------------------------------------------
// FILE *plLibOpen(fn)
//
//! Return file pointer to library file (such as a colormap palette).
//! Locations checked:
//!	PLPLOT_LIB_ENV = $(PLPLOT_LIB)
//!	current directory
//!	PLPLOT_HOME_ENV/lib = $(PLPLOT_HOME)/lib
//!	DATA_DIR
//!	PLLIBDEV
//!
//! @param fn Name of the file.
//!
//! @returns A open file pointer (if successful).
//--------------------------------------------------------------------------

FILE *
plLibOpen( const char *fn )
{
    FILE    *ret = NULL;

    PDFstrm *pdfs = plLibOpenPdfstrm( fn );
    if ( pdfs == NULL )
    {
        return NULL;
    }
    if ( pdfs->file != NULL )
    {
        ret        = pdfs->file;
        pdfs->file = NULL;
    }
    pdf_close( pdfs );
    return ret;
}

//--------------------------------------------------------------------------
// FILE *plLibOpenPdfstrm(fn)
//
//! Return file PDFstrm * to a file (originally used for loading fonts?).
//! Locations checked:
//!	PLPLOT_LIB_ENV = $(PLPLOT_LIB)
//!	current directory
//!	PLPLOT_HOME_ENV/lib = $(PLPLOT_HOME)/lib
//!	DATA_DIR
//!	PLLIBDEV
//!
//! @param fn Name of the file.
//!
//! @returns A open PDFstrm file pointer (if successful)
//--------------------------------------------------------------------------
PDFstrm *
plLibOpenPdfstrm( const char *fn )
{
    PDFstrm *file;
    char    *fs = NULL, *dn = NULL;

//***   search build tree               ***

    if ( plInBuildTree() == 1 )
    {
        plGetName( SOURCE_DIR, "data", fn, &fs );

        if ( ( file = pdf_fopen( fs, "rb" ) ) != NULL )
            goto done;
    }

//***	search PLPLOT_LIB_ENV = $(PLPLOT_LIB)	***

#if defined ( PLPLOT_LIB_ENV )
    if ( ( dn = getenv( PLPLOT_LIB_ENV ) ) != NULL )
    {
        plGetName( dn, "", fn, &fs );

        if ( ( file = pdf_fopen( fs, "rb" ) ) != NULL )
            goto done;
        fprintf( stderr, PLPLOT_LIB_ENV "=\"%s\"\n", dn ); // what IS set?
    }
#endif  // PLPLOT_LIB_ENV

//***	search current directory	***

    if ( ( file = pdf_fopen( fn, "rb" ) ) != NULL )
    {
        pldebug( "plLibOpenPdfstr", "Found file %s in current directory.\n", fn );
        free_mem( fs );
        return ( file );
    }

//***	search PLPLOT_HOME_ENV/lib = $(PLPLOT_HOME)/lib	***

#if defined ( PLPLOT_HOME_ENV )
    if ( ( dn = getenv( PLPLOT_HOME_ENV ) ) != NULL )
    {
        plGetName( dn, "lib", fn, &fs );

        if ( ( file = pdf_fopen( fs, "rb" ) ) != NULL )
            goto done;
        fprintf( stderr, PLPLOT_HOME_ENV "=\"%s\"\n", dn ); // what IS set?
    }
#endif  // PLPLOT_HOME_ENV/lib

//***   search installed location	***

#if defined ( DATA_DIR )
    plGetName( DATA_DIR, "", fn, &fs );

    if ( ( file = pdf_fopen( fs, "rb" ) ) != NULL )
        goto done;
#endif  // DATA_DIR

//***   search hardwired location	***

#ifdef PLLIBDEV
    plGetName( PLLIBDEV, "", fn, &fs );

    if ( ( file = pdf_fopen( fs, "rb" ) ) != NULL )
        goto done;
#endif  // PLLIBDEV

#ifdef macintosh
    file = plMacLibOpen( fn );
    if ( file != NULL )
        goto done;
#endif // macintosh

    if ( plplotLibDir != NULL )
    {
        plGetName( plplotLibDir, "", fn, &fs );
        if ( ( file = pdf_fopen( fs, "rb" ) ) != NULL )
            goto done;
    }

//***   not found, give up      ***
    pldebug( "plLibOpenPdfstr", "File %s not found.\n", fn );
    free_mem( fs );
    return NULL;

done:
    pldebug( "plLibOpenPdfstr", "Found file %s\n", fs );
    free_mem( fs );
    return ( file );
}

//--------------------------------------------------------------------------
// int plFindName
//
//! Authors: Paul Dubois (LLNL), others?
//! This function is in the public domain.
//!
//! Given a pathname, determine if it is a symbolic link.  If so, continue
//! searching to the ultimate terminus - there may be more than one link.
//! Use the error value to determine when the terminus is reached, and to
//! determine if the pathname really exists.  Then stat it to determine
//! whether it's executable.  Return 0 for an executable, errno otherwise.
//! Note that 'p' _must_ have at least one '/' character - it does by
//! construction in this program.  The contents of the array pointed to by
//! 'p' are changed to the actual pathname if findname is successful.
//!
//! This function is only defined under Unix for now.
//!
//! @param p Name of the executable to find.
//!
//! @returns 0 if p is found & is an executable.
//--------------------------------------------------------------------------

#ifdef __unix
int
plFindName( char *p )
{
    ssize_t     n;
    char        buf[PLPLOT_MAX_PATH], *cp;
    struct stat sbuf;

    pldebug( "plFindName", "Trying to find %s\n", p );
    while ( ( n = readlink( p, buf, PLPLOT_MAX_PATH ) ) > 0 )
    {
        pldebug( "plFindName", "Readlink read %d chars at: %s\n", n, p );
        if ( buf[0] == '/' )
        {
            // Link is an absolute path

            strncpy( p, buf, (size_t) n );
            p[n] = '\0';
            pldebug( "plFindName", "Link is absolute: %s\n", p );
        }
        else
        {
            // Link is relative to its directory; make it absolute

            cp = 1 + strrchr( p, '/' );
            strncpy( cp, buf, (size_t) n );
            cp[n] = '\0';
            pldebug( "plFindName",
                "Link is relative: %s\n\tTotal path:%s\n", cp, p );
        }
    }

// This macro not defined on the NEC SX-3

#ifdef SX
#define S_ISREG( mode )    ( mode & S_IFREG )
#endif

// SGI machines return ENXIO instead of EINVAL Dubois 11/92

    if ( errno == EINVAL || errno == ENXIO )
    {
        pldebug( "plFindName", "%s may be the one...\n", p );
        if ( ( stat( p, &sbuf ) == 0 ) && S_ISREG( sbuf.st_mode ) )
        {
            pldebug( "plFindName", "%s is a regular file\n", p );
            return ( access( p, X_OK ) );
        }
    }
    pldebug( "plFindName", "%s found but is not executable\n", p );
    return ( errno ? errno : -1 );
}

#else
int
plFindName( char *p )
{
    return 1;
}
#endif

//--------------------------------------------------------------------------
// void plGetName()
//
//! Gets search name for file by concatenating the dir, subdir, and file
//! name, allocating memory as needed.  The appropriate delimiter is added
//! after the dir specification as necessary.  The caller is responsible
//! for freeing the malloc'ed memory.
//!
//! @param dir The directory name.
//! @param subdir The sub-directory name.
//! @param filename The file name.
//! @param filespec The result of concatenating dir, subdir and filename.
//--------------------------------------------------------------------------

void
plGetName( const char *dir, const char *subdir, const char *filename, char **filespec )
{
    size_t lfilespec;

// Malloc space for filespec

    free_mem( *filespec );
    // Be slightly generous since 3 (two delimiters + NULL byte) should be
    // enough.
    lfilespec = strlen( dir ) + strlen( subdir ) + strlen( filename ) + 10;
    if ( ( *filespec = (char *) malloc( lfilespec ) ) == NULL )
    {
        plexit( "plGetName: Insufficient memory" );
    }

    strcpy( *filespec, dir );

    if ( *subdir != '\0' )
    {
        strcat_delim( *filespec );
        strcat( *filespec, subdir );
    }
    if ( *filename != '\0' )
    {
        strcat_delim( *filespec );
        strcat( *filespec, filename );
    }
    pldebug( "plGetName", "Maximum length of full pathname of file to be found is %zu\n", lfilespec - 1 );
    pldebug( "plGetName", "Full pathname of file to be found is %s\n", *filespec );
}

//--------------------------------------------------------------------------
// void strcat_delim()
//
//! Append path name deliminator if necessary (does not add one if one's
//! there already, or if dealing with a colon-terminated device name).
//!
//! @param dirspec String to add the appropriate delimiter too.
//--------------------------------------------------------------------------

void
strcat_delim( char *dirspec )
{
    size_t ldirspec = strlen( dirspec );
#if defined ( MSDOS ) || defined ( WIN32 )
    if ( dirspec[ldirspec - 1] != '\\' )
        strcat( dirspec, "\\" );
#elif defined ( macintosh )
    if ( dirspec[ldirspec - 1] != ':' )
        strcat( dirspec, ":" );
#else           // unix is the default
    if ( dirspec[ldirspec - 1] != '/' )
        strcat( dirspec, "/" );
#endif
}

//--------------------------------------------------------------------------
// plcol_interp()
//
//! Initializes device cmap 1 entry by interpolation from pls->cmap1
//! entries.  Returned PLColor is supposed to represent the i_th color
//! out of a total of ncol colors in the current color scheme.
//!
//! @param pls A plot stream structure.
//! @param newcolor A color structure to store the color in.
//! @param i Index of the desired color.
//! @param ncol Total number of colors (supported by the device?).
//--------------------------------------------------------------------------

void
plcol_interp( PLStream *pls, PLColor *newcolor, int i, int ncol )
{
    PLFLT x, delta;
    int   il, ir;

    x     = (double) ( i * ( pls->ncol1 - 1 ) ) / (double) ( ncol - 1 );
    il    = (int) x;
    ir    = il + 1;
    delta = x - il;

    if ( ir > pls->ncol1 || il < 0 )
        fprintf( stderr, "Invalid color\n" );

    else if ( ir == pls->ncol1 || ( delta == 0. ) )
    {
        newcolor->r = pls->cmap1[il].r;
        newcolor->g = pls->cmap1[il].g;
        newcolor->b = pls->cmap1[il].b;
        newcolor->a = pls->cmap1[il].a;
    }
    else
    {
        newcolor->r = (unsigned char) ( ( 1. - delta ) * pls->cmap1[il].r + delta * pls->cmap1[ir].r );
        newcolor->g = (unsigned char) ( ( 1. - delta ) * pls->cmap1[il].g + delta * pls->cmap1[ir].g );
        newcolor->b = (unsigned char) ( ( 1. - delta ) * pls->cmap1[il].b + delta * pls->cmap1[ir].b );
        newcolor->a = ( 1. - delta ) * pls->cmap1[il].a + delta * pls->cmap1[ir].a;
    }
}

//--------------------------------------------------------------------------
// plOpenFile()
//
//! Opens file for output, prompting if not set.
//! Prints extra newline at end to make output look better in batch runs.
//! A file name of "-" indicates output to stdout.
//!
//! @param pls A plot stream structure.
//--------------------------------------------------------------------------

#define MAX_NUM_TRIES    10
void
plOpenFile( PLStream *pls )
{
    int    i = 0, count = 0;
    size_t len;
    char   line[BUFFER_SIZE];

    while ( pls->OutFile == NULL )
    {
// Setting pls->FileName = NULL forces creation of a new family member
// You should also free the memory associated with it if you do this

        if ( pls->family && pls->BaseName != NULL )
            plP_getmember( pls );

// Prompt if filename still not known

        if ( pls->FileName == NULL )
        {
            do
            {
                fprintf( stdout, "Enter graphics output file name: " );
                plio_fgets( line, sizeof ( line ), stdin );
                len = strlen( line );
                if ( len )
                    len--;
                line[len] = '\0';       // strip new-line
                count++;                // count zero entries
            } while ( !len && count < MAX_NUM_TRIES );
            plP_sfnam( pls, line );
        }

// If name is "-", send to stdout

        if ( !strcmp( pls->FileName, "-" ) )
        {
            pls->OutFile     = stdout;
            pls->output_type = 1;
            break;
        }

// Need this here again, for prompted family initialization

        if ( pls->family && pls->BaseName != NULL )
            plP_getmember( pls );

        if ( ( pls->OutFile = fopen( pls->FileName, "wb+" ) ) == NULL )
			{
			if ( i++ > 10 )
				{
                fprintf( stderr, "Can't open %s.\n", pls->FileName );
				plexit( "Too many tries." );
				}
			}
        else
            pldebug( "plOpenFile", "Opened %s\n", pls->FileName );
    }
}

//--------------------------------------------------------------------------
// plCloseFile()
//
//! Closes output file unless it is associated with stdout.
//!
//! @param pls A plot stream structure.
//--------------------------------------------------------------------------

void
plCloseFile( PLStream *pls )
{
    if ( pls->OutFile != NULL )
    {
        // Don't close if the output file was stdout
        if ( pls->FileName && strcmp( pls->FileName, "-" ) == 0 )
            return;

        fclose( pls->OutFile );
        pls->OutFile = NULL;
    }
}

//--------------------------------------------------------------------------
// plP_getmember()
//
//! Sets up next file member name (in pls->FileName), but does not open it.
//!
//! @param pls A plot stream structure.
//--------------------------------------------------------------------------

void
plP_getmember( PLStream *pls )
{
    char   tmp[BUFFER_SIZE];
    char   prefix[BUFFER_SIZE];
    char   * suffix;
    char   num[BUFFER_SIZE];
    size_t maxlen;

    maxlen = strlen( pls->BaseName ) + 10;
    if ( pls->FileName == NULL )
    {
        if ( ( pls->FileName = (char *) malloc( maxlen ) ) == NULL )
        {
            plexit( "plP_getmember: Insufficient memory" );
        }
    }

    suffix = strstr( pls->BaseName, "%n" );

    snprintf( tmp, BUFFER_SIZE, "%%0%1ii", (int) pls->fflen );
    snprintf( num, BUFFER_SIZE, tmp, pls->member );

    if ( suffix == NULL )
        snprintf( pls->FileName, maxlen, "%s.%s", pls->BaseName, num );
    else
    {
        strncpy( prefix, pls->BaseName, BUFFER_SIZE - 1 );
        prefix [( suffix - pls->BaseName < BUFFER_SIZE ) ? ( suffix - pls->BaseName ) : BUFFER_SIZE - 1] = '\0';
        snprintf( pls->FileName, maxlen, "%s%s%s", prefix, num, suffix + 2 );
    }
}

//--------------------------------------------------------------------------
// plP_sfnam()
//
//! Sets up file name (with "%n" removed if present) & family stem name.
//! Reserve some extra space (10 chars) to hold an optional member number.
//!
//! @param pls A plot stream.
//! @param fnam The base file name of the plot files.
//--------------------------------------------------------------------------

void
plP_sfnam( PLStream *pls, const char *fnam )
{
    char   prefix[BUFFER_SIZE];
    char   * suffix;
    size_t maxlen;
    pls->OutFile = NULL;

    if ( pls->FileName != NULL )
        free( (void *) pls->FileName );

    maxlen = 10 + strlen( fnam );
    if ( ( pls->FileName = (char *) malloc( maxlen ) ) == NULL )
    {
        plexit( "plP_sfnam: Insufficient memory" );
    }

    suffix = (char *)strstr( fnam, "%n" );

    if ( suffix == NULL )
    {
        strncpy( pls->FileName, fnam, maxlen - 1 );
        pls->FileName[maxlen - 1] = '\0';
    }
    else
    {
        strncpy( prefix, fnam, BUFFER_SIZE - 1 );
        prefix [( suffix - fnam ) < BUFFER_SIZE ? ( suffix - fnam ) : BUFFER_SIZE - 1] = '\0';
        snprintf( pls->FileName, maxlen, "%s%s", prefix, suffix + 2 );
    }

    if ( pls->BaseName != NULL )
        free( (void *) pls->BaseName );

    if ( ( pls->BaseName = (char *) malloc( maxlen ) ) == NULL )
    {
        plexit( "plP_sfnam: Insufficient memory" );
    }

    strncpy( pls->BaseName, fnam, maxlen - 1 );
    pls->BaseName[maxlen - 1] = '\0';
}

//--------------------------------------------------------------------------
// plFamInit()
//
//! Initializes family file parameters.
//!
//! @param pls A plot stream structure.
//--------------------------------------------------------------------------

void
plFamInit( PLStream *pls )
{
    if ( pls->family )
    {
        pls->bytecnt = 0;
        if ( !pls->member )
            pls->member = 1;
        if ( !pls->finc )
            pls->finc = 1;
        if ( !pls->fflen )
            pls->fflen = 1;
        if ( !pls->bytemax )
            pls->bytemax = PL_FILESIZE_KB * 1000;
    }
}

//--------------------------------------------------------------------------
// plGetFam()
//
//! Starts new member file of family file set if necessary.
//!
//! Note each member file is a complete graphics file (can be printed
//! individually), although 'plrender' will treat a family as a single
//! logical file if given the family name instead of the member name.
//!
//! @param pls A plot stream structure.
//--------------------------------------------------------------------------

void
plGetFam( PLStream *pls )
{
    PLFLT xpmm_loc, ypmm_loc;
    if ( pls->family )
    {
        if ( pls->bytecnt > pls->bytemax || pls->famadv )
        {
            PLINT local_page_status = pls->page_status;
            plP_tidy();
            pls->member += pls->finc;
            pls->famadv  = 0;
            plP_init();
            // Restore page status (normally AT_BOP) that was changed
            // to AT_EOP by plP_init.
            pls->page_status = local_page_status;

            // Apply compensating factor to original xpmm and ypmm so that
            // character aspect ratio is preserved when overall aspect ratio
            // is changed.
            plP_gpixmm( &xpmm_loc, &ypmm_loc );
            plP_setpxl( xpmm_loc * m_plsc->caspfactor, ypmm_loc / m_plsc->caspfactor );
            return;
        }
    }
}

//--------------------------------------------------------------------------
// plRotPhy()
//
//! Rotates physical coordinates if necessary for given orientation.
//! Each time orient is incremented, the plot is rotated 90 deg clockwise.
//! Note: this is now used only to rotate by 90 degrees for devices that
//! expect portrait mode.
//!
//! @param orient New plot orientation (0-3)
//! @param xmin Current plot x minimum?
//! @param ymin Current plot y minimum?
//! @param xmax Current plot x maximum?
//! @param ymax Current plot y maximum?
//! @param px Old x coordinate mapped to new x coordinate.
//! @param py Old y coordinate mapped to new y coordinate.
//--------------------------------------------------------------------------

void
plRotPhy( PLINT orient, PLINT xmin, PLINT ymin, PLINT xmax, PLINT ymax,
          PLINT *px, PLINT *py )
{
    int x, y;

    x = *px;
    y = *py;

    switch ( orient % 4 )
    {
    case 1:
        *px = xmin + ( y - ymin );
        *py = ymin + ( xmax - x );
        break;

    case 2:
        *px = xmin + ( xmax - x );
        *py = ymin + ( ymax - y );
        break;

    case 3:
        *px = xmin + ( ymax - y );
        *py = ymin + ( x - xmin );
        break;

    default:
        break;                  // do nothing
    }
}

//--------------------------------------------------------------------------
// plAllocDev()
//
//! Allocates a standard PLDev structure for device-specific data, stores
//! the address in pls->dev, and returns the address as well.
//!
//! @param pls A plot stream structure.
//!
//! @returns A PLDev *
//--------------------------------------------------------------------------

PLDev *
plAllocDev( PLStream *pls )
{
    if ( pls->dev != NULL )
        free( (void *) pls->dev );

    pls->dev = calloc( 1, (size_t) sizeof ( PLDev ) );
    if ( pls->dev == NULL )
        plexit( "plAllocDev: cannot allocate memory\n" );

    return (PLDev *) pls->dev;
}

//--------------------------------------------------------------------------
// plGinInit()
//
//! Just fills in the PLGraphicsIn with appropriate initial values.
//!
//! @param gin A plot graphics input (i.e. keypress or mouseclick) structure.
//--------------------------------------------------------------------------

void
plGinInit( PLGraphicsIn *gin )
{
    gin->type      = 0;
    gin->state     = 0;
    gin->keysym    = 0;
    gin->button    = 0;
    gin->string[0] = '\0';
    gin->pX        = gin->pY = -1;
    gin->dX        = gin->dY = 0.;
    gin->wX        = gin->wY = 0.;
}

//--------------------------------------------------------------------------
// plGetInt()
//
//! Prompts human to input an integer in response to given message.
//!
//! @param s The prompt message.
//!
//! @returns The PLINT the human entered.
//--------------------------------------------------------------------------

PLINT
plGetInt( const char *s )
{
    int  m;
    int  i = 0;
    char line[BUFFER_SIZE];

    while ( i++ < 10 )
    {
        fputs( s, stdout );
        plio_fgets( line, sizeof ( line ), stdin );

#ifdef MSDOS
        m = atoi( line );
        return ( m );
#else
        if ( sscanf( line, "%d", &m ) == 1 )
            return ( m );
        fprintf( stdout, "No value or value out of range; please try again\n" );
#endif
    }
    plexit( "Too many tries." );
    return ( 0 );
}

//--------------------------------------------------------------------------
// plGetFlt()
//
//! Prompts human to input a float in response to given message.
//!
//! @param s The prompt message.
//!
//! @returns The PLFLT the human entered.
//--------------------------------------------------------------------------

PLFLT
plGetFlt( const char *s )
{
    PLFLT  m;
    double m1;
    int    i = 0;
    char   line[BUFFER_SIZE];

    while ( i++ < 10 )
    {
        fputs( s, stdout );
        plio_fgets( line, sizeof ( line ), stdin );

#ifdef MSDOS
        m = atof( line );
        return ( m );
#else
        if ( sscanf( line, "%lf", &m1 ) == 1 )
        {
            m = (PLFLT) m1;
            return ( m );
        }
        fprintf( stdout, "No value or value out of range; please try again\n" );
#endif
    }
    plexit( "Too many tries." );
    return ( 0. );
}

//--------------------------------------------------------------------------
// plstrdup()
//
//! A replacement for strdup(), which isn't portable.
//! Caller responsible for freeing the allocated memory.
//!
//! @param src The string to duplicate.
//!
//! @returns A copy of the string src.
//--------------------------------------------------------------------------

char *
plstrdup( const char *src )
{
    char *dest = (char *) malloc( ( strlen( src ) + 1 ) * sizeof ( char ) );
    if ( dest != NULL )
        strcpy( dest, src );
    else
        plabort( "Out of memory" );

    return dest;
}

#ifndef PL_HAVE_SNPRINTF
//--------------------------------------------------------------------------
// plsnprintf()
//
//! Dummy function for snprintf(). This function just calls
//! the unsafe function ignoring the string size. This function will
//! rarely be needed if ever.
//!
//! @param buffer String output buffer.
//! @param n Size of buffer.
//! @param format The format string.
//! @param ... The values that go in the format string (...)
//!
//! @returns The length of buffer that is actually used.
//--------------------------------------------------------------------------

int
plsnprintf( char *buffer, int n, const char *format, ... )
{
    int     ret;

    va_list args;
    va_start( args, format );
    ret = vsprintf( buffer, format, args );
    va_end( args );

    // Check if overrun occured
    if ( ret > n - 1 )
        plabort( "plsnprintf: buffer overrun" );

    return ret;
}

//--------------------------------------------------------------------------
// plsnscanf()
//
//! Dummy function for snscanf(). This function just calls
//! the unsafe function ignoring the string size. This function will
//! rarely be needed if ever.
//!
//! @param buffer String output buffer.
//! @param n Size of buffer.
//! @param format The format string.
//! @param ... The values that go in the format string (...)
//!
//! @returns The length of buffer that is actually used.
//--------------------------------------------------------------------------

int
plsnscanf( const char *buffer, int n, const char *format, ... )
{
    int     ret;

    va_list args;
    va_start( args, format );
    ret = vsscanf( buffer, format, args );
    va_end( args );

    return ret;
}

#endif // PL_HAVE_SNPRINTF

//--------------------------------------------------------------------------
// plseed()
//
//! Set the seed for the random number generator included.
//!
//! @param seed The random number generator seed value.
//--------------------------------------------------------------------------

void
c_plseed( unsigned int seed )
{
    init_genrand( seed );
}

//--------------------------------------------------------------------------
// plrandd()
//
//! @returns A random number on [0,1]-interval.
//!
//--------------------------------------------------------------------------

PLFLT
c_plrandd( void )
{
    return (PLFLT) ( genrand_real1() );
}

//--------------------------------------------------------------------------
// plsave_set_locale()
//
//! Save LC_NUMERIC locale in a string.  The pointer to that string is
//! returned. Then set LC_NUMERIC to "C" locale.
//! n.b. plsave_set_locale and plrestore_locale should always be used as
//! a pair to surround PLplot code that absolutely requires the
//! LC_NUMERIC "C" locale to be in effect.  It is one of plrestore_locale's
//! responsibilities to free the memory allocated here for the locale
//! string.
//!
//! @returns The LC_NUMERIC locale.
//--------------------------------------------------------------------------

char *
plsave_set_locale( void )
{
    char * setlocale_ptr;
    char * saved_lc_numeric_locale;

    if ( !( saved_lc_numeric_locale = (char *) malloc( 100 * sizeof ( char ) ) ) )
    {
        plexit( "plsave_set_locale: out of memory" );
    }

    //save original LC_NUMERIC locale for restore below.
    if ( !( setlocale_ptr = setlocale( LC_NUMERIC, NULL ) ) )
    {
        plexit( "plsave_set_locale: LC_NUMERIC locale could not be determined for NULL locale.\n" );
    }
    strncpy( saved_lc_numeric_locale, setlocale_ptr, 100 );
    saved_lc_numeric_locale[99] = '\0';

    // Do not use pldebug since get overflowed stack (infinite recursion)
    // if device is interactive (i.e., pls->termin is set).
    // comment out fprintf (unless there is some emergency debugging to do)
    // because output is too voluminous.
    //
    // fprintf(stderr, "plsave_set_locale: saved LC_NUMERIC locale is \"%s\"\n", saved_lc_numeric_locale);
    //

    if ( !( setlocale( LC_NUMERIC, "C" ) ) )
    {
        plexit( "plsave_set_locale: LC_NUMERIC locale could not be set to \"C\"" );
    }
    return saved_lc_numeric_locale;
}

//--------------------------------------------------------------------------
// plrestore_locale()
//
//! Restore LC_NUMERIC locale string that was determined by
//! plsave_set_locale with the pointer to that string as the argument.
//! Also, free the memory for that string.
//!
//! @param saved_lc_numeric_locale The saved numeric locale..
//--------------------------------------------------------------------------

void
plrestore_locale( char *saved_lc_numeric_locale )
{
    // Do not use pldebug since get overflowed stack (infinite recursion)
    // if device is interactive (i.e., pls->termin is set).
    // comment out fprintf (unless there is some emergency debugging to do)
    // because output is too voluminous.
    //
    // fprintf(stderr, "plrestore_locale: restored LC_NUMERIC locale is \"%s\"\n", saved_lc_numeric_locale);
    //

    if ( !( setlocale( LC_NUMERIC, saved_lc_numeric_locale ) ) )
    {
        char msgbuf[1024];
        snprintf( msgbuf, 1024, "plrestore_locale: LC_NUMERIC could not be restored to the default \"%s\" locale.\n", saved_lc_numeric_locale );
        plexit( msgbuf );
    }
    free( saved_lc_numeric_locale );
}

