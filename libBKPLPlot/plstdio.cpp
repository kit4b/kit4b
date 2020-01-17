// $Id: plstdio.c 12831 2013-12-09 14:35:42Z andrewross $
//
// Standardized I/O handler for PLplot.
//
// Copyright (C) 2006  Jim Dishaw
// Copyright (C) 2006  Hazen Babcock
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
//
#include "stdafx.h"
#include "plplotP.h"

#include <sys/types.h>
#include <fcntl.h>

// This is needed for mode flags for mkfifo as well sa for MSDOS / WIN32
#include <sys/stat.h>

// For Visual C++ 2005 and later mktemp() and open() are deprecated (see
// http://msdn.microsoft.com/en-us/library/ms235413.aspx and
// http://msdn.microsoft.com/en-us/library/ms235491.aspx). mktemp()
// is redefined to _mktemp() as well as open() to _open(). In addition
// we need to include io.h.
//
#if defined ( _MSC_VER ) && _MSC_VER >= 1400
#include <io.h>
#define mktemp    _mktemp
#define open      _open
#define fdopen    _fdopen
#endif
//
// plio_write()
//
// Writes the contents of buf to stream.  Handles any I/O error conditions
// so that the caller can "fire and forget."
//

void
plio_fwrite( void *buf, size_t size, size_t nmemb, FILE *stream )
{
    dbug_enter( "plio_fwrite" );

    // Exit if there is nothing to write
    if ( size == 0 || nmemb == 0 )
        return;

    // Clear the error flag for this steam
    clearerr( stream );

    fwrite( buf, size, nmemb, stream );

    if ( ferror( stream ) )
    {
        // Perhaps we can add a flag (global or per output stream)
        // in order to decide if we should abort or warn.  I think
        // I/O errors should generate an abort
        plabort( "Error writing to file" );
    }
}

//
// plio_fread()
//
// Read from stream into buf.  Like plio_write(), this function will
// handle any I/O error conditions.
//

void
plio_fread( void *buf, size_t size, size_t nmemb, FILE *stream )
{
    size_t bytes;

    dbug_enter( "plio_fread" );

    // If the buffer has a size of zero, we should complain
    if ( size == 0 || nmemb == 0 )
    {
        plwarn( "Zero length buffer size in plio_fread, returning" );
        return;
    }

    // Clear the error flag for this steam
    clearerr( stream );

    bytes = fread( buf, size, nmemb, stream );

    if ( ( bytes < nmemb ) && ferror( stream ) )
    {
        // The read resulted in an error
        plabort( "Error reading from file" );
    }
}

//
// plio_fgets()
//
// Read from stream into buf.  This version of fgets is designed for the occasions
// where the caller wants to ignore the return value.
//
// NOTE: If one is reading from a file until an EOF condition, fgets() is better suited
// than this function, i.e.
//
//     while(fgets(buf, size, fp) != NULL) { ... do some stuff ... }
//
// rather than
//
//     while(!feof(fp)) { plio_fgets(buf, size, fp);  ... do some stuff ... }
//
// which would require checking for an empty buffer.
//

void
plio_fgets( char *buf, int size, FILE *stream )
{
    char *s;

    dbug_enter( "plio_fgets" );

    // If the buffer has a size of zero, we should complain
    if ( size == 0 )
    {
        plwarn( "Zero length buffer size in plio_fgets, returning" );
        return;
    }

    // Clear the error flag for this steam
    clearerr( stream );

    s = fgets( buf, size, stream );

    if ( s == NULL && ferror( stream ) )
    {
        // The read resulted in an error
        plabort( "Error reading from file" );
    }
}

//
// pl_create_tempfile()
//
// Securely create a temporary file and return a file handle to it.
// This provides cross-platform compatibility and also adds some
// additional functionality over mkstemp in that it honours the TMP /
// TMPDIR / TEMP environment variables.
//
// The function returns the file handle.
//
// If the fname variable is not NULL, then on return it will contain
// a pointer to the full temporary file name. This will be allocated
// with malloc. It is the caller's responsibility to ensure this
// memory is free'd and to ensure the file is deleted after use.
// If fname is NULL then the file will be automatically deleted
// when it is closed.
//
FILE *
pl_create_tempfile( char **fname )
{
    FILE       *fd;
    const char *tmpdir;
    char       *templte;
    const char *tmpname = "plplot_XXXXXX";
#if !defined PL_HAVE_MKSTEMP
    int        flags;
#endif

#if defined ( MSDOS ) || defined ( WIN32 )
    tmpdir = getenv( "TEMP" );
#else
    tmpdir = getenv( "TMPDIR" );
#endif

// The P_TMPDIR macro is defined in stdio.h on many UNIX systems - try that
#ifdef P_TMPDIR
    if ( tmpdir == NULL )
        tmpdir = P_TMPDIR;
#endif

    if ( tmpdir == NULL )
    {
#if defined ( MSDOS ) || defined ( WIN32 )
        tmpdir = "c:\\windows\\Temp";
#else
        tmpdir = "/tmp";
#endif
    }

    // N.B. Malloc ensures templte is long enough so strcpy and strcat are safe here
    templte = (char *) malloc( sizeof ( char ) * ( strlen( tmpdir ) + strlen( tmpname ) + 2 ) );
    strcpy( templte, tmpdir );
#if defined ( MSDOS ) || defined ( WIN32 )
    strcat( templte, "\\" );
#else
    strcat( templte, "/" );
#endif
    strcat( templte, tmpname );

#ifdef PL_HAVE_MKSTEMP
    fd = fdopen( mkstemp( templte ), "wb+" );
    if ( fd == NULL )
    {
        plwarn( "pl_create_tempfile: Unable to open temporary file - returning" );
        if ( fname != NULL )
            *fname = NULL;
        free( templte );
        return NULL;
    }
    // If we are not returning the file name then unlink the file so it is
    // automatically deleted.
#ifdef PL_HAVE_UNLINK
    if ( fname == NULL )
        unlink( templte );
#endif
#else
#if !defined ( _S_IREAD )
#define _S_IREAD     256
#endif
#if !defined ( _S_IWRITE )
#define _S_IWRITE    128
#endif
    fd    = NULL;
    flags = O_RDWR | O_BINARY | O_CREAT | O_EXCL | _O_SHORT_LIVED;
    // If we are not returning the file name then add flag to automatically
    // delete file once all file handles are closed.
    if ( fname == NULL )
        flags = flags | _O_TEMPORARY;
    mktemp( templte );
    fd = fdopen( open( templte, flags, _S_IREAD | _S_IWRITE ), "wb+" );
#endif

    if ( fname != NULL )
    {
        *fname = templte;
    }
    else
    {
        free( templte );
    }

    return fd;
}

//
// pl_create_tempfifo()
//
// Securely create a temporary fifo and return the file name.
// This only works on POSIX compliant platforms at the moment.
// It creates a secure directory first using mkdtemp, then
// creates the named fifo in this directory. The combination of
// a private directory and mkfifo failing if the file name already exists
// makes this secure against race conditions / DoS attacks. This function
// includes additional functionality over mkdtemp in that it honours the
// TMP / TMPDIR / TEMP environment variables.
//
// The function returns the file name of the fifo.
//
char *
pl_create_tempfifo( const char **p_fifoname, const char **p_dirname )
{
#if !defined PL_HAVE_MKDTEMP || !defined PL_HAVE_MKFIFO
    plwarn( "Creating fifos not supported on this platform" );
    return NULL;
#else
    const char *tmpdir;
    char       *templte;
    char       *dirname;
    const char *tmpname  = "plplot_dir_XXXXXX";
    const char *fifoname = "plplot_fifo";

#if defined ( MSDOS ) || defined ( WIN32 )
    tmpdir = getenv( "TEMP" );
#else
    tmpdir = getenv( "TMPDIR" );
#endif

// The P_TMPDIR macro is defined in stdio.h on many UNIX systems - try that
#ifdef P_TMPDIR
    if ( tmpdir == NULL )
        tmpdir = P_TMPDIR;
#endif

    if ( tmpdir == NULL )
    {
#if defined ( MSDOS ) || defined ( WIN32 )
        tmpdir = "c:\\windows\\Temp";
#else
        tmpdir = "/tmp";
#endif
    }

    // N.B. Malloc ensures templte is long enough so strcpy and strcat are safe here
    dirname = (char *) malloc( sizeof ( char ) * ( strlen( tmpdir ) + strlen( tmpname ) + 2 ) );
    strcpy( dirname, tmpdir );
#if defined ( MSDOS ) || defined ( WIN32 )
    strcat( dirname, "\\" );
#else
    strcat( dirname, "/" );
#endif
    strcat( dirname, tmpname );
    // Create the temporary directory
    dirname    = mkdtemp( dirname );
    *p_dirname = dirname;

    // Now create the fifo in the directory
    templte = (char *) malloc( sizeof ( char ) * ( strlen( tmpdir ) + strlen( tmpname ) + strlen( fifoname ) + 4 ) );
    strcpy( templte, dirname );
#if defined ( MSDOS ) || defined ( WIN32 )
    strcat( templte, "\\" );
#else
    strcat( templte, "/" );
#endif
    strcat( templte, fifoname );
    *p_fifoname = templte;

    // Check that mkfifo succeeds safely
    if ( mkfifo( templte, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH ) < 0 )
    {
        plwarn( "mkfifo error" );
        free( templte );
        *p_fifoname = NULL;
        free( dirname );
        *p_dirname = NULL;
        return NULL;
    }

    return templte;
#endif
}

