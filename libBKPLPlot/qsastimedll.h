#ifndef __QSASTIME_DLL_H
#define __QSASTIME_DLL_H

#ifdef USINGDLL
#if defined ( WIN32 )
// Visual C/C++, Borland, MinGW and Watcom
#if defined ( __VISUALC__ ) || defined ( _MSC_VER ) || defined ( __BORLANDC__ ) || defined ( __GNUC__ ) || defined ( __WATCOMC__ )
#define QSASTIMEDLLEXPORT    __declspec( dllexport )
#define QSASTIMEDLLIMPORT    __declspec( dllimport )
#else
#define QSASTIMEDLLEXPORT
#define QSASTIMEDLLIMPORT
#endif
#elif defined ( __CYGWIN__ )
#define QSASTIMEDLLEXPORT    __declspec( dllexport )
#define QSASTIMEDLLIMPORT    __declspec( dllimport )
#elif defined ( __GNUC__ ) && __GNUC__ > 3
// Follow ideas in http://gcc.gnu.org/wiki/Visibility for GCC version 4.x
// The following forces exported symbols specifically designated with
// QSASTIMEDLLEXPORT to be visible.
#define QSASTIMEDLLEXPORT    __attribute__ ( ( visibility( "default" ) ) )
#define QSASTIMEDLLIMPORT
#endif
#endif

// For an unknown compiler or static built we clear the macros
#ifndef QSASTIMEDLLEXPORT
#    define QSASTIMEDLLEXPORT
#    define QSASTIMEDLLIMPORT
#endif

// The IMPEXP macros will always be set to DLLIMPORT (even for
// the static library, but DLLIMPORT is empty in this case). If
// cmake sets the corresponding macro xxxx_EXPORTS if the
// corresponding library is built DLLIMPEXP is set to DLLEXPORT
#if defined ( qsastime_EXPORTS )
#define QSASTIMEDLLIMPEXP    QSASTIMEDLLEXPORT
#define QSASTIMEDLLIMPEXP_DATA( type )    QSASTIMEDLLEXPORT type
#else
#define QSASTIMEDLLIMPEXP    QSASTIMEDLLIMPORT
#define QSASTIMEDLLIMPEXP_DATA( type )    QSASTIMEDLLIMPORT type
#endif

#endif // __QSASTIME_DLL_H
