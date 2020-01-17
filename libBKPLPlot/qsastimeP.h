#ifndef _QSASTIMEP_H_
#define _QSASTIMEP_H_

// Declaration of private structures within libqsastime which the user does
// not need to acccess.
// Also definition of functions within libqsastime that are needed
// for testing of libqsastime, but which are not normally needed for anything
// else

#include "qsastime.h"

typedef struct MJDtimeStruct
{
    //
    // MJD starts at 0h, so truncating MJD always gives the same day whatever the time (unlike JD).
    // The MJD base day is arbitrary, i.e. seconds can be greater than one day or even negative.
    //

    int    base_day; // integer part of MJD used as default
    double time_sec; // seconds from start of base_day
}MJDtime;

struct QSASConfigStruct
{
    // Values used to define the transformation between broken down time
    // and continuous time for the public API of libqsastime,
    // continuous_time_qsas, broken_down_time_qsas, and strfqsas.

    // scale multiplies the continuous time variable to convert the units to
    // days.
    double scale;

    // offset1 and offset2 (in days) specifies the amount to add to the
    // scaled continuous time to derive the MJD time value that is used
    // internally by libqsastime.  Normally epoch1 is an integral
    // value (which can be exactly stored in a double for a very wide
    // range of integers) and offset2 is normally a non-integral value
    // whose absolute value is less than 1.  This arrangement allows the
    // continuous time variable in the API to be stored as a single double
    // without compromising numerical precision if epoch1 and epoch2
    // are chosen wisely.
    double offset1, offset2;

    // The various bits of ccontrol are used as independent switches to
    // control optional additional corrections which define the
    // transformation between continuous time and broken-down time.
    //
    // If bit 0 (the lowest order bit of ccontrol) is 1 the Julian
    // proleptic calendar is used for broken-down time. Otherwise the
    // Gregorian proleptic calendar is used for broken-down time.
    //
    // If bit 1 is 1, an additional correction for the difference
    // between atomic-clock based times and UTC is applied to the broken-down
    // times.
    //
    // We reserve other bits of ccontrol for future use.
    int ccontrol;
    // index keeps track of latest bhunt_search index.
    int index;
};

// Set if the qsastime library is being tested.
/* #undef TEST_QSASTIME */

#ifdef TEST_QSASTIME
#define QSASTIME_static
#else
//#define QSASTIME_static    static
#define QSASTIME_static
#endif

QSASTIME_static void bhunt_search( const void *key, const void *base, int n, size_t size, int *low, int ( *ge )( const void *keyval, const void *datum ) );

QSASTIME_static int setFromUT( int year, int month, int day, int hour, int min, double sec, MJDtime *MJD, int forceJulian );
QSASTIME_static void breakDownMJD( int *year, int *month, int *day, int *hour, int *min, double *sec, const MJDtime *MJD, int forceJulian );
QSASTIME_static size_t strfMJD( char * buf, size_t len, const char *format, const MJDtime *MJD, int forceJulian, int if60secformat );
QSASTIME_static void normalize_MJD( MJDtime *MJD );
QSASTIME_static const char * getDayOfWeek( const MJDtime *MJD );
QSASTIME_static const char * getLongDayOfWeek( const MJDtime *MJD );
QSASTIME_static const char * getMonth( int m );
QSASTIME_static const char * getLongMonth( int m );
QSASTIME_static void getYAD( int *year, int *ifleapyear, int *doy, const MJDtime *MJD, int forceJulian );

#endif
