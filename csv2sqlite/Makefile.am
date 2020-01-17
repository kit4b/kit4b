bin_PROGRAMS = csv2sqlite
csv2sqlite_SOURCES= csv2sqlite.cpp SQLiteMarkers.cpp SQLiteMarkers.h SQLiteDE.cpp SQLiteDE.h

# set the include path found by configure
AM_CPPFLAGS= $(all_includes)

# our library
LDADD = ../libkit4b/libkit4b.a ../libzlib/libzlib.a
