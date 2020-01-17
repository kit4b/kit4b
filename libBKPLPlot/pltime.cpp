// $Id: pltime.c 11680 2011-03-27 17:57:51Z airwin $
//
//      Routines for interfacing with qsastime library routines.
//
// Copyright (C) 2009  Alan W. Irwin
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
#include "stdafx.h"
#include "plplotP.h"

// Calculate broken-down time from continuous time for current stream.
void
c_plbtime( PLINT *year, PLINT *month, PLINT *day, PLINT *hour, PLINT *min, PLFLT *sec, PLFLT ctime )
{
    btimeqsas( year, month, day, hour, min, sec, ctime, m_plsc->qsasconfig );
}

// Configure transformation between continuous and broken-down time (and
// vice versa) for current stream.
void
c_plconfigtime( PLFLT scale, PLFLT offset1, PLFLT offset2, PLINT ccontrol, PLBOOL ifbtime_offset, PLINT year, PLINT month, PLINT day, PLINT hour, PLINT min, PLFLT sec )
{
    if ( scale == 0. )
    {
        // Default transformation between continuous and broken-down time
        // (and vice versa) defined here for PLplot.
        // Note the PLplot default is not necessarily the same as the
        // libqsastime default.
        configqsas( 1. / 86400., 0., 0., 0x0, 1, 1970, 0, 1, 0, 0, 0., &( m_plsc->qsasconfig ) );
    }
    else
    {
        configqsas( scale, offset1, offset2, ccontrol, ifbtime_offset, year, month, day, hour, min, sec, &( m_plsc->qsasconfig ) );
    }
}

// Calculate continuous time from broken-down time for current stream.
void
c_plctime( PLINT year, PLINT month, PLINT day, PLINT hour, PLINT min, PLFLT sec, PLFLT *ctime )
{
    int ret;
    ret = ctimeqsas( year, month, day, hour, min, sec, ctime, m_plsc->qsasconfig );
    if ( ret )
        plabort( "plctime: ctimeqsas detected error" );
}

// Set format for date / time labels.
void
c_pltimefmt( const char *fmt )
{
    if ( m_plsc->timefmt )
        free_mem( m_plsc->timefmt );

    m_plsc->timefmt = (char *) malloc( (size_t) ( strlen( fmt ) + 1 ) );
    strcpy( m_plsc->timefmt, fmt );
}

