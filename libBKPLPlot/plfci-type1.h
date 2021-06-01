//
// This header file contains the lookup table used for converting between
// FCIs (font characterization integers) and font names for the standard
// 35 type 1 fonts.
//
// Copyright (C) 2005-2010  Alan W. Irwin
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

// This file only relevant to device drivers (currently just pdf and
// ps) that use Type1 fonts.

// There are no good choices for script fonts for Type1 so default to
// the Helvetica (sans) variants in that case.

// Default to Helvetica (sans) variants for symbol fonts to follow
// what is done for all modern unicode-aware TrueType font devices.

// N.B. if the glyph lookup comes up blank for any of the fonts below,
// then an additional search of the Type1 Symbol font glyphs is
// implemented in the Type1 device drivers as a fallback.

// N.B. When updating this table by hand be sure to keep it in
// ascending order in fci!

#define N_Type1Lookup    30
static const FCI_to_FontName_Table Type1Lookup[N_Type1Lookup] = {
    { PL_FCI_MARK | 0x000, (const uint8_t *) "Helvetica"             },
    { PL_FCI_MARK | 0x001, (const uint8_t *) "Times-Roman"           },
    { PL_FCI_MARK | 0x002, (const uint8_t *) "Courier"               },
    { PL_FCI_MARK | 0x003, (const uint8_t *) "Helvetica"             },
    { PL_FCI_MARK | 0x004, (const uint8_t *) "Helvetica"             },
    { PL_FCI_MARK | 0x010, (const uint8_t *) "Helvetica-Oblique"     },
    { PL_FCI_MARK | 0x011, (const uint8_t *) "Times-Italic"          },
    { PL_FCI_MARK | 0x012, (const uint8_t *) "Courier-Oblique"       },
    { PL_FCI_MARK | 0x013, (const uint8_t *) "Helvetica-Oblique"     },
    { PL_FCI_MARK | 0x014, (const uint8_t *) "Helvetica-Oblique"     },
    { PL_FCI_MARK | 0x020, (const uint8_t *) "Helvetica-Oblique"     },
    { PL_FCI_MARK | 0x021, (const uint8_t *) "Times-Italic"          },
    { PL_FCI_MARK | 0x022, (const uint8_t *) "Courier-Oblique"       },
    { PL_FCI_MARK | 0x023, (const uint8_t *) "Helvetica-Oblique"     },
    { PL_FCI_MARK | 0x024, (const uint8_t *) "Helvetica-Oblique"     },
    { PL_FCI_MARK | 0x100, (const uint8_t *) "Helvetica-Bold"        },
    { PL_FCI_MARK | 0x101, (const uint8_t *) "Times-Bold"            },
    { PL_FCI_MARK | 0x102, (const uint8_t *) "Courier-Bold"          },
    { PL_FCI_MARK | 0x103, (const uint8_t *) "Helvetica-Bold"        },
    { PL_FCI_MARK | 0x104, (const uint8_t *) "Helvetica-Bold"        },
    { PL_FCI_MARK | 0x110, (const uint8_t *) "Helvetica-BoldOblique" },
    { PL_FCI_MARK | 0x111, (const uint8_t *) "Times-BoldItalic"      },
    { PL_FCI_MARK | 0x112, (const uint8_t *) "Courier-BoldOblique"   },
    { PL_FCI_MARK | 0x113, (const uint8_t *) "Helvetica-BoldOblique" },
    { PL_FCI_MARK | 0x114, (const uint8_t *) "Helvetica-BoldOblique" },
    { PL_FCI_MARK | 0x120, (const uint8_t *) "Helvetica-BoldOblique" },
    { PL_FCI_MARK | 0x121, (const uint8_t *) "Times-BoldItalic"      },
    { PL_FCI_MARK | 0x122, (const uint8_t *) "Courier-BoldOblique"   },
    { PL_FCI_MARK | 0x123, (const uint8_t *) "Helvetica-BoldOblique" },
    { PL_FCI_MARK | 0x124, (const uint8_t *) "Helvetica-BoldOblique" },
};
