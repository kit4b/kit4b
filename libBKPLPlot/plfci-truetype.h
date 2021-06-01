//
// This header file contains the lookup table used for converting between
// FCIs (font characterization integers) and font names for TrueType fonts.
//
// Copyright (C) 2005  Alan W. Irwin
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

// N.B. When updating this table by hand be sure to keep it in
// ascending order in fci!
//
#define N_TrueTypeLookup    30
static const FCI_to_FontName_Table TrueTypeLookup[N_TrueTypeLookup] = {
    { PL_FCI_MARK | 0x000, (const uint8_t *) PL_FREETYPE_SANS                },
    { PL_FCI_MARK | 0x001, (const uint8_t *) PL_FREETYPE_SERIF               },
    { PL_FCI_MARK | 0x002, (const uint8_t *) PL_FREETYPE_MONO                },
    { PL_FCI_MARK | 0x003, (const uint8_t *) PL_FREETYPE_SCRIPT              },
    { PL_FCI_MARK | 0x004, (const uint8_t *) PL_FREETYPE_SYMBOL              },
    { PL_FCI_MARK | 0x010, (const uint8_t *) PL_FREETYPE_SANS_ITALIC         },
    { PL_FCI_MARK | 0x011, (const uint8_t *) PL_FREETYPE_SERIF_ITALIC        },
    { PL_FCI_MARK | 0x012, (const uint8_t *) PL_FREETYPE_MONO_ITALIC         },
    { PL_FCI_MARK | 0x013, (const uint8_t *) PL_FREETYPE_SCRIPT_ITALIC       },
    { PL_FCI_MARK | 0x014, (const uint8_t *) PL_FREETYPE_SYMBOL_ITALIC       },
    { PL_FCI_MARK | 0x020, (const uint8_t *) PL_FREETYPE_SANS_OBLIQUE        },
    { PL_FCI_MARK | 0x021, (const uint8_t *) PL_FREETYPE_SERIF_OBLIQUE       },
    { PL_FCI_MARK | 0x022, (const uint8_t *) PL_FREETYPE_MONO_OBLIQUE        },
    { PL_FCI_MARK | 0x023, (const uint8_t *) PL_FREETYPE_SCRIPT_OBLIQUE      },
    { PL_FCI_MARK | 0x024, (const uint8_t *) PL_FREETYPE_SYMBOL_OBLIQUE      },
    { PL_FCI_MARK | 0x100, (const uint8_t *) PL_FREETYPE_SANS_BOLD           },
    { PL_FCI_MARK | 0x101, (const uint8_t *) PL_FREETYPE_SERIF_BOLD          },
    { PL_FCI_MARK | 0x102, (const uint8_t *) PL_FREETYPE_MONO_BOLD           },
    { PL_FCI_MARK | 0x103, (const uint8_t *) PL_FREETYPE_SCRIPT_BOLD         },
    { PL_FCI_MARK | 0x104, (const uint8_t *) PL_FREETYPE_SYMBOL_BOLD         },
    { PL_FCI_MARK | 0x110, (const uint8_t *) PL_FREETYPE_SANS_BOLD_ITALIC    },
    { PL_FCI_MARK | 0x111, (const uint8_t *) PL_FREETYPE_SERIF_BOLD_ITALIC   },
    { PL_FCI_MARK | 0x112, (const uint8_t *) PL_FREETYPE_MONO_BOLD_ITALIC    },
    { PL_FCI_MARK | 0x113, (const uint8_t *) PL_FREETYPE_SCRIPT_BOLD_ITALIC  },
    { PL_FCI_MARK | 0x114, (const uint8_t *) PL_FREETYPE_SYMBOL_BOLD_ITALIC  },
    { PL_FCI_MARK | 0x120, (const uint8_t *) PL_FREETYPE_SANS_BOLD_OBLIQUE   },
    { PL_FCI_MARK | 0x121, (const uint8_t *) PL_FREETYPE_SERIF_BOLD_OBLIQUE  },
    { PL_FCI_MARK | 0x122, (const uint8_t *) PL_FREETYPE_MONO_BOLD_OBLIQUE   },
    { PL_FCI_MARK | 0x123, (const uint8_t *) PL_FREETYPE_SCRIPT_BOLD_OBLIQUE },
    { PL_FCI_MARK | 0x124, (const uint8_t *) PL_FREETYPE_SYMBOL_BOLD_OBLIQUE }
};
