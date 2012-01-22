// ***************************************************************************
// bamtools_global.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 19 November 2010 (DB)
// ---------------------------------------------------------------------------
// Provides the basic definitions for exporting & importing library symbols
// ***************************************************************************

#ifndef BAMTOOLS_GLOBAL_H
#define BAMTOOLS_GLOBAL_H

// BAMTOOLS_LIBRARY_EXPORT
#ifndef BAMTOOLS_LIBRARY_EXPORT
#  if defined(WIN32)
#    define BAMTOOLS_LIBRARY_EXPORT __declspec(dllexport)
#  else
#    define BAMTOOLS_LIBRARY_EXPORT __attribute__((visibility("default")))
#  endif
#endif // BAMTOOLS_LIBRARY_EXPORT

// BAMTOOLS_LIBRARY_IMPORT
#ifndef BAMTOOLS_LIBRARY_IMPORT
#  if defined(WIN32)
#    define BAMTOOLS_LIBRARY_IMPORT __declspec(dllimport)
#  else
#    define BAMTOOLS_LIBRARY_IMPORT
#  endif
#endif // BAMTOOLS_LIBRARY_IMPORT

#ifndef HAVE_FSEEK64
#define ftell64(a)     ftello(a)
#define fseek64(a,b,c) fseeko(a,b,c)
#endif

#ifndef BAMTOOLS_ASSERTS
#define BAMTOOLS_ASSERTS
# ifdef NDEBUG
# define BT_ASSERT_UNREACHABLE bamtools_noop()
# define BT_ASSERT_X( condition, message ) bamtools_noop()
# else
# include <cassert>
# include <stdexcept>
# define BT_ASSERT_UNREACHABLE assert( false )
# define BT_ASSERT_X( condition, message ) if (!( condition )) { throw std::runtime_error( message ); }
# endif
#endif // BAMTOOLS_ASSERTS

#endif // BAMTOOLS_GLOBAL_H

