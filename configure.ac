#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT(PostNonEq, 0.1, Massimo.Marchi@cea.fr)
AM_INIT_AUTOMAKE(PostNonEq, 0.1)
AC_CONFIG_SRCDIR([src/PostNonEq.cpp])

LT_INIT
AC_PROG_LIBTOOL
##AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])


AC_LANG([C++])
AC_PROG_CXX
AC_PROG_CXX_C_O 
AC_PROG_CXXCPP
AC_PROG_INSTALL

# check for maximal optimization 
AX_COMPILER_VENDOR
AX_CXX_MAXOPT

# Checks for libraries.


AC_ARG_WITH(voro-dir,
        AC_HELP_STRING(
                [--with-voro-dir=dir],
                [Provide the vor++ root directory [/usr/local/voro++] ]),
        [voro_dir=$withval],[voro_dir="/usr/local/voro++"])

AC_ARG_WITH(gmx-dir,
        AC_HELP_STRING(
                [--with-gmx-dir=dir],
                [Provide the gmx root directory [/usr/local/gmx] ]),
        [gmx_dir=$withval],[gmx_dir="/usr/local/gmx"])

AC_ARG_WITH(fftw3-dir,
        AC_HELP_STRING(
                [--with-fftw3-dir=dir],
                [Provide the fftw3 root directory [/usr/local/fftw3] ]),
        [fftw3_dir=$withval],[fftw3_dir="/usr/local/fftw3"])

BOOST_REQUIRE([1.40],,AC_MSG_ERROR([boost 1.40 and up is needed. You must provide a path to the library or use --with-boost=rootdir]))

BOOST_FIND_HEADER([boost/algorithm/string.hpp],AC_MSG_ERROR([boost is needed. You must provide a path to the library or use --with-boost=rootdir]),[])

BOOST_PROGRAM_OPTIONS()

CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS -I$voro_dir/include/voro++ -I$gmx_dir/include/gromacs -I$fftw3_dir/include"
LDFLAGS="$LDFLAGS -L$boost_ldpath -L$voro_dir/lib -L$gmx_dir/lib -L$fftw3_dir/lib"
LIBS="$LIBS -l$boost_lib "

AC_CHECK_LIB([voro++],main,,AC_MSG_ERROR([Voro++ is needed. You must provide a path to the library or use --with-voro-dir=dir]))

AC_CHECK_LIB([fftw3],main,,AC_MSG_ERROR([fftw3 is needed. You must provide a path to the libraryor use --with-fftw3-dir=dir]))


# Checks for header files.

AC_HEADER_STDC   
AC_CHECK_HEADERS([voro++.hh],,AC_MSG_ERROR([No header files found for voro++! You must provide a path to the include files or the root installation directory with --with-voro-dir=root]))

AC_CHECK_HEADERS([typedefs.h],,AC_MSG_ERROR([No header files found for GROMACS! You must provide a path to the include files or the root installation directory with --with-gmx-dir=root]))

AC_CHECK_HEADERS([fftw3.h],,AC_MSG_ERROR([No header files found for fftw3 ! You must provide a path to the include files or the root installation directory with --with-fftw3-dir=root]))


# Check whether the user wants to use a particular MPI installation

AC_ARG_ENABLE(mpi, [AC_HELP_STRING([--enable-mpi],[compile PostNonEq with MPI support])], enable_mpi=$enableval, enable_mpi=no)
AC_DEFINE([__cplusplus],[1],[Use C includes])

if test x"$enable_mpi" = x"yes"; then
   AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
   AC_CHECK_PROGS(MPICXX, mpic++ mpiCC mpicxx mpCC hcp mpxlC mpxlC_r cmpic++, "no")
   AC_SUBST(MPICXX)
   if test x"$MPICXX" = x"no"; then
      AC_MSG_ERROR([mpi is needed])
   fi
   CXX=$MPICXX
   AC_MSG_CHECKING([Searching for MPI libraries])
   AC_TRY_LINK([#include <mpi.h>],[int argc; char **argv; MPI_Init(&argc,&argv);],
   AC_MSG_RESULT([yes]),AC_MSG_ERROR([Cannot compile and link MPI code with $MPICXX]))

   AC_CHECK_PROG(MPIRUN, mpirun, mpirun)
   AC_SUBST(MPIRUN)
   AC_DEFINE([__MPI],1,[Enable MPI compilation])
   AC_CHECK_LIB([gmx_mpi],main,,AC_CHECK_LIB([gmx],main,,AC_MSG_ERROR([GROMACS 4.5.3 is needed. You must provide a path to the library or the root installation directory with --with-gmx-dir=root])))
   AC_DEFINE([GMX_LIB_MPI],1,[Enable MPI compilation])
fi
if test x"$enable_mpi" = x"no"; then
   AC_CHECK_LIB([gmx],main,,AC_MSG_ERROR([GROMACS 4.5.3 is needed. You must provide a path to the library or the root installation directory with --with-gmx-dir=root.]))
fi

# check if gmx is usable 
  AC_MSG_CHECKING([Can it link GROMACS libraries? Quite a paranoic check:])
  AC_TRY_LINK([#include <gmx_fatal.h>],[
   int fatal_errno=0;const char *file="fatal";int line=3;const char *fmt="rw";
   gmx_fatal(fatal_errno,file,line,fmt);],AC_MSG_RESULT([yes]),AC_MSG_ERROR([Cannot compile and link gmx libraries they may not be usable (e.g. the library might have been compiled with a compiler different from current)]))

AM_CONDITIONAL(MPI, test "$enable_mpi" = "yes")

# Checks for programs.


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([remove sqrt strdup mkstemp exit])

AC_CONFIG_FILES([Makefile alglib/Makefile MyLibrary/Makefile VoroElec/Makefile src/Makefile])
AC_OUTPUT

