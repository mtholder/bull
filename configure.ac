#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(BULL, 1.1, BUG-REPORT-ADDRESS)
AC_CONFIG_SRCDIR([src/bull.cpp])
AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE([1.9 foreign dist-tarZ tar-ustar filename-length-max=299])
# Checks for programs.
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX

# required because we are building a library 
AC_PROG_RANLIB
AC_PROG_INSTALL


# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h malloc.h stddef.h stdlib.h sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor memmove memset pow sqrt strchr strdup strtol])



#--------------------------------------------------------------------------------------------------#
#  Set CXXFLAGS                                                                                    #
#--------------------------------------------------------------------------------------------------#

# Initialize CXXFLAGS to an empty string to prevent it from defaulting to "-g -O2"
CXXFLAGS=""

# Check for debugging mode.
AC_ARG_ENABLE(debugging,
	AC_HELP_STRING([--enable-debugging],[build for debugging]), , [enable_debugging=no])
if test "$enable_debugging" = yes; then
	AC_MSG_NOTICE([
	   *** NOTE: debugging is enabled; optimization is suppressed!
	   ])
fi

# Initialize optimization flag in case it doesn't get set below.
CXXFLAGS_OPTIM_SPEED="-O"

# "-g" may not work with some compilers, but end users shouldn't be
if test "$enable_debugging" = yes; then
	CXXFLAGS_OPTIM_SPEED="-O0"
	CXXFLAGS="$CXXFLAGS -Wall -g -Wreturn-type -Wunused -Wredundant-decls -Wcast-align -Wcomment -Wextra"
else
	CXXFLAGS="$CXXFLAGS -DNDEBUG"
fi

if test "$CC" = "icc"; then
#	Intel C compiler for Linux
	CXXFLAGS="$CXXFLAGS -ip -fno-alias -DUNIX"
	if test "$enable_debugging" = no; then
		CXXFLAGS_OPTIM_SPEED="-O2"
	fi
elif test "$CC" = "ccc"; then
#	Compaq C compiler for Linux
	if test "x$arch" = "x"; then
		arch="host"
	fi
	if test "$enable_debugging" = no; then
		CXXFLAGS_OPTIM_SPEED="-fast -inline speed -arch $arch"
	fi
elif test "$CC" = "xlc"; then
#	IBM XL C compiler
	if test "x$arch" = "x"; then
		arch="auto"
	fi
	if test "$enable_debugging" = no; then
		CXXFLAGS_OPTIM_SPEED="-O3 -qarch=$arch"
	fi
elif test "x$GCC" = "xyes" ; then
	CXXFLAGS="$CXXFLAGS -Wno-uninitialized -fstrict-aliasing -fomit-frame-pointer -DUNIX"
	if test "$enable_debugging" = yes; then
		CXXFLAGS_OPTIM_SPEED="-O0 -Wimplicit"
	else
		CXXFLAGS_OPTIM_SPEED="-O3 -ffast-math -funroll-loops"
	fi
	case "$build_os" in
	 	darwin*) CXXFLAGS="$CXXFLAGS -Wno-long-double -include defs.h" ;;
	 	*) CXXFLAGS="$CXXFLAGS -fsigned-char";;
	esac	
fi

CXXFLAGS="$CXXFLAGS $CXXFLAGS_OPTIM_SPEED"


################################################################################
#	By default, use the bundled ncl code from source, 
#		but let the --with-ncl argument override this.
################################################################################
MAYBE_NCL=""
AC_ARG_WITH(
	[ncl],
	AC_HELP_STRING(
		[--with-ncl=DIR],
		[Specify the root directory for ncl library (parent of the include/ncl and lib directories). Omit this flag to build and use the bundled version of ncl.]
		),
	[
	if test "$withval" = "no"; then
		AC_MSG_ERROR([The ncl library is a prerequisite.  The --with-ncl option should be used to specify the parent of the include directory where ncl was installed])
	elif ! test "$withval" = "yes"; then
		NCL_INC_DIR="$withval/include"
		NCL_LIB_DIR="$withval/lib"
	fi
	],
	[
	MAYBE_NCL="ncl"
	AC_MSG_ERROR([The ncl library is a prerequisite.  The --with-ncl option should be used to specify the parent of the include directory where ncl was installed.])
	]
	)

if ! test -d "$NCL_INC_DIR/ncl" ; then
	AC_MSG_ERROR([Expecting a directory called $NCL_INC_DIR/ncl, verify that the argument to --with-ncl is correct. You may need to manually make the directory $NCL_INC_DIR/ncl copy all of the ncl headers into that directory if you are using NCL 2.0])
fi
if ! test -d "$NCL_LIB_DIR" ; then
	AC_MSG_ERROR([Expecting a directory called $NCL_LIB_DIR, verify that the argument to --with-ncl is correct.])
fi

if test -z "$MAYBE_NCL" ; 
then
	INCLUDE="$INCLUDE -I$NCL_INC_DIR -I$NCL_INC_DIR/ncl"
	if test -d "$NCL_LIB_DIR/ncl" ; then
		LDFLAGS="$LDFLAGS -L$NCL_LIB_DIR/ncl"
	else
		LDFLAGS="$LDFLAGS -L$NCL_LIB_DIR"
	fi
fi
BULL_LIBS="-lncl"


AC_SUBST([BULL_LIBS])

AC_SUBST([INCLUDE])
AC_SUBST([LDFLAGS])

AC_CONFIG_FILES([
	Makefile			\
	src/Makefile \
	test/Makefile 
	])
AC_OUTPUT
