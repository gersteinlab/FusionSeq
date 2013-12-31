#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([fusionseq], [0.7.0], [david.chen@yale.edu])
AC_CONFIG_SRCDIR([src/geneFusions.c])
AC_CONFIG_SRCDIR([src/cgi/geneFusions_cgi.c])
AC_CONFIG_SRCDIR([src/test/bp2fasta.c])
AC_CONFIG_SRCDIR([src/optional/bp2wig.c])
AC_CONFIG_SRCDIR([scripts/Confp.pm])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE([-Wall -Werror])
LT_INIT

AM_SILENT_RULES([yes])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
ROOT_PATH([], [], [AC_MSG_ERROR([Cannot find [ROOT]])])
AC_PROG_LIBTOOL

# Checks for libraries.
AC_CHECK_LIB([m], [log], [], [AC_MSG_ERROR([Cannot find standard math library])])
AC_CHECK_LIB([gslcblas], [cblas_dgemm], [], [AC_MSG_ERROR([Cannot find cblas library])])
AC_CHECK_LIB([gsl], [gsl_ran_hypergeometric_pdf], [], [AC_MSG_ERROR([Cannot find gsl library])])
AC_CHECK_LIB([bios], [needMem], [], [AC_MSG_ERROR([Cannot find bios library])])
AC_CHECK_LIB([mrf], [mrf_init], [], [AC_MSG_ERROR([Cannot find mrf library])])
AX_CHECK_GD

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])
AC_CHECK_HEADERS([gd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([getcwd pathconf rint setenv strchr strdup strstr])

AC_CONFIG_FILES([Makefile 
		 src/Makefile
		 src/optional/Makefile
		 src/cgi/Makefile
		 src/test/Makefile])
AC_OUTPUT
