AC_PREREQ([2.69])
AC_INIT([blz78suf],[1.0.0],[https://github.com/kspalaiologos/blz78suf])
AC_CONFIG_HEADERS([config.h])
AC_COPYRIGHT([Released to the public domain by Kamila Szewczyk])
AC_LANG([C++])
AC_CONFIG_MACRO_DIR([m4])

AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([-Wall color-tests foreign])

AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_CXX

AX_CHECK_COMPILE_FLAG([-std=c++23], [CXXFLAGS="$CXXFLAGS -std=c++23"],
  AC_MSG_ERROR([C++23 support is required.]))

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
