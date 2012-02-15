# Initialise
AC_INIT([The EIC smearing package], [2.0], [tpb@bnl.gov],[eic-smear])

# Runs many macros required for proper operation of the generated Makefiles.
# Report automake warnings as errors.
# Use foreign to ignore absent ChangeLog, AUTHORS etc.
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

#AC_CONFIG_AUX_DIR(config)
#AC_CONFIG_SRCDIR(toymc.cc)

# Test for the existence of this file as
# a check against incorrect --srcdir.
AC_CONFIG_SRCDIR([src/BuildTree.cxx])

AC_CONFIG_HEADER([config.h])

# Change the default behaviour for LT_INIT to enable shared libraries.
# The user may still override this default by specifying ‘--disable-shared’.
AC_ENABLE_SHARED

# Change the default behaviour for LT_INIT to disable static libraries.
# The user may still override this default by specifying ‘--enable-static’.
AC_DISABLE_STATIC

# Specify m4 as the location of additional local Autoconf macros.
#AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

# Enable support for --enable-shared and --disable-shared configure flags
AC_PROG_LIBTOOL

CFLAGS="-Wall -Wextra -O2 ${CFLAGS}"
CXXFLAGS="-Wall -Wextra -O2 ${CXXFLAGS}"
LDFLAGS="${LDFLAGS} -L`root-config --libdir`"
LIBS="${LIBS}"
CPPFLAGS="${CPPFLAGS} -I`root-config --incdir`"
FFLAGS="${FFLAGS}"

# Deal with options to disable PYTHIA
#AC_ARG_ENABLE(
#   [pythia6],
#   [AS_HELP_STRING(
#      [--enable-pythia6],
#      [include PYTHIA 6 functionality.\nMust specify the PYTHIA 6 library location with --with-pythia6-libdir]
#   )],
#   [enable_pythia6="yes"],
#   [enable_pythia6="no"]
#)
#AC_SUBST(enable_pythia6)

#AM_CONDITIONAL([ENABLE_PYTHIA6_SUPPORT], [test "X$enable_pythia6" != "Xno"])

AC_ARG_WITH(
   [pythia6-libdir],
   [AS_HELP_STRING(
      [--with-pythia6-libdir=DIR],
      [directory containing libPythia6]
   )],
   [PYTHIA6_LIBS="-L$withval -lPythia6"],
   [PYTHIA6_LIBS=""]
)
#   [PYTHIA6_LIBS="-L$withval -lPythia6"],
#   [PYTHIA6_LIBS=""]
#AM_CONDITIONAL([USE_PYTHIA], [$PYTHIA_LIBS != ""])
AM_CONDITIONAL([ENABLE_PYTHIA6_SUPPORT], [test "X$PYTHIA6_LIBS" != ""])
#AM_CONDITIONAL([USE_PYTHIA], [test ! "x$PYTHIA_LIBS" = ""])
AC_SUBST(PYTHIA6_LIBS)

AC_CONFIG_FILES([Makefile])

# And away we go...
AC_OUTPUT
