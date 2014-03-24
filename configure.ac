# Initialise
AC_INIT([The EIC smearing package], [1.0.3], [tpb@bnl.gov], [eic-smear])

# Runs many macros required for proper operation of the generated Makefiles.
# Use foreign to ignore absent ChangeLog, AUTHORS etc.
AM_INIT_AUTOMAKE([subdir-objects -Wall foreign])

# Test for the existence of this file as
# a check against incorrect --srcdir.
AC_CONFIG_SRCDIR([src/erhic/BuildTree.cxx])

AC_CONFIG_HEADER([config.h])

# Change the default behaviour for LT_INIT to enable shared libraries.
# The user may still override this default by specifying ‘--disable-shared’.
AC_ENABLE_SHARED

# Change the default behaviour for LT_INIT to disable static libraries.
# The user may still override this default by specifying ‘--enable-static’.
AC_DISABLE_STATIC

# Checks for programs.
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
m4_pattern_allow([AM_PROG_AR])
AM_PROG_AR

# Enable support for --enable-shared and --disable-shared configure flags.
AC_PROG_LIBTOOL

# Set flags, including ROOT paths via root-config.
CFLAGS="-Wall -Wextra ${CFLAGS}"
CXXFLAGS="-Wall -Wextra ${CXXFLAGS}"
LDFLAGS="${LDFLAGS} -L`root-config --libdir`"
LIBS="${LIBS}"
CPPFLAGS="${CPPFLAGS} -I`root-config --incdir`"
FFLAGS="${FFLAGS}"

# Allow compilation of PYTHIA6-related code.
# Without --with-pythia6-libdir those components are not built.
AC_ARG_WITH(
   [pythia6-libdir],
   [AS_HELP_STRING(
      [--with-pythia6-libdir=DIR],
      [directory containing libPythia6]
   )],
   [PYTHIA6_LIBS="-lEGPythia6 -L$withval -lPythia6"],
   [PYTHIA6_LIBS=""]
)
# ENABLE_PYTHIA6_SUPPORT conditional is tested in Makefile.am.
AM_CONDITIONAL([ENABLE_PYTHIA6_SUPPORT], [test x"$PYTHIA6_LIBS" != x""])
AC_SUBST(PYTHIA6_LIBS)

AC_CONFIG_FILES([Makefile])

# And away we go...
AC_OUTPUT
