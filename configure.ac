# Copyright 1997-2018, Brian Borchers.
# This file is distributed under the Eclipse Public License 2.0.

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_INIT([Csdp],[trunk],[csdp@list.coin-or.org],[csdp],[https://github.com/coin-or/csdp])

AC_COPYRIGHT([Copyright 1997-2018, Brian Borchers.
This copy of CSDP is made available under the Eclipse Public License 2.0.
See LICENCE for the details of the EPL 2.0.])

# List one file in the package so that the configure script can test whether
# the package is actually there.

AC_CONFIG_SRCDIR(lib/sdp.c)

# Do some project-level initialization work (version numbers, ...)

AC_COIN_INITIALIZE(Csdp)

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the name of the C compiler and appropriate compiler options.

AC_COIN_PROG_CC

# This is a C package, set the language accordingly.

AC_LANG_PUSH(C)

# Initialize libtool

AC_COIN_PROG_LIBTOOL

############################################################################
#                        Availability of functions, etc.                   #
############################################################################

AC_CHECK_FUNCS([signal])
AC_CHECK_FUNCS([gettimeofday])

AC_OPENMP

#############################################################################
#                           Library Dependencies                            #
#############################################################################

AC_COIN_CHK_LAPACK(CsdpLib)
if test "$coin_has_lapack" != "yes" ; then
  AC_MSG_ERROR([Lapack required for CSDP])
fi

AC_COIN_CHK_LIBM(CsdpLib)

CSDPLIB_CFLAGS="$CSDPLIB_CFLAGS -DBIT64 -DSETNUMTHREADS"

# Finalize libs and cflags for building CsdpLib.
AC_COIN_FINALIZE_FLAGS(CsdpLib)

AC_LANG_POP(C)

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile lib/Makefile solver/Makefile theta/Makefile test/Makefile example/Makefile])

# Here put the location and name of the configuration header file
#CONFIG_HEADER([src/config.h src/config_coinutils.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
