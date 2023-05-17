/*
Copyright © 2017, 2018, 2019 Andreas Enge <andreas.enge@inria.fr>

This file is part of PariTwine.

PariTwine is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

PariTwine is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PariTwine.  If not, see <http://www.gnu.org/licenses/>.
*/

AC_PREREQ(2.61)
AC_INIT(paritwine, 0.1, andreas.enge@inria.fr)
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([src/paritwine.h])
AC_CONFIG_HEADER([src/paritwine-config.h])

AM_INIT_AUTOMAKE

AC_CANONICAL_HOST

# Check for programs.
AC_PROG_CC
AC_LANG(C)

# Set up LibTool.
AC_PROG_LIBTOOL

# Extra arguments to configure to determine locations of libraries and headers.
AC_ARG_WITH([pari_include],
            [AC_HELP_STRING([--with-pari-include=DIR],
                            [PARI include directory])],
            [CPPFLAGS="-I$withval $CPPFLAGS"])
AC_ARG_WITH([pari_lib],
            [AC_HELP_STRING([--with-pari-lib=DIR],
                            [PARI lib directory])],
            [LDFLAGS="-L$withval $LDFLAGS"])
AC_ARG_WITH([pari],
            [AC_HELP_STRING([--with-pari=DIR],
                           [PARI install directory])],
            [
            if test -z "$with_pari_lib" -a -z "$with_pari_include" ; then
               CPPFLAGS="-I$withval/include $CPPFLAGS"
               LDFLAGS="-L$withval/lib $LDFLAGS"
            else
               AC_MSG_FAILURE([Do not use --with-pari and --with-pari-include/--with-pari-lib options simultaneously.])
            fi
            ])
AC_ARG_WITH([gmp_include],
            [AC_HELP_STRING([--with-gmp-include=DIR],
                            [GMP include directory])],
            [CPPFLAGS="-I$withval $CPPFLAGS"])
AC_ARG_WITH([gmp_lib],
            [AC_HELP_STRING([--with-gmp-lib=DIR],
                            [GMP lib directory])],
            [LDFLAGS="-L$withval $LDFLAGS"])
AC_ARG_WITH([gmp],
            [AC_HELP_STRING([--with-gmp=DIR],
                           [GMP install directory])],
            [
            if test -z "$with_gmp_lib" -a -z "$with_gmp_include" ; then
               CPPFLAGS="-I$withval/include $CPPFLAGS"
               LDFLAGS="-L$withval/lib $LDFLAGS"
            else
               AC_MSG_FAILURE([Do not use --with-gmp and --with-gmp-include/--with-gmp-lib options simultaneously.])
            fi
            ])
AC_ARG_WITH([mpfr_include],
            [AC_HELP_STRING([--with-mpfr-include=DIR],
                            [MPFR include directory])],
            [CPPFLAGS="-I$withval $CPPFLAGS"])
AC_ARG_WITH([mpfr_lib],
            [AC_HELP_STRING([--with-mpfr-lib=DIR],
                            [MPFR lib directory])],
            [LDFLAGS="-L$withval $LDFLAGS"])
AC_ARG_WITH([mpfr],
            [AC_HELP_STRING([--with-mpfr=DIR],
                           [MPFR install directory])],
            [
            if test -z "$with_mpfr_include" -a -z "$with_mpfr_lib" ; then
               CPPFLAGS="-I$withval/include $CPPFLAGS"
               LDFLAGS="-L$withval/lib $LDFLAGS"
            else
               AC_MSG_FAILURE([Do not use --with-mpfr and --with-mpfr-include/--with-mpfr-lib options simultaneously.])
            fi
            ])
AC_ARG_WITH([mpc_include],
            [AC_HELP_STRING([--with-mpc-include=DIR],
                            [MPC include directory])],
            [CPPFLAGS="-I$withval $CPPFLAGS"])
AC_ARG_WITH([mpc_lib],
            [AC_HELP_STRING([--with-mpc-lib=DIR],
                            [MPC lib directory])],
            [LDFLAGS="-L$withval $LDFLAGS"])
AC_ARG_WITH([mpc],
            [AC_HELP_STRING([--with-mpc=DIR],
                           [MPC install directory])],
            [
            if test -z "$with_mpc_include" -a -z "$with_mpc_lib" ; then
               CPPFLAGS="-I$withval/include $CPPFLAGS"
               LDFLAGS="-L$withval/lib $LDFLAGS"
            else
               AC_MSG_FAILURE([Do not use --with-mpc and --with-mpc-include/--with-mpc-lib options simultaneously.])
            fi
            ])
AC_ARG_WITH([flint_include],
            [AC_HELP_STRING([--with-flint-include=DIR],
                            [FLINT include directory])],
            [CPPFLAGS="-I$withval $CPPFLAGS"])
AC_ARG_WITH([flint_lib],
            [AC_HELP_STRING([--with-flint-lib=DIR],
                            [FLINT lib directory])],
            [LDFLAGS="-L$withval $LDFLAGS"])
AC_ARG_WITH([flint],
            [AC_HELP_STRING([--with-flint=DIR],
                           [FLINT install directory])],
            [
            if test -z "$with_flint_include" -a -z "$with_flint_lib" ; then
               CPPFLAGS="-I$withval/include $CPPFLAGS"
               LDFLAGS="-L$withval/lib $LDFLAGS"
            else
               AC_MSG_FAILURE([Do not use --with-flint and --with-flint-include/--with-flint-lib options simultaneously.])
            fi
            ])
AC_ARG_WITH([arb_include],
            [AC_HELP_STRING([--with-arb-include=DIR],
                            [ARB include directory])],
            [CPPFLAGS="-I$withval $CPPFLAGS"])
AC_ARG_WITH([arb_lib],
            [AC_HELP_STRING([--with-arb-lib=DIR],
                            [ARB lib directory])],
            [LDFLAGS="-L$withval $LDFLAGS"])
AC_ARG_WITH([arb],
            [AC_HELP_STRING([--with-arb=DIR],
                           [ARB install directory])],
            [
            if test -z "$with_arb_include" -a -z "$with_arb_lib" ; then
               CPPFLAGS="-I$withval/include $CPPFLAGS"
               LDFLAGS="-L$withval/lib $LDFLAGS"
            else
               AC_MSG_FAILURE([Do not use --with-arb and --with-arb-include/--with-arb-lib options simultaneously.])
            fi
            ])
AC_ARG_WITH([cmh_include],
            [AC_HELP_STRING([--with-cmh-include=DIR],
                            [CMH include directory])],
            [CPPFLAGS="-I$withval $CPPFLAGS"])
AC_ARG_WITH([cmh_lib],
            [AC_HELP_STRING([--with-cmh-lib=DIR],
                            [CMH lib directory])],
            [LDFLAGS="-L$withval $LDFLAGS"])
AC_ARG_WITH([cmh],
            [AC_HELP_STRING([--with-cmh=DIR],
                           [CMH install directory])],
            [
            if test -z "$with_cmh_include" -a -z "$with_cmh_lib" ; then
               CPPFLAGS="-I$withval/include $CPPFLAGS"
               LDFLAGS="-L$withval/lib $LDFLAGS"
            else
               AC_MSG_FAILURE([Do not use --with-cmh and --with-cmh-include/--with-cmh-lib options simultaneously.])
            fi
            ])


# Check for libraries and their headers.
AC_CHECK_HEADER([pari/pari.h], ,
   [AC_MSG_ERROR([pari headers are not found.])]
)
AC_CHECK_LIB([pari],[pari_init], ,
   [AC_MSG_ERROR([pari is not found.])]
)
AC_CHECK_HEADER([gmp.h], ,
   [AC_MSG_ERROR([gmp headers are not found.])]
)
AC_CHECK_LIB([gmp],[__gmpz_init], ,
   [AC_MSG_ERROR([gmp is not found.])]
)
GP_MPFR=$srcdir/gp/paritwine-mpfr.gp
AC_CHECK_HEADER([mpfr.h],
   [AC_CHECK_LIB([mpfr], [mpfr_init], ,
      [AC_MSG_WARN([mpfr is not found, interface will not be built.])
       GP_MPFR=/dev/null])
   ],
   [AC_MSG_WARN([mpfr headers are not found, interface will not be built.])
    GP_MPFR=/dev/null]
)
GP_MPC=$srcdir/gp/paritwine-mpc.gp
AC_CHECK_HEADER([mpc.h],
   [AC_CHECK_LIB([mpc], [mpc_init2], ,
      [AC_MSG_WARN([mpc is not found, interface will not be built.])
       GP_MPC=/dev/null])
   ],
   [AC_MSG_WARN([mpc headers are not found, interface will not be built.])
    GP_MPC=/dev/null]
)
GP_FLINT=$srcdir/gp/paritwine-flint.gp
AC_CHECK_HEADER([flint/flint.h],
   [AC_CHECK_LIB([flint], [fmpz_init], ,
      [AC_MSG_WARN([flint is not found, interface will not be built.])
       GP_FLINT=/dev/null])
   ],
   [AC_MSG_WARN([flint headers are not found, interface will not be built.])
    GP_FLINT=/dev/null]
)
GP_ARB=$srcdir/gp/paritwine-arb.gp
AC_CHECK_HEADER([arb.h],
   [AC_CHECK_LIB([arb], [arb_init], ,
      [AC_MSG_WARN([arb is not found, interface will not be built.])
       GP_ARB=/dev/null])
   ],
   [AC_MSG_WARN([arb headers are not found, interface will not be built.])
    GP_ARB=/dev/null]
)
GP_CMH=$srcdir/gp/paritwine-cmh.gp
AC_CHECK_HEADER([cmh.h],
   [AC_CHECK_LIB([cmh], [eval_10theta2_newton], ,
      [AC_MSG_WARN([cmh is not found, interface will not be built.])
       GP_CMH=/dev/null])
   ],
   [AC_MSG_WARN([cmh headers are not found, interface will not be built.])
    GP_CMH=/dev/null]
)

AC_SUBST_FILE(GP_MPFR)
AC_SUBST_FILE(GP_MPC)
AC_SUBST_FILE(GP_FLINT)
AC_SUBST_FILE(GP_ARB)
AC_SUBST_FILE(GP_CMH)
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile gp/Makefile gp/paritwine.gp doc/Makefile])
AC_OUTPUT

