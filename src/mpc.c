/*
Copyright © 2014, 2016, 2017, 2018 Andreas Enge <andreas.enge@inria.fr>

This file is part of paritwine.

Paritwine is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

Paritwine is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Paritwine.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "paritwine.h"

/* The following functions perform floating point computations using the
   mpc library. For details on the "prec" parameter, see the comments in
   the file mpfr.c.
*/

/****************************************************************************/

/* Macro for wrapping a function with one MPC input and one MPC output */
#define WRAP_C_C(name) \
GEN pari_mpc_ ## name (GEN x, long prec) \
{ \
   pari_sp ltop = avma; \
   mpfr_prec_t p = prec; \
   mpc_t z, z1; \
   pari_mpc_init2 (z, p); \
   pari_mpc_init_set_GEN (z1, x, p); \
   mpc_ ## name (z, z1, MPC_RNDNN); \
   return gerepileupto (ltop, mpc_get_GEN (z)); \
}

/* Macro for wrapping a function with two MPC inputs and one MPC output */
#define WRAP_C_CC(name) \
GEN pari_mpc_ ## name (GEN x, GEN y, long prec) \
{ \
   pari_sp ltop = avma; \
   mpfr_prec_t p = prec; \
   mpc_t z, z1, z2; \
   pari_mpc_init2 (z, p); \
   pari_mpc_init_set_GEN (z1, x, p); \
   pari_mpc_init_set_GEN (z2, y, p); \
   mpc_ ## name (z, z1, z2, MPC_RNDNN); \
   return gerepileupto (ltop, mpc_get_GEN (z)); \
}

/* Macro for wrapping a function with three MPC inputs and one MPC output */
#define WRAP_C_CCC(name) \
GEN pari_mpc_ ## name (GEN x, GEN y, GEN t, long prec) \
{ \
   pari_sp ltop = avma; \
   mpfr_prec_t p = prec; \
   mpc_t z, z1, z2, z3; \
   pari_mpc_init2 (z, p); \
   pari_mpc_init_set_GEN (z1, x, p); \
   pari_mpc_init_set_GEN (z2, y, p); \
   pari_mpc_init_set_GEN (z3, t, p); \
   mpc_ ## name (z, z1, z2, z3, MPC_RNDNN); \
   return gerepileupto (ltop, mpc_get_GEN (z)); \
}

/* Macro for wrapping a function with one MPC input and one MPFR output */
#define WRAP_F_C(name) \
GEN pari_mpc_ ## name (GEN x, long prec) \
{ \
   pari_sp ltop = avma; \
   mpfr_prec_t p = prec; \
   mpfr_t z; \
   mpc_t z1; \
   pari_mpfr_init2 (z, p); \
   pari_mpc_init_set_GEN (z1, x, p); \
   mpc_ ## name (z, z1, MPC_RNDNN); \
   return gerepileupto (ltop, mpfr_get_GEN (z)); \
}

/* Macro for wrapping a function with two unsigned long int inputs
   and one MPC output */
#define WRAP_C_uu(name) \
GEN pari_mpc_ ## name (unsigned long int i, unsigned long int j, long prec) \
{ \
   pari_sp ltop = avma; \
   mpfr_prec_t p = prec; \
   mpc_t z; \
   pari_mpc_init2 (z, p); \
   mpc_ ## name (z, i, j, MPC_RNDNN); \
   return gerepileupto (ltop, mpc_get_GEN (z)); \
}

/****************************************************************************/

#ifdef HAVE_LIBMPC

/* Basic arithmetic functions */
WRAP_C_CC(add)
WRAP_C_CC(sub)
WRAP_C_CC(mul)
WRAP_C_C(sqr)
WRAP_C_CCC(fma)
WRAP_C_CC(div)
WRAP_F_C(abs)
WRAP_F_C(norm)

/* Power functions and logarithm */
WRAP_C_C(sqrt)
WRAP_C_CC(pow)
WRAP_C_C(exp)
WRAP_C_C(log)
WRAP_C_C(log10)
/*
WRAP_C_uu(rootofunity)
*/

/* Trigonometric functions */
WRAP_C_C(sin)
WRAP_C_C(cos)
WRAP_C_C(tan)
WRAP_C_C(sinh)
WRAP_C_C(cosh)
WRAP_C_C(tanh)
WRAP_C_C(asin)
WRAP_C_C(acos)
WRAP_C_C(atan)
WRAP_C_C(asinh)
WRAP_C_C(acosh)
WRAP_C_C(atanh)

#endif

