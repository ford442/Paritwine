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
   mpfr library. They take as additional parameter a long "prec",
   coding a pari floating point precision which is used for the result.
   Floating point arguments are first transformed into variables of mpfr
   type, using their own precision; the parameter prec is also used to
   determine the floating point precision of exact arguments of type t_INT
   or t_FRAC, for which there is no intrinsic way of deriving a precision.

   If such a function is used inside GP, it should be installed with a
   parameter code of "p" for the prec parameter, which is then automagically
   replaced by the current default precision.
*/

/****************************************************************************/

/* Macro for wrapping a function with one MPFR input and one MPFR output */
#define WRAP_F_F(name) \
GEN pari_mpfr_ ## name (GEN x, long prec) \
{ \
   pari_sp ltop = avma; \
   mpfr_prec_t p = prec; \
   mpfr_t z, z1; \
   pari_mpfr_init2 (z, p); \
   pari_mpfr_init_set_GEN (z1, x, p); \
   mpfr_ ## name (z, z1, MPFR_RNDN); \
   return gerepileuptoleaf (ltop, mpfr_get_GEN (z)); \
}

/* Macro for wrapping a function with two MPFR inputs and one MPFR output */
#define WRAP_F_FF(name) \
GEN pari_mpfr_ ## name (GEN x, GEN y, long prec) \
{ \
   pari_sp ltop = avma; \
   mpfr_prec_t p = prec; \
   mpfr_t z, z1, z2; \
   pari_mpfr_init2 (z, p); \
   pari_mpfr_init_set_GEN (z1, x, p); \
   pari_mpfr_init_set_GEN (z2, y, p); \
   mpfr_ ## name (z, z1, z2, MPFR_RNDN); \
   return gerepileuptoleaf (ltop, mpfr_get_GEN (z)); \
}

/* Macro for wrapping a function with three MPFR inputs and one MPFR output */
#define WRAP_F_FFF(name) \
GEN pari_mpfr_ ## name (GEN x, GEN y, GEN t, long prec) \
{ \
   pari_sp ltop = avma; \
   mpfr_prec_t p = prec; \
   mpfr_t z, z1, z2, z3; \
   pari_mpfr_init2 (z, p); \
   pari_mpfr_init_set_GEN (z1, x, p); \
   pari_mpfr_init_set_GEN (z2, y, p); \
   pari_mpfr_init_set_GEN (z3, t, p); \
   mpfr_ ## name (z, z1, z2, z3, MPFR_RNDN); \
   return gerepileuptoleaf (ltop, mpfr_get_GEN (z)); \
}

/* Macro for wrapping a function with one unsigned long int input
   and one MPFR output */
#define WRAP_F_u(name) \
GEN pari_mpfr_ ## name (unsigned long int i, long prec) \
{ \
   pari_sp ltop = avma; \
   mpfr_prec_t p = prec; \
   mpfr_t z; \
   pari_mpfr_init2 (z, p); \
   mpfr_ ## name (z, i, MPFR_RNDN); \
   return gerepileuptoleaf (ltop, mpfr_get_GEN (z)); \
}

/* Macro for wrapping a function with one long int and one MPFR input
   and one MPFR output */
#define WRAP_F_lF(name) \
GEN pari_mpfr_ ## name (long int i, GEN x, long prec) \
{ \
   pari_sp ltop = avma; \
   mpfr_prec_t p = prec; \
   mpfr_t z, z1; \
   pari_mpfr_init2 (z, p); \
   pari_mpfr_init_set_GEN (z1, x, p); \
   mpfr_ ## name (z, i, z1, MPFR_RNDN); \
   return gerepileuptoleaf (ltop, mpfr_get_GEN (z)); \
}

/****************************************************************************/

#ifdef HAVE_LIBMPFR

/* Basic arithmetic functions */
WRAP_F_FF(add)
WRAP_F_FF(sub)
WRAP_F_FF(mul)
WRAP_F_F(sqr)
WRAP_F_FF(div)
WRAP_F_F(sqrt)
WRAP_F_F(rec_sqrt)
WRAP_F_F(cbrt)
WRAP_F_FF(pow)

/* Special functions */
WRAP_F_F(log)
WRAP_F_F(log2)
WRAP_F_F(log10)
WRAP_F_F(exp)
WRAP_F_F(exp2)
WRAP_F_F(exp10)
WRAP_F_F(cos)
WRAP_F_F(sin)
WRAP_F_F(tan)
WRAP_F_F(sec)
WRAP_F_F(csc)
WRAP_F_F(cot)
WRAP_F_F(acos)
WRAP_F_F(asin)
WRAP_F_F(atan)
WRAP_F_F(cosh)
WRAP_F_F(sinh)
WRAP_F_F(tanh)
WRAP_F_F(sech)
WRAP_F_F(csch)
WRAP_F_F(coth)
WRAP_F_F(acosh)
WRAP_F_F(asinh)
WRAP_F_F(atanh)
WRAP_F_u(fac_ui)
WRAP_F_F(log1p)
WRAP_F_F(expm1)
WRAP_F_F(eint)
WRAP_F_F(li2)
WRAP_F_F(gamma)
WRAP_F_F(lngamma)
WRAP_F_F(digamma)
WRAP_F_F(zeta)
WRAP_F_F(erf)
WRAP_F_F(erfc)
WRAP_F_F(j0)
WRAP_F_F(j1)
WRAP_F_lF(jn)
WRAP_F_F(y0)
WRAP_F_F(y1)
WRAP_F_lF(yn)
WRAP_F_FFF(fma)
WRAP_F_FFF(fms)
WRAP_F_FF(agm)
WRAP_F_FF(hypot)
WRAP_F_F(ai)

#endif

