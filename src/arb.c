/*
Copyright © 2014, 2016, 2017, 2018, 2019 Andreas Enge <andreas.enge@inria.fr>
Copyright © 2017, 2018 Fredrik Johansson <fredrik.johansson@gmail.com>

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

/****************************************************************************/

/* Macro for wrapping a function with one ACB input and one ACB output */
#define WRAP_Acb_Acb(name) \
GEN pari_acb_ ## name (GEN x, long prec) \
{ \
    acb_t z1; \
    GEN z; \
    acb_init(z1); \
    acb_set_GEN(z1, x, prec); \
    acb_ ## name(z1, z1, prec); \
    z = acb_get_GEN(z1, prec); \
    acb_clear(z1); \
    return z; \
}

/* Macro for wrapping a function with one ACB input and one ACB output,
   but not taking a prec argument*/
#define WRAP_Acb_Acb_noprec(name) \
GEN pari_acb_ ## name (GEN x, long prec) \
{ \
    acb_t z1; \
    GEN z; \
    acb_init(z1); \
    acb_set_GEN(z1, x, prec); \
    acb_ ## name(z1, z1); \
    z = acb_get_GEN(z1, prec); \
    acb_clear(z1); \
    return z; \
}

/* Macro for wrapping a function with two ACB inputs and one ACB output */
#define WRAP_Acb_AcbAcb(name) \
GEN pari_acb_ ## name (GEN x, GEN y, long prec) \
{ \
    acb_t z1, z2; \
    GEN z; \
    acb_init(z1); \
    acb_init(z2); \
    acb_set_GEN(z1, x, prec); \
    acb_set_GEN(z2, y, prec); \
    acb_ ## name(z1, z1, z2, prec); \
    z = acb_get_GEN(z1, prec); \
    acb_clear(z1); \
    acb_clear(z2); \
    return z; \
}

/****************************************************************************/

#ifdef HAVE_LIBARB

/* Misc. functions. */
GEN pari_fmpz_numbpart (GEN x) \
{ \
    fmpz_t z1; \
    GEN z; \
    fmpz_init(z1); \
    fmpz_set_GEN(z1, x); \
    partitions_fmpz_fmpz(z1, z1, 0); \
    z = fmpz_get_GEN(z1); \
    fmpz_clear(z1); \
    return z; \
}

GEN pari_arb_numbpart (GEN x, long prec) \
{ \
    fmpz_t z1; \
    arb_t z2;
    GEN z; \
    fmpz_init(z1); \
    arb_init(z2);
    fmpz_set_GEN(z1, x); \
    arb_partitions_fmpz(z2, z1, prec); \
    z = arb_get_GEN(z2, prec); \
    fmpz_clear(z1); \
    arb_clear(z2);
    return z; \
}

/****************************************************************************/

GEN pari_acb_modular_eisenstein (GEN tau, long n, long prec)
{
   acb_t tau1;
   acb_ptr r;
   GEN v;
   int i;

   acb_init(tau1);
   acb_set_GEN(tau1, tau, prec);
   r = _acb_vec_init(n);

   acb_modular_eisenstein(r, tau1, n, prec);

   v = cgetg(n+1, t_VEC);
   for (i = 1; i <= n; i++)
      gel(v, i) = acb_get_GEN(r+i-1, prec);

   acb_clear(tau1);
   _acb_vec_clear(r, n);

   return v;
}

/****************************************************************************/

GEN pari_acb_modular_theta (GEN z, GEN tau, long prec)
{
   acb_t t1, t2, t3, t4, z1, tau1;
   GEN v;

   acb_init(t1);
   acb_init(t2);
   acb_init(t3);
   acb_init(t4);
   acb_init(z1);
   acb_init(tau1);
   acb_set_GEN(z1, z, prec);
   acb_set_GEN(tau1, tau, prec);

   acb_modular_theta(t1, t2, t3, t4, z1, tau1, prec);

   v = cgetg(5, t_VEC);
   gel(v, 1) = acb_get_GEN(t1, prec);
   gel(v, 2) = acb_get_GEN(t2, prec);
   gel(v, 3) = acb_get_GEN(t3, prec);
   gel(v, 4) = acb_get_GEN(t4, prec);

   acb_clear(t1);
   acb_clear(t2);
   acb_clear(t3);
   acb_clear(t4);
   acb_clear(z1);
   acb_clear(tau1);

   return v;
}

/****************************************************************************/

WRAP_Acb_AcbAcb(add)
WRAP_Acb_AcbAcb(sub)
WRAP_Acb_AcbAcb(mul)
WRAP_Acb_AcbAcb(div)
WRAP_Acb_Acb_noprec(neg)
WRAP_Acb_Acb_noprec(conj)
WRAP_Acb_Acb(exp)
WRAP_Acb_Acb(log)
WRAP_Acb_Acb(atan)
WRAP_Acb_Acb(sin)
WRAP_Acb_Acb(cos)
WRAP_Acb_Acb(gamma)
WRAP_Acb_Acb(digamma)
WRAP_Acb_Acb(zeta)
WRAP_Acb_AcbAcb(hurwitz_zeta)
WRAP_Acb_Acb(modular_eta)
WRAP_Acb_Acb(modular_j)
WRAP_Acb_Acb(modular_delta)
WRAP_Acb_AcbAcb(elliptic_p)
WRAP_Acb_AcbAcb(elliptic_inv_p)
WRAP_Acb_AcbAcb(elliptic_zeta)
WRAP_Acb_AcbAcb(elliptic_sigma)
#endif

