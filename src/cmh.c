/*
Copyright © 2016, 2017 Andreas Enge <andreas.enge@inria.fr>

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
   cmh library. For details on the "prec" parameter, see the comments in
   the file mpfr.c.
*/

/****************************************************************************/

#ifdef HAVE_LIBCMH
GEN pari_cmh_I2I4I6I10 (GEN tau, long prec)
   /* Returns the Igusa-Clebsch invariants I2, I4, I6, I10 in the 2x2-
      period matrix tau of type t_MAT; since tau is assumed to be symmetric,
      its lower left entry is not used,
      The return value is a t_VEC with 4 components.
      Install into GP with
      install ("pari_cmh_I2I4I6I10", "Gp", "cmh_I2I4I6I10", "./libparitwine.so");
   */

{
   mpfr_prec_t p = bit_accuracy (prec);
   mpc_t t [3], theta2 [10], h [4], I [4];
   GEN res;
   int i;

   pari_sp av = avma;

   /* Copy the input. */ 
   pari_mpc_init_set_GEN (t [0], gel (gel (tau, 1), 1), p);
   pari_mpc_init_set_GEN (t [1], gel (gel (tau, 2), 1), p);
   pari_mpc_init_set_GEN (t [2], gel (gel (tau, 2), 2), p);

   /* Initialise intermediate variables. */
   for (i = 0; i < 10; i++)
      pari_mpc_init2 (theta2 [i], p);
   for (i = 0; i < 4; i++)
      pari_mpc_init2 (h [i], p);
   for (i = 0; i < 4; i++)
      pari_mpc_init2 (I [i], p);

   /* Compute the result. */
   eval_10theta2_naive (theta2, t);
   h4h10h12h16_from_th2 (h, theta2);
   I2I4I6I10_from_h4h10h12h16 (I, h);

   /* Create the output variable and copy the result. */
   res = cgetg (5, t_VEC);
   for (i = 0; i < 4; i++)
      gel (res, i+1) = mpc_get_GEN (I [i]);

   return gerepileupto (av, res);
}

/****************************************************************************/

GEN pari_cmh_4theta (GEN tau, long prec)
   /* Returns the first four theta constants in the 2x2-
      period matrix tau of type t_MAT; since tau is assumed to be symmetric,
      its lower left entry is not used,
      The return value is a t_VEC with 4 components.
      Install into GP with
      install ("pari_cmh_4theta", "Gp", "cmh_4theta", "./libparitwine.so");
   */

{
   mpfr_prec_t p = bit_accuracy (prec);
   mpc_t t [3], theta [4];
   GEN res;
   int i;

   pari_sp av = avma;

   /* Copy the input. */ 
   pari_mpc_init_set_GEN (t [0], gel (gel (tau, 1), 1), p);
   pari_mpc_init_set_GEN (t [1], gel (gel (tau, 2), 1), p);
   pari_mpc_init_set_GEN (t [2], gel (gel (tau, 2), 2), p);

   /* Initialise intermediate variables. */
   for (i = 0; i < 4; i++)
      pari_mpc_init2 (theta [i], p);

   /* Compute the result. */
   eval_4theta_naive (theta, t);

   /* Create the output variable and copy the result. */
   res = cgetg (5, t_VEC);
   for (i = 0; i < 4; i++)
      gel (res, i+1) = mpc_get_GEN (theta [i]);

   return gerepileupto (av, res);
}

/****************************************************************************/

GEN pari_cmh_10theta2 (GEN tau, long prec)
   /* Returns the squares of the ten theta constants in the 2x2-
      period matrix tau of type t_MAT; since tau is assumed to be symmetric,
      its lower left entry is not used,
      The return value is a t_VEC with 10 components.
      Install into GP with
      install ("pari_cmh_10theta2", "Gp", "cmh_10theta2", "./libparitwine.so");
   */

{
   mpfr_prec_t p = bit_accuracy (prec);
   mpc_t t [3], theta2 [10];
   GEN res;
   int i;

   pari_sp av = avma;

   /* Copy the input. */ 
   pari_mpc_init_set_GEN (t [0], gel (gel (tau, 1), 1), p);
   pari_mpc_init_set_GEN (t [1], gel (gel (tau, 2), 1), p);
   pari_mpc_init_set_GEN (t [2], gel (gel (tau, 2), 2), p);

   /* Initialise intermediate variables. */
   for (i = 0; i < 10; i++)
      pari_mpc_init2 (theta2 [i], p);

   /* Compute the result. */
   eval_10theta2_naive (theta2, t);

   /* Create the output variable and copy the result. */
   res = cgetg (11, t_VEC);
   for (i = 0; i < 10; i++)
      gel (res, i+1) = mpc_get_GEN (theta2 [i]);

   return gerepileupto (av, res);
}

/****************************************************************************/

#endif

