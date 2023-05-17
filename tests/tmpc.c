/*
Copyright © 2014 Andreas Enge <andreas.enge@inria.fr>

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
/*                                                                          */
/* Test functions for the mpc interface                                     */
/*                                                                          */
/****************************************************************************/

#ifdef HAVE_LIBMPC
int test_c_to_GEN (mpc_srcptr c)
{
   int ret;
   mpc_t c2;
   GEN x;

   mpc_init3 (c2, mpfr_get_prec (mpc_realref (c)),
                  mpfr_get_prec (mpc_imagref (c)));

   x = mpc_get_GEN (c);
   mpc_set_GEN (c2, x, MPC_RNDNN);
   ret = mpc_cmp (c, c2);

   mpc_clear (c2);

   return (ret);
}

/****************************************************************************/

int test_c_to_GEN_stack (mpc_srcptr c)
{
   int ret;
   mpc_t c2;
   GEN x;

   pari_mpc_init3 (c2, mpfr_get_prec (mpc_realref (c)),
                       mpfr_get_prec (mpc_imagref (c)));

   x = mpc_get_GEN (c);
   mpc_set_GEN (c2, x, MPC_RNDNN);
   ret = mpc_cmp (c, c2);

   return (ret);
}

/****************************************************************************/

int test_GEN_to_c (GEN x)
{
   int ret;
   mpc_t c;
   GEN x2;

   mpc_init2 (c, 100);

   mpc_set_GEN (c, x, MPC_RNDNN);
   x2 = mpc_get_GEN (c);

   mpc_clear (c);

   return (gequal (x, x2) != 1);
}

/****************************************************************************/

int test_GEN_to_c_stack (GEN x)
{
   int ret;
   mpc_t c;
   GEN x2;

   pari_mpc_init_set_GEN (c, x, 53);
   x2 = mpc_get_GEN (c);

   return (gequal (x, x2) != 1);
}
#endif

/****************************************************************************/

int main ()
{
#ifdef HAVE_LIBMPC
   mpc_t c, c2;
   gmp_randstate_t rand;
   GEN x;
   mpfr_prec_t prec;

   pari_init (1000000, 0);
   gmp_randinit_default (rand);

   mpc_init2 (c, 123);
   mpc_set_ui_ui (c, 3, 5, MPC_RNDNN);
   mpc_div_2ui (c, c, 100, MPC_RNDNN);
   assert (test_c_to_GEN (c) == 0);
   assert (test_c_to_GEN_stack (c) == 0);
   mpc_set_ui_ui (c, 0, 5, MPC_RNDNN);
   assert (test_c_to_GEN (c) == 0);
   assert (test_c_to_GEN_stack (c) == 0);
   mpc_set_ui_ui (c, 3, 0, MPC_RNDNN);
   assert (test_c_to_GEN (c) == 0);
   assert (test_c_to_GEN_stack (c) == 0);
   for (prec = 2; prec < 200; prec++) {
      mpc_set_prec (c, prec);
      mpc_urandom (c, rand);
      assert (test_c_to_GEN (c) == 0);
      assert (test_c_to_GEN_stack (c) == 0);
   }

   x = gmul2n (mkcomplex (stor (3, 6), stor (5, 6)), -200);
   assert (test_GEN_to_c (x) == 0);
   assert (test_GEN_to_c_stack (x) == 0);
   x = gmul2n (mkcomplex (gen_0, stor (5, 6)), -200);
   assert (test_GEN_to_c (x) == 0);
   assert (test_GEN_to_c_stack (x) == 0);
   x = gmul2n (mkcomplex (stor (3, 6), gen_0), -200);
   assert (test_GEN_to_c (x) == 0);
   assert (test_GEN_to_c_stack (x) == 0);
   
   mpc_clear (c);

   pari_close ();
   gmp_randclear (rand);
#else
   return (77);
#endif
}

/****************************************************************************/

