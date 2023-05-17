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
/* Test functions for the mpfr interface                                    */
/*                                                                          */
/****************************************************************************/

#ifdef HAVE_LIBMPFR
int test_f_to_GEN (mpfr_srcptr f)
{
   int ret;
   mpfr_t f2;
   GEN x;

   mpfr_init2 (f2, mpfr_get_prec (f));

   x = mpfr_get_GEN (f);
   mpfr_set_GEN (f2, x, MPFR_RNDN);
   ret = mpfr_cmp (f, f2);

   mpfr_clear (f2);

   return (ret);
}

/****************************************************************************/

int test_f_to_GEN_stack (mpfr_srcptr f)
{
   int ret;
   mpfr_t f2;
   GEN x;

   pari_mpfr_init2 (f2, mpfr_get_prec (f));

   x = mpfr_get_GEN (f);
   mpfr_set_GEN (f2, x, MPFR_RNDN);
   ret = mpfr_cmp (f, f2);

   return (ret);
}

/****************************************************************************/

int test_GEN_to_f (GEN x)
{
   int ret;
   mpfr_t f;
   GEN x2;

   mpfr_init (f);

   mpfr_set_GEN (f, x, MPFR_RNDN);
   x2 = mpfr_get_GEN (f);

   mpfr_clear (f);

   return (gcmp (x, x2));
}

/****************************************************************************/

int test_GEN_to_f_stack (GEN x)
{
   int ret;
   mpfr_t f;
   GEN x2;

   pari_mpfr_init_set_GEN (f, x, 53);
   x2 = mpfr_get_GEN (f);

   return (gcmp (x, x2));
}
#endif

/****************************************************************************/

int main ()
{
#ifdef HAVE_LIBMPFR
   mpfr_t f, f2;
   gmp_randstate_t rand;
   GEN x;
   mpfr_prec_t prec;

   pari_init (1000000, 0);
   gmp_randinit_default (rand);

   mpfr_init (f);
   mpfr_set_si (f, 3, MPFR_RNDN);
   mpfr_div_2ui (f, f, 200, MPFR_RNDN);
   assert (test_f_to_GEN (f) == 0);
   assert (test_f_to_GEN_stack (f) == 0);
   mpfr_neg (f, f, MPFR_RNDN);
   assert (test_f_to_GEN (f) == 0);
   assert (test_f_to_GEN_stack (f) == 0);
   mpfr_set_ui (f, 0, MPFR_RNDN);
   assert (test_f_to_GEN (f) == 0);
   assert (test_f_to_GEN_stack (f) == 0);
   for (prec = 2; prec < 200; prec++) {
      mpfr_set_prec (f, prec);
      mpfr_urandomb (f, rand);
      assert (test_f_to_GEN (f) == 0);
      assert (test_f_to_GEN_stack (f) == 0);
   }

   mpfr_clear (f);

   x = gmul2n (stor (3, 6), -200);
   assert (test_GEN_to_f (x) == 0);
   x = gneg (x);
   assert (test_GEN_to_f (x) == 0);
   assert (test_GEN_to_f (stor (3, 6)) == 0);
   assert (test_GEN_to_f (real_0_bit (-230)) == 0);
   assert (test_GEN_to_f (real_0_bit (230)) == 0);
   assert (test_GEN_to_f (gen_0) == 0);

   x = gmul2n (stor (3, 6), -200);
   assert (test_GEN_to_f_stack (x) == 0);
   x = gneg (x);
   assert (test_GEN_to_f_stack (x) == 0);
   assert (test_GEN_to_f_stack (stor (3, 6)) == 0);
   assert (test_GEN_to_f_stack (real_0_bit (-230)) == 0);
   assert (test_GEN_to_f_stack (real_0_bit (230)) == 0);

   pari_close ();
   gmp_randclear (rand);
#else
   return (77);
#endif
}

/****************************************************************************/

