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
/* Test functions for the gmp interface                                     */
/*                                                                          */
/****************************************************************************/

int test_z_to_GEN (mpz_srcptr z)
{
   int ret;
   mpz_t z2;
   GEN x;

   mpz_init (z2);

   x = mpz_get_GEN (z);
   mpz_set_GEN (z2, x);
   ret = mpz_cmp (z, z2);

   mpz_clear (z2);

   return (ret);
}

/****************************************************************************/

int test_GEN_to_z (GEN x)
{
   int ret;
   mpz_t z;
   GEN x2;

   mpz_init (z);

   mpz_set_GEN (z, x);
   x2 = mpz_get_GEN (z);

   mpz_clear (z);

   return (cmpii (x, x2));
}

/****************************************************************************/

int test_q_to_GEN (mpq_srcptr q)
{
   int ret;
   mpq_t q2;
   GEN x;

   mpq_init (q2);

   x = mpq_get_GEN (q);
   mpq_set_GEN (q2, x);
   ret = mpq_cmp (q, q2);

   mpq_clear (q2);

   return (ret);
}

/****************************************************************************/

int test_GEN_to_q (GEN x)
{
   int ret;
   mpq_t q;
   GEN x2;

   mpq_init (q);

   mpq_set_GEN (q, x);
   x2 = mpq_get_GEN (q);

   mpq_clear (q);

   return (gcmp (x, x2));
}

/****************************************************************************/

int main ()
{
   mpz_t z, z2;
   mpq_t q, q2;
   gmp_randstate_t rand;
   GEN x;

   pari_init (1000000, 0);
   gmp_randinit_default (rand);

   mpz_init (z);
   mpz_set_si (z, 3);
   mpz_mul_2exp (z, z, 200);
   assert (test_z_to_GEN (z) == 0);
   mpz_neg (z, z);
   assert (test_z_to_GEN (z) == 0);
   mpz_set_si (z, 0);
   assert (test_z_to_GEN (z) == 0);
   mpz_clear (z);

   x = shifti (stoi (3), 200);
   assert (test_GEN_to_z (x) == 0);
   x = negi (x);
   assert (test_GEN_to_z (x) == 0);
   assert (test_GEN_to_z (gen_0) == 0);

   mpq_init (q);
   mpq_set_si (q, 3, 1);
   mpq_div_2exp (q, q, 200);
   assert (test_q_to_GEN (q) == 0);
   mpq_neg (q, q);
   assert (test_q_to_GEN (q) == 0);
   mpq_set_ui (q, 3, 1);
   assert (test_q_to_GEN (q) == 0);
   mpq_set_ui (q, 0, 1);
   assert (test_q_to_GEN (q) == 0);
   mpq_clear (q);

   x = gmul2n (stoi (3), -200);
   assert (test_GEN_to_q (x) == 0);
   x = gneg (x);
   assert (test_GEN_to_q (x) == 0);
   assert (test_GEN_to_q (stoi (3)) == 0);
   assert (test_GEN_to_q (gen_0) == 0);

   pari_close ();
   gmp_randclear (rand);
}

/****************************************************************************/

