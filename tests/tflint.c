/*
Copyright © 2014 Andreas Enge <andreas.enge@inria.fr>
Copyright © 2017 Fredrik Johansson <fredrik.johansson@gmail.com>

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
/* Test functions for the flint interface                                   */
/*                                                                          */
/****************************************************************************/

int main ()
{
#ifdef HAVE_LIBFLINT
   fmpz_t fz, fz2;
   fmpq_t fq, fq2;

   pari_init (1000000, 0);

   fmpz_init(fz);
   fmpz_init(fz2);
   fmpq_init(fq);
   fmpq_init(fq2);

   fmpz_set_si(fz, 23);
   fmpz_set_GEN(fz2, fmpz_get_GEN(fz));
   assert (fmpz_equal(fz, fz2) != 0);
   fmpz_pow_ui(fz, fz, 100);
   fmpz_neg(fz, fz);
   fmpz_set_GEN(fz2, fmpz_get_GEN(fz));
   assert (fmpz_equal(fz, fz2) != 0);

   fmpq_set_si(fq, 23, 17);
   fmpq_set_GEN(fq2, fmpq_get_GEN(fq));
   assert (fmpq_equal(fq, fq2) != 0);
   fmpq_set_si(fq, -23, 1);
   fmpq_set_GEN(fq2, fmpq_get_GEN(fq));
   assert (fmpq_equal(fq, fq2) != 0);

   fmpz_clear(fz);
   fmpz_clear(fz2);
   fmpq_clear(fq);
   fmpq_clear(fq2);

   pari_close ();
#else
   return (77);
#endif
}

/****************************************************************************/

