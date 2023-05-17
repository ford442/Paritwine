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
/* Test functions for the arb interface                                     */
/*                                                                          */
/****************************************************************************/

int main ()
{
#ifdef HAVE_LIBARB
   arb_t ar, ar2;
   acb_t ac, ac2;

   pari_init (1000000, 0);

   arb_init(ar);
   arb_init(ar2);
   acb_init(ac);
   acb_init(ac2);

   arb_set_d(ar, 3.125);
   arb_set_GEN(ar2, arb_get_GEN(ar, 64), 64);
   assert (arb_equal(ar, ar2) != 0);
   acb_set_d_d(ac, 3.125, -0.6125);
   acb_set_GEN(ac2, acb_get_GEN(ac, 64), 64);
   assert (acb_equal(ac, ac2) != 0);

   arb_set_str(ar, "0.1e-100", 1024);
   mag_zero(arb_radref(ar));
   arb_set_GEN(ar2, arb_get_GEN(ar, 1024), 1024);
   assert (arb_equal(ar, ar2) != 0);

   arb_set_str(ar, "0.1e-100", 1024);
   arb_set_GEN(ar2, arb_get_GEN(ar, 1024), 1024);
   assert (arb_contains(ar2, ar) != 0);

   arb_set_str(acb_realref(ac), "0.1e-100", 128);
   arb_set_str(acb_imagref(ac), "-0.1e+100", 192);
   acb_set_GEN(ac2, acb_get_GEN(ac, 192), 192);
   acb_printd(ac, 30); printf("\n");
   acb_printd(ac2, 30); printf("\n");
   assert (acb_contains(ac2, ac) != 0);

   arb_clear(ar);
   arb_clear(ar2);
   acb_clear(ac);
   acb_clear(ac2);

   pari_close ();
#else
   return (77);
#endif
}

/****************************************************************************/

