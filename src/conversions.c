/*
Copyright © 2014, 2017, 2018, 2019 Andreas Enge <andreas.enge@inria.fr>
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
/* Functions converting between pari and mpz                                */
/*                                                                          */
/****************************************************************************/

void mpz_set_GEN (mpz_ptr z, GEN x)
   /* Sets z to x, which needs to be of type t_INT. */

{
   const long lx = lgefint (x) - 2;
   const long sign = signe (x);
   int i;

   assert (sizeof (long) == sizeof (mp_limb_t));

   if (typ (x) != t_INT)
      pari_err_TYPE ("mpz_set_GEN", x);
   
   if (sign == 0)
      mpz_set_ui (z, 0);
   else {
      mpz_realloc2 (z, lx * BITS_IN_LONG);
      z->_mp_size = sign * lx;
      for (i = 0; i < lx; i++)
         (z->_mp_d) [i] = *int_W (x, i);
   }
}

/****************************************************************************/

GEN mpz_get_GEN (mpz_srcptr z)
   /* Returns the GEN of type t_INT corresponding to z. */

{
   const long lz = z->_mp_size;
   const long lx = labs (lz);
   const long lx2 = lx + 2;
   int i;
   GEN x = cgeti (lx2);

   assert (sizeof (long) == sizeof (mp_limb_t));

   x [1] = evalsigne ((lz > 0 ? 1 : (lz < 0 ? -1 : 0))) | evallgefint (lx2);
   for (i = 0; i < lx; i++)
      *int_W (x, i) = (z->_mp_d) [i];

   return x;
}


/****************************************************************************/
/*                                                                          */
/* Functions converting between pari and mpq                                */
/*                                                                          */
/****************************************************************************/

void mpq_set_GEN (mpq_ptr q, GEN x)
   /* Sets q to x, which needs to be of type t_FRAC or t_INT. */

{
   if (typ (x) == t_FRAC) {
      mpz_set_GEN (&(q->_mp_num), gel (x, 1));
      mpz_set_GEN (&(q->_mp_den), gel (x, 2));
   }
   else if (typ (x) == t_INT) {
      mpz_set_GEN (&(q->_mp_num), x);
      mpz_set_GEN (&(q->_mp_den), gen_1);
   }
   else
      pari_err_TYPE ("mpq_set_GEN", x);
}

/****************************************************************************/

GEN mpq_get_GEN (mpq_srcptr q)
   /* Returns the GEN of type t_FRAC corresponding to q. */

{
   GEN x = cgetg (3, t_FRAC);

   gel (x, 1) = mpz_get_GEN (mpq_numref (q));
   gel (x, 2) = mpz_get_GEN (mpq_denref (q));

   return x;
}


/****************************************************************************/
/*                                                                          */
/* Functions converting between pari and mpfr                               */
/*                                                                          */
/****************************************************************************/

#ifdef HAVE_LIBMPFR
static void xmpn_mirrorcopy (mp_limb_t *z, mp_limb_t *x, long n)
   /* Copied from kernel/gmp/mp.c; used to revert the mantissa between
      mpfr and pari representations. */
{
  long i;
  for (i = 0; i < n; i++)
    z [i]= x [n-1-i];
}

/****************************************************************************/

void pari_mpfr_init2 (mpfr_ptr f, mpfr_prec_t prec)
   /* Works in the same way as mpfr_init2, except that the space for the
      mantissa is allocated on the pari stack. */
{
   long l;
   mp_limb_t * mant;

   assert (sizeof (long) == sizeof (mp_limb_t));
   
   l = nbits2nlong (mpfr_custom_get_size (prec) * 8);
   mant = new_chunk (l);
   mpfr_custom_init_set (f, MPFR_NAN_KIND, 0, prec, mant);
}

/****************************************************************************/

static void pari_mpfr_init_set_REAL (mpfr_ptr f, GEN x)
   /* Creates space for f on the pari stack and sets it to x, which is
      assumed to be a t_REAL.
      Unlike with the mpfr_init_set functions, the used precision is not the
      mpfr default precision, but that of x (also if x is 0). */
{
   assert (sizeof (long) == sizeof (mp_limb_t));

   if (signe (x) == 0) {
      long e = expo (x);
      mp_prec_t prec = (e > -2 ? 2 : -e);
      pari_mpfr_init2 (f, prec);
      mpfr_set_zero (f, 1);
   }
   else {
      long l = realprec (x) - 2;
      mp_limb_t * rc = new_chunk (l);

      xmpn_mirrorcopy (rc, x+2, l);
      mpfr_custom_init_set (f, signe (x) * MPFR_REGULAR_KIND,
                           expo (x) + 1, bit_prec (x), rc);
   }
}

/****************************************************************************/

static int mpfr_set_REAL (mpfr_ptr f, GEN x, mpfr_rnd_t rnd)
   /* Sets f to x, which is assumed to be of type t_REAL. */

{
   pari_sp av = avma;
   int inex;
   mpfr_t tmp;

   /* Copy x without loss. TODO: Share mantissa of x once it need not
      be reverted any more. */
   pari_mpfr_init_set_REAL (tmp, x);
   inex = mpfr_set (f, tmp, rnd);

   avma = av;

   return (inex);
}

/****************************************************************************/

static int mpfr_set_INT (mpfr_ptr f, GEN x, mpfr_rnd_t rnd)
   /* Sets f to x, which is assumed to be of type t_INT. */

{
   int inex;
   mpz_t tmp;

   mpz_init (tmp);
   mpz_set_GEN (tmp, x);
   inex = mpfr_set_z (f, tmp, rnd);
   mpz_clear (tmp);

   return (inex);
}

/****************************************************************************/

static void pari_mpfr_init_set_INT (mpfr_ptr f, GEN x, mpfr_prec_t prec)
   /* Creates space for f of mpfr precision prec on the pari stack and sets
      it to x, which is assumed to be a t_INT. */
{
   pari_mpfr_init2 (f, prec);
   mpfr_set_INT (f, x, MPFR_RNDN);
}

/****************************************************************************/

static int mpfr_set_FRAC (mpfr_ptr f, GEN x, mpfr_rnd_t rnd)
   /* Sets f to x, which is assumed to be of type t_FRAC. */

{
   int inex;
   mpq_t tmp;

   mpq_init (tmp);
   mpq_set_GEN (tmp, x);
   inex = mpfr_set_q (f, tmp, rnd);
   mpq_clear (tmp);

   return (inex);
}

/****************************************************************************/

static void pari_mpfr_init_set_FRAC (mpfr_ptr f, GEN x, mpfr_prec_t prec)
   /* Creates space for f of mpfr precision prec on the pari stack and sets
      it to x, which is assumed to be a t_FRAC. */
{
   pari_mpfr_init2 (f, prec);
   mpfr_set_FRAC (f, x, MPFR_RNDN);
}

/****************************************************************************/

int mpfr_set_GEN (mpfr_ptr f, GEN x, mpfr_rnd_t rnd)
   /* Sets f to x. */

{
   switch (typ (x)) {
      case t_REAL:
         return mpfr_set_REAL (f, x, rnd);
      case t_INT:
         return mpfr_set_INT (f, x, rnd);
      case t_FRAC:
         return mpfr_set_FRAC (f, x, rnd);
      default:
         pari_err_TYPE ("mpfr_set_GEN", x);
   }
}

/****************************************************************************/

void pari_mpfr_init_set_GEN (mpfr_ptr f, GEN x, mpfr_prec_t default_prec)
   /* Creates space for f on the pari stack and sets it to x.
      If x is of type t_REAL, its own precision is used, otherwise,
      default_prec is used for f. */

{
   switch (typ (x)) {
      case t_REAL:
         pari_mpfr_init_set_REAL (f, x);
         return;
      case t_INT:
         pari_mpfr_init_set_INT (f, x, default_prec);
         return;
      case t_FRAC:
         pari_mpfr_init_set_FRAC (f, x, default_prec);
         return;
      default:
         pari_err_TYPE ("pari_mpfr_init_set_GEN", x);
   }
}

/****************************************************************************/

GEN mpfr_get_GEN (mpfr_srcptr f)
   /* returns f as a t_REAL with the minimal precision required to represent
      it */
{
   if (!mpfr_number_p (f))
      pari_err_OVERFLOW ("mpfr_get_GEN");
   else if (mpfr_zero_p (f))
      return real_0_bit (-mpfr_get_prec (f));
   else {
      long l = nbits2nlong (mpfr_get_prec (f));
      GEN r = cgetr (l + 2);
      xmpn_mirrorcopy (r+2, mpfr_custom_get_significand (f), l);
      setsigne (r, (mpfr_signbit (f) ? -1 : 1));
      setexpo (r, mpfr_get_exp (f) - 1);
      return r;
   }
}
#endif


/****************************************************************************/
/*                                                                          */
/* Functions converting between pari and mpc                                */
/*                                                                          */
/****************************************************************************/

#ifdef HAVE_LIBMPC
void pari_mpc_init2 (mpc_ptr c, mpfr_prec_t prec)
   /* Works in the same way as mpc_init2, except that the space for the
      mantissae is allocated on the pari stack. */
{
   pari_mpfr_init2 (mpc_realref (c), prec);
   pari_mpfr_init2 (mpc_imagref (c), prec);
}

/****************************************************************************/

void pari_mpc_init3 (mpc_ptr c, mpfr_prec_t prec_re, mpfr_prec_t prec_im)
   /* Works in the same way as mpc_init3, except that the space for the
      mantissae is allocated on the pari stack. */
{
   pari_mpfr_init2 (mpc_realref (c), prec_re);
   pari_mpfr_init2 (mpc_imagref (c), prec_im);
}

/****************************************************************************/

int mpc_set_GEN (mpc_ptr c, GEN x, mpc_rnd_t rnd)
   /* Sets c to x. */

{
   int inex_re, inex_im;
   
   switch (typ (x)) {
      case t_COMPLEX:
         inex_re = mpfr_set_GEN (mpc_realref (c), greal (x),
                                 MPC_RND_RE (rnd));
         inex_im = mpfr_set_GEN (mpc_imagref (c), gimag (x),
                                 MPC_RND_IM (rnd));
         break;
      case t_REAL:
      case t_INT:
      case t_FRAC:
         inex_re = mpfr_set_GEN (mpc_realref (c), x, MPC_RND_RE (rnd));
         inex_im = mpfr_set_GEN (mpc_imagref (c), gen_0, MPC_RND_IM (rnd));
         break;
      default:
         pari_err_TYPE ("mpc_set_GEN", x);
   }

   return MPC_INEX (inex_re, inex_im);
}

/****************************************************************************/

void pari_mpc_init_set_GEN (mpc_ptr c, GEN x, mpfr_prec_t default_prec)
   /* Creates space for c on the pari stack and sets it to x.
      Components of type t_REAL get their own precision, while others
      are created with default_prec. */

{
   switch (typ (x)) {
      case t_COMPLEX:
         pari_mpfr_init_set_GEN (mpc_realref (c), greal (x), default_prec);
         pari_mpfr_init_set_GEN (mpc_imagref (c), gimag (x), default_prec);
         return;
      case t_REAL:
      case t_INT:
      case t_FRAC:
         pari_mpfr_init_set_GEN (mpc_realref (c), greal (x), default_prec);
         pari_mpfr_init_set_GEN (mpc_imagref (c), gen_0, default_prec);
         return;
      default:
         pari_err_TYPE ("mpc_init_set_GEN", x);
   }
}

/****************************************************************************/

GEN mpc_get_GEN (mpc_srcptr c)
   /* returns c as a t_COMPLEX with the minimal precision required to
      represent it */
{
   GEN x = cgetg (3, t_COMPLEX);
   
   gel (x, 1) = mpfr_get_GEN (mpc_realref (c));
   gel (x, 2) = mpfr_get_GEN (mpc_imagref (c));

   return x;
}
#endif


/****************************************************************************/
/*                                                                          */
/* Functions converting between pari and flint                              */
/*                                                                          */
/****************************************************************************/

#ifdef HAVE_LIBFLINT
GEN
fmpz_get_GEN(const fmpz_t x)
{
    if (!COEFF_IS_MPZ(*x))
    {
        return stoi(*x);
    }
    else
    {
        GEN y;
        slong ssize, size;
        __mpz_struct * z;

        z = COEFF_TO_PTR(*x);
        ssize = z->_mp_size;
        size = FLINT_ABS(ssize);

        y = cgeti(size + 2);
        y[1] = evalsigne(ssize < 0 ? -1 : 1) + evallgefint(size + 2);
        flint_mpn_copyi(&y[2], z->_mp_d, size);

        return y;
    }
}

void
fmpz_set_GEN(fmpz_t y, const GEN x)
{
    if (typ(x) == t_INT)
    {
        slong size;

        size = lgefint(x) - 2;

        if (size == 0)
        {
            fmpz_zero(y);
        }
        else if (size == 1)
        {
            if (signe(x) < 0)
                fmpz_neg_ui(y, x[2]);
            else
                fmpz_set_ui(y, x[2]);
        }
        else
        {
            __mpz_struct * z = _fmpz_promote(y);

            if (z->_mp_alloc < size)
                mpz_realloc2(z, size * FLINT_BITS);

            flint_mpn_copyi(z->_mp_d, &x[2], size);
            z->_mp_size = (signe(x) < 0) ? -size : size;
        }
    }
    else
    {
         pari_err_TYPE ("fmpz_set_GEN", x);
    }
}

GEN
fmpq_get_GEN(const fmpq_t x)
{
    if (fmpz_is_one(fmpq_denref(x)))
    {
        return fmpz_get_GEN(fmpq_numref(x));
    }
    else
    {
        GEN y = cgetg(3, t_FRAC);

        gel(y, 1) = fmpz_get_GEN(fmpq_numref(x));
        gel(y, 2) = fmpz_get_GEN(fmpq_denref(x));
        return y;
    }
}

void
fmpq_set_GEN(fmpq_t y, const GEN x)
{
    if (typ(x) == t_INT)
    {
        fmpz_set_GEN(fmpq_numref(y), x);
        fmpz_one(fmpq_denref(y));
    }
    else if (typ(x) == t_FRAC)
    {
        fmpz_set_GEN(fmpq_numref(y), gel(x, 1));
        fmpz_set_GEN(fmpq_denref(y), gel(x, 2));
    }
    else
    {
         pari_err_TYPE ("fmpq_set_GEN", x);
    }
}
#endif


/****************************************************************************/
/*                                                                          */
/* Functions converting between pari and arb                                */
/*                                                                          */
/****************************************************************************/

#ifdef HAVE_LIBARB
GEN
arf_get_GEN(const arf_t x, slong prec)
{
    if (!arf_is_finite(x))
    {
        pari_err_TYPE ("arf_get_GEN", gen_0);  /* todo: throw the right error */
        return gen_0;
        /*
        if (arf_is_pos_inf(x))
            return mkoo();
        else if (arf_is_neg_inf(x))
            return mkmoo();
        else
            ...
        */
    }
    else if (arf_is_zero(x))
    {
        return gen_0;
    }
    else
    {
        GEN y;
        slong i;
        mp_srcptr xp;
        mp_size_t xn, yn;

        ARF_GET_MPN_READONLY(xp, xn, x);

        yn = (prec + FLINT_BITS - 1) / FLINT_BITS;

        y = cgetr(yn + 2);
        y[1] = evalsigne(ARF_SGNBIT(x) ? -1 : 1) | evalexpo(ARF_EXP(x) - 1);
        for (i = 0; i < FLINT_MIN(xn, yn); i++)
            y[2 + i] = xp[xn - 1 - i];
        for (i = xn; i < yn; i++)
            y[2 + i] = 0;

        return y;
    }
}

void
arf_set_GEN(arf_t y, const GEN x)
{
    if (typ(x) == t_INT)
    {
        slong size = lgefint(x) - 2;

        if (size == 0)
            arf_zero(y);
        else
            arf_set_mpn(y, x + 2, size, signe(x) < 0);
    }
    else if (typ(x) == t_REAL)
    {
        slong expo, sgn, xn, i;
        mp_ptr yp;

        xn = lg(x) - 2;
        sgn = signe(x);
        expo = expo(x);

        while (xn > 0 && x[2 + xn - 1] == 0)
            xn--;

        if (xn == 0)
        {
            arf_zero(y);
        }
        else
        {
            ARF_GET_MPN_WRITE(yp, xn, y);

            for (i = 0; i < xn; i++)
                yp[xn - 1 - i] = x[2 + i];

            if (sgn < 0)
                ARF_NEG(y);

            fmpz_set_si(ARF_EXPREF(y), expo + 1);
        }
    }
    else
    {
         pari_err_TYPE ("arf_set_GEN", x);
    }
}

GEN
mag_get_GEN(const mag_t x)
{
    if (mag_is_zero(x))
    {
        return gen_0;
    }
    else if (mag_is_inf(x))
    {
        pari_err_TYPE ("mag_get_GEN", gen_0);  /* todo: throw the right error */
        return gen_0;
    }
    else
    {
        GEN y;
        y = cgetr(3);
        y[1] = evalsigne(1) | evalexpo(MAG_EXP(x) - 1);
        y[2] = MAG_MAN(x) << (FLINT_BITS - MAG_BITS);
        return y;
    }
}

void mag_set_GEN(mag_t y, const GEN x)
{
    arf_t t;
    arf_init(t);
    arf_set_GEN(t, x);
    arf_get_mag(y, t);
    arf_clear(t);
}

void
arb_set_GEN(arb_t y, const GEN x, slong prec)
{
    if (typ(x) == t_FRAC)
    {
        fmpq_t t;
        fmpq_init(t);
        fmpq_set_GEN(t, x);
        arb_set_fmpq(y, t, prec);
        fmpq_clear(t);
    }
    else if (typ(x) == t_VEC)
    {
        if (lg(x) != 3)
        {
            pari_err_TYPE ("arb_set_GEN", gen_0);  /* todo: throw the right error */
        }
        else
        {
            arf_set_GEN(arb_midref(y), gel(x, 1));
            mag_set_GEN(arb_radref(y), gel(x, 2));
        }
    }
    else
    {
        arf_set_GEN(arb_midref(y), x);
        mag_zero(arb_radref(y));
        arb_set_round(y, y, prec);
    }
}

GEN
arb_get_GEN(const arb_t x, slong prec)
{
    GEN y = cgetg(3, t_VEC);
    gel(y, 1) = arf_get_GEN(arb_midref(x), prec);
    gel(y, 2) = mag_get_GEN(arb_radref(x));
    return y;
}

void
acb_set_GEN(acb_t y, const GEN x, slong prec)
{
    if (typ(x) == t_COMPLEX)
    {
        arb_set_GEN(acb_realref(y), gel(x, 1), prec);
        arb_set_GEN(acb_imagref(y), gel(x, 2), prec);
    }
    else if (typ(x) == t_VEC)
    {
        if (lg(x) != 3)
        {
            pari_err_TYPE ("acb_set_GEN", gen_0);  /* todo: throw the right error */
        }
        else
        {
            /* Create two variables of type t_VEC for the real and imaginary
               parts of x, and call arb_set_GEN on them. */
            pari_sp av = avma;
            GEN re, im;
            re = cgetg (3, t_VEC);
            gel (re, 1) = greal (gel (x, 1));
            gel (re, 2) = greal (gel (x, 2));
            im = cgetg (3, t_VEC);
            gel (im, 1) = gimag (gel (x, 1));
            gel (im, 2) = gimag (gel (x, 2));

            arb_set_GEN (acb_realref (y), re, prec);
            arb_set_GEN (acb_imagref (y), im, prec);

            avma = av;
        }
    }
    else
    {
        arb_zero(acb_imagref(y));
        arb_set_GEN(acb_realref(y), x, prec);
    }
}

GEN
acb_get_GEN(const acb_t x, slong prec)
{
    GEN z, m, r;

    z = cgetg(3, t_VEC);
    gel(z, 1) = m = cgetg(3, t_COMPLEX);
    gel(z, 2) = r = cgetg(3, t_COMPLEX);

    gel(m, 1) = arf_get_GEN(arb_midref(acb_realref(x)), prec);
    gel(m, 2) = arf_get_GEN(arb_midref(acb_imagref(x)), prec);
    gel(r, 1) = mag_get_GEN(arb_radref(acb_realref(x)));
    gel(r, 2) = mag_get_GEN(arb_radref(acb_imagref(x)));

    return z;
}
#endif

