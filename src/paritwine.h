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

#include <paritwine-config.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>
#include <pari/pari.h>
#include <gmp.h>
#ifdef HAVE_LIBMPFR
#include <mpfr.h>
#endif
#ifdef HAVE_LIBMPC
#include <mpc.h>
#endif
#ifdef HAVE_LIBCMH
#include <cmh.h>
#endif
#ifdef HAVE_LIBFLINT
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#endif
#ifdef HAVE_LIBARB
#include <arb.h>
#include <acb.h>
#include <acb_modular.h>
#include <acb_elliptic.h>
#include <partitions.h>
#endif


/* Conversion functions. */

void mpz_set_GEN (mpz_ptr z, GEN x);
GEN mpz_get_GEN (mpz_srcptr z);

void mpq_set_GEN (mpq_ptr q, GEN x);
GEN mpq_get_GEN (mpq_srcptr q);

#ifdef HAVE_LIBMPFR
void pari_mpfr_init2 (mpfr_ptr f, mpfr_prec_t prec);
void pari_mpfr_init_set_GEN (mpfr_ptr f, GEN x, mpfr_prec_t default_prec);
int mpfr_set_GEN (mpfr_ptr f, GEN x, mpfr_rnd_t rnd);
GEN mpfr_get_GEN (mpfr_srcptr f);
#endif

#ifdef HAVE_LIBMPC
void pari_mpc_init2 (mpc_ptr c, mpfr_prec_t prec);
void pari_mpc_init3 (mpc_ptr c, mpfr_prec_t prec_re, mpfr_prec_t prec_im);
void pari_mpc_init_set_GEN (mpc_ptr c, GEN x, mpfr_prec_t default_prec);
int mpc_set_GEN (mpc_ptr c, GEN x, mpc_rnd_t rnd);
GEN mpc_get_GEN (mpc_srcptr c);
#endif

#ifdef HAVE_LIBFLINT
GEN fmpz_get_GEN(const fmpz_t x);
void fmpz_set_GEN(fmpz_t y, const GEN x);
GEN fmpq_get_GEN(const fmpq_t x);
void fmpq_set_GEN(fmpq_t y, const GEN x);
#endif

#ifdef HAVE_LIBARB
GEN arf_get_GEN(const arf_t x, slong prec);
void arf_set_GEN(arf_t y, const GEN x);
GEN mag_get_GEN(const mag_t x);
void mag_set_GEN(mag_t y, const GEN x);
GEN arb_get_GEN(const arb_t x, slong prec);
void arb_set_GEN(arb_t y, const GEN x, slong prec);
GEN acb_get_GEN(const acb_t x, slong prec);
void acb_set_GEN(acb_t y, const GEN x, slong prec);
#endif


/* Bindings to the different libraries. */

#ifdef HAVE_LIBMPFR
GEN pari_mpfr_add (GEN x, GEN y, long prec);
GEN pari_mpfr_sub (GEN x, GEN y, long prec);
GEN pari_mpfr_mul (GEN x, GEN y, long prec);
GEN pari_mpfr_sqr (GEN x, long prec);
GEN pari_mpfr_div (GEN x, GEN y, long prec);
GEN pari_mpfr_sqrt (GEN x, long prec);
GEN pari_mpfr_rec_sqrt (GEN x, long prec);
GEN pari_mpfr_cbrt (GEN x, long prec);
GEN pari_mpfr_log (GEN x, long prec);
GEN pari_mpfr_log2 (GEN x, long prec);
GEN pari_mpfr_log10 (GEN x, long prec);
GEN pari_mpfr_exp (GEN x, long prec);
GEN pari_mpfr_exp2 (GEN x, long prec);
GEN pari_mpfr_exp10 (GEN x, long prec);
GEN pari_mpfr_cos (GEN x, long prec);
GEN pari_mpfr_sin (GEN x, long prec);
GEN pari_mpfr_tan (GEN x, long prec);
GEN pari_mpfr_sec (GEN x, long prec);
GEN pari_mpfr_csc (GEN x, long prec);
GEN pari_mpfr_cot (GEN x, long prec);
GEN pari_mpfr_acos (GEN x, long prec);
GEN pari_mpfr_asin (GEN x, long prec);
GEN pari_mpfr_atan (GEN x, long prec);
GEN pari_mpfr_cosh (GEN x, long prec);
GEN pari_mpfr_sinh (GEN x, long prec);
GEN pari_mpfr_tanh (GEN x, long prec);
GEN pari_mpfr_sech (GEN x, long prec);
GEN pari_mpfr_csch (GEN x, long prec);
GEN pari_mpfr_coth (GEN x, long prec);
GEN pari_mpfr_acosh (GEN x, long prec);
GEN pari_mpfr_asinh (GEN x, long prec);
GEN pari_mpfr_atanh (GEN x, long prec);
GEN pari_mpfr_fac_ui (unsigned long int i, long prec);
GEN pari_mpfr_log1p (GEN x, long prec);
GEN pari_mpfr_expm1 (GEN x, long prec);
GEN pari_mpfr_eint (GEN x, long prec);
GEN pari_mpfr_li2 (GEN x, long prec);
GEN pari_mpfr_gamma (GEN x, long prec);
GEN pari_mpfr_lngamma (GEN x, long prec);
GEN pari_mpfr_digamma (GEN x, long prec);
GEN pari_mpfr_zeta (GEN x, long prec);
GEN pari_mpfr_erf (GEN x, long prec);
GEN pari_mpfr_erfc (GEN x, long prec);
GEN pari_mpfr_j0 (GEN x, long prec);
GEN pari_mpfr_j1 (GEN x, long prec);
GEN pari_mpfr_jn (long int i, GEN x, long prec);
GEN pari_mpfr_y0 (GEN x, long prec);
GEN pari_mpfr_y1 (GEN x, long prec);
GEN pari_mpfr_yn (long int i, GEN x, long prec);
GEN pari_mpfr_fma (GEN x, GEN y, GEN z, long prec);
GEN pari_mpfr_fms (GEN x, GEN y, GEN z, long prec);
GEN pari_mpfr_agm (GEN x, GEN y, long prec);
GEN pari_mpfr_hypot (GEN x, GEN y, long prec);
GEN pari_mpfr_ai (GEN x, long prec);
#endif

#ifdef HAVE_LIBMPC
GEN pari_mpc_add (GEN x, GEN y, long prec);
GEN pari_mpc_sub (GEN x, GEN y, long prec);
GEN pari_mpc_mul (GEN x, GEN y, long prec);
GEN pari_mpc_sqr (GEN x, long prec);
GEN pari_mpc_div (GEN x, GEN y, long prec);
GEN pari_mpc_fma (GEN x, GEN y, GEN z, long prec);
GEN pari_mpc_abs (GEN x, long prec);
GEN pari_mpc_norm (GEN x, long prec);
GEN pari_mpc_sqrt (GEN x, long prec);
GEN pari_mpc_pow (GEN x, GEN y, long prec);
GEN pari_mpc_exp (GEN x, long prec);
GEN pari_mpc_log (GEN x, long prec);
GEN pari_mpc_log10 (GEN x, long prec);
/*
GEN pari_mpc_rootofunity (unsigned long int n, unsigned long int k, long prec);
*/
GEN pari_mpc_sin (GEN x, long prec);
GEN pari_mpc_cos (GEN x, long prec);
GEN pari_mpc_tan (GEN x, long prec);
GEN pari_mpc_sinh (GEN x, long prec);
GEN pari_mpc_cosh (GEN x, long prec);
GEN pari_mpc_tanh (GEN x, long prec);
GEN pari_mpc_asin (GEN x, long prec);
GEN pari_mpc_acos (GEN x, long prec);
GEN pari_mpc_atan (GEN x, long prec);
GEN pari_mpc_asinh (GEN x, long prec);
GEN pari_mpc_acosh (GEN x, long prec);
GEN pari_mpc_atanh (GEN x, long prec);
#endif

#ifdef HAVE_LIBCMH
GEN pari_cmh_I2I4I6I10 (GEN tau, long prec);
GEN pari_cmh_4theta (GEN tau, long prec);
GEN pari_cmh_10theta2 (GEN tau, long prec);
#endif

#ifdef HAVE_LIBARB
GEN pari_acb_add (GEN x, GEN y, long prec);
GEN pari_acb_sub (GEN x, GEN y, long prec);
GEN pari_acb_mul (GEN x, GEN y, long prec);
GEN pari_acb_div (GEN x, GEN y, long prec);
GEN pari_acb_neg (GEN x, long prec);
GEN pari_acb_conj (GEN x, long prec);
GEN pari_acb_exp (GEN x, long prec);
GEN pari_acb_log (GEN x, long prec);
GEN pari_acb_atan (GEN x, long prec);
GEN pari_acb_sin (GEN x, long prec);
GEN pari_acb_cos (GEN x, long prec);
GEN pari_acb_gamma (GEN x, long prec);
GEN pari_acb_digamma (GEN x, long prec);
GEN pari_acb_zeta (GEN s, long prec);
GEN pari_acb_hurwitz_zeta (GEN s, GEN z, long prec);
GEN pari_acb_modular_eta (GEN tau, long prec);
GEN pari_acb_modular_j (GEN tau, long prec);
GEN pari_acb_modular_theta (GEN z, GEN tau, long prec);
GEN pari_acb_modular_delta (GEN tau, long prec);
GEN pari_acb_modular_eisenstein (GEN tau, long n, long prec);
GEN pari_acb_elliptic_p (GEN z, GEN tau, long prec);
GEN pari_acb_elliptic_inv_p (GEN z, GEN tau, long prec);
GEN pari_acb_elliptic_zeta (GEN z, GEN tau, long prec);
GEN pari_acb_elliptic_sigma (GEN z, GEN tau, long prec);
GEN pari_fmpz_numbpart (GEN x);
GEN pari_arb_numbpart (GEN x, long prec);
#endif

