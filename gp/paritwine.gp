/*
Copyright © 2014, 2017 Andreas Enge <andreas.enge@inria.fr>

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

LIBPARITWINESO="/home/enge/local/paritwine/lib/libparitwine.so";

install ("pari_mpfr_add", "GGb", "mpfr_add", LIBPARITWINESO);
install ("pari_mpfr_sub", "GGb", "mpfr_sub", LIBPARITWINESO);
install ("pari_mpfr_mul", "GGb", "mpfr_mul", LIBPARITWINESO);
install ("pari_mpfr_sqr", "Gb", "mpfr_sqr", LIBPARITWINESO);
install ("pari_mpfr_div", "GGb", "mpfr_div", LIBPARITWINESO);
install ("pari_mpfr_sqrt", "Gb", "mpfr_sqrt", LIBPARITWINESO);
install ("pari_mpfr_rec_sqrt", "Gb", "mpfr_rec_sqrt", LIBPARITWINESO);
install ("pari_mpfr_cbrt", "Gb", "mpfr_cbrt", LIBPARITWINESO);
install ("pari_mpfr_pow", "GGb", "mpfr_pow", LIBPARITWINESO);
install ("pari_mpfr_log", "Gb", "mpfr_log", LIBPARITWINESO);
install ("pari_mpfr_log2", "Gb", "mpfr_log2", LIBPARITWINESO);
install ("pari_mpfr_log10", "Gb", "mpfr_log10", LIBPARITWINESO);
install ("pari_mpfr_exp", "Gb", "mpfr_exp", LIBPARITWINESO);
install ("pari_mpfr_exp2", "Gb", "mpfr_exp2", LIBPARITWINESO);
install ("pari_mpfr_exp10", "Gb", "mpfr_exp10", LIBPARITWINESO);
install ("pari_mpfr_sin", "Gb", "mpfr_sin", LIBPARITWINESO);
install ("pari_mpfr_cos", "Gb", "mpfr_cos", LIBPARITWINESO);
install ("pari_mpfr_tan", "Gb", "mpfr_tan", LIBPARITWINESO);
install ("pari_mpfr_sec", "Gb", "mpfr_sec", LIBPARITWINESO);
install ("pari_mpfr_csc", "Gb", "mpfr_csc", LIBPARITWINESO);
install ("pari_mpfr_cot", "Gb", "mpfr_cot", LIBPARITWINESO);
install ("pari_mpfr_acos", "Gb", "mpfr_acos", LIBPARITWINESO);
install ("pari_mpfr_asin", "Gb", "mpfr_asin", LIBPARITWINESO);
install ("pari_mpfr_atan", "Gb", "mpfr_atan", LIBPARITWINESO);
install ("pari_mpfr_cosh", "Gb", "mpfr_cosh", LIBPARITWINESO);
install ("pari_mpfr_sinh", "Gb", "mpfr_sinh", LIBPARITWINESO);
install ("pari_mpfr_tanh", "Gb", "mpfr_tanh", LIBPARITWINESO);
install ("pari_mpfr_sech", "Gb", "mpfr_sech", LIBPARITWINESO);
install ("pari_mpfr_csch", "Gb", "mpfr_csch", LIBPARITWINESO);
install ("pari_mpfr_coth", "Gb", "mpfr_coth", LIBPARITWINESO);
install ("pari_mpfr_acosh", "Gb", "mpfr_acosh", LIBPARITWINESO);
install ("pari_mpfr_asinh", "Gb", "mpfr_asinh", LIBPARITWINESO);
install ("pari_mpfr_atanh", "Gb", "mpfr_atanh", LIBPARITWINESO);
install ("pari_mpfr_fac_ui", "Lb", "mpfr_fac_ui", LIBPARITWINESO);
install ("pari_mpfr_log1p", "Gb", "mpfr_log1p", LIBPARITWINESO);
install ("pari_mpfr_expm1", "Gb", "mpfr_expm1", LIBPARITWINESO);
install ("pari_mpfr_eint", "Gb", "mpfr_eint", LIBPARITWINESO);
install ("pari_mpfr_li2", "Gb", "mpfr_li2", LIBPARITWINESO);
install ("pari_mpfr_gamma", "Gb", "mpfr_gamma", LIBPARITWINESO);
install ("pari_mpfr_lngamma", "Gb", "mpfr_lngamma", LIBPARITWINESO);
install ("pari_mpfr_digamma", "Gb", "mpfr_digamma", LIBPARITWINESO);
install ("pari_mpfr_zeta", "Gb", "mpfr_zeta", LIBPARITWINESO);
install ("pari_mpfr_erf", "Gb", "mpfr_erf", LIBPARITWINESO);
install ("pari_mpfr_erfc", "Gb", "mpfr_erfc", LIBPARITWINESO);
install ("pari_mpfr_j0", "Gb", "mpfr_j0", LIBPARITWINESO);
install ("pari_mpfr_j1", "Gb", "mpfr_j1", LIBPARITWINESO);
install ("pari_mpfr_jn", "LGb", "mpfr_jn", LIBPARITWINESO);
install ("pari_mpfr_y0", "Gb", "mpfr_y0", LIBPARITWINESO);
install ("pari_mpfr_y1", "Gb", "mpfr_y1", LIBPARITWINESO);
install ("pari_mpfr_yn", "LGb", "mpfr_yn", LIBPARITWINESO);
install ("pari_mpfr_fma", "GGGb", "mpfr_fma", LIBPARITWINESO);
install ("pari_mpfr_fms", "GGGb", "mpfr_fms", LIBPARITWINESO);
install ("pari_mpfr_agm", "GGb", "mpfr_agm", LIBPARITWINESO);
install ("pari_mpfr_hypot", "GGb", "mpfr_hypot", LIBPARITWINESO);
install ("pari_mpfr_ai", "Gb", "mpfr_ai", LIBPARITWINESO);

install ("pari_mpc_add", "GGb", "mpc_add", LIBPARITWINESO);
install ("pari_mpc_sub", "GGb", "mpc_sub", LIBPARITWINESO);
install ("pari_mpc_mul", "GGb", "mpc_mul", LIBPARITWINESO);
install ("pari_mpc_sqr", "Gb", "mpc_sqr", LIBPARITWINESO);
install ("pari_mpc_div", "GGb", "mpc_div", LIBPARITWINESO);
install ("pari_mpc_fma", "GGGb", "mpc_fma", LIBPARITWINESO);
install ("pari_mpc_abs", "Gb", "mpc_abs", LIBPARITWINESO);
install ("pari_mpc_norm", "Gb", "mpc_norm", LIBPARITWINESO);
install ("pari_mpc_pow", "GGb", "mpc_pow", LIBPARITWINESO);
install ("pari_mpc_sqrt", "Gb", "mpc_sqrt", LIBPARITWINESO);
install ("pari_mpc_exp", "Gb", "mpc_exp", LIBPARITWINESO);
install ("pari_mpc_log", "Gb", "mpc_log", LIBPARITWINESO);
install ("pari_mpc_log10", "Gb", "mpc_log10", LIBPARITWINESO);
install ("pari_mpc_sin", "Gb", "mpc_sin", LIBPARITWINESO);
install ("pari_mpc_cos", "Gb", "mpc_cos", LIBPARITWINESO);
install ("pari_mpc_tan", "Gb", "mpc_tan", LIBPARITWINESO);
install ("pari_mpc_asin", "Gb", "mpc_asin", LIBPARITWINESO);
install ("pari_mpc_acos", "Gb", "mpc_acos", LIBPARITWINESO);
install ("pari_mpc_atan", "Gb", "mpc_atan", LIBPARITWINESO);
install ("pari_mpc_sinh", "Gb", "mpc_sinh", LIBPARITWINESO);
install ("pari_mpc_cosh", "Gb", "mpc_cosh", LIBPARITWINESO);
install ("pari_mpc_tanh", "Gb", "mpc_tanh", LIBPARITWINESO);
install ("pari_mpc_asinh", "Gb", "mpc_asinh", LIBPARITWINESO);
install ("pari_mpc_acosh", "Gb", "mpc_acosh", LIBPARITWINESO);
install ("pari_mpc_atanh", "Gb", "mpc_atanh", LIBPARITWINESO);

/* placeholder flint */

install ("pari_acb_add", "GGb", "acb_add", LIBPARITWINESO);
install ("pari_acb_sub", "GGb", "acb_sub", LIBPARITWINESO);
install ("pari_acb_mul", "GGb", "acb_mul", LIBPARITWINESO);
install ("pari_acb_div", "GGb", "acb_div", LIBPARITWINESO);
install ("pari_acb_neg", "Gb", "acb_neg", LIBPARITWINESO);
install ("pari_acb_conj", "Gb", "acb_conj", LIBPARITWINESO);
install ("pari_acb_log", "Gb", "acb_log", LIBPARITWINESO);
install ("pari_acb_atan", "Gb", "acb_atan", LIBPARITWINESO);
install ("pari_acb_sin", "Gb", "acb_sin", LIBPARITWINESO);
install ("pari_acb_cos", "Gb", "acb_cos", LIBPARITWINESO);
install ("pari_acb_gamma", "Gb", "acb_gamma", LIBPARITWINESO);
install ("pari_acb_digamma", "Gb", "acb_digamma", LIBPARITWINESO);
install ("pari_acb_zeta", "Gb", "acb_zeta", LIBPARITWINESO);
install ("pari_acb_hurwitz_zeta", "GGb", "acb_hurwitz_zeta", LIBPARITWINESO);
install ("pari_acb_modular_eta", "Gb", "acb_modular_eta", LIBPARITWINESO);
install ("pari_acb_modular_j", "Gb", "acb_modular_j", LIBPARITWINESO);
install ("pari_acb_modular_delta", "Gb", "acb_modular_delta", LIBPARITWINESO);
install ("pari_acb_modular_eisenstein", "GLb", "acb_modular_eisenstein", LIBPARITWINESO);
install ("pari_acb_modular_theta", "GGb", "acb_modular_theta", LIBPARITWINESO);
install ("pari_acb_elliptic_p", "GGb", "acb_elliptic_p", LIBPARITWINESO);
install ("pari_acb_elliptic_inv_p", "GGb", "acb_elliptic_inv_p", LIBPARITWINESO);
install ("pari_acb_elliptic_zeta", "GGb", "acb_elliptic_zeta", LIBPARITWINESO);
install ("pari_acb_elliptic_sigma", "GGb", "acb_elliptic_sigma", LIBPARITWINESO);

install ("pari_fmpz_numbpart", "G", "fmpz_numbpart", LIBPARITWINESO);
install ("pari_arb_numbpart", "Gb", "arb_numbpart", LIBPARITWINESO);

install ("pari_cmh_I2I4I6I10", "Gb", "cmh_I2I4I6I10", LIBPARITWINESO);
install ("pari_cmh_4theta", "Gb", "cmh_4theta", LIBPARITWINESO);
install ("pari_cmh_10theta2", "Gb", "cmh_10theta2", LIBPARITWINESO);


