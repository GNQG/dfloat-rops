//use core::ops::*;
use core::clone::Clone;
use core::marker::PhantomData;
use dfloat::DFloat;
use float_traits::IEEE754Float;
use roundops::*;
use roundops::utils::FloatSuccPred;
use safeeft::{fasttwosum, safetwosum_fma as safetwosum, safetwoproduct_fma as safetwoproduct};
use fma::{fma, Fma};

#[derive(Clone)]
pub struct RWDFloatFma<S: IEEE754Float + Fma + Clone, T: RoundOps<S>>(PhantomData<(S, T)>);

impl<S: IEEE754Float + Fma + Clone, T: RoundOps<S>> RoundingMethod for RWDFloatFma<S, T> {
    type HostMethod = T::HostMethod;
    type Num = DFloat<S>;
}

impl<S: IEEE754Float + Fma + Clone, T: RoundOps<S>> RoundAdd for RWDFloatFma<S, T> {
    fn add_up(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.into_tuple(), b.into_tuple());
        let (sh, sl) = safetwosum(a_high.clone(), b_high.clone());
        if sh.is_infinite() {
            if sh == S::infinity() {
                DFloat::infinity()
            } else {
                unsafe {
                    DFloat::from_double_components_unchecked(S::min_value().succ(),
                                                             (S::min_value() * S::eps() /
                                                              S::radix() /
                                                              S::radix())
                                                                     .pred())
                }
            }
        } else {
            #[inline(never)]
            fn tmp<S: IEEE754Float + Clone, T: RoundOps<S>>(sl: S, al: S, bl: S) -> S {
                let (w_sl, w_al, w_bl) = rnum_init!(<direction::Upward, S, T>,
                    (sl, al, bl));
                (w_sl + w_al + w_bl).extract()
            }
            let (sh, sl) = fasttwosum(sh, tmp::<S, T>(sl, a_low, b_low));
            if sh == S::neg_infinity() {
                DFloat::min_value()
            } else {
                unsafe { DFloat::from_double_components_unchecked(sh, sl) }
            }
        }
    }
    fn add_down(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.into_tuple(), b.into_tuple());
        let (sh, sl) = safetwosum(a_high.clone(), b_high.clone());
        if sh.is_infinite() {
            if sh == S::neg_infinity() {
                DFloat::neg_infinity()
            } else {
                unsafe {
                    DFloat::from_double_components_unchecked(S::max_value().pred(),
                                                             (S::max_value() * S::eps() /
                                                              S::radix() /
                                                              S::radix())
                                                                     .succ())
                }
            }
        } else {
            #[inline(never)]
            fn tmp<S: IEEE754Float + Clone, T: RoundOps<S>>(sl: S, al: S, bl: S) -> S {
                let (w_sl, w_al, w_bl) = rnum_init!(<direction::Downward, S, T>,
                    (sl, al, bl));
                (w_sl + w_al + w_bl).extract()
            }
            let (sh, sl) = fasttwosum(sh, tmp::<S, T>(sl, a_low, b_low));
            if sh == S::infinity() {
                DFloat::max_value()
            } else {
                unsafe { DFloat::from_double_components_unchecked(sh, sl) }
            }
        }
    }
}

impl<S: IEEE754Float + Fma + Clone, T: RoundOps<S>> RoundSub for RWDFloatFma<S, T> {
    fn sub_up(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.into_tuple(), b.into_tuple());
        let (sh, sl) = safetwosum(a_high.clone(), -b_high.clone());
        if sh.is_infinite() {
            if sh == S::infinity() {
                DFloat::infinity()
            } else {
                unsafe {
                    DFloat::from_double_components_unchecked(S::min_value().succ(),
                                                             (S::min_value() * S::eps() /
                                                              S::radix() /
                                                              S::radix())
                                                                     .pred())
                }
            }
        } else {
            #[inline(never)]
            fn tmp<S: IEEE754Float + Clone, T: RoundOps<S>>(sl: S, al: S, bl: S) -> S {
                let (w_sl, w_al, w_bl) = rnum_init!(<direction::Upward, S, T>,
                    (sl, al, bl));
                (w_sl + w_al - w_bl).extract()
            }
            let (sh, sl) = fasttwosum(sh, tmp::<S, T>(sl, a_low, b_low));
            if sh == S::neg_infinity() {
                DFloat::min_value()
            } else {
                unsafe { DFloat::from_double_components_unchecked(sh, sl) }
            }
        }
    }
    fn sub_down(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.into_tuple(), b.into_tuple());
        let (sh, sl) = safetwosum(a_high.clone(), -b_high.clone());
        if sh.is_infinite() {
            if sh == S::neg_infinity() {
                DFloat::neg_infinity()
            } else {
                unsafe {
                    DFloat::from_double_components_unchecked(S::max_value().pred(),
                                                             (S::max_value() * S::eps() /
                                                              S::radix() /
                                                              S::radix())
                                                                     .succ())
                }
            }
        } else {
            #[inline(never)]
            fn tmp<S: IEEE754Float + Clone, T: RoundOps<S>>(sl: S, al: S, bl: S) -> S {
                let (w_sl, w_al, w_bl) = rnum_init!(<direction::Downward, S, T>,
                    (sl, al, bl));
                (w_sl + w_al - w_bl).extract()
            }
            let (sh, sl) = fasttwosum(sh, tmp::<S, T>(sl, a_low, b_low));
            if sh == S::infinity() {
                DFloat::max_value()
            } else {
                unsafe { DFloat::from_double_components_unchecked(sh, sl) }
            }
        }
    }
}

impl<S: IEEE754Float + Fma + Clone, T: RoundOps<S>> RoundMul for RWDFloatFma<S, T> {
    fn mul_up(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.into_tuple(), b.into_tuple());
        let (mh, ml) = safetwoproduct(a_high.clone(), b_high.clone());
        let ml = ml + S::unit_underflow();

        if mh.is_infinite() {
            if mh == S::infinity() {
                DFloat::infinity()
            } else {
                unsafe {
                    DFloat::from_double_components_unchecked(S::min_value().succ().succ(),
                                                             S::min_value() * S::eps() * S::eps() /
                                                             S::radix())
                }
            }
        } else {
            let ml = {
                let (w_ml, w_ah, w_al, w_al2, w_bh, w_bl, w_bl2) =
                    rnum_init!(<direction::Upward, S, T>,
                        (ml, a_high, a_low.clone(), a_low, b_high, b_low.clone(), b_low));
                (w_ml + w_ah * w_bl + w_al * w_bh + w_al2 * w_bl2).extract()
            };
            let (mh, ml) = fasttwosum(mh, ml);
            if mh == S::neg_infinity() {
                DFloat::min_value()
            } else {
                unsafe { DFloat::from_double_components_unchecked(mh, ml) }
            }
        }
    }

    fn mul_down(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.into_tuple(), b.into_tuple());
        let (mh, ml) = safetwoproduct(a_high.clone(), b_high.clone());
        let ml = ml - S::unit_underflow();

        if mh.is_infinite() {
            if mh == S::neg_infinity() {
                DFloat::neg_infinity()
            } else {
                unsafe {
                    DFloat::from_double_components_unchecked(S::max_value().pred().pred(),
                                                             S::max_value() * S::eps() * S::eps() /
                                                             S::radix())
                }
            }
        } else {
            let ml = {
                let (w_ml, w_ah, w_al, w_al2, w_bh, w_bl, w_bl2) =
                    rnum_init!(<direction::Downward, S, T>,
                        (ml, a_high, a_low.clone(), a_low, b_high, b_low.clone(), b_low));
                (w_ml + w_ah * w_bl + w_al * w_bh + w_al2 * w_bl2).extract()
            };
            let (mh, ml) = fasttwosum(mh, ml);
            if mh == S::infinity() {
                DFloat::max_value()
            } else {
                unsafe { DFloat::from_double_components_unchecked(mh, ml) }
            }
        }
    }
}

impl<S: IEEE754Float + Fma + Clone, T: RoundOps<S>> RoundDiv for RWDFloatFma<S, T> {
    fn div_up(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.into_tuple(), b.into_tuple());
        let d_high = a_high.clone() / b_high.clone();
        if d_high.is_infinite() {
            if d_high == S::infinity() {
                DFloat::infinity()
            } else {
                unsafe {
                    DFloat::from_double_components_unchecked(S::min_value().succ(),
                                                             S::max_value().pred() * S::eps() /
                                                             S::radix() /
                                                             S::radix())
                }
            }
        } else if b_high.is_infinite() {
            DFloat::from_component(d_high)
        } else {
            let high_err = fma(b_high.clone(), -d_high.clone(), a_high);
            let low_err = fma(-b_low.clone(), d_high.clone(), a_low);
            let dl_numer = if b_high > S::zero() {
                let high_err = high_err.clone() + (S::unit_underflow() +
                               high_err.abs() * (S::eps() / S::radix() * (S::one() + S::eps())));
                let low_err = low_err.clone() + (S::unit_underflow() +
                              low_err.abs() * (S::eps() / S::radix() * (S::one() + S::eps())));
                T::add_up(high_err, low_err)
            } else {
                let high_err = high_err.clone() - (S::unit_underflow() +
                               high_err.abs() * (S::eps() / S::radix() * (S::one() + S::eps())));
                let low_err = low_err.clone() - (S::unit_underflow() +
                              low_err.abs() * (S::eps() / S::radix() * (S::one() + S::eps())));
                T::add_down(high_err, low_err)
            };

            let dl_denom = if dl_numer >= S::zero() {
                b_high.clone() - (b_high.abs() * ((S::one() + S::one() / S::radix()) * S::eps() / S::radix()))
            } else {
                b_high.clone() + (b_high.abs() * ((S::one() + S::one() / S::radix()) * S::eps() / S::radix()))
            };

            let d = fasttwosum(d_high, T::div_up(dl_numer, dl_denom));

            if d.0 == S::neg_infinity() {
                DFloat::min_value()
            } else {
                unsafe { DFloat::from_double_components_unchecked(d.0, d.1) }
            }
        }
    }
    fn div_down(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.into_tuple(), b.into_tuple());
        let d_high = a_high.clone() / b_high.clone();
        if d_high.is_infinite() {
            if d_high == S::infinity() {
                DFloat::infinity()
            } else {
                unsafe {
                    DFloat::from_double_components_unchecked(S::max_value().pred(),
                                                             S::min_value().succ() * S::eps() /
                                                             S::radix() /
                                                             S::radix())
                }
            }
        } else if b_high.is_infinite() {
            DFloat::from_component(d_high)
        } else {
            let high_err = fma(b_high.clone(), -d_high.clone(), a_high);
            let low_err = fma(-b_low.clone(), d_high.clone(), a_low);
            let dl_numer = if b_high > S::zero() {
                let high_err = high_err.clone() - (S::unit_underflow() +
                               high_err.abs() * (S::eps() / S::radix() * (S::one() + S::eps())));
                let low_err = low_err.clone() - (S::unit_underflow() +
                              low_err.abs() * (S::eps() / S::radix() * (S::one() + S::eps())));
                T::add_down(high_err, low_err)
            } else {
                let high_err = high_err.clone() + (S::unit_underflow() +
                               high_err.abs() * (S::eps() / S::radix() * (S::one() + S::eps())));
                let low_err = low_err.clone() + (S::unit_underflow() +
                              low_err.abs() * (S::eps() / S::radix() * (S::one() + S::eps())));
                T::add_up(high_err, low_err)
            };

            let dl_denom = if dl_numer >= S::zero() {
                b_high.clone() + (b_high.abs() * ((S::one() + S::one() / S::radix()) * S::eps() / S::radix()))
            } else {
                b_high.clone() - (b_high.abs() * ((S::one() + S::one() / S::radix()) * S::eps() / S::radix()))
            };

            let d = fasttwosum(d_high, T::div_down(dl_numer, dl_denom));

            if d.0 == S::infinity() {
                DFloat::max_value()
            } else {
                unsafe { DFloat::from_double_components_unchecked(d.0, d.1) }
            }
        }
    }
}

impl<S: IEEE754Float + Fma + Clone, T: RoundOps<S> + RoundSqrt> RoundSqrt for RWDFloatFma<S, T> {
    fn sqrt_up(a: DFloat<S>) -> DFloat<S> {
        let (a_h, a_l) = a.into_tuple();
        if a_h == S::infinity() {
            DFloat::infinity()
        } else if a_h == S::zero() {
            DFloat::zero()
        } else {
            let r_high = a_h.clone().sqrt();
            let ah_tester = a_h.abs() * (S::eps() / S::radix() * (S::one() + S::one() / S::radix()));
            let high_err = fma(r_high.clone(), -r_high.clone(), a_h.clone());
            let high_err = high_err.clone() + (S::unit_underflow() +
                           high_err.abs() * (S::eps() / S::radix() * (S::one() + S::eps())));
            let rl_numer = T::add_up(high_err, a_l.clone());
            let rl_denom = if rl_numer > S::zero() {
                let (w_ah_bound, w_rh) = rnum_init!(<direction::Downward, S, T>,
                    (a_h - ah_tester, r_high.clone()));
                (w_ah_bound.sqrt() + w_rh).extract()
            } else {
                let (w_ah_bound, w_rh) = rnum_init!(<direction::Upward, S, T>,
                    (a_h + ah_tester, r_high.clone()));
                (w_ah_bound.sqrt() + w_rh).extract()
            };
            let r = fasttwosum(r_high, T::div_up(rl_numer, rl_denom));
            unsafe { DFloat::from_double_components_unchecked(r.0, r.1) }
        }
    }
    fn sqrt_down(a: DFloat<S>) -> DFloat<S> {
        let (a_h, a_l) = a.into_tuple();
        if a_h == S::infinity() {
            DFloat::infinity()
        } else if a_h == S::zero() {
            DFloat::zero()
        } else {
            let r_high = a_h.clone().sqrt();
            let ah_tester = a_h.abs() * (S::eps() / S::radix() * (S::one() + S::one() / S::radix()));
            let high_err = fma(r_high.clone(), -r_high.clone(), a_h.clone());
            let high_err = high_err.clone() - (S::unit_underflow() +
                           high_err.abs() * (S::eps() / S::radix() * (S::one() + S::eps())));
            let rl_numer = T::add_down(high_err, a_l.clone());
            let rl_denom = if rl_numer > S::zero() {
                let (w_ah_bound, w_rh) = rnum_init!(<direction::Upward, S, T>,
                    (a_h + ah_tester, r_high.clone()));
                (w_ah_bound.sqrt() + w_rh).extract()
            } else {
                let (w_ah_bound, w_rh) = rnum_init!(<direction::Downward, S, T>,
                    (a_h - ah_tester, r_high.clone()));
                (w_ah_bound.sqrt() + w_rh).extract()
            };
            let r = fasttwosum(r_high, T::div_down(rl_numer, rl_denom));
            unsafe { DFloat::from_double_components_unchecked(r.0, r.1) }
        }
    }
}
