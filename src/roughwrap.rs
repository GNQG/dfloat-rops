use core::clone::Clone;
use core::marker::PhantomData;
use dfloat::DFloat;
use float_traits::IEEE754Float;
use roundops::*;
use roundops::utils::FloatSuccPred;
use safeeft::{fasttwosum, safetwosum_straight as safetwosum,
              safetwoproduct_branch as safetwoproduct};

pub struct RWDFloatRegular<S: IEEE754Float + Clone, T: RoundOps<S>>(PhantomData<(S, T)>);

impl<S: IEEE754Float + Clone, T: RoundOps<S>> RoundingMethod for RWDFloatRegular<S, T> {
    type HostMethod = T::HostMethod;
    type Num = DFloat<S>;
}

impl<S: IEEE754Float + Clone, T: RoundOps<S>> RoundAdd for RWDFloatRegular<S, T> {
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

impl<S: IEEE754Float + Clone, T: RoundOps<S>> RoundSub for RWDFloatRegular<S, T> {
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

impl<S: IEEE754Float + Clone, T: RoundOps<S>> RoundMul for RWDFloatRegular<S, T> {
    fn mul_up(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.into_tuple(), b.into_tuple());
        let (mh, ml) = safetwoproduct(a_high.clone(), b_high.clone());
        // S. Boldo, Pitfalls of a Full Floating-Point Proof: Example on the Formal Proof of
        //           the Veltkamp/Dekker Algorithms, IJCAR 2006: Automated Reasoning, 2006
        let ml = ml + (S::one() + S::one() + S::one()) * S::unit_underflow();

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
        // S. Boldo, Pitfalls of a Full Floating-Point Proof: Example on the Formal Proof of
        //           the Veltkamp/Dekker Algorithms, IJCAR 2006: Automated Reasoning, 2006
        let ml = ml - (S::one() + S::one() + S::one()) * S::unit_underflow();

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

impl<S: IEEE754Float + Clone, T: RoundOps<S>> RoundDiv for RWDFloatRegular<S, T> {
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
            let (neg_a_high, neg_a_low) = safetwoproduct(b_high.clone(), -d_high.clone());
            let dl_numer = if b_high > S::zero() {
                let (err_a_h, err_a_l) = if neg_a_high.is_infinite() {
                    // safe(b_high is not too small)
                    let (neg_a2_h, neg_a2_l) = safetwoproduct(b_high.clone() / S::radix(),
                                                              -d_high.clone());
                    // safe(-a_high/2 ~ neg_a2_h)
                    ((a_high / S::radix() + neg_a2_h) * S::radix(), neg_a2_l * S::radix())
                } else {
                    (a_high + neg_a_high,
                     neg_a_low + (S::one() + S::one() + S::one()) * S::unit_underflow())
                };
                let (w_al, w_bl, w_ah_e, w_al_e, w_dh) =
                    rnum_init!(<direction::Upward, S, T>,
                        (a_low,b_low.clone(), err_a_h, err_a_l, d_high.clone()));
                ((w_ah_e + w_al_e) + ((-w_dh) * w_bl + w_al)).extract()
            } else {
                let (err_a_h, err_a_l) = if neg_a_high.is_infinite() {
                    // safe(b_high is not too small)
                    let (neg_a2_h, neg_a2_l) = safetwoproduct(b_high.clone() / S::radix(),
                                                              -d_high.clone());
                    // safe(-a_high/2 ~ neg_a2_h)
                    ((a_high / S::radix() + neg_a2_h) * S::radix(), neg_a2_l * S::radix())
                } else {
                    // -a_high ~ neg_a_high
                    (a_high + neg_a_high,
                     neg_a_low - (S::one() + S::one() + S::one()) * S::unit_underflow())
                };
                let (w_al, w_bl, w_ah_e, w_al_e, w_dh) =
                    rnum_init!(<direction::Downward, S, T>,
                        (a_low,b_low.clone(), err_a_h, err_a_l, d_high.clone()));
                ((w_ah_e + w_al_e) + ((-w_dh) * w_bl + w_al)).extract()
            };

            let dl_denom = if dl_numer >= S::zero() {
                T::add_down(b_high, b_low)
            } else {
                T::add_up(b_high, b_low)
            };

            let (dh, dl) = fasttwosum(d_high, T::div_up(dl_numer, dl_denom));

            if dh == S::neg_infinity() {
                DFloat::min_value()
            } else {
                unsafe { DFloat::from_double_components_unchecked(dh, dl) }
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
            let (neg_a_high, neg_a_low) = safetwoproduct(b_high.clone(), -d_high.clone());
            let dl_numer = if b_high > S::zero() {
                let (err_a_h, err_a_l) = if neg_a_high.is_infinite() {
                    // safe(b_high is not too small)
                    let (neg_a2_h, neg_a2_l) = safetwoproduct(b_high.clone() / S::radix(),
                                                              -d_high.clone());
                    // safe(-a_high/2 ~ neg_a2_h)
                    ((a_high / S::radix() + neg_a2_h) * S::radix(), neg_a2_l * S::radix())
                } else {
                    (a_high + neg_a_high,
                     neg_a_low - (S::one() + S::one() + S::one()) * S::unit_underflow())
                };
                let (w_al, w_bl, w_ah_e, w_al_e, w_dh) =
                    rnum_init!(<direction::Downward, S, T>,
                        (a_low, b_low.clone(), err_a_h, err_a_l, d_high.clone()));
                ((w_ah_e + w_al_e) + ((-w_dh) * w_bl + w_al)).extract()
            } else {
                let (err_a_h, err_a_l) = if neg_a_high.is_infinite() {
                    // safe(b_high is not too small)
                    let (neg_a2_h, neg_a2_l) = safetwoproduct(b_high.clone() / S::radix(),
                                                              -d_high.clone());
                    // safe(-a_high/2 ~ neg_a2_h)
                    ((a_high / S::radix() + neg_a2_h) * S::radix(), neg_a2_l * S::radix())
                } else {
                    // -a_high ~ neg_a_high
                    (a_high + neg_a_high,
                     neg_a_low + (S::one() + S::one() + S::one()) * S::unit_underflow())
                };
                let (w_al, w_bl, w_ah_e, w_al_e, w_dh) =
                    rnum_init!(<direction::Upward, S, T>,
                        (a_low, b_low.clone(), err_a_h, err_a_l, d_high.clone()));
                ((w_ah_e + w_al_e) + ((-w_dh) * w_bl + w_al)).extract()
            };

            let dl_denom = if dl_numer >= S::zero() {
                T::add_up(b_high, b_low)
            } else {
                T::add_down(b_high, b_low)
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

impl<S: IEEE754Float + Clone, T: RoundOps<S> + RoundSqrt> RoundSqrt for RWDFloatRegular<S, T> {
    fn sqrt_up(a: DFloat<S>) -> DFloat<S> {
        let (a_high, a_low) = a.into_tuple();
        if a_high == S::infinity() {
            DFloat::infinity()
        } else if a_high == S::zero() {
            DFloat::zero()
        } else {
            let r_high = a_high.clone().sqrt();
            let (neg_a_high, neg_a_low) = safetwoproduct(-r_high.clone(), r_high.clone());
            let neg_a_low = neg_a_low + (S::one() + S::one() + S::one()) * S::unit_underflow();
            let ah_tester = a_high.abs() * (S::eps() / S::radix() * (S::one() + S::eps()));
            let (w_al, w_ah_e, w_neg_al) =
                rnum_init!(<direction::Upward, S, T>,
                    (a_low.clone(), neg_a_high + a_high.clone(), neg_a_low.clone()));
            let rl_numer = (w_ah_e + (w_al + w_neg_al)).extract();
            let rl_denom = if rl_numer > S::zero() {
                let (w_ah_bound, w_rh) = rnum_init!(<direction::Downward, S, T>,
                    (a_high + (a_low - ah_tester), r_high.clone()));
                (w_ah_bound.sqrt() + w_rh).extract()
            } else {
                let (w_ah_bound, w_rh) = rnum_init!(<direction::Upward, S, T>,
                    (a_high + (a_low + ah_tester), r_high.clone()));
                (w_ah_bound.sqrt() + w_rh).extract()
            };
            let r = fasttwosum(r_high, T::div_up(rl_numer, rl_denom));
            unsafe { DFloat::from_double_components_unchecked(r.0, r.1) }
        }
    }
    fn sqrt_down(a: DFloat<S>) -> DFloat<S> {
        let (a_high, a_low) = a.into_tuple();
        if a_high == S::infinity() {
            DFloat::infinity()
        } else if a_high == S::zero() {
            DFloat::zero()
        } else {
            let r_high = a_high.clone().sqrt();
            let (neg_a_high, neg_a_low) = safetwoproduct(-r_high.clone(), r_high.clone());
            let neg_a_low = neg_a_low - (S::one() + S::one() + S::one()) * S::unit_underflow();
            let ah_tester = a_high.abs() * (S::eps() / S::radix() * (S::one() + S::eps()));
            let (w_al, w_ah_e, w_neg_al) =
                rnum_init!(<direction::Downward, S, T>,
                    (a_low.clone(), neg_a_high + a_high.clone(), neg_a_low.clone()));
            let rl_numer = (w_ah_e + (w_al + w_neg_al)).extract();
            let rl_denom = if rl_numer > S::zero() {
                let (w_ah_bound, w_rh) = rnum_init!(<direction::Upward, S, T>,
                    (a_high + (a_low + ah_tester), r_high.clone()));
                (w_ah_bound.sqrt() + w_rh).extract()
            } else {
                let (w_ah_bound, w_rh) = rnum_init!(<direction::Downward, S, T>,
                    (a_high + (a_low - ah_tester), r_high.clone()));
                (w_ah_bound.sqrt() + w_rh).extract()
            };
            let r = fasttwosum(r_high, T::div_down(rl_numer, rl_denom));
            unsafe { DFloat::from_double_components_unchecked(r.0, r.1) }
        }
    }
}
