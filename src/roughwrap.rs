use core::clone::Clone;
use core::marker::PhantomData;
use dfloat::DFloat;
use float_traits::IEEE754Float;
use roundops::*;
use roundops::utils::FloatSuccPred;
use safeeft::{fasttwosum, safetwosum_straight as safetwosum,
              safetwoproduct_branch as safetwoproduct};

pub struct RWDFloatRegular<S: IEEE754Float + Clone, T: RoundOps<S>>(PhantomData<(S, T)>);

impl<S: IEEE754Float + Clone, T: RoundOps<S>> RoundAdd for RWDFloatRegular<S, T> {
    type Num = DFloat<S>;
    fn add_up(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.decomposite(), b.decomposite());
        let (sh, sl) = safetwosum(a_high.clone(), b_high.clone());
        if sh.is_infinite() {
            if sh == S::infinity() {
                DFloat::infinity()
            } else {
                DFloat::_from_pair_raw(S::min_value().succ(),
                                       (S::min_value() * S::eps() / S::radix() / S::radix()).pred())
            }
        } else {
            let (sh, sl) = fasttwosum(sh, T::add_up(T::add_up(a_low, b_low), sl));
            if sh == S::neg_infinity() {
                DFloat::min_value()
            } else {
                DFloat::_from_pair_raw(sh, sl)
            }
        }
    }
    fn add_down(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.decomposite(), b.decomposite());
        let (sh, sl) = safetwosum(a_high.clone(), b_high.clone());
        if sh.is_infinite() {
            if sh == S::neg_infinity() {
                DFloat::neg_infinity()
            } else {
                DFloat::_from_pair_raw(S::max_value().pred(),
                                       (S::max_value() * S::eps() / S::radix() / S::radix()).succ())
            }
        } else {
            let (sh, sl) = fasttwosum(sh, T::add_down(T::add_down(a_low, b_low), sl));
            if sh == S::infinity() {
                DFloat::max_value()
            } else {
                DFloat::_from_pair_raw(sh, sl)
            }
        }
    }
}

impl<S: IEEE754Float + Clone, T: RoundOps<S>> RoundSub for RWDFloatRegular<S, T> {
    type Num = DFloat<S>;
    fn sub_up(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.decomposite(), b.decomposite());
        let (sh, sl) = safetwosum(a_high.clone(), -b_high.clone());
        if sh.is_infinite() {
            if sh == S::infinity() {
                DFloat::infinity()
            } else {
                DFloat::_from_pair_raw(S::min_value().succ(),
                                       (S::min_value() * S::eps() / S::radix() / S::radix()).pred())
            }
        } else {
            let (sh, sl) = fasttwosum(sh, T::add_up(T::sub_up(a_low, b_low), sl));
            if sh == S::neg_infinity() {
                DFloat::min_value()
            } else {
                DFloat::_from_pair_raw(sh, sl)
            }
        }
    }
    fn sub_down(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.decomposite(), b.decomposite());
        let (sh, sl) = safetwosum(a_high.clone(), -b_high.clone());
        if sh.is_infinite() {
            if sh == S::neg_infinity() {
                DFloat::neg_infinity()
            } else {
                DFloat::_from_pair_raw(S::max_value().pred(),
                                       (S::max_value() * S::eps() / S::radix() / S::radix()).succ())
            }
        } else {
            let (sh, sl) = fasttwosum(sh, T::add_down(T::sub_down(a_low, b_low), sl));
            if sh == S::infinity() {
                DFloat::max_value()
            } else {
                DFloat::_from_pair_raw(sh, sl)
            }
        }
    }
}

impl<S: IEEE754Float + Clone, T: RoundOps<S>> RoundMul for RWDFloatRegular<S, T> {
    type Num = DFloat<S>;
    fn mul_up(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.decomposite(), b.decomposite());
        let (mh, ml) = safetwoproduct(a_high.clone(), b_high.clone());
        // S. Boldo, Pitfalls of a Full Floating-Point Proof: Example on the Formal Proof of
        //           the Veltkamp/Dekker Algorithms, IJCAR 2006: Automated Reasoning, 2006
        let ml = ml + (S::one() + S::one() + S::one()) * S::unit_underflow();

        if mh.is_infinite() {
            if mh == S::infinity() {
                DFloat::infinity()
            } else {
                DFloat::_from_pair_raw(S::min_value().succ().succ(),
                                       S::min_value() * S::eps() * S::eps() / S::radix())
            }
        } else {
            let (ahbl, albh, albl) = (T::mul_up(a_high, b_low.clone()),
                                      T::mul_up(a_low.clone(), b_high),
                                      T::mul_up(a_low, b_low));
            let (mh, ml) = fasttwosum(mh, T::add_up(T::add_up(ml, ahbl), T::add_up(albh, albl)));
            if mh == S::neg_infinity() {
                DFloat::min_value()
            } else {
                DFloat::_from_pair_raw(mh, ml)
            }
        }
    }

    fn mul_down(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.decomposite(), b.decomposite());
        let (mh, ml) = safetwoproduct(a_high.clone(), b_high.clone());
        // S. Boldo, Pitfalls of a Full Floating-Point Proof: Example on the Formal Proof of
        //           the Veltkamp/Dekker Algorithms, IJCAR 2006: Automated Reasoning, 2006
        let ml = ml - (S::one() + S::one() + S::one()) * S::unit_underflow();

        if mh.is_infinite() {
            if mh == S::neg_infinity() {
                DFloat::neg_infinity()
            } else {
                DFloat::_from_pair_raw(S::max_value().pred().pred(),
                                       S::max_value() * S::eps() * S::eps() / S::radix())
            }
        } else {
            let (ahbl, albh, albl) = (T::mul_down(a_high, b_low.clone()),
                                      T::mul_down(a_low.clone(), b_high),
                                      T::mul_down(a_low, b_low));
            let (mh, ml) = fasttwosum(mh,
                                      T::add_down(T::add_down(ml, ahbl), T::add_down(albh, albl)));
            if mh == S::infinity() {
                DFloat::max_value()
            } else {
                DFloat::_from_pair_raw(mh, ml)
            }
        }
    }
}

impl<S: IEEE754Float + Clone, T: RoundOps<S>> RoundDiv for RWDFloatRegular<S, T> {
    type Num = DFloat<S>;
    fn div_up(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.decomposite(), b.decomposite());
        let d_high = a_high.clone() / b_high.clone();
        if d_high.is_infinite() {
            if d_high == S::infinity() {
                DFloat::infinity()
            } else {
                DFloat::_from_pair_raw(S::min_value().succ(),
                                       S::max_value().pred() * S::eps() / S::radix() / S::radix())
            }
        } else if b_high.is_infinite() {
            DFloat::from_single(d_high)
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
                T::add_up(T::add_up(err_a_h, err_a_l),
                          T::add_up(T::mul_up(-d_high.clone(), b_low.clone()), a_low))
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
                T::add_down(T::add_down(err_a_h, err_a_l),
                            T::add_down(T::mul_down(-d_high.clone(), b_low.clone()), a_low))
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
                DFloat::_from_pair_raw(dh, dl)
            }
        }
    }
    fn div_down(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.decomposite(), b.decomposite());
        let d_high = a_high.clone() / b_high.clone();
        if d_high.is_infinite() {
            if d_high == S::infinity() {
                DFloat::infinity()
            } else {
                DFloat::_from_pair_raw(S::max_value().pred(),
                                       S::min_value().succ() * S::eps() / S::radix() / S::radix())
            }
        } else if b_high.is_infinite() {
            DFloat::from_single(d_high)
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
                T::add_down(T::add_down(err_a_h, err_a_l),
                            T::add_down(T::mul_down(-d_high.clone(), b_low.clone()), a_low))
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
                T::add_up(T::add_up(err_a_h, err_a_l),
                          T::add_up(T::mul_up(-d_high.clone(), b_low.clone()), a_low))
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
                DFloat::_from_pair_raw(d.0, d.1)
            }
        }
    }
}

impl<S: IEEE754Float + Clone, T: RoundOps<S> + RoundSqrt> RoundSqrt for RWDFloatRegular<S, T> {
    fn sqrt_up(a: DFloat<S>) -> DFloat<S> {
        let (a_h, a_l) = a.decomposite();
        if a_h == S::infinity() {
            DFloat::infinity()
        } else if a_h == S::zero() {
            DFloat::zero()
        } else {
            let r_high = a_h.clone().sqrt();
            let (neg_a_high, neg_a_low) = safetwoproduct(-r_high.clone(), r_high.clone());
            let neg_a_low = neg_a_low + (S::one() + S::one() + S::one()) * S::unit_underflow();
            let rl_numer = T::add_up(T::add_up(neg_a_high + a_h.clone(), a_l.clone()), // safe
                                     neg_a_low.clone());
            let rl_denom = if rl_numer > S::zero() {
                T::add_down(T::sqrt_down(T::add_down(a_h, a_l)), r_high.clone())
            } else {
                T::add_up(T::sqrt_up(T::add_up(a_h, a_l)), r_high.clone())
            };
            let r = fasttwosum(r_high, T::div_up(rl_numer, rl_denom));
            DFloat::_from_pair_raw(r.0, r.1)
        }
    }
    fn sqrt_down(a: DFloat<S>) -> DFloat<S> {
        let (a_h, a_l) = a.decomposite();
        if a_h == S::infinity() {
            DFloat::infinity()
        } else if a_h == S::zero() {
            DFloat::zero()
        } else {
            let r_high = a_h.clone().sqrt();
            let (neg_a_high, neg_a_low) = safetwoproduct(-r_high.clone(), r_high.clone());
            let neg_a_low = neg_a_low - (S::one() + S::one() + S::one()) * S::unit_underflow();
            let rl_numer = T::add_down(T::add_down(neg_a_high + a_h.clone(), a_l.clone()), // safe
                                       neg_a_low.clone());
            let rl_denom = if rl_numer > S::zero() {
                T::add_up(T::sqrt_up(T::add_up(a_h, a_l)), r_high.clone())
            } else {
                T::add_down(T::sqrt_down(T::add_down(a_h, a_l)), r_high.clone())
            };
            let r = fasttwosum(r_high, T::div_down(rl_numer, rl_denom));
            DFloat::_from_pair_raw(r.0, r.1)
        }
    }
}
