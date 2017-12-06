//use core::ops::*;
use core::clone::Clone;
use core::marker::PhantomData;
use dfloat::dfloat::DFloat;
use float_traits::IEEE754Float;
use roundops::*;
use safeeft::{fasttwosum, safetwosum_fma as safetwosum, safetwoproduct_fma as safetwoproduct};
use fma::{fma, Fma};

pub struct RWDFloatFMA<S: IEEE754Float + Fma + Clone, T: RoundOps<S>>(PhantomData<(S, T)>);

impl<S: IEEE754Float + Fma + Clone, T: RoundOps<S>> RoundAdd for RWDFloatFMA<S, T> {
    type Num = DFloat<S>;
    fn add_up(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.decomposite(), b.decomposite());
        let (sh, sl) = safetwosum(a_high.clone(), b_high.clone());
        if sh.is_infinite() {
            if sh == S::infinity() {
                DFloat::infinity()
            } else {
                DFloat::min_value()
            }
        } else {
            let (sh, sl) = fasttwosum(sh, T::add_up(T::add_up(a_low, b_low), sl));
            DFloat::_from_pair_raw(sh, sl)
        }
    }
    fn add_down(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.decomposite(), b.decomposite());
        let (sh, sl) = safetwosum(a_high.clone(), b_high.clone());
        if sh.is_infinite() {
            if sh == S::neg_infinity() {
                DFloat::neg_infinity()
            } else {
                DFloat::max_value()
            }
        } else {
            let (sh, sl) = fasttwosum(sh, T::add_down(T::add_down(a_low, b_low), sl));
            DFloat::_from_pair_raw(sh, sl)
        }
    }
}

impl<S: IEEE754Float + Fma + Clone, T: RoundOps<S>> RoundSub for RWDFloatFMA<S, T> {
    type Num = DFloat<S>;
    fn sub_up(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.decomposite(), b.decomposite());
        let (sh, sl) = safetwosum(a_high.clone(), -b_high.clone());
        if sh.is_infinite() {
            if sh == S::infinity() {
                DFloat::infinity()
            } else {
                DFloat::min_value()
            }
        } else {
            let (sh, sl) = fasttwosum(sh, T::add_up(T::sub_up(a_low, b_low), sl));
            DFloat::_from_pair_raw(sh, sl)
        }
    }
    fn sub_down(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.decomposite(), b.decomposite());
        let (sh, sl) = safetwosum(a_high.clone(), -b_high.clone());
        if sh.is_infinite() {
            if sh == S::neg_infinity() {
                DFloat::neg_infinity()
            } else {
                DFloat::max_value()
            }
        } else {
            let (sh, sl) = fasttwosum(sh, T::add_down(T::sub_down(a_low, b_low), sl));
            DFloat::_from_pair_raw(sh, sl)
        }
    }
}

impl<S: IEEE754Float + Fma + Clone, T: RoundOps<S>> RoundMul for RWDFloatFMA<S, T> {
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
                DFloat::min_value()
            }
        } else {
            let tmp_prod = (T::mul_up(a_low.clone(), b_high),
                            T::mul_up(a_high.clone() +
                                      (a_high.abs() *
                                       (S::eps() / S::radix() * (S::one() + S::eps())) +
                                       a_low),
                                      b_low));
            let (mh, ml) = fasttwosum(mh, T::add_up(T::add_up(ml, tmp_prod.0), tmp_prod.1));
            DFloat::_from_pair_raw(mh, ml)
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
                DFloat::max_value()
            }
        } else {
            let tmp_prod = (T::mul_down(a_low.clone(), b_high),
                            T::mul_down(a_high.clone() -
                                        (a_high.abs() *
                                         (S::eps() / S::radix() * (S::one() + S::eps())) +
                                         a_low),
                                        b_low));
            let (mh, ml) = fasttwosum(mh, T::add_down(T::add_down(ml, tmp_prod.0), tmp_prod.1));
            DFloat::_from_pair_raw(mh, ml)
        }
    }
}

impl<S: IEEE754Float + Fma + Clone, T: RoundOps<S>> RoundDiv for RWDFloatFMA<S, T> {
    type Num = DFloat<S>;
    fn div_up(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.decomposite(), b.decomposite());
        let d_high = a_high.clone() / b_high.clone();
        if d_high.is_infinite() {
            if d_high == S::infinity() {
                DFloat::infinity()
            } else {
                DFloat::min_value()
            }
        } else if b_high.is_infinite() {
            DFloat::from_single(d_high)
        } else {
            let sign_b_high = b_high.clone() / b_high.clone().abs();
            // let sign_b_high = b_high.sign();
            let (add_offset, mul_offset) = (sign_b_high.clone() *
                                            ((S::one() + S::eps()) * (S::eps() / S::radix())),
                                            sign_b_high.clone() * S::unit_underflow());
            //let tmp = mul_offset.clone();
            let (adder, muler) = (|a: &S, b: &S, p: &S| {
                                      let s = a.clone() + b.clone();
                                      s.clone() + s.abs() * p.clone()
                                  },
                                  |a: &S, b: &S, p1: &S, p2: &S| {
                                      let p = a.clone() * b.clone() + p1.clone();
                                      p.clone() + p.abs() * p2.clone()
                                  });
            let (near_ma_h, near_ma_l) = safetwoproduct(b_high.clone(), -d_high.clone());
            if near_ma_h.is_infinite() {
                let near_ma = safetwoproduct(b_high.clone() / S::radix(), -d_high.clone());
                let (near_ma_h, near_ma_l) =
                    (near_ma.0, near_ma.1 + (S::one() + S::one() + S::one()) * mul_offset.clone());
                let (half_a_h, half_b_l, half_a_l) = (a_high / S::radix() + mul_offset.clone(),
                                                      b_low.clone() / S::radix() +
                                                      mul_offset.clone(),
                                                      a_low / S::radix() + mul_offset.clone());
                let d = adder(&adder(&near_ma_h, &half_a_h, &add_offset),
                              &muler(&(-d_high.clone()), &half_b_l, &mul_offset, &add_offset),
                              &add_offset);
                let d = adder(&d, &half_a_l, &add_offset);
                let d = adder(&d, &near_ma_l, &add_offset);
                let e = if d > S::zero() {
                    T::add_down(b_high, b_low) / S::radix() - S::unit_underflow()
                } else {
                    T::add_up(b_high, b_low) / S::radix() + S::unit_underflow()
                };
                let (d_high, d_low) = fasttwosum(d_high, T::div_up(d, e));
                DFloat::_from_pair_raw(d_high, d_low)
            } else {
                let (near_ma_h, near_ma_l) =
                    (near_ma_h, near_ma_l + (S::one() + S::one() + S::one()) * mul_offset.clone());
                let d = adder(&adder(&near_ma_h, &a_high, &add_offset),
                              &muler(&(-d_high.clone()), &b_low, &mul_offset, &add_offset),
                              &add_offset);
                let d = adder(&d, &a_low, &add_offset);
                let d = adder(&d, &near_ma_l, &add_offset);
                let e = if d > S::zero() {
                    T::add_down(b_high, b_low) / S::radix() - S::unit_underflow()
                } else {
                    T::add_up(b_high, b_low) / S::radix() + S::unit_underflow()
                };
                let (d_high, d_low) = fasttwosum(d_high, T::div_up(d, e));
                DFloat::_from_pair_raw(d_high, d_low)
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
                DFloat::min_value()
            }
        } else if b_high.is_infinite() {
            DFloat::from_single(d_high)
        } else {
            let sign_b_high = b_high.clone() / b_high.clone().abs();
            // let sign_b_high = b_high.sign();
            let (add_offset, mul_offset) = (-sign_b_high.clone() *
                                            ((S::one() + S::eps()) * (S::eps() / S::radix())),
                                            -sign_b_high.clone() * S::unit_underflow());
            //let tmp = mul_offset.clone();
            let (adder, muler) = (|a: &S, b: &S, p: &S| {
                                      let s = a.clone() + b.clone();
                                      s.clone() + s.abs() * p.clone()
                                  },
                                  |a: &S, b: &S, p1: &S, p2: &S| {
                                      let p = a.clone() * b.clone() + p1.clone();
                                      p.clone() + p.abs() * p2.clone()
                                  });
            let (near_ma_h, near_ma_l) = safetwoproduct(b_high.clone(), -d_high.clone());
            if near_ma_h.is_infinite() {
                let near_ma = safetwoproduct(b_high.clone() / S::radix(), -d_high.clone());
                let (near_ma_h, near_ma_l) =
                    (near_ma.0, near_ma.1 + (S::one() + S::one() + S::one()) * mul_offset.clone());
                let (half_a_h, half_b_l, half_a_l) = (a_high / S::radix() + mul_offset.clone(),
                                                      b_low.clone() / S::radix() +
                                                      mul_offset.clone(),
                                                      a_low / S::radix() + mul_offset.clone());
                let d = adder(&adder(&near_ma_h, &half_a_h, &add_offset),
                              &muler(&(-d_high.clone()), &half_b_l, &mul_offset, &add_offset),
                              &add_offset);
                let d = adder(&d, &half_a_l, &add_offset);
                let d = adder(&d, &near_ma_l, &add_offset);
                let e = if d > S::zero() {
                    T::add_up(b_high, b_low) / S::radix() + S::unit_underflow()
                } else {
                    T::add_down(b_high, b_low) / S::radix() - S::unit_underflow()
                };
                let (d_high, d_low) = fasttwosum(d_high, T::div_down(d, e));
                DFloat::_from_pair_raw(d_high, d_low)
            } else {
                let (near_ma_h, near_ma_l) =
                    (near_ma_h, near_ma_l + (S::one() + S::one() + S::one()) * mul_offset.clone());
                let d = adder(&adder(&near_ma_h, &a_high, &add_offset),
                              &muler(&(-d_high.clone()), &b_low, &mul_offset, &add_offset),
                              &add_offset);
                let d = adder(&d, &a_low, &add_offset);
                let d = adder(&d, &near_ma_l, &add_offset);
                let e = if d > S::zero() {
                    T::add_up(b_high, b_low)
                } else {
                    T::add_down(b_high, b_low)
                };
                let (d_high, d_low) = fasttwosum(d_high, T::div_down(d, e));
                DFloat::_from_pair_raw(d_high, d_low)
            }
        }
    }
}

impl<S: IEEE754Float + Fma + Clone, T: RoundOps<S>> RoundSqrt for RWDFloatFMA<S, T> {
    fn sqrt_up(a: DFloat<S>) -> DFloat<S> {
        unimplemented!()
    }
    fn sqrt_down(a: DFloat<S>) -> DFloat<S> {
        unimplemented!()
    }
}
