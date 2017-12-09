//use core::ops::*;
use core::clone::Clone;
use core::marker::PhantomData;
use dfloat::DFloat;
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
        let ml = ml + S::unit_underflow();

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
        let ml = ml - S::unit_underflow();

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
            let high_err = fma(b_high.clone(), -d_high.clone(), a_high);
            let low_err = fma(-b_low.clone(), d_high.clone(), a_low);
            if b_high > S::zero() {
                let high_err = (high_err.clone() + S::unit_underflow()) +
                               high_err.abs() * (S::eps() / S::radix() * (S::one() + S::eps()));
                let low_err = (low_err.clone() + S::unit_underflow()) +
                              low_err.abs() * (S::eps() / S::radix() * (S::one() + S::eps()));
                let d_low_tmp = T::add_up(high_err, low_err);
                let tmp = if d_low_tmp >= S::zero() {
                    T::add_down(b_high, b_low)
                } else {
                    T::add_up(b_high, b_low)
                };
                let d = safetwosum(d_high, T::div_up(d_low_tmp, tmp));
                DFloat::_from_pair_raw(d.0, d.1)
            } else {
                let high_err = (high_err.clone() - S::unit_underflow()) -
                               high_err.abs() * (S::eps() / S::radix() * (S::one() + S::eps()));
                let low_err = (low_err.clone() - S::unit_underflow()) -
                              low_err.abs() * (S::eps() / S::radix() * (S::one() + S::eps()));
                let d_low_tmp = T::add_down(high_err, low_err);
                let tmp = if d_low_tmp >= S::zero() {
                    T::add_down(b_high, b_low)
                } else {
                    T::add_up(b_high, b_low)
                };
                let d = safetwosum(d_high, T::div_up(d_low_tmp, tmp));
                DFloat::_from_pair_raw(d.0, d.1)
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
            let high_err = fma(b_high.clone(), -d_high.clone(), a_high);
            let low_err = fma(-b_low.clone(), d_high.clone(), a_low);
            if b_high > S::zero() {
                let high_err = (high_err.clone() - S::unit_underflow()) -
                               high_err.abs() * (S::eps() / S::radix() * (S::one() + S::eps()));
                let low_err = (low_err.clone() - S::unit_underflow()) -
                              low_err.abs() * (S::eps() / S::radix() * (S::one() + S::eps()));
                let d_low_tmp = T::add_down(high_err, low_err);
                let tmp = if d_low_tmp >= S::zero() {
                    T::add_up(b_high, b_low)
                } else {
                    T::add_down(b_high, b_low)
                };
                let d = safetwosum(d_high, T::div_down(d_low_tmp, tmp));
                DFloat::_from_pair_raw(d.0, d.1)
            } else {
                let high_err = (high_err.clone() + S::unit_underflow()) +
                               high_err.abs() * (S::eps() / S::radix() * (S::one() + S::eps()));
                let low_err = (low_err.clone() + S::unit_underflow()) +
                              low_err.abs() * (S::eps() / S::radix() * (S::one() + S::eps()));
                let d_low_tmp = T::add_up(high_err, low_err);
                let tmp = if d_low_tmp >= S::zero() {
                    T::add_up(b_high, b_low)
                } else {
                    T::add_down(b_high, b_low)
                };
                let d = safetwosum(d_high, T::div_down(d_low_tmp, tmp));
                DFloat::_from_pair_raw(d.0, d.1)
            }
        }
    }
}

impl<S: IEEE754Float + Fma + Clone, T: RoundOps<S> + RoundSqrt> RoundSqrt for RWDFloatFMA<S, T> {
    fn sqrt_up(a: DFloat<S>) -> DFloat<S> {
        let (a_h, a_l) = a.decomposite();
        if a_h == S::infinity() {
            DFloat::infinity()
        } else if a_h == S::zero() {
            DFloat::zero()
        } else {
            let r_high = a_h.clone().sqrt();
            let high_err = fma(r_high.clone(), -r_high.clone(), a_h.clone());
            let high_err = (high_err.clone() + S::unit_underflow()) +
                           high_err.abs() * (S::eps() / S::radix() * (S::one() + S::eps()));
            let r_low_tmp = T::add_up(high_err, a_l.clone());
            let tmp = if r_low_tmp > S::zero() {
                T::add_down(T::sqrt_down(T::add_down(a_h, a_l)), r_high.clone())
            } else {
                T::add_up(T::sqrt_up(T::add_up(a_h, a_l)), r_high.clone())
            };
            let r = safetwosum(r_high, T::div_up(r_low_tmp, tmp));
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
            let high_err = fma(r_high.clone(), -r_high.clone(), a_h.clone());
            let high_err = (high_err.clone() - S::unit_underflow()) -
                           high_err.abs() * (S::eps() / S::radix() * (S::one() + S::eps()));
            let r_low_tmp = T::add_down(high_err, a_l.clone());
            let tmp = if r_low_tmp > S::zero() {
                T::add_up(T::sqrt_up(T::add_up(a_h, a_l)), r_high.clone())
            } else {
                T::add_down(T::sqrt_down(T::add_down(a_h, a_l)), r_high.clone())
            };
            let r = safetwosum(r_high, T::div_down(r_low_tmp, tmp));
            DFloat::_from_pair_raw(r.0, r.1)
        }
    }
}
