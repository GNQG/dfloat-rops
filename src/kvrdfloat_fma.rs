//use core::ops::*;
use core::clone::Clone;
use core::marker::PhantomData;
use dfloat::dfloat::DFloat;
use float_traits::IEEE754Float;
use roundops::*;
use safeeft::{fasttwosum, safetwosum_branch as safetwosum,
              safetwoproduct_fma as safetwoproduct};
use fma::{fma,Fma};

pub struct KVRDFloatFMA<S: IEEE754Float + Fma + Clone, T: RoundOps<S>>(PhantomData<(S, T)>);

impl<S: IEEE754Float + Fma + Clone, T: RoundOps<S>> RoundAdd for KVRDFloatFMA<S, T> {
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
            let (sh, sl) = safetwosum(sh, T::add_up(T::add_up(a_low, b_low), sl));
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
            let (sh, sl) = safetwosum(sh, T::add_down(T::add_down(a_low, b_low), sl));
            DFloat::_from_pair_raw(sh, sl)
        }
    }
}

impl<S: IEEE754Float + Fma + Clone, T: RoundOps<S>> RoundSub for KVRDFloatFMA<S, T> {
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
            let (sh, sl) = safetwosum(sh, T::add_up(T::sub_up(a_low, b_low), sl));
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
            let (sh, sl) = safetwosum(sh, T::add_down(T::sub_down(a_low, b_low), sl));
            DFloat::_from_pair_raw(sh, sl)
        }
    }
}

impl<S: IEEE754Float + Fma + Clone, T: RoundOps<S>> RoundMul for KVRDFloatFMA<S, T> {
    type Num = DFloat<S>;
    fn mul_up(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.decomposite(), b.decomposite());
        let (mh, ml) = safetwoproduct(a_high.clone(), b_high.clone());
        if mh.is_infinite() {
            if mh == S::infinity() {
                DFloat::infinity()
            } else {
                DFloat::min_value()
            }
        } else {
            let ml = ml + (S::one() + S::one() + S::one()) * S::unit_underflow();
            let (ahbl, albh, albl) = (T::mul_up(a_high, b_low.clone()),
                                      T::mul_up(a_low.clone(), b_high),
                                      T::mul_up(a_low, b_low));
            let (mh, ml) = safetwosum(mh, T::add_up(ml, T::add_up(T::add_up(albl, albh), ahbl)));
            DFloat::_from_pair_raw(mh, ml)
        }
    }

    fn mul_down(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((a_high, a_low), (b_high, b_low)) = (a.decomposite(), b.decomposite());
        let (mh, ml) = safetwoproduct(a_high.clone(), b_high.clone());
        if mh.is_infinite() {
            if mh == S::neg_infinity() {
                DFloat::neg_infinity()
            } else {
                DFloat::max_value()
            }
        } else {
            let ml = ml + (S::one() + S::one() + S::one()) * S::unit_underflow();
            let (ahbl, albh, albl) = (T::mul_down(a_high, b_low.clone()),
                                      T::mul_down(a_low.clone(), b_high),
                                      T::mul_down(a_low, b_low));
            let (mh, ml) = safetwosum(mh,
                                      T::add_down(ml, T::add_down(T::add_down(albl, albh), ahbl)));
            DFloat::_from_pair_raw(mh, ml)
        }
    }
}

impl<S: IEEE754Float + Fma + Clone, T: RoundOps<S>> RoundDiv for KVRDFloatFMA<S, T> {
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
            if b_high > S::zero() {
                let (near_ma_h, near_ma_l) = safetwoproduct(b_high.clone(), -d_high.clone());
                if near_ma_h.is_infinite() {
                    // b_high is not too small. -> b_high / S::radix() is safe.
                    let (near_ma_h, near_ma_l) = safetwoproduct(b_high.clone() / S::radix(),
                                                                -d_high.clone());
                    let d_low_tmp =                                           //  __| safe |__
                        T::add_up(T::add_up(T::add_up(T::add_up(near_ma_h, a_high / S::radix()),
                                                      T::mul_up(-d_high.clone(), b_low.clone() / S::radix())),
                                            a_low / S::radix()), // <-    provably unsafe ->     """"""""""
                                  near_ma_l);
                    let tmp = if d_low_tmp >= S::zero() {
                        T::add_down(b_high / S::radix(), b_low / S::radix())
                    } else {
                        //                                     ^ provably unsafe
                        T::add_up(b_high / S::radix(), b_low / S::radix())
                    }; //                                      ^ provably unsafe
                    let d = safetwosum(d_high, T::div_up(d_low_tmp, tmp));
                    DFloat::_from_pair_raw(d.0, d.1)
                } else {
                    let d_low_tmp = T::add_up(T::add_up(T::add_up(T::add_up(near_ma_h, a_high),
                                                                  T::mul_up(-d_high.clone(),
                                                                            b_low.clone())),
                                                        a_low),
                                              near_ma_l);
                    let tmp = if d_low_tmp >= S::zero() {
                        T::add_down(b_high, b_low)
                    } else {
                        T::add_up(b_high, b_low)
                    };
                    let d = safetwosum(d_high, T::div_up(d_low_tmp, tmp));
                    DFloat::_from_pair_raw(d.0, d.1)
                }
            } else {
                let (near_ma_h, near_ma_l) = safetwoproduct(b_high.clone(), -d_high.clone());
                if near_ma_h.is_infinite() {
                    let (near_ma_h, near_ma_l) = safetwoproduct(b_high.clone() / S::radix(),
                                                                -d_high.clone());
                    let d_low_tmp =
                        T::add_down(T::add_down(T::add_down(T::add_down(near_ma_h,
                                                                        a_high / S::radix()),
                                                            T::mul_down(-d_high.clone(),
                                                                        b_low.clone() /
                                                                        S::radix())),
                                                a_low / S::radix()),
                                    near_ma_l);
                    let tmp = if d_low_tmp >= S::zero() {
                        T::add_down(b_high / S::radix(), b_low / S::radix())
                    } else {
                        T::add_up(b_high / S::radix(), b_low / S::radix())
                    };
                    let d = safetwosum(d_high, T::div_up(d_low_tmp, tmp));
                    DFloat::_from_pair_raw(d.0, d.1)
                } else {
                    let d_low_tmp =
                        T::add_down(T::add_down(T::add_down(T::add_down(near_ma_h, a_high),
                                                            T::mul_down(-d_high.clone(),
                                                                        b_low.clone())),
                                                a_low),
                                    near_ma_l);
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
            if b_high > S::zero() {
                let (near_ma_h, near_ma_l) = safetwoproduct(b_high.clone(), -d_high.clone());
                if near_ma_h.is_infinite() {
                    // b_high is not too small. -> b_high / S::radix() is safe.
                    let (near_ma_h, near_ma_l) = safetwoproduct(b_high.clone() / S::radix(),
                                                                -d_high.clone());
                    let d_low_tmp =                                           //  __| safe |__
                        T::add_down(T::add_down(T::add_down(T::add_down(near_ma_h, a_high / S::radix()),
                                                      T::mul_down(-d_high.clone(), b_low.clone() / S::radix())),
                                            a_low / S::radix()), // <-    provably unsafe ->     """"""""""
                                  near_ma_l);
                    let tmp = if d_low_tmp >= S::zero() {
                        T::add_up(b_high / S::radix(), b_low / S::radix())
                    } else {
                        //                                     ^ provably unsafe
                        T::add_down(b_high / S::radix(), b_low / S::radix())
                    }; //                                      ^ provably unsafe
                    let d = safetwosum(d_high, T::div_down(d_low_tmp, tmp));
                    DFloat::_from_pair_raw(d.0, d.1)
                } else {
                    let d_low_tmp =
                        T::add_down(T::add_down(T::add_down(T::add_down(near_ma_h, a_high),
                                                            T::mul_down(-d_high.clone(),
                                                                        b_low.clone())),
                                                a_low),
                                    near_ma_l);
                    let tmp = if d_low_tmp >= S::zero() {
                        T::add_up(b_high, b_low)
                    } else {
                        T::add_down(b_high, b_low)
                    };
                    let d = safetwosum(d_high, T::div_down(d_low_tmp, tmp));
                    DFloat::_from_pair_raw(d.0, d.1)
                }
            } else {
                let (near_ma_h, near_ma_l) = safetwoproduct(b_high.clone(), -d_high.clone());
                if near_ma_h.is_infinite() {
                    let (near_ma_h, near_ma_l) = safetwoproduct(b_high.clone() / S::radix(),
                                                                -d_high.clone());
                    let d_low_tmp =
                        T::add_up(T::add_up(T::add_up(T::add_up(near_ma_h, a_high / S::radix()),
                                                      T::mul_up(-d_high.clone(),
                                                                b_low.clone() / S::radix())),
                                            a_low / S::radix()),
                                  near_ma_l);
                    let tmp = if d_low_tmp >= S::zero() {
                        T::add_up(b_high / S::radix(), b_low / S::radix())
                    } else {
                        T::add_down(b_high / S::radix(), b_low / S::radix())
                    };
                    let d = safetwosum(d_high, T::div_down(d_low_tmp, tmp));
                    DFloat::_from_pair_raw(d.0, d.1)
                } else {
                    let d_low_tmp = T::add_up(T::add_up(T::add_up(T::add_up(near_ma_h, a_high),
                                                                  T::mul_up(-d_high.clone(),
                                                                            b_low.clone())),
                                                        a_low),
                                              near_ma_l);
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
}

impl<S: IEEE754Float + Fma + Clone, T: RoundOps<S> + RoundSqrt> RoundSqrt for KVRDFloatFMA<S, T> {
    fn sqrt_up(a: DFloat<S>) -> DFloat<S> {
        let (a_h, a_l) = a.decomposite();
        if a_h == S::infinity() {
            DFloat::infinity()
        } else if a_h == S::zero() {
            DFloat::zero()
        } else {
            let r_high = a_h.clone().sqrt();
            let (near_ma_h, near_ma_l) = safetwoproduct(-r_high.clone(), r_high.clone());
            let r_low_tmp = T::add_up(T::add_up(T::add_up(near_ma_h, a_h.clone()), a_l.clone()),
                                      near_ma_l.clone());
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
            let (near_ma_h, near_ma_l) = safetwoproduct(-r_high.clone(), r_high.clone());
            let r_low_tmp = T::add_down(T::add_down(T::add_down(near_ma_h, a_h.clone()),
                                                    a_l.clone()),
                                        near_ma_l.clone());
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
