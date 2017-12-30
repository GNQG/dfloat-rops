use core::clone::Clone;
use core::marker::PhantomData;
use dfloat::DFloat;
use float_traits::IEEE754Float;
use roundops::*;
use roundops::rmode::*;
use safeeft::safetwosum_branch as safetwosum;

#[derive(Clone)]
pub struct KVRDFloat<S: IEEE754Float + Clone, T: RoundingMethod<HostMethod = U, Num = S> + RoundOps<S> + Clone, U: NativeRoundingMode>(PhantomData<(S, T)>);

impl<S, T, U> RoundingMethod for KVRDFloat<S, T, U>
    where S: IEEE754Float + Clone,
          T: RoundingMethod<HostMethod = U, Num = S> + RoundOps<S> + Clone,
          U: NativeRoundingMode
{
    type HostMethod = U;
    type Num = DFloat<S>;
}

trait LocalTool<M: NativeRoundingMode>: RoundingMethod {
    type InternalNum;
    fn safetwoproduct_up(a: Self::InternalNum,
                         b: Self::InternalNum)
                         -> (Self::InternalNum, Self::InternalNum);
    fn safetwoproduct_down(a: Self::InternalNum,
                           b: Self::InternalNum)
                           -> (Self::InternalNum, Self::InternalNum);
    fn add_up_lower(a: (Self::InternalNum, Self::InternalNum),
                    b: (Self::InternalNum, Self::InternalNum),
                    rh: (Self::InternalNum, Self::InternalNum))
                    -> Self::InternalNum;
    fn add_down_lower(a: (Self::InternalNum, Self::InternalNum),
                      b: (Self::InternalNum, Self::InternalNum),
                      rh: (Self::InternalNum, Self::InternalNum))
                      -> Self::InternalNum;
    fn sub_up_lower(a: (Self::InternalNum, Self::InternalNum),
                    b: (Self::InternalNum, Self::InternalNum),
                    rh: (Self::InternalNum, Self::InternalNum))
                    -> Self::InternalNum;
    fn sub_down_lower(a: (Self::InternalNum, Self::InternalNum),
                      b: (Self::InternalNum, Self::InternalNum),
                      rh: (Self::InternalNum, Self::InternalNum))
                      -> Self::InternalNum;
    fn mul_up_lower(a: (Self::InternalNum, Self::InternalNum),
                    b: (Self::InternalNum, Self::InternalNum),
                    rh: (Self::InternalNum, Self::InternalNum))
                    -> Self::InternalNum;
    fn mul_down_lower(a: (Self::InternalNum, Self::InternalNum),
                      b: (Self::InternalNum, Self::InternalNum),
                      rh: (Self::InternalNum, Self::InternalNum))
                      -> Self::InternalNum;
    fn div_up_lower(a: (Self::InternalNum, Self::InternalNum),
                    b: (Self::InternalNum, Self::InternalNum),
                    rh: (Self::InternalNum, Self::InternalNum))
                    -> Self::InternalNum;
    fn div_down_lower(a: (Self::InternalNum, Self::InternalNum),
                      b: (Self::InternalNum, Self::InternalNum),
                      rh: (Self::InternalNum, Self::InternalNum))
                      -> Self::InternalNum;
    fn sqrt_up_lower(a: (Self::InternalNum, Self::InternalNum),
                     rh: (Self::InternalNum, Self::InternalNum))
                     -> Self::InternalNum;
    fn sqrt_down_lower(a: (Self::InternalNum, Self::InternalNum),
                       rh: (Self::InternalNum, Self::InternalNum))
                       -> Self::InternalNum;
}

impl<S, T> LocalTool<Switchable> for KVRDFloat<S, T, Switchable>
    where S: EditRoundingMode + IEEE754Float + Clone,
          T: RoundingMethod<HostMethod = Switchable, Num = S> + RoundOps<S> + RoundSqrt + Clone
{
    type InternalNum = S;
    #[inline]
    fn safetwoproduct_up(a: S, b: S) -> (S, S) {
        use safeeft::split;
        let prod = a.clone() * b.clone();
        let ((a1, a2), (b1, b2)) =
            if a.clone().abs() >= S::one() / (S::min_positive() / S::epsilon()) {
                (split(a * S::epsilon()), split(b / S::epsilon()))
            } else if b.clone().abs() >= (S::min_positive() / S::epsilon()) {
                (split(a / S::epsilon()), split(b * S::epsilon()))
            } else {
                (split(a), split(b))
            };
        unsafe {
            S::rmode_controler()
                .unwrap()
                .upward_session(|| {
                    let tmp = if prod.abs() > S::radix() / S::min_positive() {
                        ((a1.clone() / S::radix()) * b1.clone() - (prod.clone() / S::radix())) *
                        S::radix()
                    } else {
                        a1.clone() * b1.clone() - prod.clone()
                    };
                    (prod, (tmp + a2.clone() * b1.clone()) + a1 * b2.clone() + a2 * b2)
                })
        }
    }
    #[inline]
    fn safetwoproduct_down(a: S, b: S) -> (S, S) {
        use safeeft::split;
        let prod = a.clone() * b.clone();
        let ((a1, a2), (b1, b2)) =
            if a.clone().abs() >= S::one() / (S::min_positive() / S::epsilon()) {
                (split(a * S::epsilon()), split(b / S::epsilon()))
            } else if b.clone().abs() >= (S::min_positive() / S::epsilon()) {
                (split(a / S::epsilon()), split(b * S::epsilon()))
            } else {
                (split(a), split(b))
            };
        unsafe {
            S::rmode_controler()
                .unwrap()
                .downward_session(|| {
                    let tmp = if prod.abs() > S::radix() / S::min_positive() {
                        ((a1.clone() / S::radix()) * b1.clone() - (prod.clone() / S::radix())) *
                        S::radix()
                    } else {
                        a1.clone() * b1.clone() - prod.clone()
                    };
                    (prod, (tmp + a2.clone() * b1.clone()) + a1 * b2.clone() + a2 * b2)
                })
        }
    }
    #[inline]
    fn add_up_lower(a: (S, S), b: (S, S), tmp: (S, S)) -> S {
        let (sl, al, bl) = (tmp.1, a.1, b.1);
        unsafe {
            S::rmode_controler()
                .unwrap()
                .upward_session(|| sl + al + bl)
        }
    }
    #[inline]
    fn add_down_lower(a: (S, S), b: (S, S), tmp: (S, S)) -> S {
        let (sl, al, bl) = (tmp.1, a.1, b.1);
        unsafe {
            S::rmode_controler()
                .unwrap()
                .downward_session(|| sl + al + bl)
        }
    }
    #[inline]
    fn sub_up_lower(a: (S, S), b: (S, S), tmp: (S, S)) -> S {
        let (sl, al, bl) = (tmp.1, a.1, b.1);
        unsafe {
            S::rmode_controler()
                .unwrap()
                .upward_session(|| sl + al - bl)
        }
    }
    #[inline]
    fn sub_down_lower(a: (S, S), b: (S, S), tmp: (S, S)) -> S {
        let (sl, al, bl) = (tmp.1, a.1, b.1);
        unsafe {
            S::rmode_controler()
                .unwrap()
                .downward_session(|| sl + al - bl)
        }
    }
    #[inline]
    fn mul_up_lower(a: (S, S), b: (S, S), tmp: (S, S)) -> S {
        let ((_, ml), (ah, al), (bh, bl)) = (tmp, a, b);
        unsafe {
            S::rmode_controler()
                .unwrap()
                .upward_session(|| ml + ah * bl.clone() + al.clone() * bh + al * bl)
        }
    }
    #[inline]
    fn mul_down_lower(a: (S, S), b: (S, S), tmp: (S, S)) -> S {
        let ((_, ml), (ah, al), (bh, bl)) = (tmp, a, b);
        unsafe {
            S::rmode_controler()
                .unwrap()
                .downward_session(|| ml + ah * bl.clone() + al.clone() * bh + al * bl)
        }
    }
    #[inline]
    fn div_up_lower(a: (S, S), b: (S, S), tmp: (S, S)) -> S {
        let ((dh, _), (ah, al), (bh, bl)) = (tmp, a, b);
        let mut c = S::rmode_controler().unwrap();
        if bh > S::zero() {
            let (neg_a_h, neg_a_l) = Self::safetwoproduct_up(bh.clone(), -dh.clone());
            if neg_a_h.is_infinite() {
                // bh is not too small. -> bh / S::radix() is safe.
                let (neg_a_h, neg_a_l) = Self::safetwoproduct_up(bh.clone() / S::radix(),
                                                                 -dh.clone());
                unsafe {
                    let d_low_tmp = c.upward_then(|| {
                                                      (neg_a_h / S::radix() + ah / S::radix()) +
                                                      ((-dh) * bl.clone() / S::radix()) +
                                                      al / S::radix() +
                                                      neg_a_l / S::radix()
                                                  });
                    let tmp = if d_low_tmp >= S::zero() {
                        // down -> up
                        let r = c.downward_then(|| bh / S::radix() + bl / S::radix());
                        S::upward();
                        r
                    } else {
                        // up
                        bh / S::radix() + bl / S::radix()
                    };
                    d_low_tmp / tmp
                }
            } else {
                unsafe {
                    let d_low_tmp =
                        c.upward_then(|| (neg_a_h + ah) + ((-dh) * bl.clone()) + al + neg_a_l);
                    let tmp = if d_low_tmp >= S::zero() {
                        c.downward_session(|| bh + bl)
                    } else {
                        bh + bl
                    };
                    d_low_tmp / tmp
                }
            }
        } else {
            let (neg_a_h, neg_a_l) = Self::safetwoproduct_down(bh.clone(), -dh.clone());
            if neg_a_h.is_infinite() {
                let (neg_a_h, neg_a_l) = Self::safetwoproduct_down(bh.clone() / S::radix(),
                                                                   -dh.clone());
                unsafe {
                    let d_low_tmp = c.downward_then(|| {
                                                        (neg_a_h / S::radix() + ah / S::radix()) +
                                                        ((-dh) * bl.clone() / S::radix()) +
                                                        al / S::radix() +
                                                        neg_a_l / S::radix()
                                                    });
                    let tmp = if d_low_tmp >= S::zero() {
                        // down -> up
                        let t = bh / S::radix() + bl / S::radix();
                        S::upward();
                        t
                    } else {
                        // up
                        c.upward_then(|| bh / S::radix() + bl / S::radix())
                    };
                    // up
                    d_low_tmp / tmp
                }
            } else {
                unsafe {
                    let d_low_tmp =
                        c.downward_then(|| (neg_a_h + ah) + ((-dh) * bl.clone()) + al + neg_a_l);
                    let tmp = if d_low_tmp >= S::zero() {
                        // down -> up
                        let t = bh + bl;
                        S::upward();
                        t
                    } else {
                        // up
                        c.upward_then(|| bh + bl)
                    };
                    // up
                    d_low_tmp / tmp
                }
            }
        }
    }
    #[inline]
    fn div_down_lower(a: (S, S), b: (S, S), tmp: (S, S)) -> S {
        let ((dh, _), (ah, al), (bh, bl)) = (tmp, a, b);
        let mut c = S::rmode_controler().unwrap();
        if bh > S::zero() {
            let (neg_a_h, neg_a_l) = Self::safetwoproduct_down(bh.clone(), -dh.clone());
            if neg_a_h.is_infinite() {
                // bh is not too small. -> bh / S::radix() is safe.
                let (neg_a_h, neg_a_l) = Self::safetwoproduct_down(bh.clone() / S::radix(),
                                                                   -dh.clone());
                unsafe {
                    let d_low_tmp = c.downward_then(|| {
                                                        (neg_a_h / S::radix() + ah / S::radix()) +
                                                        ((-dh) * bl.clone() / S::radix()) +
                                                        al / S::radix() +
                                                        neg_a_l / S::radix()
                                                    });
                    let tmp = if d_low_tmp >= S::zero() {
                        // up -> down
                        let r = c.upward_then(|| bh / S::radix() + bl / S::radix());
                        S::downward();
                        r
                    } else {
                        // down
                        bh / S::radix() + bl / S::radix()
                    };
                    // down
                    d_low_tmp / tmp
                }
            } else {
                unsafe {
                    let d_low_tmp =
                        c.downward_then(|| (neg_a_h + ah) + ((-dh) * bl.clone()) + al + neg_a_l);
                    let tmp = if d_low_tmp >= S::zero() {
                        let r = c.upward_then(|| bh + bl);
                        S::downward();
                        r
                    } else {
                        bh + bl
                    };
                    d_low_tmp / tmp
                }
            }
        } else {
            let (neg_a_h, neg_a_l) = Self::safetwoproduct_up(bh.clone(), -dh.clone());
            if neg_a_h.is_infinite() {
                let (neg_a_h, neg_a_l) = Self::safetwoproduct_up(bh.clone() / S::radix(),
                                                                 -dh.clone());
                unsafe {
                    let d_low_tmp = c.upward_then(|| {
                                                      (neg_a_h / S::radix() + ah / S::radix()) +
                                                      ((-dh) * bl.clone() / S::radix()) +
                                                      al / S::radix() +
                                                      neg_a_l / S::radix()
                                                  });
                    let tmp = if d_low_tmp >= S::zero() {
                        // up -> down
                        let t = bh / S::radix() + bl / S::radix();
                        S::downward();
                        t
                    } else {
                        // down
                        c.downward_then(|| bh / S::radix() + bl / S::radix())
                    };
                    // down
                    d_low_tmp / tmp
                }
            } else {
                unsafe {
                    let d_low_tmp =
                        c.upward_then(|| (neg_a_h + ah) + ((-dh) * bl.clone()) + al + neg_a_l);
                    let tmp = if d_low_tmp >= S::zero() {
                        // up -> down
                        let t = bh + bl;
                        S::downward();
                        t
                    } else {
                        // down
                        c.downward_then(|| bh + bl)
                    };
                    // down
                    d_low_tmp / tmp
                }
            }
        }
    }
    #[inline]
    fn sqrt_up_lower(a: (S, S), tmp: (S, S)) -> S {
        let ((rh, _), (ah, al)) = (tmp, a);
        let mut c = S::rmode_controler().unwrap();
        let (neg_a_h, neg_a_l) = Self::safetwoproduct_up(-rh.clone(), rh.clone());
        unsafe {
            let rl_tmp = c.upward_then(|| neg_a_h + ah.clone() + al.clone() + neg_a_l);
            let tmp = if rl_tmp > S::zero() {
                // down -> up
                let r = c.downward_then(|| (ah + al).sqrt() + rh);
                S::upward();
                r
            } else {
                // up
                (ah + al).sqrt() + rh
            };
            rl_tmp / tmp
        }
    }
    #[inline]
    fn sqrt_down_lower(a: (S, S), tmp: (S, S)) -> S {
        let ((rh, _), (ah, al)) = (tmp, a);
        let mut c = S::rmode_controler().unwrap();
        let (neg_a_h, neg_a_l) = Self::safetwoproduct_down(-rh.clone(), rh.clone());
        unsafe {
            let rl_tmp = c.downward_then(|| neg_a_h + ah.clone() + al.clone() + neg_a_l);
            let tmp = if rl_tmp > S::zero() {
                // up -> down
                let r = c.upward_then(|| (ah + al).sqrt() + rh);
                S::downward();
                r
            } else {
                // down
                (ah + al).sqrt() + rh
            };
            rl_tmp / tmp
        }
    }
}

impl<S, T> LocalTool<DefaultRounding> for KVRDFloat<S, T, DefaultRounding>
    where S: IEEE754Float + Clone,
          T: RoundingMethod<HostMethod = DefaultRounding, Num = S> + RoundOps<S> + RoundSqrt + Clone
{
    type InternalNum = S;
    #[inline]
    fn safetwoproduct_up(a: S, b: S) -> (S, S) {
        use safeeft::split;
        let prod = a.clone() * b.clone();
        let ((a1, a2), (b1, b2)) =
            if a.clone().abs() >= S::one() / (S::min_positive() / S::epsilon()) {
                (split(a * S::epsilon()), split(b / S::epsilon()))
            } else if b.clone().abs() >= (S::min_positive() / S::epsilon()) {
                (split(a / S::epsilon()), split(b * S::epsilon()))
            } else {
                (split(a), split(b))
            };
        let tmp = if prod.clone().abs() > S::radix() / S::min_positive() {
            let (wprod, wa1, wb1) =
                rnum_init!(<direction::Upward, S, T>,
                (prod.clone() / S::radix(), a1.clone() / S::radix(), b1.clone()));
            (wa1 * wb1 - wprod).extract() * S::radix()
        } else {
            let (wprod, wa1, wb1) = rnum_init!(<direction::Upward, S, T>,
                (prod.clone(), a1.clone(), b1.clone()));
            (wa1 * wb1 - wprod).extract()
        };
        let (wtmp, wa1, wa2, wb1, wb2) = rnum_init!(<direction::Upward, S, T>,
                (tmp, a1, a2, b1, b2));
        (prod, (wtmp + wa2.clone() * wb1 + wa1 * wb2.clone() + wa2 * wb2).extract())
    }
    #[inline]
    fn safetwoproduct_down(a: S, b: S) -> (S, S) {
        use safeeft::split;
        let prod = a.clone() * b.clone();
        let ((a1, a2), (b1, b2)) =
            if a.clone().abs() >= S::one() / (S::min_positive() / S::epsilon()) {
                (split(a * S::epsilon()), split(b / S::epsilon()))
            } else if b.clone().abs() >= (S::min_positive() / S::epsilon()) {
                (split(a / S::epsilon()), split(b * S::epsilon()))
            } else {
                (split(a), split(b))
            };
        let tmp = if prod.clone().abs() > S::radix() / S::min_positive() {
            let (wprod, wa1, wb1) =
                rnum_init!(<direction::Downward, S, T>,
                (prod.clone() / S::radix(), a1.clone() / S::radix(), b1.clone()));
            (wa1 * wb1 - wprod).extract() * S::radix()
        } else {
            let (wprod, wa1, wb1) = rnum_init!(<direction::Downward, S, T>,
                (prod.clone(), a1.clone(), b1.clone()));
            (wa1 * wb1 - wprod).extract()
        };
        let (wtmp, wa1, wa2, wb1, wb2) = rnum_init!(<direction::Downward, S, T>,
                (tmp, a1, a2, b1, b2));
        (prod, (wtmp + wa2.clone() * wb1 + wa1 * wb2.clone() + wa2 * wb2).extract())
    }
    #[inline]
    fn add_up_lower(a: (S, S), b: (S, S), tmp: (S, S)) -> S {
        let (wsl, wal, wbl) = rnum_init!(<direction::Upward, S, T>, (tmp.1, a.1, b.1));
        (wsl + wal + wbl).extract()
    }
    #[inline]
    fn add_down_lower(a: (S, S), b: (S, S), tmp: (S, S)) -> S {
        let (wsl, wal, wbl) = rnum_init!(<direction::Downward, S, T>, (tmp.1, a.1, b.1));
        (wsl + wal + wbl).extract()
    }
    #[inline]
    fn sub_up_lower(a: (S, S), b: (S, S), tmp: (S, S)) -> S {
        let (wsl, wal, wbl) = rnum_init!(<direction::Upward, S, T>, (tmp.1, a.1, b.1));
        (wsl + wal - wbl).extract()
    }
    #[inline]
    fn sub_down_lower(a: (S, S), b: (S, S), tmp: (S, S)) -> S {
        let (wsl, wal, wbl) = rnum_init!(<direction::Downward, S, T>, (tmp.1, a.1, b.1));
        (wsl + wal - wbl).extract()
    }
    #[inline]
    fn mul_up_lower(a: (S, S), b: (S, S), tmp: (S, S)) -> S {
        let (wml, wah, wal, wbh, wbl) =
            rnum_init!(<direction::Upward, S, T>, (tmp.1, a.0, a.1, b.0, b.1));
        (wml + wal.clone() * wbh + wah * wbl.clone() + wal * wbl).extract()
    }
    #[inline]
    fn mul_down_lower(a: (S, S), b: (S, S), tmp: (S, S)) -> S {
        let (wml, wah, wal, wbh, wbl) =
            rnum_init!(<direction::Downward, S, T>, (tmp.1, a.0, a.1, b.0, b.1));
        (wml + wal.clone() * wbh + wah * wbl.clone() + wal * wbl).extract()
    }
    #[inline]
    fn div_up_lower(a: (S, S), b: (S, S), tmp: (S, S)) -> S {
        let ((dh, _), (ah, al), (bh, bl)) = (tmp, a, b);
        if bh > S::zero() {
            let (neg_a_h, neg_a_l) = Self::safetwoproduct_up(bh.clone(), -dh.clone());
            if neg_a_h.is_infinite() {
                // bh is not too small. -> bh / S::radix() is safe.
                let (neg_a_h, neg_a_l) = Self::safetwoproduct_up(bh.clone() / S::radix(),
                                                                 -dh.clone());
                let (w_ah2, w_al2, w_bl2, w_neg_ah2, w_neg_al2, w_dh) =
                    rnum_init!(<direction::Upward, S, T>,
                            (ah / S::radix(), al / S::radix(), bl.clone() / S::radix(),
                             neg_a_h, neg_a_l, dh.clone()));
                let d_low_tmp = ((w_neg_ah2 + w_ah2) + ((-w_dh) * w_bl2) + w_al2 + w_neg_al2)
                    .extract();
                let tmp = if d_low_tmp >= S::zero() {
                    T::add_down(bh / S::radix(), bl / S::radix())
                } else {
                    //                                     ^ provably unsafe
                    T::add_up(bh / S::radix(), bl / S::radix())
                }; //                                      ^ provably unsafe
                T::div_up(d_low_tmp, tmp)
            } else {
                let (w_ah, w_al, w_bl, w_neg_ah, w_neg_al, w_dh) =
                    rnum_init!(<direction::Upward, S, T>,
                            (ah, al, bl.clone(), neg_a_h, neg_a_l, dh.clone()));
                let d_low_tmp = ((w_neg_ah + w_ah) + ((-w_dh) * w_bl) + w_al + w_neg_al).extract();
                let tmp = if d_low_tmp >= S::zero() {
                    T::add_down(bh, bl)
                } else {
                    T::add_up(bh, bl)
                };
                T::div_up(d_low_tmp, tmp)
            }
        } else {
            let (neg_a_h, neg_a_l) = Self::safetwoproduct_down(bh.clone(), -dh.clone());
            if neg_a_h.is_infinite() {
                let (neg_a_h, neg_a_l) = Self::safetwoproduct_down(bh.clone() / S::radix(),
                                                                   -dh.clone());
                let (w_ah2, w_al2, w_bl2, w_neg_ah2, w_neg_al2, w_dh) =
                    rnum_init!(<direction::Downward, S, T>,
                            (ah / S::radix(), al / S::radix(), bl.clone() / S::radix(),
                             neg_a_h, neg_a_l, dh.clone()));
                let d_low_tmp = ((w_neg_ah2 + w_ah2) + ((-w_dh) * w_bl2) + w_al2 + w_neg_al2)
                    .extract();
                let tmp = if d_low_tmp >= S::zero() {
                    T::add_down(bh / S::radix(), bl / S::radix())
                } else {
                    T::add_up(bh / S::radix(), bl / S::radix())
                };
                T::div_up(d_low_tmp, tmp)
            } else {
                let (w_ah, w_al, w_bl, w_neg_ah, w_neg_al, w_dh) =
                    rnum_init!(<direction::Downward, S, T>,
                            (ah, al, bl.clone(), neg_a_h, neg_a_l, dh.clone()));
                let d_low_tmp = ((w_neg_ah + w_ah) + ((-w_dh) * w_bl) + w_al + w_neg_al).extract();
                let tmp = if d_low_tmp >= S::zero() {
                    T::add_down(bh, bl)
                } else {
                    T::add_up(bh, bl)
                };
                T::div_up(d_low_tmp, tmp)
            }
        }
    }
    #[inline]
    fn div_down_lower(a: (S, S), b: (S, S), tmp: (S, S)) -> S {
        let ((dh, _), (ah, al), (bh, bl)) = (tmp, a, b);
        if bh > S::zero() {
            let (neg_a_h, neg_a_l) = Self::safetwoproduct_down(bh.clone(), -dh.clone());
            if neg_a_h.is_infinite() {
                // bh is not too small. -> bh / S::radix() is safe.
                let (neg_a_h, neg_a_l) = Self::safetwoproduct_down(bh.clone() / S::radix(),
                                                                   -dh.clone());
                let (w_ah2, w_al2, w_bl2, w_neg_ah2, w_neg_al2, w_dh) =
                    rnum_init!(<direction::Downward, S, T>,
                        (ah / S::radix(), al / S::radix(), bl.clone() / S::radix(),
                         neg_a_h, neg_a_l, dh.clone()));
                let d_low_tmp = ((w_neg_ah2 + w_ah2) + ((-w_dh) * w_bl2) + w_al2 + w_neg_al2)
                    .extract();
                let tmp = if d_low_tmp >= S::zero() {
                    T::add_up(bh / S::radix(), bl / S::radix())
                } else {
                    //                                     ^ provably unsafe
                    T::add_down(bh / S::radix(), bl / S::radix())
                }; //                                      ^ provably unsafe
                T::div_down(d_low_tmp, tmp)
            } else {
                let (w_ah, w_al, w_bl, w_neg_ah, w_neg_al, w_dh) =
                    rnum_init!(<direction::Downward, S, T>,
                        (ah,al, bl.clone(), neg_a_h, neg_a_l, dh.clone()));
                let d_low_tmp = ((w_neg_ah + w_ah) + ((-w_dh) * w_bl) + w_al + w_neg_al).extract();
                let tmp = if d_low_tmp >= S::zero() {
                    T::add_up(bh, bl)
                } else {
                    T::add_down(bh, bl)
                };
                T::div_down(d_low_tmp, tmp)
            }
        } else {
            let (neg_a_h, neg_a_l) = Self::safetwoproduct_up(bh.clone(), -dh.clone());
            if neg_a_h.is_infinite() {
                let (neg_a_h, neg_a_l) = Self::safetwoproduct_up(bh.clone() / S::radix(),
                                                                 -dh.clone());
                let (w_ah2, w_al2, w_bl2, w_neg_ah2, w_neg_al2, w_dh) =
                    rnum_init!(<direction::Upward, S, T>,
                            (ah / S::radix(), al / S::radix(), bl.clone() / S::radix(),
                             neg_a_h, neg_a_l, dh.clone()));
                let d_low_tmp = ((w_neg_ah2 + w_ah2) + ((-w_dh) * w_bl2) + w_al2 + w_neg_al2)
                    .extract();
                let tmp = if d_low_tmp >= S::zero() {
                    T::add_up(bh / S::radix(), bl / S::radix())
                } else {
                    T::add_down(bh / S::radix(), bl / S::radix())
                };
                T::div_down(d_low_tmp, tmp)
            } else {
                let (w_ah, w_al, w_bl, w_neg_ah, w_neg_al, w_dh) =
                    rnum_init!(<direction::Upward, S, T>,
                            (ah, al, bl.clone(), neg_a_h, neg_a_l, dh.clone()));
                let d_low_tmp = ((w_neg_ah + w_ah) + ((-w_dh) * w_bl) + w_al + w_neg_al).extract();
                let tmp = if d_low_tmp >= S::zero() {
                    T::add_up(bh, bl)
                } else {
                    T::add_down(bh, bl)
                };
                T::div_down(d_low_tmp, tmp)
            }
        }
    }
    #[inline]
    fn sqrt_up_lower(a: (S, S), tmp: (S, S)) -> S {
        let ((rh, _), (ah, al)) = (tmp, a);
        let (neg_a_h, neg_a_l) = Self::safetwoproduct_up(-rh.clone(), rh.clone());
        let (w_ah, w_al, w_neg_ah, w_neg_al) =
            rnum_init!(<direction::Upward,S,T>,
                    (ah.clone(), al.clone(), neg_a_h, neg_a_l));
        let r_low_tmp = (w_neg_ah + w_ah + w_al + w_neg_al).extract();
        let tmp = if r_low_tmp > S::zero() {
            let (w_ah, w_al, w_rh) = rnum_init!(<direction::Downward, S, T>, (ah, al, rh.clone()));
            ((w_ah + w_al).sqrt() + w_rh).extract()
        } else {
            let (w_ah, w_al, w_rh) = rnum_init!(<direction::Upward, S, T>, (ah, al, rh.clone()));
            ((w_ah + w_al).sqrt() + w_rh).extract()
        };
        T::div_up(r_low_tmp, tmp)
    }
    #[inline]
    fn sqrt_down_lower(a: (S, S), tmp: (S, S)) -> S {
        let ((rh, _), (ah, al)) = (tmp, a);
        let (neg_a_h, neg_a_l) = Self::safetwoproduct_down(-rh.clone(), rh.clone());
        let (w_ah, w_al, w_neg_ah, w_neg_al) =
            rnum_init!(<direction::Downward, S, T>, 
                    (ah.clone(), al.clone(), neg_a_h, neg_a_l));
        let r_low_tmp = (w_neg_ah + w_ah + w_al + w_neg_al).extract();
        let tmp = if r_low_tmp > S::zero() {
            let (w_ah, w_al, w_rh) = rnum_init!(<direction::Upward, S, T>, (ah, al, rh.clone()));
            ((w_ah + w_al).sqrt() + w_rh).extract()
        } else {
            let (w_ah, w_al, w_rh) = rnum_init!(<direction::Downward, S, T>, (ah, al, rh.clone()));
            ((w_ah + w_al).sqrt() + w_rh).extract()
        };
        T::div_down(r_low_tmp, tmp)
    }
}

macro_rules! impl_for_native {
    ($U:ident $(,$apdbound:ident)*) => (
impl<S, T> RoundAdd for KVRDFloat<S, T, $U> 
    where S: $($apdbound+)* IEEE754Float + Clone,
          T: RoundingMethod<HostMethod = $U, Num = S> + RoundOps<S> + RoundSqrt + Clone
{
    fn add_up(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((ah, al), (bh, bl)) = (a.decomposite(), b.decomposite());
        let (sh, sl) = safetwosum(ah.clone(), bh.clone());
        if sh.is_infinite() {
            if sh == S::infinity() {
                DFloat::infinity()
            } else {
                DFloat::min_value()
            }
        } else {
            let (sh, sl) = safetwosum(sh.clone(),
                                      Self::add_up_lower((ah, al), (bh, bl), (sh, sl)));
            DFloat::_from_pair_raw(sh, sl)
        }
    }
    fn add_down(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((ah, al), (bh, bl)) = (a.decomposite(), b.decomposite());
        let (sh, sl) = safetwosum(ah.clone(), bh.clone());
        if sh.is_infinite() {
            if sh == S::neg_infinity() {
                DFloat::neg_infinity()
            } else {
                DFloat::max_value()
            }
        } else {
            let (sh, sl) = safetwosum(sh.clone(),
                                      Self::add_down_lower((ah, al), (bh, bl), (sh, sl)));
            DFloat::_from_pair_raw(sh, sl)
        }
    }
}

impl<S, T> RoundSub for KVRDFloat<S, T, $U> 
    where S: $($apdbound+)* IEEE754Float + Clone,
          T: RoundingMethod<HostMethod = $U, Num = S> + RoundOps<S> + RoundSqrt + Clone
{
    fn sub_up(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((ah, al), (bh, bl)) = (a.decomposite(), b.decomposite());
        let (sh, sl) = safetwosum(ah.clone(), -bh.clone());
        if sh.is_infinite() {
            if sh == S::infinity() {
                DFloat::infinity()
            } else {
                DFloat::min_value()
            }
        } else {
            let (sh, sl) = safetwosum(sh.clone(),
                                      Self::sub_up_lower((ah, al), (bh, bl), (sh, sl)));
            DFloat::_from_pair_raw(sh, sl)
        }
    }
    fn sub_down(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((ah, al), (bh, bl)) = (a.decomposite(), b.decomposite());
        let (sh, sl) = safetwosum(ah.clone(), -bh.clone());
        if sh.is_infinite() {
            if sh == S::neg_infinity() {
                DFloat::neg_infinity()
            } else {
                DFloat::max_value()
            }
        } else {
            let (sh, sl) = safetwosum(sh.clone(),
                                      Self::sub_down_lower((ah, al), (bh, bl), (sh, sl)));
            DFloat::_from_pair_raw(sh, sl)
        }
    }
}

impl<S, T> RoundMul for KVRDFloat<S, T, $U>
    where S: $($apdbound+)* IEEE754Float + Clone,
          T: RoundingMethod<HostMethod = $U, Num = S> + RoundOps<S> + RoundSqrt + Clone
{
    fn mul_up(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((ah, al), (bh, bl)) = (a.decomposite(), b.decomposite());
        let (mh, ml) = Self::safetwoproduct_up(ah.clone(), bh.clone());
        if mh.is_infinite() {
            if mh == S::infinity() {
                DFloat::infinity()
            } else {
                DFloat::min_value()
            }
        } else {
            let (mh, ml) = safetwosum(mh.clone(),
                                      Self::mul_up_lower((ah, al), (bh, bl), (mh, ml)));
            DFloat::_from_pair_raw(mh, ml)
        }
    }

    fn mul_down(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((ah, al), (bh, bl)) = (a.decomposite(), b.decomposite());
        let (mh, ml) = Self::safetwoproduct_down(ah.clone(), bh.clone());
        if mh.is_infinite() {
            if mh == S::neg_infinity() {
                DFloat::neg_infinity()
            } else {
                DFloat::max_value()
            }
        } else {
            let (mh, ml) = safetwosum(mh.clone(),
                                      Self::mul_down_lower((ah, al), (bh, bl), (mh, ml)));
            DFloat::_from_pair_raw(mh, ml)
        }
    }
}

impl<S, T> RoundDiv for KVRDFloat<S, T, $U>
    where S: $($apdbound+)* IEEE754Float + Clone,
          T: RoundingMethod<HostMethod = $U, Num = S> + RoundOps<S> + RoundSqrt + Clone
{
    fn div_up(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((ah, al), (bh, bl)) = (a.decomposite(), b.decomposite());
        let dh = ah.clone() / bh.clone();
        if dh.is_infinite() {
            if dh == S::infinity() {
                DFloat::infinity()
            } else {
                DFloat::min_value()
            }
        } else if bh.is_infinite() {
            DFloat::from_single(dh)
        } else {
            let d = safetwosum(dh.clone(),
                               Self::div_down_lower((ah, al), (bh, bl), (dh, S::zero())));
            DFloat::_from_pair_raw(d.0, d.1)
        }
    }
    fn div_down(a: DFloat<S>, b: DFloat<S>) -> DFloat<S> {
        let ((ah, al), (bh, bl)) = (a.decomposite(), b.decomposite());
        let dh = ah.clone() / bh.clone();
        if dh.is_infinite() {
            if dh == S::infinity() {
                DFloat::infinity()
            } else {
                DFloat::min_value()
            }
        } else if bh.is_infinite() {
            DFloat::from_single(dh)
        } else {
            let d = safetwosum(dh.clone(),
                               Self::div_down_lower((ah, al), (bh, bl), (dh, S::zero())));
            DFloat::_from_pair_raw(d.0, d.1)
        }
    }
}

impl<S, T> RoundSqrt for KVRDFloat<S, T, $U>
    where S: $($apdbound+)* IEEE754Float + Clone,
          T: RoundingMethod<HostMethod = $U, Num = S> + RoundOps<S> + RoundSqrt + Clone
{
    fn sqrt_up(a: DFloat<S>) -> DFloat<S> {
        let (ah, al) = a.decomposite();
        if ah == S::infinity() {
            DFloat::infinity()
        } else if ah == S::zero() {
            DFloat::zero()
        } else {
            let rh = ah.clone().sqrt();
            let r = safetwosum(rh.clone(),
                               Self::sqrt_up_lower((ah, al), (rh, S::zero())));
            DFloat::_from_pair_raw(r.0, r.1)
        }
    }
    fn sqrt_down(a: DFloat<S>) -> DFloat<S> {
        let (ah, al) = a.decomposite();
        if ah == S::infinity() {
            DFloat::infinity()
        } else if ah == S::zero() {
            DFloat::zero()
        } else {
            let rh = ah.clone().sqrt();
            let r = safetwosum(rh.clone(),
                               Self::sqrt_down_lower((ah, al), (rh, S::zero())));
            DFloat::_from_pair_raw(r.0, r.1)
        }
    }
}
)}

impl_for_native!(Switchable, EditRoundingMode);
impl_for_native!(DefaultRounding);