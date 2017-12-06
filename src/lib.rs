#![cfg_attr(feature = "use-fma", feature(cfg_target_feature))]
extern crate core;
extern crate dfloat;
extern crate float_traits;
extern crate safeeft;
extern crate roundops;
#[cfg(feature = "use-fma")]
extern crate fma;

mod kvrdfloat;
pub use kvrdfloat::KVRDFloatRegular;

#[cfg(feature = "use-fma")]
mod kvrdfloat_fma;
#[cfg(feature = "use-fma")]
pub use kvrdfloat_fma::KVRDFloatFMA;

#[cfg(not(feature = "use-fma"))]
pub use KVRDFloatRegular as KVRDFloat;
#[cfg(feature = "use-fma")]
pub use KVRDFloatFMA as KVRDFloat;

mod roughwrap;
pub use roughwrap::RWDFloatRegular;

#[cfg(feature = "use-fma")]
mod roughwrap_fma;
#[cfg(feature = "use-fma")]
pub use roughwrap_fma::RWDFloatFMA;

#[cfg(not(feature = "use-fma"))]
pub use RWDFloatRegular as RWDFloat;
#[cfg(feature = "use-fma")]
pub use RWDFloatFMA as RWDFloat;


#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
