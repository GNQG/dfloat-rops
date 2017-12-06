#[cfg_attr(feature = "use-fma", feature(cfg_target_feature))]
extern crate core;
extern crate dfloat;
extern crate float_traits;
extern crate safeeft;
extern crate roundops;
#[cfg(feature = "use-fma")]
extern crate fma;

mod kvrdfloat;
pub use kvrdfloat::KVRDFloatRegular;
mod roughwrap;
pub use roughwrap::RWDFloatRegular;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
