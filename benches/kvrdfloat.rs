#![feature(test)]

extern crate test;
extern crate dfloat;
extern crate roundops;
extern crate dfloat_rops;

use dfloat::dfloat::DFloat;
use roundops::*;
use roundops::methods::Emulation;
use dfloat_rops::KVRDFloatRegular;

type KVRDF64 = KVRDFloatRegular<f64, Emulation<f64>>;

#[bench]
fn bench_div(b: &mut test::Bencher) {
    let mut v1 = vec![DFloat::<f64>::one()];
    v1.push(DFloat::from_pair(234.65, 231.7));
    v1.push(DFloat::from_pair(-64.65, -0.14357));
    v1.push(DFloat::from_pair(23.65, 234631.7));
    v1.push(DFloat::from_pair(-2.644525, 1.7));
    v1.push(DFloat::from_pair(145e-20, 0.));
    let mut v2 = vec![DFloat::from_pair(234.65, 231.7)];
    v2.push(DFloat::from_pair(-64.65, -0.14357));
    v2.push(DFloat::from_pair(23.65, 234631.7));
    v2.push(DFloat::from_pair(-2.644525, 1.7));
    v2.push(DFloat::from_pair(145e-20, 0.));
    v2.push(DFloat::<f64>::one());

    b.iter(|| for (df1, df2) in v1.iter().zip(v2.iter()) {
               test::black_box(KVRDF64::div_up(df1.clone(), df2.clone()));
           })
}
