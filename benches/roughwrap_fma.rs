#![cfg(feature = "use-fma")]
#![feature(test)]

extern crate test;
extern crate rand;
extern crate dfloat;
extern crate roundops;
extern crate dfloat_rops;

use rand::Rng;
use dfloat::DFloat;
use roundops::*;
use roundops::methods::RoughWrappingUnchecked;
use dfloat_rops::RWDFloatFma;

type RWDF64 = RWDFloatFma<f64, RoughWrappingUnchecked<f64>>;

#[bench]
fn bench_add(b: &mut test::Bencher) {
    let mut rng = rand::thread_rng();
    let mut v1 =  Vec::<DFloat<f64>>::new();
    let mut v2 =  Vec::<DFloat<f64>>::new();
    for _ in 0..1000000 {
        v1.push(DFloat::from_pair(rng.next_f64(),rng.next_f64()));
        v2.push(DFloat::from_pair(rng.next_f64(),rng.next_f64()));
    }

    b.iter(|| for (df1, df2) in v1.iter().zip(v2.iter()) {
               test::black_box(RWDF64::add_up(df1.clone(), df2.clone()));
           })
}

#[bench]
fn bench_sub(b: &mut test::Bencher) {
    let mut rng = rand::thread_rng();
    let mut v1 = Vec::<DFloat<f64>>::new();
    let mut v2 = Vec::<DFloat<f64>>::new();
    for _ in 0..1000000 {
        v1.push(DFloat::from_pair(rng.next_f64(),rng.next_f64()));
        v2.push(DFloat::from_pair(rng.next_f64(),rng.next_f64()));
    }

    b.iter(|| for (df1, df2) in v1.iter().zip(v2.iter()) {
               test::black_box(RWDF64::sub_up(df1.clone(), df2.clone()));
           })
}

#[bench]
fn bench_mul(b: &mut test::Bencher) {
    let mut rng = rand::thread_rng();
    let mut v1 = Vec::<DFloat<f64>>::new();
    let mut v2 = Vec::<DFloat<f64>>::new();
    for _ in 0..1000000 {
        v1.push(DFloat::from_pair(rng.next_f64(),rng.next_f64()));
        v2.push(DFloat::from_pair(rng.next_f64(),rng.next_f64()));
    }

    b.iter(|| for (df1, df2) in v1.iter().zip(v2.iter()) {
               test::black_box(RWDF64::mul_up(df1.clone(), df2.clone()));
           })
}

#[bench]
fn bench_div(b: &mut test::Bencher) {
    let mut rng = rand::thread_rng();
    let mut v1 = Vec::<DFloat<f64>>::new();
    let mut v2 = Vec::<DFloat<f64>>::new();
    for _ in 0..1000000 {
        v1.push(DFloat::from_pair(rng.next_f64(),rng.next_f64()));
        v2.push(DFloat::from_pair(rng.next_f64(),rng.next_f64()));
    }

    b.iter(|| for (df1, df2) in v1.iter().zip(v2.iter()) {
               test::black_box(RWDF64::div_up(df1.clone(), df2.clone()));
           })
}

#[bench]
fn bench_sqrt(b: &mut test::Bencher) {
    let mut rng = rand::thread_rng();
    let mut v1 = Vec::<DFloat<f64>>::new();
    for _ in 0..1000000 {
        v1.push(DFloat::from_pair(rng.next_f64(),rng.next_f64()));
    }

    b.iter(|| for df1 in v1.iter() {
               test::black_box(RWDF64::sqrt_up(df1.clone()));
           })
}
