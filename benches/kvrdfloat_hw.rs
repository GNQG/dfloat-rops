#![cfg(feature = "hwrm")]
#![feature(test)]

extern crate test;
extern crate rand;
extern crate dfloat;
extern crate roundops;
extern crate dfloat_rops;

use rand::Rng;
use dfloat::DFloat;
use roundops::*;
use roundops::methods::Hardware;
use dfloat_rops::KVRDFloat;

type KVRDF64 = KVRDFloat<f64, Hardware<f64>, rmode::Switchable>;

#[bench]
fn bench_add(b: &mut test::Bencher) {
    let mut rng = rand::thread_rng();
    let mut v1 = Vec::<DFloat<f64>>::new();
    let mut v2 = Vec::<DFloat<f64>>::new();
    for _ in 0..1000000 {
        v1.push(DFloat::from_two_components(rng.next_f64(),rng.next_f64()));
        v2.push(DFloat::from_two_components(rng.next_f64(),rng.next_f64()));
    }

    b.iter(|| for (df1, df2) in v1.iter().zip(v2.iter()) {
               test::black_box(KVRDF64::add_up(df1.clone(), df2.clone()));
           })
}

#[bench]
fn bench_sub(b: &mut test::Bencher) {
    let mut rng = rand::thread_rng();
    let mut v1 = Vec::<DFloat<f64>>::new();
    let mut v2 = Vec::<DFloat<f64>>::new();
    for _ in 0..1000000 {
        v1.push(DFloat::from_two_components(rng.next_f64(),rng.next_f64()));
        v2.push(DFloat::from_two_components(rng.next_f64(),rng.next_f64()));
    }

    b.iter(|| for (df1, df2) in v1.iter().zip(v2.iter()) {
               test::black_box(KVRDF64::sub_up(df1.clone(), df2.clone()));
           })
}

#[bench]
fn bench_mul(b: &mut test::Bencher) {
    let mut rng = rand::thread_rng();
    let mut v1 = Vec::<DFloat<f64>>::new();
    let mut v2 = Vec::<DFloat<f64>>::new();
    for _ in 0..1000000 {
        v1.push(DFloat::from_two_components(rng.next_f64(),rng.next_f64()));
        v2.push(DFloat::from_two_components(rng.next_f64(),rng.next_f64()));
    }

    b.iter(|| for (df1, df2) in v1.iter().zip(v2.iter()) {
               test::black_box(KVRDF64::mul_up(df1.clone(), df2.clone()));
           })
}

#[bench]
fn bench_div(b: &mut test::Bencher) {
    let mut rng = rand::thread_rng();
    let mut v1 = Vec::<DFloat<f64>>::new();
    let mut v2 = Vec::<DFloat<f64>>::new();
    for _ in 0..1000000 {
        v1.push(DFloat::from_two_components(rng.next_f64(),rng.next_f64()));
        v2.push(DFloat::from_two_components(rng.next_f64(),rng.next_f64()));
    }

    b.iter(|| for (df1, df2) in v1.iter().zip(v2.iter()) {
               test::black_box(KVRDF64::div_up(df1.clone(), df2.clone()));
           })
}

#[bench]
fn bench_sqrt(b: &mut test::Bencher) {
    let mut rng = rand::thread_rng();
    let mut v1 = Vec::<DFloat<f64>>::new();
    for _ in 0..1000000 {
        v1.push(DFloat::from_two_components(rng.next_f64(),rng.next_f64()));
    }

    b.iter(|| for df1 in v1.iter() {
               test::black_box(KVRDF64::sqrt_up(df1.clone()));
           })
}
