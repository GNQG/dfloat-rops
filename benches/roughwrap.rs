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
use dfloat_rops::RWDFloatRegular;

fn gen_df64(rng: &mut rand::ThreadRng) -> DFloat<f64>{
    let h = (rng.next_f64()+1.) * 2f64.powi(rng.gen_range(-25, 25)) *
                 ((rng.gen_range::<i16>(0, 2) * 2 - 1) as f64);
    let l = h * (rng.next_f64() + 1.) * 2f64.powi(-53) * ((rng.gen_range::<i16>(0, 2) * 2 - 1) as f64);
    DFloat::from_two_components(h, l)
}

type RWDF64 = RWDFloatRegular<f64, RoughWrappingUnchecked<f64>>;

#[bench]
fn bench_add(b: &mut test::Bencher) {
    let mut rng = rand::thread_rng();
    let mut v1 =  Vec::<DFloat<f64>>::new();
    let mut v2 =  Vec::<DFloat<f64>>::new();
    for _ in 0..1000 {
        v1.push(gen_df64(&mut rng));
        v2.push(gen_df64(&mut rng));
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
    for _ in 0..1000 {
        v1.push(gen_df64(&mut rng));
        v2.push(gen_df64(&mut rng));
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
    for _ in 0..1000 {
        v1.push(gen_df64(&mut rng));
        v2.push(gen_df64(&mut rng));
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
    for _ in 0..1000 {
        v1.push(gen_df64(&mut rng));
        v2.push(gen_df64(&mut rng));
    }

    b.iter(|| for (df1, df2) in v1.iter().zip(v2.iter()) {
               test::black_box(RWDF64::div_up(df1.clone(), df2.clone()));
           })
}

#[bench]
fn bench_sqrt(b: &mut test::Bencher) {
    let mut rng = rand::thread_rng();
    let mut v1 = Vec::<DFloat<f64>>::new();
    for _ in 0..1000 {
        v1.push(gen_df64(&mut rng));
    }

    b.iter(|| for df1 in v1.iter() {
               test::black_box(RWDF64::sqrt_up(df1.clone()));
           })
}
