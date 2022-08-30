use integrals::elliptic::{ellipe, ellpk};
use once_cell::sync::Lazy;

pub const SAMPLES: usize = 100_000;

pub static E_K: Lazy<Vec<f64>> = Lazy::new(|| {
    (0..SAMPLES)
        .map(|f| ellpk(f as f64 / SAMPLES as f64))
        .collect()
});

pub static E_E: Lazy<Vec<f64>> = Lazy::new(|| {
    (0..SAMPLES)
        .map(|f| ellipe(f as f64 / SAMPLES as f64))
        .collect()
});
