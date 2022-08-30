use std::f64::consts::PI;

use nalgebra::Vector3;

use crate::{tables, Precision};

const PERMEABILITY: f64 = 4e-7 * PI;

#[derive(Debug, Clone, Copy)]
pub struct Field {
    pub radius: Precision,
    pub position: Vector3<Precision>,
    pub current: Precision,
}

impl Field {
    // https://tiggerntatie.github.io/emagnet-py/offaxis/off_axis_loop.html
    pub fn b_field(&self, particle_pos: &Vector3<Precision>) -> Vector3<Precision> {
        // todo: adjust with orientation
        let a = self.radius;
        let diff_vector = particle_pos - self.position;
        let x = diff_vector.z;
        let r = Vector3::new(diff_vector.x, diff_vector.y, 0.0).magnitude();
        let i = self.current;

        let alpha = r / a;
        let beta = x / a;
        let gamma = x / r;

        let q = (1.0 + alpha).powi(2) + beta.powi(2);
        let k = (4.0 * alpha / q).sqrt();
        let b_0 = (i * PERMEABILITY) / (2.0 * a);

        let index = (k * tables::SAMPLES as f64) as usize;
        let ek = tables::E_E[index];
        let kk = tables::E_K[index];

        let b_x = b_0
            * (1.0 / (PI * q.sqrt()))
            * ((ek * (1.0 - alpha.powi(2) - beta.powi(2)) / (q - 4.0 * alpha)) + kk);

        let b_r = b_0
            * (gamma / (PI * q.sqrt()))
            * ((ek * (1.0 + alpha.powi(2) + beta.powi(2)) / (q - 4.0 * alpha)) - kk);

        Vector3::new(
            particle_pos.x / self.radius * b_r,
            particle_pos.y / self.radius * b_r,
            b_x,
        )
    }
}

// https://tiggerntatie.github.io/emagnet-py/solenoids/current_loop.html
#[test]
fn on_axis_field_works() {
    let field = Field {
        radius: 1.0,
        position: Vector3::new(0.0, 0.0, 0.0),
        current: 1.0,
    };

    use approx::assert_relative_eq as apeq;

    let vals = [
        (0.0, 1.0),
        (0.5555555555555556, 0.6679881072249486),
        (1.1111111111111112, 0.2993709572454057),
        (1.6666666666666667, 0.13619005290728642),
        (2.2222222222222223, 0.06910507040256261),
        (2.7777777777777777, 0.03886158205395477),
        (3.333333333333333, 0.023725972202725775),
        (3.8888888888888884, 0.015445883786110361),
        (4.444444444444444, 0.010577327665007765),
        (4.999999999999999, 0.007542928274545547),
    ];

    for (x, y) in vals {
        let b = field.b_field(&Vector3::new(0.0, 0.0, x));
        apeq!(y, b.z / 6.283185307179586e-7);
    }

    // for z in crate::util::arange(0.0, 5.0, 10) {
    //     let b = field.b_field(&Vector3::new(0.0, 0.0, z));
    //     println!("{} {}", z, b.z / 6.283185307179586e-7);
    // }
}

// const PERMEABILITY: f64 = 1.25663706e-6;

// let b_r = particle_pos.z * (particle_pos.z - self.position.x);

// complete elliptic integral of the first kind
// let k = (int)(4.0f * a / (Q*Q) * 10000000);
// let k = 4.0 * alpha / q.powi(2) * 10_000_000.0;
// Q = sqrt((1.0f + a) * (1.0f + a) + B*B);
// q = sqrt((1 + a)^2 + b^2)

// pz * ((pz - px)/r)

// let B_r = self.position.z * ((particle_pos.z - coils[i].x) / r) / (PI * Q)
//     * ((E_k * (1.0 + alpha * alpha + B * B) / ((Q * Q) - 4.0 * alpha)) - K_k);

// let B_r = self.position.z * ((particle_pos.z - coils[i].x) / r) / (PI * Q)
//     * ((E_k * (1.0 + alpha * alpha + B * B) / ((Q * Q) - 4.0 * alpha)) - K_k);

// b.x += pos.x / r * B_r;
// b.y += pos.y / r * B_r;
// b.z += coils[i].z * (1.0f)/( M_PI * Q) * (( E_k * (1.0f - a*a - B*B)/((Q*Q)-4.0f*a)) + K_k);
