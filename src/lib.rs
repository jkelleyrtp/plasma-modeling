use std::f64::consts::PI;

mod integrals;
mod polynomial;
use nalgebra::Vector3;

pub struct Config {
    pub dt: Precision,
    pub steps: usize,
    pub num_particles: usize,
    pub micro_step: usize,
}

type Precision = f64;

#[derive(Clone)]
pub struct Electron {
    pub position: Vector3<Precision>,
    pub velocity: Vector3<Precision>,
}

impl Electron {
    fn new(energy: Precision) -> Self {
        let mut rng = rand::thread_rng();

        let mut position = Vector3::default();
        let mut velocity = Vector3::default();

        Electron { position, velocity }
    }
}

pub struct Field {
    pub radius: Precision,
    pub position: Vector3<Precision>,
    pub orientation: Vector3<Precision>,
}

pub struct Simulation {
    /// Scale up calculations with tiny numbers to prevent roundoff errors
    /// Important to do it when things are added to the simulation
    pub scale_factor: Precision,

    pub electron_mass: Precision,
    pub electron_charge: Precision,

    pub fields: Vec<Field>,
    pub particles: Vec<Electron>,

    pub history: Vec<Vec<Electron>>,
}

impl Simulation {
    pub fn new() -> Self {
        let scale_factor = 1e20;

        let electron_mass = 9.10938291e-31;
        let electron_charge = 1.6021766208e-19;

        Self {
            scale_factor,
            electron_mass,
            electron_charge,
            fields: Vec::new(),
            particles: Vec::new(),
            history: Vec::new(),
        }
    }

    // Rescale everything in the simulation to prevent float errors
    pub fn rescale(&mut self, scale_factor: Precision) {
        self.scale_factor = scale_factor;
        self.electron_charge *= scale_factor;
        self.electron_mass *= scale_factor;
    }

    pub fn populate(&mut self) {}

    pub fn tick_adpative(&mut self, dt: Precision) {
        //
    }

    pub fn tick_rk4(&mut self, dt: Precision) {
        for x in 0..self.particles.len() {
            let particle = &self.particles[x];

            let k1 = dt * self.net_accel(particle, Vector3::default());
            let l1 = dt * particle.velocity;

            let k2 = dt * self.net_accel(particle, 0.5 * l1);
            let l2 = dt * particle.velocity + (0.5 * k1);

            let k3 = dt * self.net_accel(particle, 0.5 * l2);
            let l3 = dt * particle.velocity + (0.5 * k2);

            let k4 = self.net_accel(particle, l3);
            let l4 = dt * (particle.velocity + k3);

            self.particles[x].velocity += (k1 + (2.0 * k2) + (2.0 * k3) + k4) / 6.0;
            self.particles[x].position += (l1 + (2.0 * l2) + (2.0 * l3) + l4) / 6.0;
        }
    }

    pub fn tick_eulers(&mut self, dt: Precision) {
        for x in 0..self.particles.len() {
            let acceleration = self.net_accel(&self.particles[x], Vector3::default());
            let velocity = self.particles[x].velocity;

            // v = v + a*dt
            self.particles[x].velocity += acceleration * dt;

            // p = p + v*dt
            self.particles[x].position += velocity * dt;
        }
    }

    fn net_accel(&self, particle: &Electron, offset: Vector3<Precision>) -> Vector3<Precision> {
        let position = particle.position + offset;

        // net b
        let net_b = self
            .fields
            .iter()
            .fold(Vector3::new(0.0, 0.0, 0.0), |acc, field| {
                acc + field.b_field(&position)
            });

        // f = qv x b
        let force = self.electron_charge * particle.velocity.cross(&net_b);

        // f = ma
        force / self.electron_mass
    }
}

impl Field {
    // https://tiggerntatie.github.io/emagnet-py/offaxis/off_axis_loop.html
    fn b_field(&self, particle_pos: &Vector3<Precision>) -> Vector3<Precision> {
        let mut total_beta = Vector3::default();

        let alpha = self.radius / self.position.y;
        let beta = (particle_pos.z - self.position.z) / self.position.y;
        let gamma = (particle_pos.z - self.position.z);

        // Q = sqrt((1.0f + a) * (1.0f + a) + B*B);
        let q = ((1.0 + alpha).powi(2) + beta.powi(2)).sqrt();

        // complete elliptic integral of the first kind
        // k = (int)(4.0f * a / (Q*Q) *10000000);
        // let k = 4.0 * alpha / q.powi(2) * 10_000_000.0;

        // E_k = ee[k];
        // K_k = ek[k];

        // let B_r = self.position.z * ((particle_pos.z - coils[i].x) / r) / (PI * Q)
        //     * ((E_k * (1.0 + alpha * alpha + B * B) / ((Q * Q) - 4.0 * alpha)) - K_k);

        // b.z += coils[i].z * (1.0f)/( M_PI * Q) * (( E_k * (1.0f - a*a - B*B)/((Q*Q)-4.0f*a)) + K_k);
        // b.x += pos.x / r * B_r;
        // b.y += pos.y / r * B_r;

        total_beta
    }
}
