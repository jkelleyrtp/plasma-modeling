mod fields;
mod matrix;
pub mod polynomial;
mod tables;
mod util;
pub use fields::Field;
use nalgebra::Vector3;

#[derive(Debug)]
pub struct Config {
    pub dt: Precision,
    pub steps: usize,
    pub num_particles: usize,
    pub micro_step: usize,
}

pub type Precision = f64;

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
