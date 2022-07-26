use fusion::*;
use nalgebra::Vector3;

static CFG: Config = Config {
    dt: 0.000001,
    num_particles: 10,
    steps: 100_000,
    micro_step: 50,
};

fn main() {
    let mut simulation = Simulation::new();

    simulation.history.reserve(CFG.steps);

    // A biconic cusp, made of two fields
    simulation.fields.push(Field {
        radius: 0.50,
        position: [-0.50, 0.0, 0.0],
        orientation: [1.0, 0.0, 0.0],
    });

    simulation.fields.push(Field {
        radius: 0.50,
        position: [0.50, 0.0, 0.0],
        orientation: [-1.0, 0.0, 0.0],
    });

    // draw a circle around the cusp
    for n in 0..CFG.num_particles {
        let angle = n as f64 * 2.0 * std::f64::consts::PI / CFG.num_particles as f64;
        simulation.particles.push(Electron {
            position: Vector3::new(-0.60, angle.cos(), angle.sin()),
            velocity: Vector3::new(0.0, 0.0, 0.0),
        });
    }

    // tick the simulation
    for _ in 0..CFG.steps {
        for _ in 0..CFG.micro_step {
            simulation.tick_eulers(CFG.dt);
        }

        // save the current state
        simulation.history.push(simulation.particles.clone());
    }

    // draw the simulation history
}
