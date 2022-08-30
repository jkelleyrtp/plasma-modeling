use fusion::*;
use nalgebra::Vector3;

static CFG: Config = Config {
    dt: 0.0000000000002,
    num_particles: 1,
    steps: 100,
    micro_step: 5000,
};

static COIL_1: Field = Field {
    radius: 0.50,
    position: Vector3::new(0.0, 0.0, -0.50),
    current: 1000.0,
};

static COIL_2: Field = Field {
    radius: 0.50,
    position: Vector3::new(0.0, 0.0, 0.50),
    current: 1000.0,
};

fn main() {
    let mut simulation = Simulation::new();

    simulation.history.reserve(CFG.steps);

    // A biconic cusp, made of two fields
    simulation.fields.push(COIL_1);
    simulation.fields.push(COIL_2);

    // todo: draw a circle around the cusp
    // currently places particle near axis but slightly off
    for _n in 0..CFG.num_particles {
        // let angle = n as f64 * 2.0 * std::f64::consts::PI / CFG.num_particles as f64;

        simulation.particles.push(Electron {
            position: Vector3::new(0.0001, 0.0001, -0.60),
            velocity: Vector3::new(0.0, 0.0, 1e3),
        });
    }

    println!("----------------------------");
    println!("Simulating fusion experiment");
    println!("----------------------------");
    println!("- num fields: {}", simulation.fields.len());
    println!("- num_particles: {}", CFG.num_particles);
    println!("- dt: {:?}", CFG.dt);
    println!("- steps: {}", CFG.steps);
    println!("- micro_steps: {:?}", CFG.micro_step);
    println!("- total steps: {:?}", CFG.steps * CFG.micro_step);
    println!();

    // tick the simulation
    for _ in 0..CFG.steps {
        for _ in 0..CFG.micro_step {
            simulation.tick_eulers(CFG.dt);
        }

        // save the current state
        simulation.history.push(simulation.particles.clone());
    }

    // draw the simulation history
    let mut table = prettytable::Table::new();
    table.set_titles(prettytable::row!["Time", "Pos X", "Pos Y", "Pos Z"]);

    for (id, row) in simulation
        .history
        .into_iter()
        .flatten()
        .map(|f| f.position)
        .enumerate()
    {
        table.add_row(prettytable::row![id, row.x, row.y, row.z]);
    }

    table.set_format(*prettytable::format::consts::FORMAT_NO_LINESEP_WITH_TITLE);

    table.printstd();
}
