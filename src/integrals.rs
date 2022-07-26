use crate::polynomial::Polynomial;

mod constants {
    pub const MACHEP: f64 = 1.11022302462515654042E-16;
}

pub fn ellpe(x: f64) -> f64 {
    static P: Polynomial = Polynomial::new([
        1.53552577301013293365E-4,
        2.50888492163602060990E-3,
        8.68786816565889628429E-3,
        1.07350949056076193403E-2,
        7.77395492516787092951E-3,
        7.58395289413514708519E-3,
        1.15688436810574127319E-2,
        2.18317996015557253103E-2,
        5.68051945617860553470E-2,
        4.43147180560990850618E-1,
        1.00000000000000000299E0,
    ]);

    static Q: Polynomial = Polynomial::new([
        3.27954898576485872656E-5,
        1.00962792679356715133E-3,
        6.50609489976927491433E-3,
        1.68862163993311317300E-2,
        2.61769742454493659583E-2,
        3.34833904888224918614E-2,
        4.27180926518931511717E-2,
        5.85936634471101055642E-2,
        9.37499997197644278445E-2,
        2.49999999999888314361E-1,
        1.00000000000000000299E0,
    ]);

    if x <= 0.0 || x > 1.0 {
        if x == 0.0 {
            return 1.0;
        } else {
            return 0.0;
        }
    }

    P.eval(x) - x.ln() * (x * Q.eval(x))
}

/// Returns the incomplete elliptic integral of the first kind.
fn ellpk(x: f64) -> f64 {
    static P: Polynomial = Polynomial::new([
        1.37982864606273237150E-4,
        2.28025724005875567385E-3,
        7.97404013220415179367E-3,
        9.85821379021226008714E-3,
        6.87489687449949877925E-3,
        6.18901033637687613229E-3,
        8.79078273952743772254E-3,
        1.49380448916805252718E-2,
        3.08851465246711995998E-2,
        9.65735902811690126535E-2,
        1.38629436111989062502E0,
    ]);

    static Q: Polynomial = Polynomial::new([
        2.94078955048598507511E-5,
        9.14184723865917226571E-4,
        5.94058303753167793257E-3,
        1.54850516649762399335E-2,
        2.39089602715924892727E-2,
        3.01204715227604046988E-2,
        3.73774314173823228969E-2,
        4.88280347570998239232E-2,
        7.03124996963957469739E-2,
        1.24999999999870820058E-1,
        4.99999999999999999821E-1,
    ]);

    /* log(4) */
    const C1: f64 = 1.386_294_361_119_890_6;

    if (x < 0.0) || (x > 1.0) {
        return 0.0;
    }

    if x > constants::MACHEP {
        P.eval(x) - x.ln() * (x * Q.eval(x))
    } else if x == 0.0 {
        f64::MAX
    } else {
        C1 - 0.5 * x.ln()
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn ellpe_in_range() {
        dbg!(ellpk(0.1));
    }
}