use std::ops::{Add, Mul, Sub};

fn main() {
    println!("Hello, world!");
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Vector3 {
    x: f64,
    y: f64,
    z: f64,
}

impl Add for Vector3 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl Sub for Vector3 {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl Mul<f64> for Vector3 {
    type Output = Self;
    fn mul(self, other: f64) -> Self {
        Self {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}

impl Vector3 {
    fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    fn norm(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }
}

#[derive(Debug, Clone, Copy)]
struct UnitCell {
    a: f64,
    b: f64,
    c: f64,
}

impl UnitCell {
    fn minimum_image(&self, diff: Vector3) -> Vector3 {
        fn wrap_to_half(coord: f64, lattice: f64) -> f64 {
            let mut c = coord;
            while c > lattice / 2.0 {
                c -= lattice;
            }
            while c < -lattice / 2.0 {
                c += lattice;
            }
            c
        }

        Vector3 {
            x: wrap_to_half(diff.x, self.a),
            y: wrap_to_half(diff.y, self.b),
            z: wrap_to_half(diff.z, self.c),
        }
    }
}

#[derive(Debug, Clone)]
struct System {
    pos: Vec<Vector3>,
    cell: UnitCell,
}

impl System {
    // Return the minimum image distance between atoms i and j
    fn distance(&self, i: usize, j: usize) -> f64 {
        let diff = self.pos[j] - self.pos[i];
        let wrapped = self.cell.minimum_image(diff);
        wrapped.norm()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vector3_norm() {
        let vec = Vector3::new(3.0, 4.0, 5.0);
        let result = vec.norm();
        let expected = 50.0_f64.sqrt();
        assert!((result - expected).abs() < 1e-10);
    }

    #[test]
    fn test_minimum_image() {
        let cell = UnitCell {
            a: 2.0,
            b: 4.0,
            c: 8.0,
        };
        let diff = Vector3::new(0.5, -2.0, 6.0);
        let result = cell.minimum_image(diff);
        let expected = Vector3::new(0.5, -2.0, -2.0);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_distance_with_pbc() {
        let pos = vec![Vector3::new(1.0, 2.0, 4.0), Vector3::new(4.5, 3.0, 1.0)];
        let cell = UnitCell {
            a: 5.0,
            b: 5.0,
            c: 5.0,
        };
        let result = System { pos, cell }.distance(0, 1);
        let expected = (1.5 * 1.5 + 1.0 * 1.0 + 2.0 * 2.0_f64).sqrt();
        assert!((result - expected).abs() < 1e-10);
    }

    #[test]
    fn test_distance_with_pbc_large_diff() {
        // Atom B is far outside the cell (diff > L)
        let pos = vec![Vector3::new(1.0, 1.0, 1.0), Vector3::new(14.0, -8.0, 22.0)];
        let cell = UnitCell {
            a: 5.0,
            b: 5.0,
            c: 5.0,
        };
        // diff = (13, -9, 21) → wrap → (-2, 1, 1)
        let result = System { pos, cell }.distance(0, 1);
        let expected = (2.0 * 2.0 + 1.0 * 1.0 + 1.0 * 1.0_f64).sqrt();
        assert!((result - expected).abs() < 1e-10);
    }

    #[test]
    fn test_vector3_add() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 5.0, 6.0);
        assert_eq!(a + b, Vector3::new(5.0, 7.0, 9.0))
    }

    #[test]
    fn test_vector3_sub() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 8.0, 16.0);
        assert_eq!(b - a, Vector3::new(3.0, 6.0, 13.0));
    }

    #[test]
    fn test_vector3_mul() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        assert_eq!(a * 5.5, Vector3::new(5.5, 11.0, 16.5));
    }
}
