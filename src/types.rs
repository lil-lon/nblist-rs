use std::ops::{Add, Mul, Sub};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vector3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
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
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn norm(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }
}

#[derive(Debug, Clone, Copy)]
pub struct UnitCell {
    pub a: f64,
    pub b: f64,
    pub c: f64,
}

/// Positions are stored in Cartesian coordinates.
#[derive(Debug, Clone)]
pub struct System {
    pub pos: Vec<Vector3>,
    pub cell: UnitCell,
}

#[derive(Debug)]
pub struct Neighbor {
    pub i: usize,
    pub j: usize,
    pub offset: [i32; 3],
    pub distance: f64,
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
    fn test_vector3_add() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 5.0, 6.0);
        assert_eq!(a + b, Vector3::new(5.0, 7.0, 9.0));
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
