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

    pub fn cross(&self, other: &Vector3) -> Vector3 {
        Vector3::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }
}

#[derive(Debug, Clone, Copy)]
pub struct UnitCell {
    pub a: Vector3,
    pub b: Vector3,
    pub c: Vector3,
    // Inverse of the lattice matrix H⁻¹, stored in row-major order
    inv: [f64; 9],
    volume: f64,
    // Perpendicular heights: d_a = V / |b × c|, etc.
    heights: [f64; 3],
}

impl UnitCell {
    pub fn new(a: Vector3, b: Vector3, c: Vector3) -> Result<Self, String> {
        let (inv, volume) = Self::calc_inverse(a, b, c)?;
        let heights = [
            volume / b.cross(&c).norm(),
            volume / a.cross(&c).norm(),
            volume / a.cross(&b).norm(),
        ];
        Ok(Self {
            a,
            b,
            c,
            inv,
            volume,
            heights,
        })
    }

    fn calc_inverse(a: Vector3, b: Vector3, c: Vector3) -> Result<([f64; 9], f64), String> {
        // row-major order
        let h = [a.x, b.x, c.x, a.y, b.y, c.y, a.z, b.z, c.z];
        let det = h[0] * (h[4] * h[8] - h[5] * h[7]) - h[1] * (h[3] * h[8] - h[5] * h[6])
            + h[2] * (h[3] * h[7] - h[4] * h[6]);

        if det.abs() < 1e-12 {
            return Err("Lattice vectors are illegal: linearly dependent (det ≈ 0)".to_string());
        }

        let inv_det = 1.0 / det;
        let inv_mat = [
            (h[4] * h[8] - h[5] * h[7]) * inv_det,
            (h[2] * h[7] - h[1] * h[8]) * inv_det,
            (h[1] * h[5] - h[2] * h[4]) * inv_det,
            (h[5] * h[6] - h[3] * h[8]) * inv_det,
            (h[0] * h[8] - h[2] * h[6]) * inv_det,
            (h[2] * h[3] - h[0] * h[5]) * inv_det,
            (h[3] * h[7] - h[4] * h[6]) * inv_det,
            (h[1] * h[6] - h[0] * h[7]) * inv_det,
            (h[0] * h[4] - h[1] * h[3]) * inv_det,
        ];
        Ok((inv_mat, det.abs()))
    }

    pub fn get_cartesian(&self, frac: &Vector3) -> Vector3 {
        self.a * frac.x + self.b * frac.y + self.c * frac.z
    }

    pub fn get_fractional(&self, cart: &Vector3) -> Vector3 {
        Vector3::new(
            cart.x * self.inv[0] + cart.y * self.inv[1] + cart.z * self.inv[2],
            cart.x * self.inv[3] + cart.y * self.inv[4] + cart.z * self.inv[5],
            cart.x * self.inv[6] + cart.y * self.inv[7] + cart.z * self.inv[8],
        )
    }

    pub fn get_volume(&self) -> f64 {
        self.volume
    }

    /// Perpendicular heights [d_a, d_b, d_c]
    pub fn get_heights(&self) -> [f64; 3] {
        self.heights
    }
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

    #[test]
    fn test_vector3_cross_basis() {
        // i × j = k
        let i = Vector3::new(1.0, 0.0, 0.0);
        let j = Vector3::new(0.0, 1.0, 0.0);
        let k = i.cross(&j);
        assert_eq!(k, Vector3::new(0.0, 0.0, 1.0));

        // Anti-commutativity: j × i = -k
        let neg_k = j.cross(&i);
        assert_eq!(neg_k, Vector3::new(0.0, 0.0, -1.0));
    }

    #[test]
    fn test_vector3_cross_general() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 5.0, 6.0);
        let c = a.cross(&b);
        // (2*6 - 3*5, 3*4 - 1*6, 1*5 - 2*4) = (-3, 6, -3)
        assert_eq!(c, Vector3::new(-3.0, 6.0, -3.0));
    }

    #[test]
    fn test_unitcell_new_triclinic() {
        let cell = UnitCell::new(
            Vector3::new(5.0, 0.0, 0.0),
            Vector3::new(1.0, 6.0, 0.0),
            Vector3::new(0.5, 0.5, 7.0),
        );
        assert!(cell.is_ok());
    }

    #[test]
    fn test_unitcell_new_linearly_dependent() {
        let cell = UnitCell::new(
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
        );
        assert!(cell.is_err());
    }

    #[test]
    fn test_unitcell_inverse_correctness() {
        // Triclinic cell: verify H * H⁻¹ = I
        let cell = UnitCell::new(
            Vector3::new(5.0, 0.0, 0.0),
            Vector3::new(1.0, 6.0, 0.0),
            Vector3::new(0.5, 0.5, 7.0),
        )
        .unwrap();

        let h = [
            cell.a.x, cell.b.x, cell.c.x, cell.a.y, cell.b.y, cell.c.y, cell.a.z, cell.b.z,
            cell.c.z,
        ];

        for i in 0..3 {
            for j in 0..3 {
                let mut sum = 0.0;
                for k in 0..3 {
                    sum += h[i * 3 + k] * cell.inv[k * 3 + j];
                }
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (sum - expected).abs() < 1e-10,
                    "H * H⁻¹ [{},{}] = {}, expected {}",
                    i,
                    j,
                    sum,
                    expected
                );
            }
        }
    }

    #[test]
    fn test_get_cartesian_triclinic() {
        let cell = UnitCell::new(
            Vector3::new(5.0, 0.0, 0.0),
            Vector3::new(1.0, 6.0, 0.0),
            Vector3::new(0.5, 0.5, 7.0),
        )
        .unwrap();

        // x = 0.5*5 + 0.3*1 + 0.2*0.5 = 2.5 + 0.3 + 0.1 = 2.9
        // y = 0.5*0 + 0.3*6 + 0.2*0.5 = 0 + 1.8 + 0.1 = 1.9
        // z = 0.5*0 + 0.3*0 + 0.2*7   = 0 + 0 + 1.4   = 1.4
        let cart = cell.get_cartesian(&Vector3::new(0.5, 0.3, 0.2));
        assert!((cart.x - 2.9).abs() < 1e-10);
        assert!((cart.y - 1.9).abs() < 1e-10);
        assert!((cart.z - 1.4).abs() < 1e-10);
    }

    #[test]
    fn test_get_fractional_triclinic() {
        // Same cell, verify the reverse: cart=(2.9, 1.9, 1.4) → frac=(0.5, 0.3, 0.2)
        let cell = UnitCell::new(
            Vector3::new(5.0, 0.0, 0.0),
            Vector3::new(1.0, 6.0, 0.0),
            Vector3::new(0.5, 0.5, 7.0),
        )
        .unwrap();
        let frac = cell.get_fractional(&Vector3::new(2.9, 1.9, 1.4));
        assert!((frac.x - 0.5).abs() < 1e-10);
        assert!((frac.y - 0.3).abs() < 1e-10);
        assert!((frac.z - 0.2).abs() < 1e-10);
    }

    #[test]
    fn test_roundtrip_cart_frac_cart() {
        // Verify cart -> frac -> cart recovers the original vector
        let cell = UnitCell::new(
            Vector3::new(6.0, 1.0, 0.5),
            Vector3::new(0.5, 8.0, 1.5),
            Vector3::new(1.0, 0.5, 5.0),
        )
        .unwrap();
        let original = Vector3::new(3.7, 2.1, 5.5);
        let frac = cell.get_fractional(&original);
        let recovered = cell.get_cartesian(&frac);
        assert!((recovered.x - original.x).abs() < 1e-10);
        assert!((recovered.y - original.y).abs() < 1e-10);
        assert!((recovered.z - original.z).abs() < 1e-10);
    }

    #[test]
    fn test_roundtrip_frac_cart_frac() {
        // Verify frac -> cart -> frac recovers the original vector
        let cell = UnitCell::new(
            Vector3::new(6.0, 1.0, 0.5),
            Vector3::new(0.5, 8.0, 1.5),
            Vector3::new(1.0, 0.5, 5.0),
        )
        .unwrap();
        let original = Vector3::new(0.7, 0.3, 0.9);
        let cart = cell.get_cartesian(&original);
        let recovered = cell.get_fractional(&cart);
        assert!((recovered.x - original.x).abs() < 1e-10);
        assert!((recovered.y - original.y).abs() < 1e-10);
        assert!((recovered.z - original.z).abs() < 1e-10);
    }

    #[test]
    fn test_volume_triclinic() {
        // a=(5,0,0), b=(1,6,0), c=(0.5,0.5,7)
        // V = a · (b × c)
        // b × c = (6*7 - 0*0.5, 0*0.5 - 1*7, 1*0.5 - 6*0.5)
        //       = (42, -7, -2.5)
        // a · (b × c) = 5*42 + 0*(-7) + 0*(-2.5) = 210
        let cell = UnitCell::new(
            Vector3::new(5.0, 0.0, 0.0),
            Vector3::new(1.0, 6.0, 0.0),
            Vector3::new(0.5, 0.5, 7.0),
        )
        .unwrap();
        assert!((cell.get_volume() - 210.0).abs() < 1e-10);
    }

    #[test]
    fn test_heights_triclinic() {
        // a=(5,0,0), b=(1,6,0), c=(0.5,0.5,7), V=210
        // d_a = V / |b × c| = 210 / |(42, -7, -2.5)| = 210 / sqrt(1764+49+6.25)
        //     = 210 / sqrt(1819.25)
        // d_b = V / |a × c| = 210 / |(0*7-0*0.5, 0*0.5-5*7, 5*0.5-0*0.5)|
        //     = 210 / |(0, -35, 2.5)| = 210 / sqrt(1231.25)
        // d_c = V / |a × b| = 210 / |(0*0-0*6, 0*1-5*0, 5*6-0*1)|
        //     = 210 / |(0, 0, 30)| = 210 / 30 = 7
        let cell = UnitCell::new(
            Vector3::new(5.0, 0.0, 0.0),
            Vector3::new(1.0, 6.0, 0.0),
            Vector3::new(0.5, 0.5, 7.0),
        )
        .unwrap();
        let h = cell.get_heights();
        assert!((h[0] - 210.0 / 1819.25_f64.sqrt()).abs() < 1e-10);
        assert!((h[1] - 210.0 / 1231.25_f64.sqrt()).abs() < 1e-10);
        assert!((h[2] - 7.0).abs() < 1e-10);
    }
}
