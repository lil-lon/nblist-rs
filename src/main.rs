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

// Positions are stored in Cartesian coordinates.
#[derive(Debug, Clone)]
struct System {
    pos: Vec<Vector3>,
    cell: UnitCell,
}

impl System {
    fn build_neighbor_list(&self, cutoff: f64) -> Vec<Neighbor> {
        let mut result = Vec::new();
        let a_replicas = (cutoff / self.cell.a).ceil() as i32;
        let b_replicas = (cutoff / self.cell.b).ceil() as i32;
        let c_replicas = (cutoff / self.cell.c).ceil() as i32;
        let n = self.pos.len();
        for i in 0..n {
            for j in 0..n {
                for a_ind in -a_replicas..=a_replicas {
                    for b_ind in -b_replicas..=b_replicas {
                        for c_ind in -c_replicas..=c_replicas {
                            if i == j && (a_ind, b_ind, c_ind) == (0, 0, 0) {
                                continue;
                            }
                            let pos_j = self.pos[j]
                                + Vector3::new(
                                    a_ind as f64 * self.cell.a,
                                    b_ind as f64 * self.cell.b,
                                    c_ind as f64 * self.cell.c,
                                );
                            let distance = (pos_j - self.pos[i]).norm();
                            if distance < cutoff {
                                result.push(Neighbor {
                                    i,
                                    j,
                                    offset: [a_ind, b_ind, c_ind],
                                    distance,
                                })
                            }
                        }
                    }
                }
            }
        }
        result
    }
}

#[derive(Debug)]
struct Neighbor {
    i: usize,
    j: usize,
    offset: [i32; 3],
    distance: f64,
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
    fn test_neighbor_list_basic() {
        // 4 atoms in a 10x10x10 box, cutoff = 3.0
        // A(0,0,0) B(2,0,0) C(5,5,5) D(8,0,0)
        // A-B: 2.0 ✓, B-A: 2.0 ✓
        // A-D via PBC: 2.0 ✓, D-A via PBC: 2.0 ✓
        let sys = System {
            pos: vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(2.0, 0.0, 0.0),
                Vector3::new(5.0, 5.0, 5.0),
                Vector3::new(8.0, 0.0, 0.0),
            ],
            cell: UnitCell {
                a: 10.0,
                b: 10.0,
                c: 10.0,
            },
        };
        let result = sys.build_neighbor_list(3.0);
        // Both directions: (0,1), (1,0), (0,3), (3,0)
        let pairs: Vec<(usize, usize)> = result.iter().map(|n| (n.i, n.j)).collect();
        assert!(pairs.contains(&(0, 1)));
        assert!(pairs.contains(&(1, 0)));
        assert!(pairs.contains(&(0, 3)));
        assert!(pairs.contains(&(3, 0)));
        assert_eq!(result.len(), 4);
    }

    #[test]
    fn test_neighbor_list_no_neighbors() {
        let sys = System {
            pos: vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(5.0, 5.0, 5.0)],
            cell: UnitCell {
                a: 20.0,
                b: 20.0,
                c: 20.0,
            },
        };
        let result = sys.build_neighbor_list(1.0);
        assert!(result.is_empty());
    }

    #[test]
    fn test_neighbor_list_pbc_corner() {
        // A(0.5, 0.5, 0.5) B(9.5, 9.5, 9.5) in 10x10x10 box
        // dist via PBC = sqrt(3) ≈ 1.73
        let sys = System {
            pos: vec![Vector3::new(0.5, 0.5, 0.5), Vector3::new(9.5, 9.5, 9.5)],
            cell: UnitCell {
                a: 10.0,
                b: 10.0,
                c: 10.0,
            },
        };
        let result = sys.build_neighbor_list(2.0);
        // Both directions
        assert_eq!(result.len(), 2);
        assert!(result.iter().any(|n| n.i == 0 && n.j == 1));
        assert!(result.iter().any(|n| n.i == 1 && n.j == 0));
    }
}
