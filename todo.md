# Step 1: Data Structures & Minimum Image

## Data Structures
- [ ] `Vector3`: struct with x, y, z (f64)
  - [ ] `new(x, y, z) -> Self`
  - [ ] `norm(&self) -> f64`
  - [ ] `impl Add<Vector3>` for `+`
  - [ ] `impl Sub<Vector3>` for `-`
  - [ ] `impl Mul<f64>` for scalar `*`
- [ ] `UnitCell`: struct with a, b, c (f64, orthogonal only)
- [ ] `System`: struct with `positions: Vec<Vector3>` and `cell: UnitCell`

## Minimum Image Convention
- [ ] `UnitCell::minimum_image(&self, diff: Vector3) -> Vector3`
  - Apply PBC: each component wraps to [-L/2, L/2)
- [ ] `System::distance(&self, i: usize, j: usize) -> f64`
  - Compute `positions[j] - positions[i]`, apply minimum image, return norm

## Validation
- [ ] `#[test]` for Vector3 operations (add, sub, norm)
- [ ] `#[test]` for minimum image (e.g. atoms near cell boundary)
- [ ] `#[test]` for System::distance with PBC

# Step 2: Naive O(N^2) Neighbor List

- [ ] `fn build_neighbor_list(system: &System, cutoff: f64) -> Vec<(usize, usize)>`
  - Double loop over all pairs (i < j)
  - Use System::distance, collect pairs within cutoff
- [ ] `#[test]` for small known system (e.g. 4 atoms in a box)
- [ ] Benchmark with N = 1000, 5000, 10000

# Step 3: Cell List (O(N))

- [ ] Divide box into cells of size >= cutoff
- [ ] Assign atoms to cells
- [ ] Search only neighboring cells (27 cells in 3D)
- [ ] `#[test]` verify same result as naive
- [ ] Benchmark and compare with naive

# Step 4: Triclinic Cell Support

- [ ] `UnitCell` with 3 lattice vectors (Vector3 x3)
- [ ] Fractional <-> Cartesian coordinate conversion
- [ ] Minimum image in fractional coordinates
- [ ] Update cell list for triclinic

# Step 5: Further Optimization

- [ ] Verlet list (skin distance, rebuild check)
- [ ] Data layout optimization (AoS vs SoA)
- [ ] Parallelization with rayon