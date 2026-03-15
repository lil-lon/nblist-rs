# nblist-rs

CPU-based neighbor list construction in Rust for atomistic simulations.

Given a set of atomic positions and a periodic unit cell, find all atom pairs within a cutoff distance — a core primitive for molecular dynamics and machine learning interatomic potentials.

## Features

- **Full triclinic cell support** — fractional-coordinate binning handles arbitrary lattice vectors, not just orthorhombic boxes
- **Per-axis PBC** — `pbc: [bool; 3]` enables bulk, slab, wire, and isolated geometries in a single interface
- **Zero dependencies**

## Algorithms

| Algorithm | Complexity | Description |
|-----------|-----------|-------------|
| `naive` | O(N²) | All-pairs with image enumeration. Used as ground truth. |
| `cell_list` | O(N) | Fractional-space binning with perpendicular-height-based bin counts. |

Both algorithms share the same test suite to guarantee equivalence.
