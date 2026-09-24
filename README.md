# Orientifold Database Generator for Calabi-Yau Threefolds in Toric Ambient Spaces

This repository contains tools for generating databases of orientifold configurations of Calabi-Yau threefolds. It uses `cytools` or Hugging Face's dataset for retrieving Calabi-Yau threefolds, and systematically enumerates the discrete symmetries to compute geometric data for string theory models.

## Prerequisites

You need the `cytools` environment.

Note: To compute prime toric divisor cohomologies, you also need to have `cohomcalg` installed on your system and accessible via the `PATH` environment variable, or by setting `COHOMCALG_PATH` explicitly.

## Scripts

### 1. `generate_complete_orientifold_db_mixed.py`

This script generates the orientifold configurations and saves them in a JSON database in the `databases/` folder.

#### Usage:

You can use the built-in help flag to see all options:
```bash
conda run -n cytools python generate_complete_orientifold_db_mixed.py --help
```

**Output:**
```
usage: generate_complete_orientifold_db_mixed.py [-h] --h11 H11
                                                 [--compute_cohomologies]
                                                 [--workers WORKERS]
                                                 [--source {cytools,hf}]

Generate a complete orientifold database for a specified Hodge number h1,1.
This script queries a dataset of polytopes (from TU Wien via cytools or Hugging Face),
computes valid orientifold involutions, evaluates O-planes, and extracts
tadpole constraints (N_D3, chi_Y4) alongside geometrical properties.

options:
  -h, --help            show this help message and exit
  --h11 H11             h11 value to process
  --compute_cohomologies
                        Run cohomcalg for O7 planes to compute divisor cohomologies.
  --workers WORKERS     Number of workers for multiprocessing.
  --source {cytools,hf}
                        Source of the polytopes database. 'cytools' uses the TU Wien database, 'hf' uses Hugging Face 'calabi-yau-data/polytopes-4d' dataset.

Examples:
  python generate_complete_orientifold_db_mixed.py --h11 1 --workers 4
  python generate_complete_orientifold_db_mixed.py --h11 2 --source hf --compute_cohomologies
```

#### Demo Command:

To test the generation quickly, you can run it for $h^{1,1} = 1$ using the `hf` source:

```bash
conda run -n cytools python generate_complete_orientifold_db_mixed.py --h11 1 --source hf
```

This will run relatively quickly and produce an output file at `databases/cy_orient_h11_1.json`.

### 2. `generate_single_polytope_orientifold_mixed.py`

Similar to the main database script, but specifically computes orientifold data for a single user-specified polytope. It shares the exact same JSON output format as the complete database.

#### Usage:
```bash
python generate_single_polytope_orientifold_mixed.py -h
```

**Output:**
```
usage: generate_single_polytope_orientifold_mixed.py [-h]
                                                     [--compute_cohomologies]
                                                     [--out OUT]
                                                     [--points POINTS]

Generate orientifolds and tadpole constraints for a single specified polytope.
Computes valid orientifold involutions, evaluates O-planes, and extracts
tadpole constraints alongside geometrical properties.

options:
  -h, --help            show this help message and exit
  --compute_cohomologies
                        Run cohomcalg for O7 planes to compute divisor cohomologies.
  --out OUT             Output JSON file.
  --points POINTS       JSON string or Python list string of the polytope points.
                        Example: '[[-1, -6, -8, 1], [0, 0, 1, 1], [1, 0, 1, 1], [0, 1, 1, 1], [0, -1, -1, 1], [0, 0, 0, -1]]'

Examples:
  python generate_single_polytope_orientifold_mixed.py
  python generate_single_polytope_orientifold_mixed.py --compute_cohomologies
  python generate_single_polytope_orientifold_mixed.py --points '[[-1,-1,-1,1], [1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]]'
```

### 3. `read_orientifold_db.py`

A helper script designed to easily parse the generated JSON database files and provide a structural summary. It serves as both a CLI utility and an importable module for your own scripts.

#### Command Line Usage

Run the script directly on a generated database file to get a statistical summary:

```bash
conda run -n cytools python read_orientifold_db.py databases/cy_orient_h11_1.json
```

**Example Output:**
```
Loaded database with 5 Calabi-Yau threefolds.
Total triangulations: 5
Total orientifolds configurations: 271

--- Example CY (First Entry) ---
POLYID: 1
h11: 1
h21: 11
Number of orientifolds: 2
```

#### Programmatic Usage (Importing)

You can import functions from `read_orientifold_db.py` into your own Python code to load and analyze the database programmatically. 

```python
from read_orientifold_db import load_orientifold_db, print_summary

# 1. Load the database (returns a list of dictionaries, one for each Calabi-Yau geometry)
db = load_orientifold_db("databases/cy_orient_h11_1.json")

# 2. Iterate through the Calabi-Yau configurations
for cy in db:
    poly_id = cy.get("POLYID")
    h11 = cy.get("h11")
    
    # 3. Iterate through all orientifold entries for this specific geometry
    for orientifold in cy.get("orientifolds", []):
        invol = orientifold["INVOL"]
        if orientifold["tadpole_data"]:
            chi_D7 = orientifold["tadpole_data"][0]["chi_D7"]
            print(f"POLYID {poly_id} | Involution {invol} | D7 Euler Char: {chi_D7}")
```

- `load_orientifold_db(filepath)`: Reads the target JSON and parses it into memory as a standard Python `list` of `dict`s.
- `print_summary(db)`: Parses the loaded list and prints the total count of CY threefolds, triangulations, and total orientifold constraints across the entire dataset, followed by a sample breakdown of the very first CY geometry.

### 4. `generate_trilayer_orientifolds.py`

This script computes the orientifold data specifically for the favorable geometries contained within the `databases/trilayer.json` database. The trilayer polytopes define Calabi--Yau hypersurfaces within toric ambient spaces constructed via the Kreuzer--Skarke algorithm, providing the geometric data required for the F-theory fourfold uplift. 

Under the Sen limit, the resolution of the F-theory Calabi--Yau fourfold involves refining the normal fan to resolve $I_0^*$ singularities located over the prime toric divisors hosting non-Higgsable clusters (NHCs). The script parses the specified geometries, evaluates the valid involutions, and computes the required tadpole conditions, effectively mirroring the single-polytope pipeline across the trilayer dataset. The output is stored in `databases/trilayer_orientifolds.json`.

#### Usage:
```bash
python generate_trilayer_orientifolds.py -h
```

**Output:**
```
usage: generate_trilayer_orientifolds.py [-h] [--compute_cohomologies]
                                         [--workers WORKERS]
                                         [--h11_max H11_MAX]

options:
  -h, --help            show this help message and exit
  --compute_cohomologies
                        Run cohomcalg for O7 planes
  --workers WORKERS     Number of workers
  --h11_max H11_MAX     Maximum h11 (Picard_number_N) to process
```

### 5. `verify_all_intB_relations.py`

This utility script explicitly verifies the topological intersection relations between the Type IIB double cover Calabi--Yau threefold $\IX_3$ and its orientifold quotient base $B_3$. The continuous topological quantities on the quotient base uplift to the Calabi--Yau threefold with an overall factor of two, arising from the double cover branched over the O7-plane locus.

Using the `cytools` package, the script builds the corresponding toric varieties and systematically evaluates three topological consistency relations:
1. **Triple Intersections**:
   $$\int_{\IX_3} k_a \wedge k_b \wedge k_c = 2 \int_{B_3} j_a \wedge j_b \wedge j_c$$
2. **Second Chern Class Contractions**:
   $$\int_{\IX_3} c_2(T\IX_3) \wedge k_a = 2 \int_{B_3} j_a \wedge \left(c_1^2(TB_3) + c_2(TB_3)\right)$$
3. **Euler Characteristics**:
   $$\chi(\IX_3) = 2 \int_{B_3} \left(c_3(TB_3) - c_1(TB_3)c_2(TB_3) - 2c_1(TB_3)^3\right)$$
where $k_a$ correspond to the K\"ahler forms in $\IX_3$, and $j_a$ denote the divisors in the base $B_3$.

#### Usage:
```bash
python verify_all_intB_relations.py
```

### Geometric admissibility (base locus check)

The script checks geometric admissibility to prevent base locus singularities on the Calabi-Yau hypersurface. It computes the monomial exponents $E_{m, j} = \langle m, r_j \rangle + 1$ for all dual polytope points $m \in \Delta^\circ$ and rays $r_j \in \Delta$. For each sign configuration, it finds invariant monomials with a sign parity matching the canonical sign $S_{\text{canon}}$. It evaluates these monomials against the coordinate intersections permitted by the Stanley-Reisner ideal. If all invariant monomials vanish on an allowed toric stratum, the script rejects the orientifold configuration. Rejected configurations do not appear in the JSON database.

## Understanding the JSON Output

The generated JSON file is a list of dictionaries. Each dictionary represents a favorable Calabi--Yau threefold and contains:

- `POLYID`: A unique identifier for the polytope.
- `h11`: Hodge number $h^{1,1}$.
- `h21`: Hodge number $h^{2,1}$.
- `polytope_points`: Points in the $N$-lattice defining the dual polytope.
- `triangulations`: Geometric data for triangulations, including intersections, divisor basis, and the second Chern class.
- `orientifolds`: A list of valid involutions alongside fixed planes (O3, O5, O7, O9) and D3/D7 tadpole cancellation properties.

Each orientifold entry explicitly contains the following fields:

### Basic Data & Involution
- `TRIANGN`: Index of the triangulation used for this orientifold configuration.
- `INVOLN`: Index of the involution.
- `INVOL`: String representation of the coordinate permutation $x_i \to x_{P(i)}$.
- `SCANON`: Canonical sign $S_{\text{canon}} \in \{1, -1\}$ chosen for the orientifold action.
- `BCANON`: List of signs $b_{\text{canon}} \in \{1, -1\}$ for each homogeneous coordinate, representing the diagonal sign shift of the canonical representative involution orbit.
- `SOMEGA`: The overall sign of the orientifold action $S_\Omega$. Mathematically, this is defined as:
  $$S_\Omega = S_{\text{canon}} \cdot \text{sign}(P) \cdot \prod b_{\text{canon}}$$
  In code, this is computed as `S_Omega = S_canon * sign_perm * prod_b_canon`, where `sign_perm` is the parity of the coordinate permutation $P$ and `prod_b_canon` is the product of signs acquired by the defining polynomials.
- `is_reducible`: A boolean flag determining if the CY vanishes identically under the orientifold equations: `S_canon * parity_lam == -1`, where `parity_lam = (-1) ** ((h11 - trace) // 2)` is the determinant of the induced action on the $H^{1,1}$ cohomology.
- `integer_hodge`: Boolean indicating whether the equivariant Hodge numbers are strictly integers.

### Equivariant Hodge Numbers
The orientifold projection splits the cohomology groups into even and odd eigenspaces, decomposing the moduli spaces. The dimensions of the eigenspaces are determined by:
$$h^{1,1}_{\pm} = \frac{1}{2} \left(h^{1,1} \pm \text{Tr}(\Lambda)\right)$$
$$h^{2,1}_+ = \frac{1}{2} \left( h^{2,1} + 1 - S_\Omega + h^{1,1}_+ - h^{1,1}_- - \frac{1}{2} \chi_{\text{fix}} \right)$$
$$h^{2,1}_- = h^{2,1} - h^{2,1}_+$$
where $S_\Omega = +1$ corresponds to O3/O7-plane configurations, and $S_\Omega = -1$ defines O5/O9-planes. The total Euler characteristic of the fixed loci, $\chi_{\text{fix}}$, decomposes into the sum of the Euler characteristics of the fixed surfaces and the number of isolated fixed points.

In code, these are implemented via the permutation matrix $P$ acting on the divisor basis $Q$, with $\Lambda = Q P Q^{-1}$ and $\text{trace} = \text{Tr}(\Lambda)$:
```python
# h11 and h21 are the Hodge numbers of the original Calabi--Yau
h11_plus = round((h11 + trace) / 2)
h11_minus = round((h11 - trace) / 2)
h21_diff = (2 - 2 * S_Omega + 2 * (h11_plus - h11_minus) - chi_fix) / 2
h21_plus = (h21 + h21_diff) / 2
h21_minus = (h21 - h21_diff) / 2
```
The dictionary stores `h11_plus`, `h11_minus`, `h21_plus`, and `h21_minus`.

### Tadpole Data
The `tadpole_data` field contains the data necessary to satisfy the D3/D7 tadpole cancellation conditions. It branches into all combinations of weak-coupling (Whitney umbrella, `WU`) and F-theory $\text{SO}(8)$ limits for the O7-planes.

Each entry corresponds to a specific combination of choices for the O7-planes:
- `choices`: Dictionary mapping O7-plane ideals to their choice (`WU` or `SO8`) and the resulting D7-brane Euler characteristic $\chi(S_{D7})$.
- `N_O3`: Total number of O3-planes.
- `chi_O7_total`: Total Euler characteristic of O7-planes.
- `chi_D7`: D7-brane Euler characteristic for this combination.

**D3-Brane Charges:**
The net D3-brane charge $N_{D3}$ induced by the localized sources is evaluated in the double cover $\IX_3$ as:
$$N_{D3} = \frac{N_{O3}}{2} + \frac{\chi(S)}{6} + \frac{\chi(S_{D7})}{24}$$
- `N_D3`: Exact fractional D3-charge.
- `N_flux_frac`: Fractional part of the flux.
- `N_D3_mobile_max`: Maximum integer mobile D3 branes.

For a generic smooth divisor $S$, the Euler characteristic is computed by integrating its top Chern class via the adjunction formula:
$$\chi(S) = \int_{\IX_3} \left( [S]^3 + c_2(T\IX_3) \wedge [S] \right)$$

If an O7-plane wraps a rigid divisor, the tadpole is typically canceled by placing four D7-branes and their orientifold images on top of the O7-plane, generating an $\text{SO}(8)$ gauge group. Their contribution to the Euler characteristic is:
$$\chi(S_{D7_{\mathrm{SO}(8)}}) = 8 \chi(S) = \int_{\IX_3} \left( 8 [S]^3 + 8 c_2(T\IX_3) \wedge [S] \right)$$

Alternatively, if the divisor admits complex structure deformations, the branes can recombine into a single Whitney umbrella wrapping the homology class $8[S]$. Due to its non-trivial singular locus, the Euler characteristic is corrected to:
$$\chi(S_{D7_{\text{WU}}}) = \int_{\IX_3} \left( 344 [S]^3 + 8 c_2(T\IX_3) \wedge [S] \right)$$

In code, the $\text{SO}(8)$ limit contributes $8 \cdot (\text{c2D\_op} + \text{D3\_op})$, whereas the Whitney umbrella limit contributes $8 \cdot \text{c2D\_op} + 344 \cdot \text{D3\_op}$.

**Fourfold Euler Characteristics:**
- `chi_Y4_corrected`: Total $\chi(Y_4)$ computed via $12 \cdot N_{D3}$.
- `chi_Y4_correction_O3`: The O3-plane nodal shift correction $6 N_{O3}$.
- `chi_Y4_naive`: Uncorrected $\chi(Y_4)$ limit, `chi_Y4_corrected - 6 N_{O3}`.
- `nodal_singularities`: Boolean indicating if $n^S_{df=0} \neq 0$.
- `n_S_df0`: Singularities count $n^S_{df=0} = \int_{O7} (c_2(S) - c_1(S)^2)$.

**ED3 Instanton Divisors:**
- `ed3_instantons`: List containing data on rigid divisors wrapping ED3-instantons, including naive and corrected vertical fluxes, and cohomologies (`h10`, `h20`, `h11`). Complete rigidity ($h^{1,0} = h^{2,0} = 0$) is the simplest prerequisite for instanton contributions to the non-perturbative superpotential.

### O-Planes (`OPLANES`)
List of fixed O-planes resulting from the orientifold equations. Each entry contains:
- `OIDEAL`: List of coordinate constraints setting the fixed locus.
- `ODIM`: Dimension of the O-plane. Derived from the number of constraint equations relative to the ambient dimension.
- `chi`: The computed Euler characteristic of the fixed plane.

For an O7-plane located at divisor $S \equiv D_{O7}$, the geometric quantities are mathematically computed via sparse intersection numbers ($\kappa_{ijk}$) and the second Chern class ($c_{2,i}$):
```python
D3 = np.einsum("ijk,i,j,k->", kappa, D, D, D) # Self intersection
c2D = np.dot(c2, D)                           # Chern intersection
chi = round(D3 + c2D)
```
The individual terms are stored as `D3` and `c2D`.

### Example Orientifold Entry

Below is an explicit example of an entry in the database `orientifolds` list:

```json
{
  "POLYID": 11,
  "h11": 2,
  "h21": 272,
  "polytope_points": [
    [0, 0, 0, 0],
    [-1, -1, 1, 0],
    [-1, -1, 1, 1],
    [-1, -1, 4, -1],
    [-1, 2, -1, 0],
    [1, -1, 0, 0],
    [-1, -1, 2, 0],
    [-1, 0, 1, 0],
    [-1, 1, 0, 0],
    [0, -1, 1, 0]
  ],
  "triangulations": [
    "... (geometric data) ..."
  ],
  "orientifolds": [
    {
      "TRIANGN": 1,
      "INVOLN": 17,
      "INVOL": "{x1->x2, x2->x1}",
      "SCANON": 1,
      "SOMEGA": -1,
      "BCANON": [-1, -1, -1, 1, -1, 1],
      "h11_plus": 2,
      "h11_minus": 0,
      "h21_plus": 128,
      "h21_minus": 144,
      "integer_hodge": true,
      "is_reducible": false,
      "tadpole_data": [
        {
          "choices": {
            "['x1-x2']": [
              "WU",
              288
            ]
          },
          "N_O3": 4,
          "chi_O7_total": 36,
          "chi_D7": 288,
          "N_D3": 20.0,
          "N_flux_frac": 0.0,
          "N_D3_mobile_max": 20,
          "chi_Y4_correction_O3": 24,
          "chi_Y4_naive": 216,
          "chi_Y4_corrected": 240,
          "nodal_singularities": true,
          "n_S_df0": 30,
          "ed3_instantons": [
            {
              "divisor": "x6",
              "h10": 0,
              "h20": 0,
              "chi": 3,
              "h11": 1,
              "n_O3": 1,
              "chi_correction_vertical_divisor": 2,
              "chi_vertical_naive": -12,
              "chi_vertical_corrected": -10
            }
          ]
        },
        {
          "choices": {
            "['x1-x2']": [
              "SO8",
              288
            ]
          },
          "N_O3": 4,
          "chi_O7_total": 36,
          "chi_D7": 288,
          "N_D3": 20.0,
          "N_flux_frac": 0.0,
          "N_D3_mobile_max": 20,
          "chi_Y4_correction_O3": 24,
          "chi_Y4_naive": 216,
          "chi_Y4_corrected": 240,
          "nodal_singularities": true,
          "n_S_df0": 30,
          "ed3_instantons": [
            {
              "divisor": "x6",
              "h10": 0,
              "h20": 0,
              "chi": 3,
              "h11": 1,
              "n_O3": 1,
              "chi_correction_vertical_divisor": 2,
              "chi_vertical_naive": 12,
              "chi_vertical_corrected": 14
            }
          ]
        }
      ],
      "OPLANES": [
        {
          "OIDEAL": [
            "x1+x2",
            "x3",
            "x5"
          ],
          "ODIM": 3,
          "chi": 3
        },
        {
          "OIDEAL": [
            "x1+x2",
            "x3",
            "x6"
          ],
          "ODIM": 3,
          "chi": 1
        },
        {
          "OIDEAL": [
            "x1-x2"
          ],
          "ODIM": 7,
          "chi": 36,
          "D3": 0,
          "c2D": 36,
          "cohomologies": {
            "h10": 0,
            "h20": 2,
            "h11": 30
          },
          "rigid": false,
          "completely_rigid": false
        }
      ]
    },
    "... (other orientifolds omitted) ..."
  ]
}
```
