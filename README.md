# Orientifold Database Generator for Calabi-Yau Threefolds in Toric Ambient Spaces

This repository provides tools to generate databases of orientifold configurations for Calabi--Yau threefolds. Using `cytools` or datasets from Hugging Face to retrieve Calabi--Yau threefolds, we systematically enumerate discrete symmetries and compute the geometric data required for string theory models.

This repository accompanies the work presented in [arXiv:2609.35958](https://arxiv.org/abs/2609.35958). The computations rely on [Solver Agent](https://github.com/starrfree/solver-agent) (AGPL-3.0), a large language model-based framework for proofs and computations in mathematics and theoretical physics. The system splits the solution process, delegating tasks to specialized sub-agents with access to external symbolic, numerical, and domain-specific algorithms. The framework tracks this process through a persistent ledger that records assumptions, derivations, and computations, while independent agents verify intermediate steps and final results. To ensure the reproducibility of the computer-assisted calculations, the system exports exact sessions as a raw `ledger.json`, a human-readable `report.html`, and a printable `report.tex`. Note that executing Solver Agent requires an external API subject to usage fees.

## Prerequisites

We require the `cytools` environment.

Computing prime toric divisor cohomologies requires `cohomcalg` to be installed and accessible via the `PATH` environment variable, or by setting `COHOMCALG_PATH` explicitly.

## Input Databases and Theoretical Background

The algebraic torus $T$ defining the ambient toric variety specifies a lattice $N = \mathrm{Hom}(\mathbb{C}^*, T) \simeq \mathbb{Z}^d$ and a dual character lattice $M = \mathrm{Hom}(T, \mathbb{C}^*) \simeq \mathbb{Z}^d$, alongside their real vector spaces $N_\mathbb{R} = N\otimes_\mathbb{Z} \mathbb{R}$ and $M_\mathbb{R} = M\otimes_\mathbb{Z} \mathbb{R}$. The Calabi--Yau geometry is determined by a polytope $\Delta\subset M_\mathbb{R}$ and its dual $\Delta^\circ\subset N_\mathbb{R}$. The geometric data of the ambient space is captured by the normal fan $\Sigma_\Delta \subset N_\mathbb{R}$, constructed from the cones over the proper faces of $\Delta^\circ$. This database contains the fundamental three-dimensional lattice polytopes $\Delta_3$.

Calabi--Yau threefold hypersurfaces in toric ambient spaces are constructed via the correspondence between toric varieties and reflexive lattice polytopes. For a reflexive polytope $\Delta$ and its dual $\Delta^\circ$, the normal fan of $\Delta^\circ$ defines a toric fourfold containing a Calabi--Yau hypersurface. This hypersurface is resolved by refining the fan through a Fine, Regular, Star Triangulation (FRST) of $\Delta^\circ$. The Calabi--Yau is defined as the zero locus of a generic anticanonical polynomial, with its topological invariants—such as Hodge numbers and intersection numbers—extracted directly from the polytope combinatorics and the gauged linear sigma model (GLSM) charges of the toric divisors.

To uplift a given Calabi--Yau threefold $X_3$ to an elliptically fibered Calabi--Yau fourfold $\pi: Y_4\to B_3$, a Weierstrass model is constructed over the orientifold quotient base $B_3= X_3/\iota$. A controlled class of such Calabi--Yau orientifolds is obtained through the trilayer polytope construction. Here, $X_3$ is realized as a hypersurface in a four-dimensional ambient toric variety determined by a triangulation of a reflexive lattice polytope $\Delta_\text{tr}$. This four-dimensional polytope is systematically constructed from a three-dimensional lattice polytope $\Delta_3$ via 

$$\Delta_\text{tr} = \text{conv}\left(\{(v,1)\mid v\in\text{Vert}(\Delta_3)\} \cup\{(0,-1)\}\right)\,.$$


### `databases/3dpoly.json`

This database is obtained from the original [Kreuzer-Skarke database](https://hep.itp.tuwien.ac.at/~kreuzer/CY/CYcy.html) of 3d reflexive polytopes and is provided as pre-compiled JSON files.

**Example entry (first component):**
```json
{
  "id": 1,
  "dimension": 3,
  "num_vertices": 4,
  "M_lattice": {
    "points": 5,
    "vertices": 4
  },
  "N_lattice": {
    "points": 35,
    "vertices": 4
  },
  "Picard_number": 19,
  "Correction_term": 0,
  "coordinates": [
    [
      1,
      0,
      0,
      -1
    ],
    [
      0,
      1,
      0,
      -1
    ],
    [
      0,
      0,
      1,
      -1
    ]
  ]
}
```

### `databases/trilayer.json`

This database contains the four-dimensional reflexive polytopes $\Delta_\text{tr}$. These are obtained by taking the three-dimensional polytopes $\Delta_3$ from the `3dpoly.json` database and applying the trilayer construction.

**Example entry (first component):**
```json
{
  "id": 1,
  "dimension": 4,
  "num_vertices": 5,
  "M_lattice": {
    "points": 201,
    "vertices": 5
  },
  "N_lattice": {
    "points": 7,
    "vertices": 5
  },
  "Picard_number_N": 1,
  "Picard_number_M": 149,
  "coordinates": [
    [
      -1,
      0,
      0,
      0,
      1
    ],
    [
      -1,
      0,
      0,
      1,
      0
    ],
    [
      -1,
      0,
      1,
      0,
      0
    ],
    [
      1,
      -1,
      1,
      1,
      1
    ]
  ],
  "cy_is_favorable": true
}
```

## Scripts

The following scripts to classify and store orientifold configurations for Calabi--Yau threefolds in toric ambient spaces have been written mostly based on the work by:

1. **Andres Collinucci, Frederik Denef, Mboyo Esole**, *"D-brane Deconstructions in IIB Orientifolds"*. [arXiv:0805.1573](https://arxiv.org/abs/0805.1573)
2. **Ralph Blumenhagen, Volker Braun, Thomas W. Grimm, Timo Weigand**, *"GUTs in Type IIB Orientifold Compactifications"*. [arXiv:0811.2936](https://arxiv.org/abs/0811.2936)
3. **Andres Collinucci**, *"New F-theory lifts"*. [arXiv:0812.0175](https://arxiv.org/abs/0812.0175)
4. **Andres Collinucci**, *"New F-theory lifts II: Permutation orientifolds and enhanced singularities"*. [arXiv:0906.0003](https://arxiv.org/abs/0906.0003)
5. **Chiara Crinò, Fernando Quevedo, Andreas Schachner, Roberto Valandro**, *"A Database of Calabi-Yau Orientifolds and the Size of D3-Tadpoles"*. [arXiv:2204.13115](https://arxiv.org/abs/2204.13115) | [GitHub Repository](https://github.com/AndreasSchachner/CY_Orientifold_database)
6. **Patrick Jefferson, Manki Kim**, *"On the intermediate Jacobian of M5-branes"*. [arXiv:2211.00210](https://arxiv.org/abs/2211.00210)
7. **Jakob Moritz**, *"Orientifolding Kreuzer-Skarke"*. [arXiv:2305.06363](https://arxiv.org/abs/2305.06363)
8. **Bjoern Hassfeld, Jakob Moritz**, *"Calabi-Yau Orientifold Hypersurfaces and their F-theory Uplifts"*. [arXiv:2606.19423](https://arxiv.org/abs/2606.19423) | [GitHub Repository](https://github.com/B-Hassfeld/CYorientifolds_FTheoryUplifts)

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

### 3. `generate_trilayer_orientifolds.py`

This script computes the orientifold data for the geometries contained in the `databases/trilayer.json` database. These trilayer polytopes define Calabi--Yau threefold hypersurfaces embedded in a four-dimensional ambient toric variety. For cases where the origin is the unique lattice point with a vanishing fourth coordinate in the four-dimensional polytope, the orientifold quotient $B_3 = X_3/\iota$ is a toric threefold. 

The F-theory uplift is constructed directly as a generic Weierstrass model over $B_3$. When $B_3$ contains rigid toric divisors intersecting the O7-plane, non-Higgsable clusters (NHCs) arise, forcing an $I_0^*$ degeneration of the elliptic fiber. Resolving these singularities requires crepant toric blowups in the ambient space. The script parses the input geometries, extracts the allowed orientifold involutions, and evaluates the corresponding D3-brane tadpole cancellation conditions. The output is stored in `databases/trilayer_orientifolds.json`.

#### Usage:

You can use the built-in help flag to see all options:
```bash
conda run -n cytools python generate_trilayer_orientifolds.py -h
```

**Output:**
```
usage: generate_trilayer_orientifolds.py [-h] [--compute_cohomologies]
                                         [--workers WORKERS]
                                         [--h11_max H11_MAX]

options:
  -h, --help            show this help message and exit
  --compute_cohomologies
                        Run cohomcalg for O7 planes to compute divisor cohomologies.
  --workers WORKERS     Number of workers for multiprocessing.
  --h11_max H11_MAX     Maximum h11 (Picard_number_N) to process.
```

### 4. `read_orientifold_db.py`

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

### 5. `verify_all_intB_relations.py`

This script verifies the topological intersection relations between the Calabi--Yau threefold $X_3$ and its orientifold quotient base $B_3$. In the trilayer construction, the Calabi--Yau threefold $X_3$ embeds as a bisection of an auxiliary genus one fibered Calabi--Yau fourfold $\widetilde{Y}_4$ defined over $B_3$. Applying the adjunction formula to this geometry evaluates the topological invariants of $X_3$ directly from the intersection data of $B_3$. This script requires the `databases/3dpoly.json` and `databases/trilayer.json` databases to run.

Using `cytools`, the script constructs the corresponding toric varieties and verifies the following identities:

1. **Triple Intersections**:
   $$\int_{X_3} k_a \wedge k_b \wedge k_c = 2 \kappa_{abc} = 2 \int_{B_3} j_a \wedge j_b \wedge j_c$$
2. **Second Chern Class Contractions**:
   $$\int_{X_3} c_2(TX_3) \wedge k_a = 2 \int_{B_3} j_a \wedge \left(c_1^2(TB_3) + c_2(TB_3)\right)$$
3. **Euler Characteristics**:
   $$\chi(X_3) = 2 \int_{B_3} \left(c_3(TB_3) - c_1(TB_3)c_2(TB_3) - 2c_1(TB_3)^3\right)$$

Here, $k_a$ are the K\"ahler cone generators on $X_3$ induced by its embedding into $\widetilde{Y}_4$, $\kappa_{abc}$ are the intersection numbers on $B_3$, and $j_a$ are the K\"ahler cone generators on the base $B_3$.

#### Usage:

This script runs automatically across the entire configured database without requiring any command-line arguments.

```bash
conda run -n cytools python verify_all_intB_relations.py
```

### Geometric admissibility (base locus check)

The script enforces geometric admissibility to prevent base locus singularities on the Calabi--Yau hypersurface. We compute the monomial exponents $E_{m, j} = \langle m, r_j \rangle + 1$ for all points $m \in \Delta^\circ$ and rays $r_j \in \Delta$. For each sign configuration, the script identifies invariant monomials matching the canonical sign $S_{\text{canon}}$, and evaluates these monomials against the coordinate intersections permitted by the Stanley--Reisner ideal. If all invariant monomials vanish on an allowed toric stratum, the orientifold configuration is rejected and excluded from the database.

## Understanding the JSON Output

The generated JSON file is a list of dictionaries, where each entry represents a favorable Calabi--Yau threefold and contains:

- `POLYID`: The unique identifier for the polytope.
- `h11`: The Hodge number $h^{1,1}$.
- `h21`: The Hodge number $h^{2,1}$.
- `polytope_points`: Points in the $N$-lattice defining the dual polytope.
- `triangulations`: Geometric data for triangulations, including intersections, divisor basis, and the second Chern class.
- `orientifolds`: A list of valid involutions alongside fixed planes (O3, O5, O7, O9) and D3/D7 tadpole cancellation properties.

Each orientifold entry explicitly contains the following fields:

### Basic Data & Involution
- `TRIANGN`: The index of the triangulation.
- `INVOLN`: The index of the involution.
- `INVOL`: A string representation of the coordinate permutation $x_i \to x_{P(i)}$.
- `SCANON`: The canonical sign $S_{\text{canon}} \in \{1, -1\}$ chosen for the orientifold action.
- `BCANON`: A list of signs $b_{\text{canon}} \in \{1, -1\}$ for each homogeneous coordinate, defining the diagonal sign shift of the canonical representative involution orbit.
- `SOMEGA`: The orientifold projection choice $S_O$. Mathematically, this is defined as:
  $$S_O = S_{\text{canon}} \cdot \text{sign}(P) \cdot \prod b_{\text{canon}}$$
  
  Setting $S_O = +1$ corresponds to the O3/O7-plane configurations, while $S_O = -1$ defines O5/O9-planes. In the code, this is computed as `S_Omega = S_canon * sign_perm * prod_b_canon`, where `sign_perm` is the parity of the coordinate permutation $P$ and `prod_b_canon` is the product of signs acquired by the defining polynomials.
- `is_reducible`: A boolean flag indicating whether the Calabi--Yau vanishes identically under the orientifold equations: `S_canon * parity_lam == -1`, where `parity_lam = (-1) ** ((h11 - trace) // 2)` is the determinant of the induced action on the $H^{1,1}$ cohomology.
- `integer_hodge`: A boolean indicating whether the equivariant Hodge numbers are strictly integers.

### Equivariant Hodge Numbers
The orientifold projection splits the cohomology groups $H^{p,q}(X_3)$ into even and odd eigenspaces under the pull-back of the involution, decomposing the K\"ahler and complex structure moduli spaces into invariant and anti-invariant sectors. The dimensions of these eigenspaces are determined by:

$$h^{1,1}_{\pm} = \frac{1}{2} \left(h^{1,1} \pm \text{Tr}(\Lambda)\right)$$

$$h^{2,1}_+ = \frac{1}{2} \left( h^{2,1} + 1 - S_O + h^{1,1}_+ - h^{1,1}_- - \frac{1}{2} \chi_{\text{fix}} \right)$$

$$h^{2,1}_- = h^{2,1} - h^{2,1}_+$$
The parameter $S_O$ (labeled `SOMEGA` in the JSON) encodes the choice of the orientifold projection; setting $S_O = +1$ corresponds to O3/O7-plane configurations, while $S_O = -1$ defines O5/O9-planes. The total Euler characteristic of the fixed loci, $\chi_{\text{fix}}$, decomposes into the sum of the Euler characteristics of the codimension-one fixed surfaces and the number of isolated fixed points.

In the code, this decomposition is evaluated by projecting the coordinate permutation matrix $P$ onto the Picard lattice using the GLSM charge matrix $Q$. This defines $\Lambda = Q P Q^{-1}$, whose trace determines the K\"ahler moduli splitting:
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
The `tadpole_data` array contains the physical quantities governing the D3/D7 tadpole cancellation conditions. The evaluation branches depending on the geometric resolution of the O7-plane tadpole: recombination into a single Whitney umbrella (`WU`) or placement of an $\text{SO}(8)$ stack.

For a specific combination of O7-plane configurations, the JSON logs:
- `choices`: A dictionary mapping the O7-plane ideals to their specified configuration (`WU` or `SO8`) and the resulting D7-brane Euler characteristic $\chi(S_{D7})$.
- `N_O3`: The number of isolated O3-planes.
- `chi_O7_total`: The sum of the Euler characteristics of the O7-planes.
- `chi_D7`: The total D7-brane Euler characteristic for this configuration.

**D3-Brane Charges:**
The net D3-brane charge $N_{D3}$ induced by the localized sources is evaluated in the double cover $X_3$:
$$N_{D3} = \frac{N_{O3}}{2} + \frac{\chi(S)}{6} + \frac{\chi(S_{D7})}{24}$$
- `N_D3`: The fractional D3-charge budget.
- `N_flux_frac`: The fractional part of the gauge flux contribution.
- `N_D3_mobile_max`: The maximum integer number of mobile D3-branes.

For a generic smooth divisor $S_{D7}$, the Euler characteristic is computed by integrating the top Chern class via the adjunction formula:
$$\chi(S_{D7}) = \int_{X_3} \left( [S_{D7}]^3 + c_2(TX_3) \wedge [S_{D7}] \right)$$

If an O7-plane wraps a rigid divisor, the tadpole is canceled locally by an $\text{SO}(8)$ stack—four D7-branes and their orientifold images placed exactly on the O7-plane locus $S$. Their total Euler characteristic is eight times that of $S$:
$$\chi(S_{D7_{\mathrm{SO}(8)}}) = 8 \chi(S) = \int_{X_3} \left( 8 [S]^3 + 8 c_2(TX_3) \wedge [S] \right)$$

If the wrapped divisor admits complex structure deformations, the branes can recombine into a Whitney umbrella wrapping the homology class $8[S]$. The singular locus of the umbrella shifts the Euler characteristic:
$$\chi(S_{D7_{\text{WU}}}) = \int_{X_3} \left( 344 [S]^3 + 8 c_2(TX_3) \wedge [S] \right)$$

In the code, the $\text{SO}(8)$ stack contributes $8 \cdot (\text{c2D\_ op} + \text{D3\_ op})$, whereas the Whitney umbrella contributes $8 \cdot \text{c2D\_ op} + 344 \cdot \text{D3\_ op}$.

**Fourfold Euler Characteristics:**
The D3-brane tadpole evaluated in the F-theory uplift requires the "stringy" Euler characteristic of the Calabi--Yau fourfold, $\chi_{\text{st}}(Y_4)$, which includes stringy corrections from terminal $\mathbb{Z}_2$ quotient singularities induced by O3-planes:
- `chi_Y4_corrected`: The stringy Euler characteristic $\chi_{\text{st}}(Y_4) = 24 \cdot \frac{N_{D3}}{2}$.
- `chi_Y4_correction_O3`: The stringy correction $6 N_{O3}$.
- `chi_Y4_naive`: The naive integral of the top Chern class, derived by subtracting $6 N_{O3}$ from $\chi_{\text{st}}(Y_4)$.
- `nodal_singularities`: A boolean evaluating whether $n^S_{df=0} \neq 0$.
- `n_S_df0`: The singularity count $n^S_{df=0} = \int_{O7} (c_2(S) - c_1(S)^2)$.

**ED3 Instanton Divisors:**
To generate a non-perturbative superpotential to stabilize K\"ahler moduli, an ED3-instanton must carry exactly two neutral fermionic zero modes. Complete geometric rigidity ($h^{1,0}(D) = h^{2,0}(D) = 0$) guarantees the absence of additional deformation moduli or Wilson lines.
- `ed3_instantons`: Evaluates completely rigid divisors wrapping ED3-instantons. This includes their Hodge numbers (`h10`, `h20`, `h11`) and the stringy vertical divisor Euler characteristic $\chi(\overline{D}) = \chi_{\text{naive}}(\overline{D}) + 6N_{O3}(\widehat{D})$.

### O-Planes (`OPLANES`)
The list of fixed O-planes resulting from the orientifold equations. Each entry contains:
- `OIDEAL`: A list of coordinate constraints setting the fixed locus.
- `ODIM`: The dimension of the O-plane, determined by the number of constraint equations relative to the ambient dimension.
- `chi`: The computed Euler characteristic of the fixed plane.

For an O7-plane located at divisor $S \equiv D_{O7}$, the geometric quantities are evaluated via sparse intersection numbers ($\kappa_{ijk}$) and the second Chern class ($c_{2,i}$):
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
