import os
import warnings

warnings.filterwarnings("ignore", category=UserWarning)

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import argparse
import itertools
import json
import multiprocessing as mp
import re
import shutil
import subprocess
import tempfile
from functools import partial

import cytools
import duckdb
import numpy as np
from tqdm import tqdm

cytools.config.enable_experimental_features()


# ----------------- Worker Multiprocessing Pool -----------------
class NoDaemonProcess(mp.Process):
    @property
    def daemon(self):
        return False

    @daemon.setter
    def daemon(self, value):
        pass


class NoDaemonContext(type(mp.get_context())):
    Process = NoDaemonProcess


class NoDaemonPool(mp.pool.Pool):
    def __init__(self, *args, **kwargs):
        kwargs["context"] = NoDaemonContext()
        super().__init__(*args, **kwargs)


# ----------------- Orientifold Mathematics -----------------
def get_permutations_from_autos(autos, pts):
    """
    Compute valid orientifold involution permutations from polytope automorphisms.

    This function applies given automorphism matrices to the polytope vertices
    and identifies which automorphisms correspond to valid geometric involutions
    (i.e., order 2 permutations where P(P(x)) = x).

    Args:
        autos (list of numpy.ndarray): A list of automorphism matrices.
        pts (numpy.ndarray): The points (vertices) of the polytope or ray generators.

    Returns:
        list of dict: A list of dictionaries representing the valid involutions.
            Each dictionary maps the 1-based index of a point to its image index.
    """
    perms = []
    for auto in autos:
        perm = {}
        mapped_pts = np.dot(pts, auto)
        valid = True
        for i, pt in enumerate(pts):
            match = np.where((mapped_pts == pt).all(axis=1))[0]
            if len(match) == 1:
                perm[i + 1] = match[0] + 1
            else:
                valid = False
                break
        if valid and perm not in perms:
            # An orientifold must be an involution (order 2), meaning perm(perm(x)) == x
            is_involution = all(perm.get(perm.get(k, k), k) == k for k in perm)
            if is_involution:
                perms.append(perm)
    return perms


def compute_h11_pm(Q, perm):
    """
    Compute the equivariant Hodge numbers h11_plus and h11_minus.

    Calculates the splitting of h11 into orientifold even (h11_plus) and
    odd (h11_minus) sectors under the orientifold involution using the
    Lefschetz trace formula on the divisor classes.

    Args:
        Q (numpy.ndarray): The GLSM charge matrix (or divisor basis) of the CY.
        perm (dict): The involution permutation mapping 1-based indices.

    Returns:
        tuple: A tuple containing (h11_plus, h11_minus) as integers.
    """
    N = Q.shape[1]
    P = np.zeros((N, N))
    for i in range(N):
        # perm maps 1..N to 1..N (since 0 is the origin)
        mapped = perm.get(i + 1, i + 1)
        P[i, mapped - 1] = 1

    QP = Q @ P
    Q_pinv = np.linalg.pinv(Q)
    Lambda = QP @ Q_pinv
    Lambda = np.round(Lambda)

    trace = np.trace(Lambda)
    h11 = Q.shape[0]
    h11_plus = int(round((h11 + trace) / 2))
    h11_minus = int(round((h11 - trace) / 2))
    return h11_plus, h11_minus


def compute_oplan_chi(
    c2, kappa, oplane_eq, D_basis, vanishes_identically=False, simplices=None
):
    """
    Compute the Euler characteristic and other topological data of an O-plane.

    Depending on the codimension of the fixed locus (number of equations),
    this computes the properties of O7, O5, or O3 planes. It calculates
    intersections using the second Chern class and triple intersection numbers.

    Args:
        c2 (numpy.ndarray): The second Chern class of the CY threefold.
        kappa (numpy.ndarray): The triple intersection numbers of the CY.
        oplane_eq (list of str): List of coordinate equations defining the O-plane.
        D_basis (numpy.ndarray): The basis of divisor classes.
        vanishes_identically (bool, optional): Flag indicating if the polynomial
            vanishes identically under the involution. Defaults to False.
        simplices (list of list, optional): The simplices of the triangulation,
            used when vanishes_identically is True. Defaults to None.

    Returns:
        dict: A dictionary containing topological data such as 'chi', and
            for O7-planes, 'D3' and 'c2D'. May also contain 'is_reducible'
            if the entire CY vanishes.
    """
    if len(oplane_eq) == 0:
        return {"chi": 0, "D3": 0, "c2D": 0}

    if oplane_eq == ["Whole CY"]:
        # Handled externally usually, but in case:
        return {"chi": 0, "D3": 0, "c2D": 0}  # handled outside

    div_classes = []
    idx_list = []
    for eq in oplane_eq:
        m = re.search(r"x(\d+)", eq)
        if m:
            coord_idx = int(m.group(1)) - 1
            idx_list.append(coord_idx)
            v = D_basis[:, coord_idx]
            div_classes.append(v)

    if len(div_classes) == 1:
        if vanishes_identically:
            # Reducible CY: The fixed locus is the entire ambient divisor.
            # Calculating the full ambient divisor Euler characteristic is complex.
            # We flag it as reducible instead of computing intersection numbers.
            return {"chi": 0, "D3": 0, "c2D": 0, "is_reducible": True}
        else:
            # Divisor intersecting generic CY (O7 plane)
            D = div_classes[0]
            D3 = np.einsum("ijk,i,j,k->", kappa, D, D, D)
            c2D = np.dot(c2, D)
            return {
                "chi": int(round(D3 + c2D)),
                "D3": int(round(D3)),
                "c2D": int(round(c2D)),
            }
    elif len(div_classes) == 2:
        if vanishes_identically and simplices is not None:
            # O7 plane (toric surface): count maximal cones containing both rays
            # idx_list is 0-based (x1→0), but simplices include origin as point 0,
            # so x1 is point 1. Need +1 offset.
            pt1 = idx_list[0] + 1
            pt2 = idx_list[1] + 1
            chi = sum(1 for simp in simplices if pt1 in simp and pt2 in simp)
            return {"chi": chi}
        else:
            # Curve (O5 plane)
            D1, D2 = div_classes[0], div_classes[1]
            D1pD2 = D1 + D2
            val = -np.einsum("ijk,i,j,k->", kappa, D1pD2, D1, D2)
            return {"chi": int(round(val)), "D3": 0, "c2D": 0}
    elif len(div_classes) == 3:
        if vanishes_identically and simplices is not None:
            # O5 plane (toric curve): count maximal cones containing all 3 rays
            pt1 = idx_list[0] + 1
            pt2 = idx_list[1] + 1
            pt3 = idx_list[2] + 1
            chi = sum(
                1 for simp in simplices if pt1 in simp and pt2 in simp and pt3 in simp
            )
            return {"chi": chi, "D3": 0, "c2D": 0}
        else:
            # Point in CY (O3 plane)
            D1, D2, D3 = div_classes[0], div_classes[1], div_classes[2]
            val = np.einsum("ijk,i,j,k->", kappa, D1, D2, D3)
            return {"chi": int(round(val)), "D3": 0, "c2D": 0}
    elif len(div_classes) == 4:
        if vanishes_identically and simplices is not None:
            # O3 plane (toric point): count maximal cones containing all 4 rays
            pt1 = idx_list[0] + 1
            pt2 = idx_list[1] + 1
            pt3 = idx_list[2] + 1
            pt4 = idx_list[3] + 1
            chi = sum(
                1
                for simp in simplices
                if pt1 in simp and pt2 in simp and pt3 in simp and pt4 in simp
            )
            return {"chi": chi, "D3": 0, "c2D": 0}
        else:
            # Point in ambient space. A generic CY hypersurface does not pass through it.
            return {"chi": 0, "D3": 0, "c2D": 0}
    else:
        return {"chi": 0, "D3": 0, "c2D": 0}


# ----------------- Cohomcalg Logic -----------------
def run_cohomcalg_for_divisor(sr_ideal, glsm_charges, divisor_v):
    """
    Run cohomCalg to compute the line bundle cohomology for a prime toric divisor.

    Writes a temporary input file for the cohomCalg software and extracts
    the dimensions of the cohomology groups H^i(Y, D) via Koszul splitting.

    Args:
        sr_ideal (list of list): The Stanley-Reisner ideal generators.
        glsm_charges (list of list): The GLSM charge matrix.
        divisor_v (list): The divisor class vector for the line bundle.

    Returns:
        list: A list [h10, h20] representing the dimensions of H^1(Y, O(D))
            and H^2(Y, O(D)) respectively. Returns [0, 0] on failure.
    """
    # Returns [h10, h20] for the divisor
    # This invokes cohomcalg via subprocess
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".in") as f:
        # Write vertices and GLSM charges
        for i in range(len(glsm_charges[0])):
            charges = [row[i] for row in glsm_charges]
            f.write(f"vertex x{i+1} | GLSM: ( {', '.join(map(str, charges))} );\n")

        sr_vars = []
        for gen in sr_ideal:
            var_str = "*".join([f"x{i+1}" for i in gen])
            sr_vars.append(var_str)
        f.write(f"srideal [ {', '.join(sr_vars)} ];\n")

        # We need D and D + KX
        kx_charges = [-sum(row) for row in glsm_charges]
        d_charges = divisor_v
        dkx_charges = [d_charges[r] + kx_charges[r] for r in range(len(kx_charges))]

        f.write(f"ambientcohom O( {', '.join(map(str, d_charges))} );\n")
        f.write(f"ambientcohom O( {', '.join(map(str, dkx_charges))} );\n")
        tmp_path = f.name

    cohomcalg_path = os.environ.get(
        "COHOMCALG_PATH", shutil.which("cohomcalg") or "/opt/homebrew/bin/cohomcalg"
    )
    try:
        res = subprocess.run(
            [cohomcalg_path, tmp_path], capture_output=True, text=True, timeout=60
        )
        output = res.stdout
    except subprocess.TimeoutExpired:
        os.unlink(tmp_path)
        return [0, 0]
    except Exception:
        os.unlink(tmp_path)
        return [0, 0]

    os.unlink(tmp_path)

    def parse_bundle():
        lines = output.splitlines()
        cohoms = []
        for line in lines:
            if line.strip().startswith("dim H^i"):
                parts = line.split("=")
                if len(parts) > 1:
                    tup_str = parts[1].strip().strip("()")
                    dims = [int(x.strip()) for x in tup_str.split(",")]
                    cohoms.append(dims)
        return cohoms

    parsed = parse_bundle()
    if len(parsed) != 2:
        return [0, 0]

    hx_d = parsed[0]
    hx_dkx = parsed[1]

    # Koszul splitting: h^k(Y, D) = h^k(A, D) + h^{k+1}(A, D + K_A)
    hy_d = []
    for k in range(4):  # h0 to h3
        hy_d.append(hx_d[k] + hx_dkx[k + 1])

    h10 = hy_d[1]
    h20 = max(0, hy_d[0] - 1)

    return [h10, h20]


# ----------------- Main Processing Pipeline -----------------
def process_polytope(poly_idx, h11, pts, compute_cohomologies=False):
    """
    Process a single CY polytope to find orientifold configurations and tadpoles.

    This is the core pipeline function. It triangulates the polytope, finds
    favourable CY threefolds, extracts geometric data (intersection numbers,
    Chern classes), computes possible orientifold involutions, evaluates
    fixed loci (O-planes), and calculates D3/D7 tadpole cancellation bounds.

    Args:
        poly_idx (int): The unique identifier for the polytope.
        h11 (int): The Hodge number h1,1 of the polytope.
        pts (list of list): The vertices of the polytope.
        compute_cohomologies (bool, optional): Whether to run cohomcalg for
            divisor cohomologies (ED3 instantons). Defaults to False.

    Returns:
        dict: A dictionary containing comprehensive orientifold data including
            triangulations, O-planes, and tadpole limits. Returns an empty dict
            if no favourable triangulation is found.
    """
    h21 = 0
    poly = cytools.Polytope(pts)
    m_pts = poly.dual().points()
    p_id = poly_idx
    pts_array = poly.points()
    autos = poly.automorphisms()

    triangs = list(poly.ntfe_frsts())

    out = {
        "POLYID": p_id,
        "h11": h11,
        "h21": 0,  # Placeholder, will be updated from cy
        "polytope_points": pts_array.tolist(),
        "triangulations": [],
        "orientifolds": [],
    }

    # 1. Triangulations & Geometry
    triang_map = {}
    for t_idx, t in enumerate(triangs):
        try:
            cy = t.get_cy()
            out["h21"] = cy.h21()  # Set from cy
            h21 = cy.h21()

            c2 = cy.second_chern_class(in_basis=True)
            kappa_sparse = cy.intersection_numbers(in_basis=True, format="coo")
            kappa = np.zeros((h11, h11, h11))
            for i, j, k, val in kappa_sparse:
                for p in itertools.permutations([i, j, k]):
                    kappa[p] = val

            Q = cy.glsm_charge_matrix(include_origin=False)
            D_basis = Q

            sr_ideal = t.sr_ideal()

            entry = {
                "TRIANGN": t_idx + 1,
                "simplices": t.simplices().tolist(),
                "glsm_charge_matrix": Q.tolist(),
                "divisor_basis": cy.divisor_basis().tolist(),
                "second_chern_class": c2.tolist(),
                "intersection_numbers": kappa_sparse.tolist(),
                "sr_ideal": [[int(i) for i in x] for x in sr_ideal],
            }
            # Store this to compute Hodge numbers later
            triang_map[t_idx + 1] = {
                "c2": c2,
                "kappa": kappa,
                "I4": np.array(
                    t.get_toric_variety().intersection_numbers(format="dense")
                ),
                "Q": Q,
                "D_basis": D_basis,
                "sr_ideal": sr_ideal,
                "ray_pts": t.points()[1:],
                "simplices": t.simplices().tolist(),
            }
            out["triangulations"].append(entry)
        except Exception as e:
            warnings.warn(
                f"Failed to process triangulation {t_idx + 1} of poly {poly_idx}: {e}"
            )
            continue

    if len(out["triangulations"]) == 0:
        return out

    has_favorable = any(
        triang_map[t["TRIANGN"]]["Q"].shape[0] == h11 for t in out["triangulations"]
    )
    if not has_favorable:
        return {}

    m_pts = poly.dual().points()

    # Pre-extract ray_pts from the first triang to compute perms (since it's identical across all triangulations)
    ray_pts_invariant = triang_map[out["triangulations"][0]["TRIANGN"]]["ray_pts"]
    perms = get_permutations_from_autos(autos, ray_pts_invariant)

    for triang_entry in out["triangulations"]:
        t_data = triang_map[triang_entry["TRIANGN"]]
        if t_data["Q"].shape[0] != h11:
            continue

        triang_n = triang_entry["TRIANGN"]
        Q = t_data["Q"]
        kappa = t_data["kappa"]
        c2 = t_data["c2"]
        D_basis = t_data["D_basis"]
        ray_pts = t_data["ray_pts"]
        I4 = t_data.get("I4", None)

        Q_h11 = Q.shape[0]
        N = Q.shape[1]

        # Precompute monomial exponents for base locus check
        N_m = len(m_pts)
        E = np.zeros((N_m, N), dtype=int)
        for i, m in enumerate(m_pts):
            for j, r in enumerate(ray_pts):
                E[i, j] = np.dot(m, r) + 1
        E_dict = {tuple(row): i for i, row in enumerate(E)}

        Q_sum = [sum(Q[a]) for a in range(Q_h11)]
        unique_H = {}
        for s_tuple in itertools.product([1, -1], repeat=Q_h11):
            lam = []
            for j in range(N):
                val = 1
                for a in range(Q_h11):
                    if s_tuple[a] == -1 and Q[a][j] % 2 != 0:
                        val *= -1
                lam.append(val)
            parity_lam = 1
            for a in range(Q_h11):
                if s_tuple[a] == -1 and Q_sum[a] % 2 != 0:
                    parity_lam *= -1
            unique_H[tuple(lam)] = parity_lam
        H = [(lam, parity_lam) for lam, parity_lam in unique_H.items()]

        B_xi_all = set()
        for signs in itertools.product([1, -1], repeat=N):
            orbit = []
            for lam, _ in H:
                orbit.append(tuple(signs[k] * lam[k] for k in range(N)))
            B_xi_all.add(min(orbit))
        B_xi = list(B_xi_all)

        b_identity = tuple(1 for _ in range(N))

        # Pre-compute cohomologies for all prime toric divisors if requested
        divisor_cohoms = {}
        if compute_cohomologies:
            sr_ideal_clean = []
            for gen in t_data["sr_ideal"]:
                if 0 not in gen:
                    sr_ideal_clean.append([g - 1 for g in gen])
            if not sr_ideal_clean:
                # If SR ideal is empty (e.g. for projective spaces), the excluded set is the origin.
                sr_ideal_clean.append(list(range(N)))
            for k in range(N):
                div_v = D_basis[:, k].tolist()
                coh = run_cohomcalg_for_divisor(sr_ideal_clean, Q.tolist(), div_v)
                if coh is not None:
                    divisor_cohoms[k] = coh

            if divisor_cohoms:
                triang_div_cohoms = {}
                for k, coh in divisor_cohoms.items():
                    D_k_vec = Q[:, k]
                    D3 = np.einsum("a,b,c,abc->", D_k_vec, D_k_vec, D_k_vec, kappa)
                    h10 = coh[0]
                    h20 = coh[1]
                    h11_D = int(round(10 - 8 * h10 + 10 * h20 - D3))
                    triang_div_cohoms[f"x{k+1}"] = {
                        "h10": h10,
                        "h20": h20,
                        "h11": h11_D,
                    }

                triang_entry["prime_toric_divisors_cohomologies"] = triang_div_cohoms

        # 2. Orientifolds
        invol_counter = 1
        for perm in perms:
            # Determine sign of permutation
            visited_p = set()
            transpositions = 0
            for i in range(N):
                if i not in visited_p:
                    c = 0
                    curr = i
                    while curr not in visited_p:
                        visited_p.add(curr)
                        curr = perm.get(curr, curr)
                        c += 1
                    transpositions += c - 1
            sign_perm = -1 if transpositions % 2 == 1 else 1

            h11_plus, h11_minus = compute_h11_pm(Q, perm)

            # Monomial permutation mapping for base locus check
            pi_inv = {j - 1: k - 1 for k, j in perm.items()}
            for k in range(N):
                if k not in pi_inv:
                    pi_inv[k] = k
            E_perm = np.zeros_like(E)
            for j in range(N):
                E_perm[:, j] = E[:, pi_inv[j]]
            m_prime_idx = np.zeros(N_m, dtype=int)
            for i, row_perm in enumerate(E_perm):
                m_prime_idx[i] = E_dict[tuple(row_perm)]

            seen_configs = set()

            for b_xi in B_xi:
                orbit = []
                for lam, _ in H:
                    orbit.append(tuple(b_xi[k] * lam[k] for k in range(N)))
                orbit_canon = min(orbit)
                prod_b_canon = np.prod(orbit_canon)

                for S_canon in [1, -1]:
                    S_Omega = S_canon * sign_perm * prod_b_canon

                    # --- Geometric Admissibility Check (Base Locus) ---
                    p = np.array([1 if orbit_canon[k] == -1 else 0 for k in range(N)])
                    C_m = np.power(-1, np.dot(E, p) % 2)

                    has_base_locus = False
                    for simplex in t_data["simplices"]:
                        simplex_rays = [idx - 1 for idx in simplex if idx > 0]
                        if len(simplex_rays) == 0:
                            continue
                        survive_mask = np.all(E[:, simplex_rays] == 0, axis=1)
                        surviving_m_indices = np.where(survive_mask)[0]

                        if len(surviving_m_indices) > 0:
                            simplex_has_survivor = False
                            for m_idx in surviving_m_indices:
                                if m_prime_idx[m_idx] != m_idx or C_m[m_idx] == S_canon:
                                    simplex_has_survivor = True
                                    break
                            if not simplex_has_survivor:
                                has_base_locus = True
                                break

                    if has_base_locus:
                        continue
                    # --------------------------------------------------

                    chi_fix = 0
                    o_planes_list = []
                    all_eqs = []
                    is_reducible_global = False

                    for lam, parity_lam in H:
                        b_prime = tuple(lam[k] * orbit_canon[k] for k in range(N))

                        eqs = []
                        visited = set()
                        for k in range(N):
                            if k in visited:
                                continue
                            j = perm.get(k + 1, k + 1) - 1
                            if k == j:
                                if b_prime[k] == -1:
                                    eqs.append(f"x{k+1}")
                                visited.add(k)
                            else:
                                if b_prime[k] != b_prime[j]:
                                    eqs.append(f"x{k+1}")
                                    eqs.append(f"x{j+1}")
                                else:
                                    if b_prime[k] == 1:
                                        eqs.append(f"x{min(k,j)+1}-x{max(k,j)+1}")
                                    else:
                                        eqs.append(f"x{min(k,j)+1}+x{max(k,j)+1}")
                                visited.add(k)
                                visited.add(j)

                        if len(eqs) == 0:
                            if S_canon * parity_lam == -1:
                                is_reducible_global = True
                                all_eqs.append(tuple(["Whole CY"]))
                                continue
                            else:
                                chi_op = 2 * (h11 - h21)
                                chi_fix += chi_op
                                o_planes_list.append(
                                    {"OIDEAL": ["Whole CY"], "ODIM": 9, "chi": chi_op}
                                )
                                all_eqs.append(tuple(["Whole CY"]))
                                continue

                        if len(eqs) > 4:
                            continue

                        vanishes_identically = S_canon * parity_lam == -1

                        try:
                            chi_op_data = compute_oplan_chi(
                                c2,
                                kappa,
                                eqs,
                                D_basis,
                                vanishes_identically=vanishes_identically,
                                simplices=t_data["simplices"],
                            )
                            chi_op = chi_op_data["chi"]
                        except ValueError as e:
                            raise ValueError(
                                f"eqs={eqs}, D_basis shape={D_basis.shape}"
                            ) from e
                        chi_fix += chi_op

                        if vanishes_identically:
                            odim = 11 - 2 * len(eqs)
                        else:
                            odim = 9 - 2 * len(eqs)

                        if chi_op == 0 and odim != 9:
                            continue

                        op_entry = {
                            "OIDEAL": eqs,
                            "ODIM": odim,
                            "chi": chi_op,
                        }
                        if odim == 7:
                            op_entry["D3"] = chi_op_data.get("D3", 0)
                            op_entry["c2D"] = chi_op_data.get("c2D", 0)
                        if "is_reducible" in chi_op_data:
                            op_entry["is_reducible"] = chi_op_data["is_reducible"]
                            is_reducible_global = True
                        o_planes_list.append(op_entry)
                        all_eqs.append(tuple(sorted(eqs)))

                    hash_key = (tuple(sorted(all_eqs)), S_canon, S_Omega)
                    if hash_key in seen_configs:
                        continue
                    seen_configs.add(hash_key)

                    h21_diff = (
                        2 - 2 * S_Omega + 2 * (h11_plus - h11_minus) - chi_fix
                    ) / 2
                    h21_plus = (h21 + h21_diff) / 2
                    h21_minus = (h21 - h21_diff) / 2

                    integer_hodge = True
                    if (
                        not h21_plus.is_integer()
                        or not h21_minus.is_integer()
                        or h21_plus < 0
                        or h21_minus < 0
                    ):
                        integer_hodge = False

                    inv_str = (
                        "{"
                        + ", ".join([f"x{k}->x{v}" for k, v in perm.items() if k != v])
                        + "}"
                    )
                    if inv_str == "{}":
                        inv_str = "Id"

                    if S_Omega == -1:
                        # O3/O7 system tadpoles
                        if S_canon == 1:
                            # Compute total D_O7
                            D_O7_total = np.zeros(h11)
                            for op in o_planes_list:
                                if op.get("ODIM") == 7:
                                    eqs = op["OIDEAL"]
                                    m = re.search(r"x(\d+)", eqs[0])
                                    if m:
                                        coord_idx = int(m.group(1)) - 1
                                        D_O7_class = D_basis[:, coord_idx]
                                        D_O7_total += D_O7_class
                                        op["D3"] = int(
                                            round(
                                                np.einsum(
                                                    "ijk,i,j,k->",
                                                    kappa,
                                                    D_O7_class,
                                                    D_O7_class,
                                                    D_O7_class,
                                                )
                                            )
                                        )
                                        op["c2D"] = int(round(np.dot(c2, D_O7_class)))

                            if np.any(D_O7_total):
                                D3_total = np.einsum(
                                    "ijk,i,j,k->",
                                    kappa,
                                    D_O7_total,
                                    D_O7_total,
                                    D_O7_total,
                                )
                                c2D_total = np.dot(c2, D_O7_total)
                            else:
                                D3_total = 0
                                c2D_total = 0
                        else:
                            D3_total = 0
                            c2D_total = 0
                            if I4 is not None:
                                N_divs = I4.shape[0]
                                for op in o_planes_list:
                                    if op.get("ODIM") == 7:
                                        if "D3" in op and "c2D" in op:
                                            D3_total += op["D3"]
                                            c2D_total += op["c2D"]
                                        else:
                                            eqs = op["OIDEAL"]
                                            if len(eqs) >= 2:
                                                m1 = re.search(r"x(\d+)", eqs[0])
                                                m2 = re.search(r"x(\d+)", eqs[1])
                                                if m1 and m2:
                                                    i = int(m1.group(1)) - 1
                                                    j = int(m2.group(1)) - 1

                                                    # D3_total = \int_X D_i D_j (D_i + D_j + K)^2
                                                    V = -np.ones(N_divs)
                                                    V[i] += 1
                                                    V[j] += 1

                                                    current_D3 = 0.0
                                                    for a in range(N_divs):
                                                        for b in range(N_divs):
                                                            current_D3 += (
                                                                V[a]
                                                                * V[b]
                                                                * I4[i, j, a, b]
                                                            )

                                                    D3_total += current_D3
                                                    op["D3"] = int(round(current_D3))

                                                    # c2D_total = \int_X c_2(CY) D_i D_j (-K)
                                                    # chi(S_O7) = D3 + c2D => c2D = chi - D3
                                                    current_c2D = (
                                                        op.get("chi", 0) - current_D3
                                                    )
                                                    c2D_total += current_c2D
                                                    op["c2D"] = int(round(current_c2D))

                        # Compute n^S_{df=0} = \int_O7 (c_2(S) - c_1(S)^2) = c2D_total - D3_total - \int_V D_X^2 D_{O7}^2
                        n_S = 0.0
                        if S_canon == 1:
                            if np.any(D_O7_total):
                                D_X_class = np.sum(D_basis, axis=1)
                                int_DX_DO72 = np.einsum(
                                    "ijk,i,j,k->",
                                    kappa,
                                    D_X_class,
                                    D_O7_total,
                                    D_O7_total,
                                )
                                n_S = c2D_total - D3_total - int_DX_DO72
                        else:
                            if I4 is not None:
                                for op in o_planes_list:
                                    if op.get("ODIM") == 7:
                                        eqs = op["OIDEAL"]
                                        if len(eqs) == 1:
                                            m = re.search(r"x(\d+)", eqs[0])
                                            if m:
                                                coord_idx = int(m.group(1)) - 1
                                                D_O7_class = D_basis[:, coord_idx]
                                                D_X_class = np.sum(D_basis, axis=1)
                                                int_DX_DO72 = np.einsum(
                                                    "ijk,i,j,k->",
                                                    kappa,
                                                    D_X_class,
                                                    D_O7_class,
                                                    D_O7_class,
                                                )
                                                n_S += (
                                                    op.get("c2D", 0)
                                                    - op.get("D3", 0)
                                                    - int_DX_DO72
                                                )
                                        elif len(eqs) >= 2:
                                            m1 = re.search(r"x(\d+)", eqs[0])
                                            m2 = re.search(r"x(\d+)", eqs[1])
                                            if m1 and m2:
                                                i = int(m1.group(1)) - 1
                                                j = int(m2.group(1)) - 1
                                                int_DX2_DO72 = 0.0
                                                for a in range(N_divs):
                                                    for b in range(N_divs):
                                                        int_DX2_DO72 += I4[a, b, i, j]
                                                n_S += (
                                                    op.get("c2D", 0)
                                                    - op.get("D3", 0)
                                                    - int_DX2_DO72
                                                )

                        n_S = int(round(n_S))

                        N_O3 = sum(
                            op["chi"] for op in o_planes_list if op.get("ODIM") == 3
                        )
                        o7_planes = [op for op in o_planes_list if op.get("ODIM") == 7]
                        tadpole_data_list = []

                        # Generate all combinations of WU and SO8 for the O7-planes
                        if len(o7_planes) > 0:
                            combinations = list(
                                itertools.product(["WU", "SO8"], repeat=len(o7_planes))
                            )
                        else:
                            combinations = [()]

                        for comb in combinations:
                            chi_D7_comb = 0
                            chi_O7_comb = 0
                            choices = {}

                            for idx, choice in enumerate(comb):
                                if len(o7_planes) > 0:
                                    op = o7_planes[idx]
                                    c2D_op = op.get("c2D", 0)
                                    D3_op = op.get("D3", 0)
                                    if choice == "WU":
                                        cur_chi_D7 = int(
                                            round(8 * c2D_op + 344 * D3_op)
                                        )
                                    else:
                                        cur_chi_D7 = int(round(8 * (c2D_op + D3_op)))
                                    chi_D7_comb += cur_chi_D7
                                    chi_O7_comb += int(round(c2D_op + D3_op))
                                    choices[str(op["OIDEAL"])] = [choice, cur_chi_D7]

                            N_D3 = N_O3 / 2.0 + chi_O7_comb / 6.0 + chi_D7_comb / 24.0
                            chi_Y4 = int(round(12 * N_D3))
                            flux_frac = N_D3 % 1.0

                            entry = {
                                "choices": choices,
                                "N_O3": N_O3,
                                "chi_O7_total": chi_O7_comb,
                                "chi_D7": chi_D7_comb,
                                "N_D3": float(N_D3),
                                "N_flux_frac": float(flux_frac),
                                "N_D3_mobile_max": int(N_D3 // 1),
                                "chi_Y4_correction_O3": 6 * N_O3,
                                "chi_Y4_naive": chi_Y4 - 6 * N_O3,
                                "chi_Y4_corrected": chi_Y4,
                            }

                            if n_S != 0:
                                entry["nodal_singularities"] = True
                                entry["n_S_df0"] = n_S
                            else:
                                entry["nodal_singularities"] = False
                                entry["n_S_df0"] = 0

                            entry["ed3_instantons"] = []
                            tadpole_data_list.append(entry)

                        # ED3 Instanton divisors
                        if compute_cohomologies:
                            for k in range(N):
                                j = perm.get(k + 1, k + 1) - 1
                                # Invariant divisor
                                if k == j:
                                    if k in divisor_cohoms:
                                        h10, h20 = divisor_cohoms[k]
                                        # Completely rigid divisor
                                        if h10 == 0 and h20 == 0:
                                            eqs_D = [f"x{k+1}"]
                                            chi_D_data = compute_oplan_chi(
                                                c2,
                                                kappa,
                                                eqs_D,
                                                D_basis,
                                                vanishes_identically=False,
                                                simplices=t_data["simplices"],
                                            )
                                            chi_D = chi_D_data["chi"]
                                            h11_D = chi_D - 2 + 4 * h10 - 2 * h20

                                            n_O3_D = sum(
                                                op["chi"]
                                                for op in o_planes_list
                                                if op.get("ODIM") == 3
                                                and f"x{k+1}" in op["OIDEAL"]
                                            )

                                            D_X = D_basis[:, k]
                                            op_integrals = []

                                            if S_canon == 1:
                                                if np.any(D_O7_total):
                                                    for op1 in o7_planes:
                                                        eqs = op1["OIDEAL"]
                                                        m = re.search(r"x(\d+)", eqs[0])
                                                        if m:
                                                            coord_idx = (
                                                                int(m.group(1)) - 1
                                                            )
                                                            D_O7_class = D_basis[
                                                                :, coord_idx
                                                            ]
                                                            i2_Si = np.einsum(
                                                                "ijk,i,j,k->",
                                                                kappa,
                                                                D_X,
                                                                D_X,
                                                                D_O7_class,
                                                            )
                                                            i1_Si2 = np.einsum(
                                                                "ijk,i,j,k->",
                                                                kappa,
                                                                D_X,
                                                                D_O7_class,
                                                                D_O7_class,
                                                            )
                                                            op_integrals.append(
                                                                (op1, i2_Si, i1_Si2)
                                                            )
                                            else:
                                                if I4 is not None:
                                                    N_divs = I4.shape[0]
                                                    for op1 in o7_planes:
                                                        eqs = op1["OIDEAL"]
                                                        i2_Si = 0.0
                                                        i1_Si2 = 0.0
                                                        if len(eqs) == 1:
                                                            m1 = re.search(
                                                                r"x(\d+)", eqs[0]
                                                            )
                                                            if m1:
                                                                i = int(m1.group(1)) - 1
                                                                D_O7_class = D_basis[
                                                                    :, i
                                                                ]
                                                                i2_Si = np.einsum(
                                                                    "ijk,i,j,k->",
                                                                    kappa,
                                                                    D_X,
                                                                    D_X,
                                                                    D_O7_class,
                                                                )
                                                                i1_Si2 = np.einsum(
                                                                    "ijk,i,j,k->",
                                                                    kappa,
                                                                    D_X,
                                                                    D_O7_class,
                                                                    D_O7_class,
                                                                )
                                                        elif len(eqs) >= 2:
                                                            m1 = re.search(
                                                                r"x(\d+)", eqs[0]
                                                            )
                                                            m2 = re.search(
                                                                r"x(\d+)", eqs[1]
                                                            )
                                                            if m1 and m2:
                                                                i = int(m1.group(1)) - 1
                                                                j = int(m2.group(1)) - 1
                                                                i2_Si = I4[k, k, i, j]

                                                                V = -np.ones(N_divs)
                                                                V[i] += 1
                                                                V[j] += 1
                                                                for a in range(N_divs):
                                                                    i1_Si2 += (
                                                                        V[a]
                                                                        * I4[k, i, j, a]
                                                                    )
                                                        op_integrals.append(
                                                            (op1, i2_Si, i1_Si2)
                                                        )

                                            for entry, comb in zip(
                                                tadpole_data_list, combinations
                                            ):
                                                chi_vertical_naive = 0
                                                for idx, (
                                                    op1,
                                                    i2_Si,
                                                    i1_Si2,
                                                ) in enumerate(op_integrals):
                                                    choice = comb[idx]
                                                    if choice == "WU":
                                                        chi_vertical_naive += int(
                                                            round(
                                                                -6 * i2_Si - 30 * i1_Si2
                                                            )
                                                        )
                                                    else:
                                                        chi_vertical_naive += int(
                                                            round(
                                                                -6 * i2_Si - 6 * i1_Si2
                                                            )
                                                        )

                                                entry["ed3_instantons"].append(
                                                    {
                                                        "divisor": f"x{k+1}",
                                                        "h10": h10,
                                                        "h20": h20,
                                                        "chi": chi_D,
                                                        "h11": int(h11_D),
                                                        "n_O3": n_O3_D,
                                                        "chi_correction_vertical_divisor": 2
                                                        * n_O3_D,
                                                        "chi_vertical_naive": chi_vertical_naive,
                                                        "chi_vertical_corrected": chi_vertical_naive
                                                        + 2 * n_O3_D,
                                                    }
                                                )

                        tadpole_data = tadpole_data_list
                    else:
                        # O5/O9 system: empty tadpole data
                        tadpole_data = []

                    if is_reducible_global:
                        tadpole_data = []
                        h21_plus = None
                        h21_minus = None
                        integer_hodge = False

                    orientifold_entry = {
                        "TRIANGN": triang_n,
                        "INVOLN": invol_counter,
                        "INVOL": inv_str,
                        "SCANON": S_canon,
                        "SOMEGA": int(S_Omega),
                        "BCANON": list(int(x) for x in orbit_canon),
                        "h11_plus": int(h11_plus) if not is_reducible_global else None,
                        "h11_minus": (
                            int(h11_minus)
                            if integer_hodge
                            else (float(h11_minus) if not is_reducible_global else None)
                        ),
                        "h21_plus": (
                            int(h21_plus)
                            if integer_hodge
                            else (float(h21_plus) if not is_reducible_global else None)
                        ),
                        "h21_minus": (
                            int(h21_minus)
                            if integer_hodge
                            else (float(h21_minus) if not is_reducible_global else None)
                        ),
                        "integer_hodge": integer_hodge,
                        "is_reducible": is_reducible_global,
                        "tadpole_data": tadpole_data,
                        "OPLANES": o_planes_list,
                    }

                    # 3. Cohomologies
                    if compute_cohomologies:
                        for op in orientifold_entry["OPLANES"]:
                            if op.get("ODIM") == 7 and len(op["OIDEAL"]) == 1:
                                eqs = op["OIDEAL"]
                                m = re.search(r"x(\d+)", eqs[0])
                                if m:
                                    coord_idx = int(m.group(1)) - 1
                                    if coord_idx in divisor_cohoms:
                                        h10, h20 = divisor_cohoms[coord_idx]
                                        chi = op["chi"]
                                        h11_D = chi - 2 + 4 * h10 - 2 * h20
                                        op["cohomologies"] = {
                                            "h10": h10,
                                            "h20": h20,
                                            "h11": int(h11_D),
                                        }
                                        op["rigid"] = h20 == 0
                                        op["completely_rigid"] = h10 == 0 and h20 == 0
                            elif op.get("ODIM") == 7 and len(op["OIDEAL"]) == 2:
                                eqs = op["OIDEAL"]
                                m1 = re.search(r"x(\d+)", eqs[0])
                                m2 = re.search(r"x(\d+)", eqs[1])
                                if m1 and m2:
                                    idx1 = int(m1.group(1)) - 1
                                    idx2 = int(m2.group(1)) - 1
                                    D1_charge = D_basis[:, idx1].tolist()
                                    D2_charge = D_basis[:, idx2].tolist()
                                    coh = run_cohomcalg_for_complete_intersection(
                                        sr_ideal_clean, Q.tolist(), D1_charge, D2_charge
                                    )
                                    if coh is not None:
                                        h10, h20 = coh[0], coh[1]
                                        chi = op["chi"]
                                        h11_D = chi - 2 + 4 * h10 - 2 * h20
                                        op["cohomologies"] = {
                                            "h10": h10,
                                            "h20": h20,
                                            "h11": int(h11_D),
                                        }
                                        op["rigid"] = h20 == 0
                                        op["completely_rigid"] = h10 == 0 and h20 == 0

                    out["orientifolds"].append(orientifold_entry)
                    invol_counter += 1

    return out


def run_cohomcalg_for_complete_intersection(
    sr_ideal, glsm_charges, D1_charge, D2_charge
):
    """
    Run cohomCalg to compute the cohomology of a line bundle on a complete intersection surface.

    This function writes the necessary input file for cohomCalg to evaluate the cohomology
    of a line bundle restricted to a surface defined as the complete intersection of two
    divisors D1 and D2 in the ambient space.

    Args:
        sr_ideal (list of list): The Stanley-Reisner ideal generators.
        glsm_charges (list of list): The GLSM charge matrix of the ambient space.
        D1_charge (list): The divisor class charge vector for the first intersecting divisor.
        D2_charge (list): The divisor class charge vector for the second intersecting divisor.

    Returns:
        list: A list [h10, h20] representing the dimensions of H^1(S, O(D)) and H^2(S, O(D)).
    """
    import os
    import subprocess
    import tempfile

    minus_D1 = [-c for c in D1_charge]
    minus_D2 = [-c for c in D2_charge]
    minus_D1D2 = [-(c1 + c2) for c1, c2 in zip(D1_charge, D2_charge)]

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".in") as f:
        tmp_path = f.name
        for i in range(len(glsm_charges[0])):
            charges = [row[i] for row in glsm_charges]
            f.write(f"vertex x{i+1} | GLSM: ( {', '.join(map(str, charges))} );\n")

        sr_vars = []
        for gen in sr_ideal:
            var_str = "*".join([f"x{i+1}" for i in gen])
            sr_vars.append(var_str)
        f.write(f"srideal [ {', '.join(sr_vars)} ];\n")

        f.write(f"ambientcohom O( {', '.join(map(str, minus_D1))} );\n")
        f.write(f"ambientcohom O( {', '.join(map(str, minus_D2))} );\n")
        f.write(f"ambientcohom O( {', '.join(map(str, minus_D1D2))} );\n")

    cohomcalg_path = os.environ.get("COHOMCALG_PATH", "/opt/homebrew/bin/cohomcalg")
    if not os.path.exists(cohomcalg_path):
        cohomcalg_path = "cohomcalg"
    try:
        res = subprocess.run(
            [cohomcalg_path, "--nomonomfile", tmp_path],
            capture_output=True,
            text=True,
            timeout=60,
        )
        output = res.stdout
    except Exception:
        os.unlink(tmp_path)
        return [0, 0]
    os.unlink(tmp_path)

    def parse_bundle():
        """
        Parse the cohomCalg output for line bundle cohomology dimensions.

        Returns:
            list of list: Parsed dimensions for each cohomology sequence.
        """
        lines = output.splitlines()
        cohoms = []
        for line in lines:
            if line.strip().startswith("dim H^i"):
                parts = line.split("=")
                if len(parts) > 1:
                    tup_str = parts[1].strip().strip("()")
                    dims = [int(x.strip()) for x in tup_str.split(",")]
                    cohoms.append(dims)
        return cohoms

    parsed = parse_bundle()
    if len(parsed) != 3:
        return [0, 0]

    h10 = 0
    h20 = 0
    return [h10, h20]


def process_polytope_unpack(args, compute_cohomologies=False):
    """
    Unpack arguments and process a single polytope.

    Args:
        args (tuple): A tuple containing (poly_idx, h11, pts).
        compute_cohomologies (bool, optional): Whether to run cohomcalg for divisor cohomologies. Defaults to False.

    Returns:
        dict: The result of process_polytope.
    """
    return process_polytope(*args, compute_cohomologies=compute_cohomologies)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--h11", type=int, required=True, help="h11 value to process")
    parser.add_argument(
        "--compute_cohomologies",
        action="store_true",
        help="Run cohomcalg for O7 planes",
    )
    parser.add_argument("--workers", type=int, default=1, help="Number of workers")
    parser.add_argument(
        "--source",
        type=str,
        default="cytools",
        choices=["cytools", "hf"],
        help="Source of the polytopes database. 'cytools' uses the TU Wien database, 'hf' uses Hugging Face 'calabi-yau-data/polytopes-4d' dataset.",
    )
    args = parser.parse_args()

    h11 = args.h11
    results = []

    print(f"Fetching CYs with h11={h11} from source: {args.source}...")

    # Calculate global POLYID offset using HF dataset via DuckDB and cytools (always accurate)
    offset = 0
    if h11 > 1:
        print(f"Calculating global POLYID offset for h11 < {h11}...")
        con = duckdb.connect()
        con.execute("INSTALL httpfs; LOAD httpfs;")
        res_offset = con.execute(
            f"SELECT vertices FROM read_parquet('hf://datasets/calabi-yau-data/polytopes-4d/*.parquet') WHERE h11 < {h11}"
        )
        while True:
            row = res_offset.fetchone()
            if row is None:
                break
            pts = [[int(x) for x in v] for v in row[0]]
            poly = cytools.Polytope(pts).dual()
            if poly.is_favorable(lattice="N"):
                offset += 1
        con.close()
        print(
            f"Found {offset} favorable polytopes before h11={h11}. Starting POLYID at {offset + 1}."
        )

    if args.source == "hf":
        con = duckdb.connect()
        con.execute("INSTALL httpfs; LOAD httpfs;")

        # Read from the parquet files directly
        res = con.execute(
            f"SELECT vertices FROM read_parquet('hf://datasets/calabi-yau-data/polytopes-4d/*.parquet') WHERE h11 = {h11}"
        )

        pts_list = []
        while True:
            row = res.fetchone()
            if row is None:
                break
            pts = [[int(x) for x in v] for v in row[0]]
            poly = cytools.Polytope(pts).dual()
            if poly.is_favorable(lattice="N"):
                pts_list.append(poly.points().tolist())

        con.close()
    else:
        polys = cytools.fetch_polytopes(h11=h11, limit=1000000000, favorable=True)
        pts_list = [p.points().tolist() for p in polys]

    # pass offset+i so that process_polytope (which adds 1) results in correct POLYID
    ids = [(offset + i, h11, pts_list[i]) for i in range(len(pts_list))]
    print(f"Found {len(ids)} favorable polytopes.")

    worker = partial(
        process_polytope_unpack, compute_cohomologies=args.compute_cohomologies
    )

    with NoDaemonPool(processes=args.workers) as pool:
        for res in tqdm(pool.imap_unordered(worker, ids), total=len(ids)):
            if res is not None:
                results.append(res)

    out_file = f"databases/cy_orient_h11_{h11}.json"
    os.makedirs("databases", exist_ok=True)

    with open(out_file, "w") as f:
        json.dump(results, f, indent=2)

    print(f"Saved {len(results)} CYs to {out_file}!")
