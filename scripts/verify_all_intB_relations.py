"""
Verification script to check intersection relations between orientifold bases (B3)
and their Calabi-Yau threefold (X3) double covers.
"""

import os
import warnings

warnings.filterwarnings("ignore", category=UserWarning)

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import json

import cytools
import numpy as np
from tqdm import tqdm


def verify_all_polytopes(limit=None):
    """
    Verify intersection relations for all polytopes between a Calabi-Yau threefold (X3)
    and its orientifold base (B3).

    This function computes the properties of the base B3 and the Calabi-Yau X3.
    It checks three consistency relations:
    1. Triple intersection: \int_{X_3} D'_{m(a)} D'_{m(b)} D'_{m(c)} = 2 \int_{B_3} D_a D_b D_c
    2. c2 relation: \int_{X_3} c_2(X_3) D'_{m(a)} = 2 \int_{B_3} D_a (c_1^2 + c_2)
    3. chi relation: \chi(X_3) = 2 \int_{B_3} (c_3 - c_1 c_2 - 2 c_1^3)

    Args:
        limit (int, optional): The maximum number of polytopes to verify. Defaults to None.
    """
    with open("databases/3dpoly.json", "r") as f:
        b3_data = json.load(f)

    with open("databases/trilayer.json", "r") as f:
        x3_data = json.load(f)

    x3_dict = {entry["id"]: entry for entry in x3_data}
    if limit:
        b3_data = b3_data[:limit]

    print(f"Verifying {len(b3_data)} polytopes using prime divisor basis...")

    # Statistics
    count_success = 0
    count_fail = 0

    for b3_entry in tqdm(b3_data):
        poly_id = b3_entry["id"]

        pts_B3 = list(zip(*b3_entry["coordinates"]))
        p3 = cytools.Polytope(pts_B3)
        t3 = p3.triangulate()
        b3_tv = t3.get_toric_variety()

        if poly_id not in x3_dict:
            continue

        x3_entry = x3_dict[poly_id]
        if not x3_entry.get("cy_is_favorable", True):
            continue

        pts_X3 = list(zip(*x3_entry["coordinates"]))
        p4 = cytools.Polytope(pts_X3)
        t4 = p4.triangulate()
        cy3 = t4.get_cy()

        # Build mapping from B3 prime divisors to X3 prime divisors
        # p3.points() has origin at index 0, and prime divisors at indices 1 to N_b
        p3_pts = p3.points()
        p4_pts = p4.points()

        mapping = {}
        for i in range(1, len(p3_pts)):
            target = tuple(list(p3_pts[i]) + [1])
            for j in range(1, len(p4_pts)):
                if tuple(p4_pts[j]) == target:
                    mapping[i] = j
                    break

        # Compute B3 properties
        # kappa_raw has shape (N_b + 1, N_b + 1, N_b + 1) where index 0 is Canonical divisor
        kappa_raw = b3_tv.intersection_numbers(format="dense")
        kappa = np.array(kappa_raw)
        if kappa.ndim == 0:
            kappa = np.array([[[kappa]]])

        N_b = len(p3_pts) - 1

        # Compute X3 properties
        lhs_chi = cy3.chi()

        try:
            lhs_k3_raw = cy3.intersection_numbers(format="dense")
            if isinstance(lhs_k3_raw, dict):
                # If sparse format is returned, convert to dense
                # N_v is the number of points in p4 minus 1 (origin)
                N_v = len(p4_pts) - 1
                # The shape is (N_v + 1, N_v + 1, N_v + 1)
                lhs_k3 = np.zeros((N_v + 1, N_v + 1, N_v + 1))
                for (a, b, c), val in lhs_k3_raw.items():
                    lhs_k3[a, b, c] = val
                    lhs_k3[a, c, b] = val
                    lhs_k3[b, a, c] = val
                    lhs_k3[b, c, a] = val
                    lhs_k3[c, a, b] = val
                    lhs_k3[c, b, a] = val
            else:
                lhs_k3 = np.array(lhs_k3_raw)
                if lhs_k3.ndim == 0:
                    lhs_k3 = np.array([[[lhs_k3]]])

            # second_chern_class also has size (N_v + 1), index 0 is canonical divisor
            lhs_c2_raw = cy3.second_chern_class()
            if isinstance(lhs_c2_raw, dict):
                N_v = len(p4_pts) - 1
                lhs_c2 = np.zeros(N_v + 1)
                for a, val in lhs_c2_raw.items():
                    lhs_c2[a] = val
            else:
                lhs_c2 = np.array(lhs_c2_raw)

        except Exception as e:
            print(f"Error computing X3 properties for ID {poly_id}: {e}")
            count_fail += 1
            continue

        failed = False

        # 3. Check triple intersection: \int_{X_3} D'_{m(a)} D'_{m(b)} D'_{m(c)} = 2 \int_{B_3} D_a D_b D_c
        for a in range(1, N_b + 1):
            for b in range(1, N_b + 1):
                for c in range(1, N_b + 1):
                    lhs_val = lhs_k3[mapping[a], mapping[b], mapping[c]]
                    rhs_val = 2 * kappa[a, b, c]
                    if lhs_val != rhs_val:
                        print(
                            f"Failed triple intersection for ID {poly_id} at ({a},{b},{c}): LHS={lhs_val}, RHS={rhs_val}"
                        )
                        failed = True
                        break
                if failed:
                    break
            if failed:
                break

        if failed:
            count_fail += 1
            continue

        # 4. Check c2 relation: \int_{X_3} c_2(X_3) D'_{m(a)} = 2 \int_{B_3} D_a (c_1^2 + c_2)
        # Note: c_1(B_3) = \sum_{i=1}^{N_b} D_i.
        for a in range(1, N_b + 1):
            lhs_val = lhs_c2[mapping[a]]

            int_Da_c1_sq = np.sum(kappa[a, 1:, 1:])

            int_Da_c2 = 0
            for i in range(1, N_b + 1):
                for j in range(i + 1, N_b + 1):
                    int_Da_c2 += kappa[a, i, j]

            rhs_val = 2 * (int_Da_c1_sq + int_Da_c2)

            if lhs_val != rhs_val:
                print(
                    f"Failed c2 relation for ID {poly_id} at divisor {a}: LHS={lhs_val}, RHS={rhs_val}"
                )
                failed = True
                break

        if failed:
            count_fail += 1
            continue

        # 5. Check chi relation: \chi(X_3) = 2 \int_{B_3} (c_3 - c_1 c_2 - 2 c_1^3)
        # c_3 is the Euler characteristic of the base, which is the number of 3-simplices in the triangulation
        int_c3 = len(t3.simplices())

        int_c1_c2 = 0
        for a in range(1, N_b + 1):
            for i in range(1, N_b + 1):
                for j in range(i + 1, N_b + 1):
                    int_c1_c2 += kappa[a, i, j]

        int_c1_3 = np.sum(kappa[1:, 1:, 1:])

        rhs_chi = 2 * (int_c3 - int_c1_c2 - 2 * int_c1_3)
        if lhs_chi != rhs_chi:
            print(f"Failed chi relation for ID {poly_id}: LHS={lhs_chi}, RHS={rhs_chi}")
            failed = True

        if not failed:
            count_success += 1
        else:
            count_fail += 1

    print(
        f"\nVerification finished: {count_success} successful, {count_fail} failed out of {count_success+count_fail} favourable polytopes processed."
    )


if __name__ == "__main__":
    verify_all_polytopes(limit=None)
