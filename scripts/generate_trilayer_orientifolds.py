import argparse
import json
import multiprocessing as mp
import os
from functools import partial

import numpy as np
from generate_single_polytope_orientifold_mixed import process_polytope
from tqdm import tqdm


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


def main():
    """
    Main entry point for generating orientifolds for trilayer polytopes.

    Parses arguments, filters favorable polytopes, sets up a multiprocessing pool,
    and runs process_polytope on each favorable polytope. The results are saved
    to a JSON file.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--compute_cohomologies",
        action="store_true",
        help="Run cohomcalg for O7 planes",
    )
    parser.add_argument("--workers", type=int, default=1, help="Number of workers")
    parser.add_argument(
        "--h11_max",
        type=int,
        default=None,
        help="Maximum h11 (Picard_number_N) to process",
    )
    args = parser.parse_args()

    input_file = "databases/trilayer.json"
    print(f"Loading data from {input_file}...")
    with open(input_file, "r") as f:
        data = json.load(f)

    # Filter favorable polytopes
    if args.h11_max is not None:
        favorable_polys = [
            d
            for d in data
            if d.get("cy_is_favorable") == True
            and d.get("Picard_number_N", 0) <= args.h11_max
        ]
        print(
            f"Found {len(favorable_polys)} favorable polytopes with h11 <= {args.h11_max}."
        )
    else:
        favorable_polys = [d for d in data if d.get("cy_is_favorable") == True]
        print(f"Found {len(favorable_polys)} favorable polytopes.")

    # Prepare inputs for multiprocessing pool
    ids = []
    for poly in favorable_polys:
        poly_id = poly["id"]
        h11 = poly["Picard_number_N"]
        # Coordinates in trilayer.json are in columns (4xN). We need them as rows (Nx4).
        pts = np.array(poly["coordinates"]).T.tolist()
        ids.append((poly_id, h11, pts))

    worker = partial(
        process_polytope_unpack, compute_cohomologies=args.compute_cohomologies
    )

    results = []
    with NoDaemonPool(processes=args.workers) as pool:
        for res in tqdm(pool.imap_unordered(worker, ids), total=len(ids)):
            if res is not None and len(res) > 0:
                results.append(res)

    out_file = "databases/trilayer_orientifolds.json"
    os.makedirs("databases", exist_ok=True)

    # Sort results by POLYID to keep it organized
    results.sort(key=lambda x: x.get("POLYID", 0))

    with open(out_file, "w") as f:
        json.dump(results, f, indent=2)

    print(f"Saved {len(results)} CYs to {out_file}!")


if __name__ == "__main__":
    main()
