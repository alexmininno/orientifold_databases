"""
Utility script to read and summarize an orientifold database JSON file.
"""

import os
import warnings

warnings.filterwarnings("ignore", category=UserWarning)

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import argparse
import json


def load_orientifold_db(filepath):
    """
    Load an orientifold database JSON file and return its contents as a list of dictionaries.

    Args:
        filepath (str): Path to the database JSON file.

    Returns:
        list: A list of Calabi-Yau configurations and their orientifolds.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Database file not found: {filepath}")

    with open(filepath, "r") as f:
        db = json.load(f)

    return db


def print_summary(db):
    """
    Print a basic summary of the loaded database.

    Calculates and displays the total number of Calabi-Yau threefolds,
    triangulations, and orientifold configurations in the database.
    It also prints an example sample from the first Calabi-Yau entry.

    Args:
        db (list): A list of dictionaries representing the database.
    """
    print(f"Loaded database with {len(db)} Calabi-Yau threefolds.")

    total_triangulations = sum(len(cy.get("triangulations", [])) for cy in db)
    total_orientifolds = sum(len(cy.get("orientifolds", [])) for cy in db)

    print(f"Total triangulations: {total_triangulations}")
    print(f"Total orientifolds configurations: {total_orientifolds}")

    # Print a sample from the first CY
    if db:
        first_cy = db[0]
        print("\n--- Example CY (First Entry) ---")
        print(f"POLYID: {first_cy.get('POLYID')}")
        print(f"h11: {first_cy.get('h11')}")
        print(f"h21: {first_cy.get('h21')}")
        print(f"Number of orientifolds: {len(first_cy.get('orientifolds', []))}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Read and summarize an orientifold database JSON file."
    )
    parser.add_argument(
        "db_file",
        type=str,
        help="Path to the JSON database file (e.g., Database/cy_orient_h11_1.json)",
    )

    args = parser.parse_args()

    try:
        db = load_orientifold_db(args.db_file)
        print_summary(db)
    except Exception as e:
        print(f"Error: {e}")
