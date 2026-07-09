from __future__ import annotations

import argparse
import os
import sys
import glob
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
import pymatching
import stim


def build_uniform_rotated_surface_code(distance: int, error_rate: float, round: int) -> stim.Circuit:
    circuit = stim.Circuit.generated(
        "surface_code:rotated_memory_x",
        distance=distance,
        rounds=round,
        after_clifford_depolarization=error_rate,
        before_round_data_depolarization=error_rate,
        before_measure_flip_probability=error_rate,
        after_reset_flip_probability=error_rate,
    )
    return circuit


def use_observable_as_detector(
    dem: stim.DetectorErrorModel,
    *,
    observable_id: int = 0,
) -> stim.DetectorErrorModel:

    flattened = dem.flattened()
    extra_detector = flattened.num_detectors
    result = stim.DetectorErrorModel()

    for instruction in flattened:
        if instruction.type != "error":
            result.append(instruction)
            continue

        targets = []
        for target in instruction.targets_copy():
            if target.is_logical_observable_id() and target.val == observable_id:
                targets.append(stim.DemTarget.relative_detector_id(extra_detector))
            else:
                targets.append(target)
        result.append(
            "error",
            instruction.args_copy(),
            targets,
            tag=instruction.tag,
        )
    return result


def constrained_matching_weights(
    matching: pymatching.Matching,
    det: np.ndarray,
) -> tuple[float, float]:
    weights = []
    for forced_observable in [0, 1]:
        constrained_syndrome = np.concatenate(
            [det.astype(np.uint8, copy=False), np.array([forced_observable], dtype=np.uint8)]
        )
        _, weight = matching.decode(
            constrained_syndrome,
            return_weight=True,
            enable_correlations=True,
        )
        weights.append(float(weight))
    return weights[0], weights[1]


def weight_gap_to_db(gap: float) -> float:
    return gap * 10 / np.log(10)


def bin_value(value: float, bin_width: float) -> int:
    return int(round(value / bin_width))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute a complementary-gap distribution for a 10d-round rotated surface code.",
    )
    parser.add_argument("distance", type=int, help="Surface-code distance d.")
    parser.add_argument("shots", type=int, help="Number of Monte Carlo shots.")
    parser.add_argument("--seed", type=int, default=42, help="Stim sampler seed.")
    return parser.parse_args()

def run(distance: int, error_rate: float, shots: int, seed: int, bin_width: float, rounds: int) -> dict[int, int]:
    # create circuit, matching, and sampler
    circuit = build_uniform_rotated_surface_code(distance, error_rate, rounds)
    dem = circuit.detector_error_model(decompose_errors=True)
    constrained_dem = use_observable_as_detector(dem)
    matching = pymatching.Matching.from_detector_error_model(
        constrained_dem,
        enable_correlations=True,
    )
    sampler = circuit.compile_detector_sampler(seed=seed)

    dets, observables = sampler.sample(shots, separate_observables=True)
    counter: dict[int, int] = {}
    failures = 0
    signed_gaps: list[float] = []

    for det, observable in zip(dets, observables):
        weight0, weight1 = constrained_matching_weights(matching, det)
        prediction = int(weight1 < weight0)
        gap_db = weight_gap_to_db(abs(weight1 - weight0))
        correct = prediction == int(observable[0])
        signed_gap = gap_db if correct else -gap_db
        failures += int(not correct)
        signed_gaps.append(signed_gap)
        binned_gap = int(bin_value(signed_gap, bin_width))
        counter[binned_gap] = counter.get(binned_gap, 0) + 1
    return counter

def eval_logical_error_rate(distance: int, error_rate: float, shots: int, seed: int, rounds: int) -> float:
    circuit = build_uniform_rotated_surface_code(distance, error_rate, rounds)
    dem = circuit.detector_error_model(decompose_errors=True)
    matching = pymatching.Matching.from_detector_error_model(
        dem,
        enable_correlations=True,
    )
    sampler = circuit.compile_detector_sampler(seed=seed)

    dets, observables = sampler.sample(shots, separate_observables=True)
    preds = matching.decode_batch(dets, enable_correlations=True)
    num_errors = 0
    for shot in range(shots):
        actual_for_shot = observables[shot]
        predicted_for_shot = preds[shot]
        if not np.array_equal(actual_for_shot, predicted_for_shot):
            num_errors += 1
    logical_error_rate = num_errors / shots
    return logical_error_rate
    
def main() -> None:
    if len(sys.argv) < 5:
        print("Usage: python gap_eval.py run <num_patch> <num_round_coef> <distance> <shots>")
        sys.exit(1)
    num_patch = int(sys.argv[2])
    num_round_coef = int(sys.argv[3])
    distance = int(sys.argv[4])
    shots = int(sys.argv[5])

    seed = np.random.randint(0, 2**32 - 1)
    error_rate = 1e-3
    bin_width = 1
    if distance < 2:
        raise ValueError("distance must be at least 2")
    if shots < 1:
        raise ValueError("shots must be positive")

    rounds = num_round_coef * distance
    counter = run(distance, error_rate, shots, seed, bin_width, rounds)
    filename = f"./data/gap_stat_{num_round_coef}_{distance}.txt"

    if os.path.exists(filename):
        DBs, counts = np.loadtxt(filename, unpack=True, dtype=int)
        for DB, count in zip(DBs, counts):
            counter[DB] = counter.get(DB, 0) + int(count)

    with open(filename, "w") as f:
        for db, count in sorted(counter.items()):
            f.write(f"{int(db)} {int(count)}\n")

    logical_error_rate_per_patch = sum(count for db, count in counter.items() if db < 0) / sum(list(counter.values()))
    logical_error_rate = 1 - pow(1-logical_error_rate_per_patch, num_patch)
    logical_error_rate_per_round = 1-(1-logical_error_rate)**(1./(rounds))
    logical_error_rate_per_round *= 2 # account pL = px+pz
    print(f"0D Yoke: n={num_patch} r={num_round_coef}d d={distance}: logical error rate per round = {logical_error_rate_per_round:.6e}")

    logical_error_rate_per_patch_debug = eval_logical_error_rate(distance, error_rate, shots, seed, rounds)
    print(logical_error_rate_per_patch, logical_error_rate_per_patch_debug)



def plot() -> None:
    for num_round_coef in [4,6,8,10,12]:
        plt.figure(figsize=(10, 10))
        data_exist = False
        for distance in range(2, 12):
            filename = f"./data/gap_stat_{num_round_coef}_{distance}.txt"
            if not os.path.exists(filename):
                continue
            data_exist = True
            DBs, counts = np.loadtxt(filename, unpack=True, dtype=int)
            num_shot = sum(counts)
            plt.plot(DBs, counts/num_shot, marker="o", label=f"r={num_round_coef}d d={distance}")
        if data_exist is False:
            continue
        plt.xlabel("Complementary gap [dB]")
        plt.ylabel("Frequency")
        plt.xticks(np.arange(-100, 201, 10), rotation=90)
        plt.xlim(-100,200)
        plt.ylim(1e-9,1.1e0)
        plt.yscale("log")
        plt.title(f"Complementary gap distribution for r={num_round_coef}d")
        plt.grid()
        plt.legend()
        plt.savefig(f"./figure/gap_stat_r={num_round_coef}d.png")

if __name__ == "__main__":
    if sys.argv[1] == "run":
        main()
    else:
        plot()
