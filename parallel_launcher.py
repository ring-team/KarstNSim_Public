#!/usr/bin/env python3
"""
Parallel launcher for KarstNSim simulation sequences.

Each worker is assigned one seed. Within a worker, several simulations can be
executed sequentially, for example:
    seed 1 -> simulation 1 -> simulation 2 -> simulation 3
    seed 2 -> simulation 1 -> simulation 2 -> simulation 3

The script can be used either through command-line arguments or by editing the
USER CONFIGURATION block in main().
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import dataclasses
import datetime as dt
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, Sequence


@dataclasses.dataclass(frozen=True)
class SimulationStep:
    """
    Description of one simulation step in a sequential workflow.

    instruction_file:
        Path to the reference instruction file to be patched for each seed.
    name_suffix:
        Optional suffix added to the generated temporary instruction filename.
    """

    instruction_file: Path
    name_suffix: str = ""


@dataclasses.dataclass(frozen=True)
class RunConfig:
    """
    Global run configuration.

    executable:
        Path to karstnsim.exe.
    repository_root:
        Root directory containing Input_files/ and outputs/.
    steps:
        Ordered list of simulation steps. Steps are run sequentially for each seed.
    seeds:
        Seeds distributed across parallel workers.
    max_workers:
        Maximum number of simultaneous KarstNSim processes.
    run_root:
        Directory used for generated instruction files and logs.
    stop_on_error:
        If true, a failed step stops the sequence for the corresponding seed.
    keep_generated_instructions:
        If false, generated temporary instruction files are removed after success.
    """

    executable: Path
    repository_root: Path
    steps: Sequence[SimulationStep]
    seeds: Sequence[int]
    max_workers: int
    run_root: Path
    stop_on_error: bool = True
    keep_generated_instructions: bool = True


def _replace_or_append_parameter(text: str, key: str, value: str) -> str:
    """
    Replace an active 'key: value' line, or append the parameter if absent.

    Commented lines are intentionally ignored. This avoids modifying commented
    documentation examples inside instruction files.
    """

    pattern = re.compile(rf"^(?!\s*(//|#|%|::))\s*{re.escape(key)}\s*:.*$", re.MULTILINE)
    replacement = f"{key}: {value}"

    if pattern.search(text):
        return pattern.sub(replacement, text, count=1)

    if not text.endswith("\n"):
        text += "\n"

    return text + replacement + "\n"


def _patch_instruction_file(
        source_instruction: Path,
        destination_instruction: Path,
        *,
        seed: int,
        repository_root: Path,
        simulation_input_dir: Path,
        network_name_suffix: str | None = None,
) -> None:
    """
    Create a patched instruction file for one seed and one simulation step.

    The original instruction file is not modified. The generated file sets:
      - selected_seed to the worker seed;
      - number_of_iterations to 1;
      - vary_seed to false;
      - main_repository to the absolute repository root;
      - simulation_input_dir to the original instruction directory.

    simulation_input_dir is explicitly set because KarstNSim otherwise uses the
    directory of the temporary instruction file. Keeping it pointed to the
    original instruction folder preserves auxiliary files such as connectivity
    matrices.
    """

    text = source_instruction.read_text(encoding="utf-8")

    if network_name_suffix:
        if match := re.search(r"^(?!\s*(//|#|%|::))\s*karstic_network_name\s*:\s*(\S+).*$", text, flags=re.MULTILINE, ):
            original_name = match[2]
            text = _replace_or_append_parameter(
                text,
                "karstic_network_name",
                f"{original_name}{network_name_suffix}",
            )

    text = _replace_or_append_parameter(text, "main_repository", _as_posix(repository_root))
    text = _replace_or_append_parameter(text, "simulation_input_dir", _as_posix(simulation_input_dir))
    text = _replace_or_append_parameter(text, "selected_seed", str(seed))
    text = _replace_or_append_parameter(text, "number_of_iterations", "1")
    text = _replace_or_append_parameter(text, "vary_seed", "false")

    destination_instruction.parent.mkdir(parents=True, exist_ok=True)
    destination_instruction.write_text(text, encoding="utf-8")


def _as_posix(path: Path) -> str:
    """
    Return a normalized path string accepted by KarstNSim on Windows and Linux.
    """

    return str(path.resolve()).replace("\\", "/")


def _run_one_process(
        executable: Path,
        instruction_file: Path,
        *,
        cwd: Path,
        log_file: Path,
) -> int:
    """
    Launch one KarstNSim process and write stdout/stderr to a log file.

    A newline is provided on stdin because the current C++ main() waits for
    Enter before exiting, even in command-line mode.
    """

    log_file.parent.mkdir(parents=True, exist_ok=True)

    command = [str(executable.resolve()), _as_posix(instruction_file)]

    with log_file.open("w", encoding="utf-8", errors="replace") as log:
        completed = write_log(log, command, cwd)
    return completed.returncode


def write_log(log, command, cwd):
    log.write(f"Command: {' '.join(command)}\n")
    log.write(f"Working directory: {cwd.resolve()}\n")
    log.write(f"Started: {dt.datetime.now().isoformat(timespec='seconds')}\n\n")
    log.flush()

    result = subprocess.run(command, cwd=str(cwd.resolve()), input="\n", text=True, stdout=log,
                            stderr=subprocess.STDOUT, check=False, )

    log.write(f"\nFinished: {dt.datetime.now().isoformat(timespec='seconds')}\n")
    log.write(f"Return code: {result.returncode}\n")

    return result


def _run_sequence_for_seed(config: RunConfig, seed: int) -> tuple[int, bool, list[tuple[Path, int]]]:
    """
    Run all simulation steps sequentially for one seed.

    Returns:
        seed, success flag, and the list of generated instruction files with
        their corresponding return codes.
    """

    seed_dir = config.run_root / f"seed_{seed}"
    results: list[tuple[Path, int]] = []

    for step_index, step in enumerate(config.steps, start=1):
        source_instruction = step.instruction_file.resolve()

        if not source_instruction.exists():
            raise FileNotFoundError(f"Instruction file not found: {source_instruction}")

        suffix = f"_{step.name_suffix}" if step.name_suffix else ""
        generated_instruction = seed_dir / f"step_{step_index:02d}{suffix}_instructions_seed_{seed}.txt"
        log_file = seed_dir / f"step_{step_index:02d}{suffix}_seed_{seed}.log"

        _patch_instruction_file(
            source_instruction=source_instruction,
            destination_instruction=generated_instruction,
            seed=seed,
            repository_root=config.repository_root,
            simulation_input_dir=source_instruction.parent,
            network_name_suffix=None,
        )

        return_code = _run_one_process(
            executable=config.executable,
            instruction_file=generated_instruction,
            cwd=config.executable.parent,
            log_file=log_file,
        )

        results.append((generated_instruction, return_code))

        if return_code != 0 and config.stop_on_error:
            return seed, False, results

        if return_code == 0 and not config.keep_generated_instructions:
            generated_instruction.unlink(missing_ok=True)

    success = all(return_code == 0 for _, return_code in results)
    return seed, success, results


def run_parallel_sequences(config: RunConfig) -> int:
    """
    Execute seed-parallel, step-sequential KarstNSim workflows.

    Parallelism is applied only between seeds. For a given seed, steps are
    executed in the exact order provided in config.steps.
    """

    config.run_root.mkdir(parents=True, exist_ok=True)

    if not config.executable.exists():
        raise FileNotFoundError(f"Executable not found: {config.executable}")

    if not config.steps:
        raise ValueError("At least one simulation step must be provided.")

    if not config.seeds:
        raise ValueError("At least one seed must be provided.")

    max_workers = min(config.max_workers, len(config.seeds))
    max_workers = max(1, max_workers)

    print(f"Executable: {config.executable.resolve()}")
    print(f"Repository root: {config.repository_root.resolve()}")
    print(f"Run directory: {config.run_root.resolve()}")
    print(f"Seeds: {list(config.seeds)}")
    print(f"Parallel workers: {max_workers}")
    print(f"Sequential steps per seed: {len(config.steps)}")

    failures: list[int] = []

    with cf.ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(_run_sequence_for_seed, config, seed): seed
            for seed in config.seeds
        }

        for future in cf.as_completed(futures):
            seed = futures[future]
            try:
                completed_seed, success, results = future.result()
            except Exception as exc:
                failures.append(seed)
                print(f"[seed {seed}] failed before completion: {exc}", file=sys.stderr)
                continue

            status = "OK" if success else "FAILED"
            print(f"[seed {completed_seed}] {status}")

            for instruction_file, return_code in results:
                print(f"  return_code={return_code} | {instruction_file}")

            if not success:
                failures.append(completed_seed)

    if failures:
        print(f"Failed seeds: {failures}", file=sys.stderr)
        return 1

    return 0


def _parse_seed_list(raw: str) -> list[int]:
    """
    Parse seeds from either a comma-separated list or an inclusive range.

    Examples:
        1,2,3,10
        1-20
    """

    raw = raw.strip()

    if "-" in raw and "," not in raw:
        start_text, end_text = raw.split("-", maxsplit=1)
        start = int(start_text)
        end = int(end_text)
        if end < start:
            raise ValueError("Seed range end must be greater than or equal to start.")
        return list(range(start, end + 1))

    return [int(item.strip()) for item in raw.split(",") if item.strip()]


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run KarstNSim simulations in parallel over seeds and sequentially over instruction files."
    )

    parser.add_argument(
        "--exe",
        type=Path,
        required=False,
        help="Path to karstnsim.exe.",
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        required=False,
        help="Path to the KarstNSim repository root.",
    )
    parser.add_argument(
        "--instructions",
        type=Path,
        nargs="+",
        required=False,
        help="Ordered list of instruction files. They are run sequentially for each seed.",
    )
    parser.add_argument(
        "--seeds",
        type=str,
        required=False,
        help="Seed list, e.g. '1,2,3,4', or inclusive range, e.g. '1-20'.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, (os.cpu_count() or 2) - 1),
        help="Maximum number of parallel KarstNSim processes.",
    )
    parser.add_argument(
        "--run-root",
        type=Path,
        default=Path("batch_runs"),
        help="Directory for generated instruction files and logs.",
    )
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Continue the sequence for a seed even if one step fails.",
    )
    parser.add_argument(
        "--clean-generated-instructions",
        action="store_true",
        help="Remove generated temporary instruction files after successful runs.",
    )

    return parser


def _config_from_args(args: argparse.Namespace) -> RunConfig:
    if missing := [name for name, value in
                   {"--exe": args.exe, "--repo-root": args.repo_root, "--instructions": args.instructions,
                    "--seeds": args.seeds, }.items() if value is None]:
        raise ValueError(
            "Missing command-line arguments: "
            + ", ".join(missing)
            + ". Either provide all required arguments or edit the USER CONFIGURATION block in main()."
        )

    steps = [
        SimulationStep(instruction_file=instruction, name_suffix=instruction.parent.name)
        for instruction in args.instructions
    ]

    return RunConfig(
        executable=args.exe,
        repository_root=args.repo_root,
        steps=steps,
        seeds=_parse_seed_list(args.seeds),
        max_workers=args.workers,
        run_root=args.run_root,
        stop_on_error=not args.continue_on_error,
        keep_generated_instructions=not args.clean_generated_instructions,
    )


def main() -> int:
    """

    Two usage modes are supported:
      1. Command-line mode with argparse.
      2. Direct configuration mode by editing the USER CONFIGURATION block below.

    If at least one command-line argument is provided, argparse mode is used.
    If no command-line argument is provided, the direct configuration block is used.

    Usage example for argparse :

    python scripts/run_parallel_sequences.py ^
      --exe KarstNSim/build/release/karstnsim.exe ^
      --repo-root . ^
      --instructions ^
        Input_files/1_base/instructions.txt ^
        Input_files/2_polygenic_karst/instructions.txt ^
        Input_files/3_amplification/instructions.txt ^
        Input_files/4_karst_section_generation/instructions.txt ^
      --seeds 1-12 ^
      --workers 6

    """

    if len(sys.argv) > 1:
        parser = _build_parser()
        args = parser.parse_args()
        config = _config_from_args(args)
        return run_parallel_sequences(config)

    # ---------------------------------------------------------------------
    # USER CONFIGURATION BLOCK
    # Edit this block for direct execution from an IDE or by double-clicking.
    # ---------------------------------------------------------------------

    repository_root = Path(__file__).resolve().parents[0]

    config = RunConfig(
        executable=repository_root / "KarstNSim" / "build" / "release" / "karstnsim.exe",
        repository_root=repository_root,
        steps=[
            SimulationStep(repository_root / "Input_files" / "1_base" / "instructions.txt", "base"),
            SimulationStep(repository_root / "Input_files" / "2_polygenic_karst" / "instructions.txt", "polygenic"),
            SimulationStep(repository_root / "Input_files" / "3_amplification" / "instructions.txt", "amplification"),
            SimulationStep(repository_root / "Input_files" / "4_karst_section_generation" / "instructions.txt", "sections"),
        ],
        seeds=list(range(4)),
        max_workers=4,
        run_root=repository_root / "batch_runs",
        stop_on_error=True,
        keep_generated_instructions=True,
    )

    # ---------------------------------------------------------------------
    # END OF USER CONFIGURATION BLOCK
    # ---------------------------------------------------------------------

    return run_parallel_sequences(config)


if __name__ == "__main__":
    raise SystemExit(main())
