#!/usr/bin/env python3
"""Create matched fixed- and random-effects Manhattan plots from GWAMA output."""

from __future__ import annotations

import argparse
import math
from collections import defaultdict
from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as error:
    raise SystemExit(
        "ERROR: matplotlib is required for plotting. "
        "Install it with: python3 -m pip install matplotlib"
    ) from error


MARKER_COLUMNS = ("rs_number", "MARKERNAME", "markername", "SNP", "ID")
LOG_P_COLUMNS = ("_-log10_p-value", "-log10_p-value", "LOG10P", "log10p")
P_COLUMNS = ("p-value", "P", "p", "p_value", "pvalue")
SPECIAL_CHROMOSOMES = {"X": 23, "Y": 24, "XY": 25, "MT": 26}
COLORS = ("#2F6B9A", "#D88032")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create separate fixed- and random-effects Manhattan plots using "
            "the same chromosome layout and Y-axis scale. MARKERNAME must begin "
            "with CHR:POS, for example 1:12345:A:C."
        )
    )
    parser.add_argument("--fixed", required=True, type=Path, help="FFX GWAMA .out file")
    parser.add_argument("--random", required=True, type=Path, help="RFX GWAMA .out file")
    parser.add_argument(
        "--output-prefix",
        required=True,
        type=Path,
        help="Output prefix; .fixed.manhattan.png and .random.manhattan.png are added",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=5e-8,
        help="Genome-wide significance threshold (default: 5e-8)",
    )
    parser.add_argument("--dpi", type=int, default=180, help="PNG resolution (default: 180)")
    return parser.parse_args()


def first_present(fieldnames: list[str], candidates: tuple[str, ...]) -> str | None:
    names = set(fieldnames)
    return next((candidate for candidate in candidates if candidate in names), None)


def normalize_chromosome(raw: str) -> str:
    chromosome = raw.strip().upper()
    if chromosome.startswith("CHR"):
        chromosome = chromosome[3:]
    if chromosome == "M":
        chromosome = "MT"
    if chromosome.isdigit():
        chromosome = str(int(chromosome))
    if not chromosome:
        raise ValueError("empty chromosome")
    return chromosome


def chromosome_sort_key(chromosome: str) -> tuple[int, int | str]:
    if chromosome.isdigit():
        return (0, int(chromosome))
    if chromosome in SPECIAL_CHROMOSOMES:
        return (0, SPECIAL_CHROMOSOMES[chromosome])
    return (1, chromosome)


def parse_marker(marker: str) -> tuple[str, int]:
    fields = marker.split(":")
    if len(fields) < 2:
        raise ValueError("marker does not begin with CHR:POS")
    chromosome = normalize_chromosome(fields[0])
    position = int(fields[1].replace(",", ""))
    if position <= 0:
        raise ValueError("position must be positive")
    return chromosome, position


def read_gwama(path: Path) -> tuple[list[tuple[str, int, float]], int]:
    variants: list[tuple[str, int, float]] = []
    skipped = 0

    with path.open("r", encoding="utf-8", newline="") as handle:
        header_line = handle.readline()
        if not header_line:
            raise ValueError(f"empty GWAMA output: {path}")

        fieldnames = header_line.split()
        marker_column = first_present(fieldnames, MARKER_COLUMNS)
        log_p_column = first_present(fieldnames, LOG_P_COLUMNS)
        p_column = first_present(fieldnames, P_COLUMNS)

        if marker_column is None:
            raise ValueError(f"{path} has no recognized marker column")
        if log_p_column is None and p_column is None:
            raise ValueError(f"{path} has no recognized P-value column")

        marker_index = fieldnames.index(marker_column)
        log_p_index = fieldnames.index(log_p_column) if log_p_column else None
        p_index = fieldnames.index(p_column) if p_column else None

        for line_number, line in enumerate(handle, start=2):
            if not line.strip():
                continue
            values = line.split()
            if len(values) != len(fieldnames):
                raise ValueError(
                    f"{path}:{line_number}: expected {len(fieldnames)} fields, "
                    f"found {len(values)}"
                )

            try:
                chromosome, position = parse_marker(values[marker_index])
                if log_p_index is not None:
                    log_p = float(values[log_p_index])
                else:
                    p_value = float(values[p_index])  # type: ignore[index]
                    if not 0 < p_value <= 1:
                        raise ValueError("P-value outside (0, 1]")
                    log_p = -math.log10(max(p_value, 1e-300))
                if not math.isfinite(log_p) or log_p < 0:
                    raise ValueError("invalid -log10(P)")
            except (ValueError, OverflowError):
                skipped += 1
                continue

            variants.append((chromosome, position, log_p))

    if not variants:
        raise ValueError(
            f"{path} has no plottable variants. Marker IDs must begin with "
            "CHR:POS, for example 1:12345:A:C."
        )
    return variants, skipped


def build_layout(
    fixed: list[tuple[str, int, float]],
    random: list[tuple[str, int, float]],
) -> tuple[list[str], dict[str, int], list[float]]:
    maximum_position: dict[str, int] = defaultdict(int)
    for chromosome, position, _ in fixed + random:
        maximum_position[chromosome] = max(maximum_position[chromosome], position)

    chromosomes = sorted(maximum_position, key=chromosome_sort_key)
    total_length = sum(maximum_position.values())
    gap = max(1, round(total_length * 0.005))
    offsets: dict[str, int] = {}
    ticks: list[float] = []
    current_offset = 0

    for chromosome in chromosomes:
        offsets[chromosome] = current_offset
        ticks.append(current_offset + maximum_position[chromosome] / 2)
        current_offset += maximum_position[chromosome] + gap

    return chromosomes, offsets, ticks


def plot_manhattan(
    variants: list[tuple[str, int, float]],
    chromosomes: list[str],
    offsets: dict[str, int],
    ticks: list[float],
    y_limit: float,
    threshold: float,
    title: str,
    output: Path,
    dpi: int,
) -> None:
    grouped: dict[str, list[tuple[int, float]]] = defaultdict(list)
    for chromosome, position, log_p in variants:
        grouped[chromosome].append((position, log_p))

    figure, axis = plt.subplots(figsize=(16, 5.5))
    for index, chromosome in enumerate(chromosomes):
        points = grouped.get(chromosome, [])
        if not points:
            continue
        x_values = [offsets[chromosome] + position for position, _ in points]
        y_values = [log_p for _, log_p in points]
        axis.scatter(
            x_values,
            y_values,
            s=5,
            alpha=0.7,
            color=COLORS[index % len(COLORS)],
            linewidths=0,
            rasterized=True,
        )

    axis.axhline(
        -math.log10(threshold),
        color="#B22222",
        linestyle="--",
        linewidth=1,
        label=f"P = {threshold:g}",
    )
    axis.set_xticks(ticks, chromosomes)
    axis.set_xlim(left=0)
    axis.set_ylim(0, y_limit)
    axis.set_xlabel("Chromosome")
    axis.set_ylabel(r"$-\log_{10}(P)$")
    axis.set_title(title)
    axis.grid(axis="y", color="#D9D9D9", linewidth=0.5, alpha=0.7)
    axis.legend(loc="upper right", frameon=False)
    figure.tight_layout()

    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def create_plots(
    fixed_path: Path,
    random_path: Path,
    output_prefix: Path,
    threshold: float = 5e-8,
    dpi: int = 180,
) -> tuple[Path, Path, int, int]:
    if not 0 < threshold <= 1:
        raise ValueError("threshold must be in (0, 1]")
    if dpi <= 0:
        raise ValueError("dpi must be positive")

    fixed, fixed_skipped = read_gwama(fixed_path)
    random, random_skipped = read_gwama(random_path)
    chromosomes, offsets, ticks = build_layout(fixed, random)
    maximum_log_p = max(value for *_, value in fixed + random)
    y_limit = max(-math.log10(threshold) + 1, math.ceil(maximum_log_p * 1.05))

    fixed_output = Path(f"{output_prefix}.fixed.manhattan.png")
    random_output = Path(f"{output_prefix}.random.manhattan.png")
    plot_manhattan(
        fixed,
        chromosomes,
        offsets,
        ticks,
        y_limit,
        threshold,
        "GWAMA fixed-effect meta-analysis",
        fixed_output,
        dpi,
    )
    plot_manhattan(
        random,
        chromosomes,
        offsets,
        ticks,
        y_limit,
        threshold,
        "GWAMA random-effects meta-analysis",
        random_output,
        dpi,
    )
    return fixed_output, random_output, fixed_skipped, random_skipped


def main() -> int:
    args = parse_args()
    try:
        fixed_output, random_output, fixed_skipped, random_skipped = create_plots(
            args.fixed,
            args.random,
            args.output_prefix,
            args.threshold,
            args.dpi,
        )
    except (OSError, ValueError) as error:
        raise SystemExit(f"ERROR: {error}") from error

    print(f"Fixed-effect Manhattan plot:  {fixed_output}")
    print(f"Random-effects Manhattan plot: {random_output}")
    if fixed_skipped or random_skipped:
        print(
            "Skipped variants with invalid coordinates/P-values: "
            f"fixed={fixed_skipped}, random={random_skipped}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
