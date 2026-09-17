#!/usr/bin/env python3
"""Create a marker-level comparison of fixed- and random-effects GWAMA output."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Iterable


MARKER_COLUMNS = ("rs_number", "MARKERNAME", "markername", "SNP", "ID")
EFFECT_COLUMNS = ("OR", "BETA", "beta", "effect")
SE_COLUMNS = ("OR_se", "BETA_se", "beta_se", "SE", "se")
P_COLUMNS = ("p-value", "P", "p", "p_value", "pvalue")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare fixed- and random-effects GWAMA result files."
    )
    parser.add_argument("--fixed", required=True, type=Path)
    parser.add_argument("--random", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def first_present(fieldnames: Iterable[str], candidates: Iterable[str]) -> str | None:
    names = set(fieldnames)
    return next((candidate for candidate in candidates if candidate in names), None)


def read_gwama(path: Path) -> tuple[dict[str, dict[str, str]], dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        header_line = handle.readline()
        if not header_line:
            raise ValueError(f"empty GWAMA output: {path}")

        fieldnames = header_line.split()
        marker_column = first_present(fieldnames, MARKER_COLUMNS)
        effect_column = first_present(fieldnames, EFFECT_COLUMNS)
        se_column = first_present(fieldnames, SE_COLUMNS)
        p_column = first_present(fieldnames, P_COLUMNS)

        required = {
            "marker": marker_column,
            "effect": effect_column,
            "se": se_column,
            "p": p_column,
        }
        missing = [name for name, column in required.items() if column is None]
        if missing:
            raise ValueError(
                f"{path} is missing recognizable {', '.join(missing)} column(s); "
                f"found: {', '.join(fieldnames)}"
            )

        rows: dict[str, dict[str, str]] = {}
        for line_number, line in enumerate(handle, start=2):
            if not line.strip():
                continue
            values = line.split()
            if len(values) != len(fieldnames):
                raise ValueError(
                    f"{path}:{line_number}: expected {len(fieldnames)} fields, "
                    f"found {len(values)}"
                )
            row = dict(zip(fieldnames, values))
            marker = row[marker_column]  # type: ignore[index]
            if marker in rows:
                raise ValueError(f"{path}:{line_number}: duplicate marker {marker}")
            rows[marker] = row

    columns = {
        "marker": marker_column,
        "effect": effect_column,
        "se": se_column,
        "p": p_column,
    }
    return rows, columns  # type: ignore[return-value]


def optional_value(row: dict[str, str] | None, *columns: str) -> str:
    if row is None:
        return ""
    for column in columns:
        if column in row:
            return row[column]
    return ""


def absolute_difference(left: str, right: str) -> str:
    if not left or not right:
        return ""
    try:
        difference = abs(float(left) - float(right))
    except ValueError:
        return ""
    if not math.isfinite(difference):
        return ""
    return f"{difference:.12g}"


def compare(fixed_path: Path, random_path: Path, output_path: Path) -> int:
    fixed_rows, fixed_columns = read_gwama(fixed_path)
    random_rows, random_columns = read_gwama(random_path)

    if fixed_columns["effect"].lower() != random_columns["effect"].lower():
        raise ValueError(
            "fixed and random outputs use different effect types: "
            f"{fixed_columns['effect']} versus {random_columns['effect']}"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "MARKERNAME",
        "FFX_EFFECT",
        "FFX_SE",
        "FFX_P",
        "RFX_EFFECT",
        "RFX_SE",
        "RFX_P",
        "EFFECT_ABS_DIFF",
        "P_ABS_DIFF",
        "Q_STATISTIC",
        "Q_P_VALUE",
        "I2",
        "N_STUDIES",
        "STATUS",
    ]

    markers = sorted(set(fixed_rows) | set(random_rows))
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for marker in markers:
            fixed = fixed_rows.get(marker)
            random = random_rows.get(marker)
            fixed_effect = optional_value(fixed, fixed_columns["effect"])
            random_effect = optional_value(random, random_columns["effect"])
            fixed_p = optional_value(fixed, fixed_columns["p"])
            random_p = optional_value(random, random_columns["p"])
            source = random or fixed

            if fixed is not None and random is not None:
                status = "both"
            elif fixed is not None:
                status = "fixed_only"
            else:
                status = "random_only"

            writer.writerow(
                {
                    "MARKERNAME": marker,
                    "FFX_EFFECT": fixed_effect,
                    "FFX_SE": optional_value(fixed, fixed_columns["se"]),
                    "FFX_P": fixed_p,
                    "RFX_EFFECT": random_effect,
                    "RFX_SE": optional_value(random, random_columns["se"]),
                    "RFX_P": random_p,
                    "EFFECT_ABS_DIFF": absolute_difference(fixed_effect, random_effect),
                    "P_ABS_DIFF": absolute_difference(fixed_p, random_p),
                    "Q_STATISTIC": optional_value(source, "q_statistic", "Q"),
                    "Q_P_VALUE": optional_value(source, "q_p-value", "Q_P"),
                    "I2": optional_value(source, "i2", "I2"),
                    "N_STUDIES": optional_value(source, "n_studies", "N_STUDIES"),
                    "STATUS": status,
                }
            )

    return len(markers)


def main() -> int:
    args = parse_args()
    try:
        count = compare(args.fixed, args.random, args.output)
    except (OSError, ValueError) as error:
        raise SystemExit(f"ERROR: {error}") from error
    print(f"Compared {count} markers: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
