"""Workflow mode guard helpers."""

from __future__ import annotations

import argparse
from pathlib import Path

from .utils import read_yaml


def check_reference_inputs(config_path: Path, output: Path) -> None:
    config = read_yaml(config_path)
    reference_inputs = config.get("reference_inputs", {})
    missing = []
    for name, value in reference_inputs.items():
        if str(value) in {"", "TODO", "NOT_ASSESSED"}:
            missing.append(f"{name}=TODO")
        elif not Path(value).exists():
            missing.append(f"{name}={value}")
    if missing:
        raise SystemExit(
            "reference_only mode requires user-provided local files before execution. "
            f"Missing or unresolved inputs: {', '.join(missing)}"
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("reference_only input check passed\n", encoding="utf-8")


def fail_full_mode() -> None:
    raise SystemExit(
        "full mode is intentionally not implemented in the MVP. "
        "Use demo, metadata_only, or reference_only with local files."
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Workflow mode guards.")
    subparsers = parser.add_subparsers(dest="command", required=True)
    reference = subparsers.add_parser("check-reference-inputs")
    reference.add_argument("--config", type=Path, required=True)
    reference.add_argument("--output", type=Path, required=True)
    subparsers.add_parser("fail-full")
    args = parser.parse_args()
    if args.command == "check-reference-inputs":
        check_reference_inputs(args.config, args.output)
    elif args.command == "fail-full":
        fail_full_mode()


if __name__ == "__main__":
    main()

