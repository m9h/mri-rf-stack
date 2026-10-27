import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

from .fields import analyze_hertzon, load_npz
from .inverse import minimum_norm_operator, point_spread, reconstruct
from .jemris import export_jemris
from .model import Birdcage, write_birdcage


def main(argv=None):
    parser = argparse.ArgumentParser(description="MRI RF research benchmarks and coil encoding")
    commands = parser.add_subparsers(dest="command", required=True)
    model = commands.add_parser("model", help="generate native Hertzon birdcage project")
    model.add_argument("output", type=Path)
    model.add_argument("--config", type=Path, help="JSON object overriding Birdcage defaults")
    model.add_argument("--passes", type=int)
    model.add_argument("--capacitance-pf", type=float)
    solve = commands.add_parser("solve", help="run the pinned Hertzon solver in a subprocess")
    solve.add_argument("project", type=Path)
    solve.add_argument("output", type=Path)
    solve.add_argument("--design", default="Birdcage")
    solve.add_argument("--setup", default="RF")
    analyze = commands.add_parser("analyze", help="derive power-normalized B1 and local SAR")
    analyze.add_argument("results", type=Path)
    analyze.add_argument("materials", type=Path)
    analyze.add_argument("output", type=Path)
    analyze.add_argument("--label", default="LastAdaptive")
    analyze.add_argument("--frequency-index", type=int, default=0)
    analyze.add_argument(
        "--drive", nargs="+", type=complex, help="complex incident sqrt(W) per port"
    )
    analyze.add_argument("--accepted-w", type=float, default=1.0)
    analyze.add_argument("--duty-cycle", type=float, default=1.0)
    jemris = commands.add_parser("jemris", help="export relative B1 to JEMRIS ExternalCoil")
    jemris.add_argument("fields", type=Path)
    jemris.add_argument("output", type=Path)
    jemris.add_argument("--points", type=int, default=24)
    jemris.add_argument("--extent-mm", type=float, default=200)
    jemris.add_argument("--reference-b1-ut", type=float, required=True)
    jemris.add_argument("--component", choices=("plus", "minus"), default="plus")
    jemris.add_argument("--conjugate", action="store_true")
    inverse = commands.add_parser(
        "inverse", help="analyze a supplied complex coil forward operator"
    )
    inverse.add_argument(
        "input", type=Path, help="NPZ containing forward, covariance; optional data"
    )
    inverse.add_argument("output", type=Path)
    inverse.add_argument("--alpha", type=float, default=0.1)
    inverse.add_argument("--source", type=int, default=0)
    args = parser.parse_args(argv)
    try:
        if args.command == "model":
            config = json.loads(args.config.read_text()) if args.config else {}
            if args.passes is not None:
                config["passes"] = args.passes
            if args.capacitance_pf is not None:
                config["capacitance_f"] = args.capacitance_pf * 1e-12
            print(write_birdcage(args.output, Birdcage(**config)))
        elif args.command == "solve":
            from importlib.metadata import version

            if version("hertzon") != "0.1.25" or version("emerge") != "3.0.0a25":
                raise ValueError("solver must be Hertzon 0.1.25 and EMerge 3.0.0a25")
            if args.output.exists() and any(args.output.iterdir()):
                raise ValueError("solver output directory must be new or empty")
            args.output.mkdir(parents=True, exist_ok=True)
            env = os.environ.copy()
            for key, sub in (
                ("MPLCONFIGDIR", "mpl"),
                ("NUMBA_CACHE_DIR", "numba"),
                ("XDG_CACHE_HOME", "cache"),
            ):
                env.setdefault(key, str((args.output / sub).resolve()))
            command = [
                sys.executable,
                "-m",
                "hertzon.solver.run",
                str(args.project.resolve()),
                args.design,
                args.setup,
                str(args.output.resolve()),
            ]
            with (args.output / "solver.log").open("w") as log:
                result = subprocess.run(
                    command, check=False, env=env, stdout=log, stderr=subprocess.STDOUT
                )
            if result.returncode:
                raise ValueError(f"Hertzon failed; see {args.output / 'solver.log'}")
            status = json.loads((args.output / "status.json").read_text())
            if status.get("state") != "done":
                raise ValueError("Hertzon exited without a completed result")
            print(args.output)
        elif args.command == "analyze":
            report = analyze_hertzon(
                args.results,
                json.loads(args.materials.read_text()),
                args.output,
                label=args.label,
                frequency_index=args.frequency_index,
                drive=args.drive,
                target_accepted_w=args.accepted_w,
                duty_cycle=args.duty_cycle,
            )
            print(json.dumps(report, indent=2))
        elif args.command == "jemris":
            print(
                json.dumps(
                    export_jemris(
                        args.fields,
                        args.output,
                        points=args.points,
                        extent_mm=args.extent_mm,
                        reference_b1_t=args.reference_b1_ut * 1e-6,
                        component=args.component,
                        conjugate=args.conjugate,
                    ),
                    indent=2,
                )
            )
        else:
            data = load_npz(args.input)
            w, info = minimum_norm_operator(data["forward"], data["covariance"], args.alpha)
            result = dict(info, point_spread=point_spread(w, data["forward"], args.source))
            if "data" in data:
                y = data["data"]
                if y.ndim not in (1, 2) or y.shape[0] != w.shape[1] or not np.isfinite(y).all():
                    raise ValueError("data must be finite (channels,) or (channels,time)")
                result["reconstruction"] = reconstruct(w, y)
            args.output.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(args.output, **result)
            print(f"effective degrees of freedom: {info['effective_dof']:.4g}")
    except (ValueError, OSError, KeyError, ImportError, TypeError) as exc:
        parser.exit(2, f"mri-rf: {exc}\n")


if __name__ == "__main__":
    main()
