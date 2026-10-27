"""Generate a native Hertzon project: an eight-rung low-pass birdcage.

All dimensions in the public configuration are SI. PEC sheets and shield are an
intentional first benchmark approximation; capacitor value is a tuning input.
"""

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from .physics import PROTON_HZ_PER_T


@dataclass
class Birdcage:
    b0_t: float = 3.0
    rungs: int = 8
    radius_m: float = 0.15
    length_m: float = 0.26
    strip_width_m: float = 0.012
    gap_m: float = 0.008
    capacitance_f: float = 30e-12
    shield_radius_m: float = 0.20
    shield_length_m: float = 0.42
    phantom_radius_m: float = 0.09
    phantom_length_m: float = 0.18
    epsilon_r: float = 60.0
    sigma_s_m: float = 0.5
    density_kg_m3: float = 1000.0
    air_margin_m: float = 0.12
    mesh_m: float = 0.035
    passes: int = 12

    def validate(self):
        vals = asdict(self)
        if any(not np.isfinite(v) or v <= 0 for v in vals.values()):
            raise ValueError("all benchmark parameters must be finite and positive")
        if self.rungs < 4 or self.rungs % 4 or int(self.rungs) != self.rungs:
            raise ValueError("rungs must be a multiple of four, at least four")
        if self.passes < 1 or int(self.passes) != self.passes:
            raise ValueError("passes must be a positive integer")
        if not self.phantom_radius_m < self.radius_m - self.strip_width_m:
            raise ValueError("phantom must fit inside coil")
        if not self.radius_m + self.strip_width_m < self.shield_radius_m:
            raise ValueError("coil must fit inside shield")
        if not self.phantom_length_m < self.length_m < self.shield_length_m:
            raise ValueError("expected phantom length < coil length < shield length")
        if self.gap_m >= self.length_m / 8:
            raise ValueError("capacitor/feed gaps overlap")
        if self.strip_width_m >= 2 * self.radius_m * np.sin(np.pi / self.rungs):
            raise ValueError("rung strips overlap")


def build_birdcage(config: Birdcage):
    from hertzon.model import schema
    from hertzon.model.design import Assignment, Command, Project, Setup, Sweep
    from hertzon.model.materials import Material
    from hertzon.model.validate import validate

    config.validate()
    c = config
    p = Project("MRI_Birdcage_3T")
    d = p.designs[0]
    d.name = "Birdcage"
    d.units = "mm"
    d.settings.update(AutoOpen="Auto Open - Filter", FilterAirMargin=f"{c.air_margin_m}m")
    p.materials.append(
        Material(
            "phantom",
            permittivity=str(c.epsilon_r),
            conductivity=str(c.sigma_s_m),
            color="#d789ae",
            description="Synthetic homogeneous phantom, not anatomical tissue",
        )
    )
    d.add_object(
        "CreateCylinder",
        name="Phantom",
        material="phantom",
        CenterPosition=f"0,0,{-c.phantom_length_m / 2}m",
        Axis="Z",
        Radius=f"{c.phantom_radius_m}m",
        Height=f"{c.phantom_length_m}m",
    )
    shield = d.add_object(
        "CreateCylinder",
        name="Shield",
        material="vacuum",
        CenterPosition=f"0,0,{-c.shield_length_m / 2}m",
        Radius=f"{c.shield_radius_m}m",
        Height=f"{c.shield_length_m}m",
        Axis="Z",
    )
    # OCC's primitive cylinder has lateral face 0, followed by its two caps.
    # Keep only that sheet: an open PEC RF shield, with no thin volume to mesh.
    shield.history.append(Command.new("KeepFaces", Faces="0"))

    def assign(kind, obj, name, excitation=False, **kwargs):
        registry = schema.EXCITATIONS if excitation else schema.BOUNDARIES
        props = schema.defaults(registry[kind])
        props.update({key: str(val) for key, val in kwargs.items()})
        assignment = Assignment(name, kind, props, objects=[obj.name])
        (d.excitations if excitation else d.boundaries).append(assignment)

    assign("Perfect E", shield, "PEC_Shield")

    # Planar annular end rings; vertical strips intersect them along shared edges.
    for end, z in enumerate((-c.length_m / 2, c.length_m / 2)):
        ring = d.add_object(
            "CreateCircle",
            name=f"Ring{end}",
            CenterPosition=f"0,0,{z}m",
            Axis="Z",
            Radius=f"{c.radius_m + c.strip_width_m / 2}m",
        )
        hole = d.add_object(
            "CreateCircle",
            name=f"RingHole{end}",
            CenterPosition=f"0,0,{z}m",
            Axis="Z",
            Radius=f"{c.radius_m - c.strip_width_m / 2}m",
        )
        ring.history.append(Command.new("Subtract", Tools=hole.name))
        hole.consumed_by = ring.name
        assign("Perfect E", ring, f"PEC_Ring{end}")

    def strip(name, theta, z1, z2):
        radial = np.array([np.cos(theta), np.sin(theta)]) * c.radius_m
        tangent = np.array([-np.sin(theta), np.cos(theta)]) * c.strip_width_m / 2
        corners = [
            (*xy, z)
            for xy, z in (
                (radial - tangent, z1),
                (radial + tangent, z1),
                (radial + tangent, z2),
                (radial - tangent, z2),
            )
        ]
        points = ";".join(",".join(f"{v:.12g}m" for v in point) for point in corners)
        obj = d.add_object(
            "CreatePolyline", name=name, Points=points, Closed="True", Covered="True"
        )
        line = ";".join(f"{radial[0]:.12g}m,{radial[1]:.12g}m,{z:.12g}m" for z in (z1, z2))
        return obj, line

    for rung in range(c.rungs):
        theta = 2 * np.pi * rung / c.rungs
        gaps = [(0.0, "cap")]
        if rung == 0:
            gaps.insert(0, (-c.length_m / 4, "feed"))
        last = -c.length_m / 2
        for index, (z, kind) in enumerate(gaps):
            lo, hi = z - c.gap_m / 2, z + c.gap_m / 2
            conductor, _ = strip(f"Rung{rung}_{index}", theta, last, lo)
            assign("Perfect E", conductor, f"PEC_{conductor.name}")
            gap, line = strip(f"{kind}_{rung}", theta, lo, hi)
            if kind == "cap":
                assign(
                    "Lumped RLC",
                    gap,
                    f"C{rung}",
                    UseResistance="False",
                    UseCapacitance="True",
                    Capacitance=f"{c.capacitance_f}F",
                    IntegrationLine=line,
                )
            else:
                assign(
                    "Lumped Port",
                    gap,
                    "Feed",
                    excitation=True,
                    PortNumber=1,
                    Impedance="50ohm",
                    IntegrationLine=line,
                )
            last = hi
        conductor, _ = strip(f"Rung{rung}_end", theta, last, c.length_m / 2)
        assign("Perfect E", conductor, f"PEC_{conductor.name}")

    mesh_props = schema.defaults(schema.MESH_OPS["Length Based Inside Selection"])
    mesh_props["MaxLength"] = f"{c.mesh_m}m"
    d.mesh_ops.append(
        Assignment("PhantomMesh", "Length Based Inside Selection", mesh_props, objects=["Phantom"])
    )
    setup = Setup("RF")
    setup.props.update(
        Frequency=f"{c.b0_t * PROTON_HZ_PER_T}Hz",
        MaximumPasses=str(c.passes),
        MaximumRefinement="60000",
        MinimumPasses="2",
        MaximumDeltaS="0.02",
        Solver="SuperLU",
        Threads="2",
        CurvedElements="False",
        SaveFields="True",
        BasisOrder="Zero Order",
        LambdaTarget="0.1",
    )
    sweep = Sweep("Resonance")
    sweep.props.update(
        Enabled="False",
        SweepType="Discrete",
        Distribution="Linear Count",
        Start="90MHz",
        Stop="170MHz",
        Count="17",
        SaveFields="False",
    )
    setup.sweeps.append(sweep)
    d.setups.append(setup)
    errors = [check.message for check in validate(p, d) if check.level == "error"]
    if errors:
        raise ValueError("; ".join(errors))
    return p


def write_birdcage(folder, config):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    project = build_birdcage(config)
    project.save(folder / "birdcage.hzn")
    (folder / "config.json").write_text(json.dumps(asdict(config), indent=2) + "\n")
    materials = {"Phantom": {"sigma_s_m": config.sigma_s_m, "density_kg_m3": config.density_kg_m3}}
    (folder / "materials.json").write_text(json.dumps(materials, indent=2) + "\n")
    return folder / "birdcage.hzn"
