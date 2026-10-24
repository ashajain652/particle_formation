"""Apply sub-plan 18 Task 8 (flags, run names, columns, frame fields) to the scratch copy."""
import sys

root = sys.argv[1]


def patch(path, pairs):
    s = open(path).read()
    for old, new in pairs:
        assert s.count(old) == 1, (path, old[:80], s.count(old))
        s = s.replace(old, new)
    open(path, "w").write(s)


patch(root + "/reentry_model/cli.py", [
    ("from . import __version__, aero, atmosphere, body, compare, coupled, dispersion, fap, heating, material, mesh, sesam_io, spray, surface_flow, thermal, viz\n",
     "from . import __version__, aero, atmosphere, body, compare, coupled, dispersion, fap, heating, material, mesh, sesam_io, skin, spray, surface_flow, thermal, viz\n"),
    ("""                   deep_runoff=True, molten_cascade=True, seed=DEFAULT_SEED, rigid_substrate=True):
    \"\"\"The run name encodes the configuration, so runs of different settings can share an output directory. A melting
    run with the deep runoff off (amendment of 2026-10-02) ends in `_deeprunoff-off`, one with the molten cascade off
    (amendment of 2026-10-03) in `_moltencascade-off`, one with the rigid substrate off (amendment of 2026-10-06) in
    `_rigidsubstrate-off`, and any run with a seed other than DEFAULT_SEED (amendment of 2026-10-05) in `_seed-<n>`;
    the defaults keep the old name.\"\"\"
""", """                   deep_runoff=True, molten_cascade=True, seed=DEFAULT_SEED, rigid_substrate=True, surface_model="elements",
                   skin_suffix=""):
    \"\"\"The run name encodes the configuration, so runs of different settings can share an output directory. A melting
    run with the deep runoff off (amendment of 2026-10-02) ends in `_deeprunoff-off`, one with the molten cascade off
    (amendment of 2026-10-03) in `_moltencascade-off`, one with the rigid substrate off (amendment of 2026-10-06) in
    `_rigidsubstrate-off`, one with skins (sub-plan 18) in `_surface-skin` followed by `skin_suffix` (the skins'
    non-default settings, `skin_name_suffix`), and any run with a seed other than DEFAULT_SEED (amendment of
    2026-10-05) in `_seed-<n>`; the defaults keep the old name.\"\"\"
"""),
    ("""            + ("_rigidsubstrate-off" if melt and not rigid_substrate else "")
            + ("_seed-{}".format(seed) if seed != DEFAULT_SEED else ""))
""", """            + ("_rigidsubstrate-off" if melt and not rigid_substrate else "")
            + ("_surface-skin" + skin_suffix if melt and surface_model == "skin" else "")
            + ("_seed-{}".format(seed) if seed != DEFAULT_SEED else ""))


def skin_name_suffix(thickness_mm, cell_um, substep_ms, preset):
    \"\"\"The run-name suffix of a skin run's non-default settings (sub-plan 18): `_skin-dev` for the dev preset,
    `_skin<mm>mm` for another thickness, `_cell<um>um` for an explicit cell other than the preset's, `_sub<ms>ms` for
    another sub-step; empty for the defaults.\"\"\"
    out = "_skin-dev" if preset == "dev" else ""
    if thickness_mm != 0.4:
        out += "_skin{:g}mm".format(thickness_mm)
    if cell_um is not None and abs(cell_um * 1e-6 - skin.SKIN_PRESETS[preset]) > 1e-12:
        out += "_cell{:g}um".format(cell_um)
    if substep_ms != 10.0:
        out += "_sub{:g}ms".format(substep_ms)
    return out
"""),
    ("""    me.add_argument("--rt-spray", choices=("on", "off"), default="on",
""", """    me.add_argument("--surface-model", choices=body.SURFACE_MODEL_NAMES, default="elements",
                    help="skin: a one-dimensional skin of fine cells under every melting patch resolves the melt layer and "
                         "feeds the film continuously, coupled to the 3D conduction through an implicit interface (Step 3 "
                         "sub-plan 18). elements: the surface elements melt and feed the film, as in every run before it (default)")
    me.add_argument("--skin-thickness", type=float, default=0.4, help="skin target thickness [mm] (default %(default)s)")
    me.add_argument("--skin-cell", type=float, default=None, help="skin cell size [um] (default: the preset's)")
    me.add_argument("--skin-substep", type=float, default=10.0, help="skin sub-step [ms] (default %(default)s)")
    me.add_argument("--skin-preset", choices=tuple(skin.SKIN_PRESETS), default="production",
                    help="skin cell-size preset: production 10 um, dev 20 um (default %(default)s)")
    me.add_argument("--rt-spray", choices=("on", "off"), default="on",
"""),
    ("""                                          rigid_substrate=args.rigid_substrate == "on")
        the_body = body.MeltingBody(""", """                                          rigid_substrate=args.rigid_substrate == "on", surface_model=args.surface_model,
                                          skin=skin_settings(args))
        the_body = body.MeltingBody("""),
    ("""                     "molten_cascade": args.molten_cascade, "rigid_substrate": args.rigid_substrate, "rt_spray": args.rt_spray,
""", """                     "molten_cascade": args.molten_cascade, "rigid_substrate": args.rigid_substrate, "rt_spray": args.rt_spray,
                     "surface_model": args.surface_model,
                     "skin": {"thickness_mm": args.skin_thickness, "cell_um": skin_settings(args).cell * 1e6,
                              "substep_ms": args.skin_substep, "preset": args.skin_preset},
"""),
    ("""    return the_body, heating_model, info


def cmd_run(args, parser):
""", """    return the_body, heating_model, info


def skin_settings(args):
    \"\"\"The skins' settings from the command line (sub-plan 18).\"\"\"
    cell = args.skin_cell * 1e-6 if args.skin_cell is not None else skin.SKIN_PRESETS[args.skin_preset]
    return skin.SkinSettings(thickness=args.skin_thickness * 1e-3, cell=cell, substep=args.skin_substep * 1e-3)


def cmd_run(args, parser):
"""),
    ("""    if not 0.0 < args.demise_fraction < 1.0:
        parser.error("--demise-fraction must be within (0, 1)")
""", """    if not 0.0 < args.demise_fraction < 1.0:
        parser.error("--demise-fraction must be within (0, 1)")
    for label, value in (("--skin-thickness", args.skin_thickness), ("--skin-substep", args.skin_substep),
                         ("--skin-cell", 1.0 if args.skin_cell is None else args.skin_cell)):
        if value <= 0.0:
            parser.error("{} must be > 0".format(label))
    if 2.0 * skin_settings(args).cell > args.skin_thickness * 1e-3 * (1.0 + 1e-12):
        parser.error("--skin-cell must be at most half of --skin-thickness")
"""),
    ("""                                       rigid_substrate=args.rigid_substrate == "on")
    run_dir = os.path.join(args.outdir, name)
""", """                                       rigid_substrate=args.rigid_substrate == "on", surface_model=args.surface_model,
                                       skin_suffix=skin_name_suffix(args.skin_thickness, args.skin_cell, args.skin_substep,
                                                                    args.skin_preset))
    run_dir = os.path.join(args.outdir, name)
"""),
])

patch(root + "/reentry_model/coupled.py", [
    ("""                "nonrigid_depth_mean_mm", "slurry_thick_fraction", "rigid_thin_fraction", "slurry_held_mass_kg"]
""", """                "nonrigid_depth_mean_mm", "slurry_thick_fraction", "rigid_thin_fraction", "slurry_held_mass_kg"]
SKIN_COLUMNS = ["skin_count", "skin_mass_kg", "skin_thickness_min_mm", "skin_thickness_median_mm", "skin_drawn_mass_kg",
                "skin_handover_mass_kg", "skin_short_events", "interface_mismatch_J", "interface_repeats",
                "liquid_depth_mean_um", "liquid_depth_max_um", "nonrigid_depth_skin_mean_um", "skin_wall_s"]   # sub-plan 18
MELT_COLUMNS = MELT_COLUMNS + SKIN_COLUMNS
"""),
    ("""        poly.cell_data["nonrigid_depth"] = carry(body.last_nonrigid, np.nan)
""", """        poly.cell_data["nonrigid_depth"] = carry(body.last_nonrigid, np.nan)
        # the skins (sub-plan 18): top temperature, resolved liquid and non-rigid depths, thickness; nan without a skin
        sk = getattr(body, "skins", None)
        if sk is not None:
            on = body.skin_of_patch >= 0
            r = body.skin_of_patch[on]
            d_liq, _ = sk.depth_above(body.material.T_feed)
            d_nr, _ = sk.depth_above(body.material.T_rigid)
            for key, v in (("skin_T_top", sk.T[:, 0] if sk.n else None), ("skin_liquid_depth", d_liq),
                           ("skin_nonrigid_depth", d_nr), ("skin_thickness", sk.thickness())):
                out = np.full(surface.n_patches, np.nan)
                if sk.n:
                    out[on] = v[r]
                poly.cell_data[key] = out
"""),
    ("""                "removed_enthalpy_J": body.removed_enthalpy, "frozen_mass_kg": body.frozen_mass,
                "melt_energy_balance_residual": body.energy_balance_residual()}
""", """                "removed_enthalpy_J": body.removed_enthalpy, "frozen_mass_kg": body.frozen_mass,
                "skins_created": body.skins_created, "skin_drawn_mass_kg": body.skin_drawn_mass,
                "skin_handover_mass_kg": body.skin_handover_mass, "interface_repeats": body.interface_repeats,
                "skin_wall_s": body.skin_wall,
                "melt_energy_balance_residual": body.energy_balance_residual()}
"""),
])

patch(root + "/reentry_model/body.py", [
    ("""                "n_dead_elements": float(self.mesh.n_elements - self.mesh.n_active)}
""", """                "n_dead_elements": float(self.mesh.n_elements - self.mesh.n_active), **self._skin_stats()}

    def _skin_stats(self):
        \"\"\"The skins' history columns (sub-plan 18); zero or nan without skins.\"\"\"
        sk = self.skins
        th = sk.thickness() if sk is not None and sk.n else np.zeros(0)
        on = self.skin_of_patch >= 0
        liq = self.last_skin_liquid if self.last_skin_liquid is not None and len(self.last_skin_liquid) == len(on) else None
        nr = self.last_skin_nonrigid if self.last_skin_nonrigid is not None and len(self.last_skin_nonrigid) == len(on) else None
        some = bool(on.any())
        return {"skin_count": float(sk.n) if sk is not None else 0.0, "skin_mass_kg": sk.mass() if sk is not None else 0.0,
                "skin_thickness_min_mm": float(th.min() * 1e3) if th.size else float("nan"),
                "skin_thickness_median_mm": float(np.median(th) * 1e3) if th.size else float("nan"),
                "skin_drawn_mass_kg": self.skin_drawn_mass, "skin_handover_mass_kg": self.skin_handover_mass,
                "skin_short_events": float(self.skin_short_events), "interface_mismatch_J": self.last_mismatch,
                "interface_repeats": float(self.interface_repeats),
                "liquid_depth_mean_um": float(liq[on].mean() * 1e6) if liq is not None and some else float("nan"),
                "liquid_depth_max_um": float(liq[on].max() * 1e6) if liq is not None and some else float("nan"),
                "nonrigid_depth_skin_mean_um": float(nr[on].mean() * 1e6) if nr is not None and some else float("nan"),
                "skin_wall_s": self.skin_wall}
"""),
])
print("task 8 patched")
