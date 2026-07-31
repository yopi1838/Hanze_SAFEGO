# -*- coding: ascii -*-
"""
Central path resolution for the Hanze_SAFEGO post-processing scripts.

Everything resolves relative to the repository root (the folder containing
this file), so a fresh clone works with no editing. Each location can be
overridden by an environment variable if you keep data outside the repo.

  SAFEGO_EXP_DATA   experimental data folder   (default: <repo>/EXP_DATA)
  SAFEGO_SIM_DIR    primary simulation folder  (default: <repo>/stratC_results_NODAMP_v6_NEW)
  SAFEGO_CACHE      staged-script JSON cache   (default: <repo>/.cache)

Derived experimental CSVs are resolved with exp_derived(), which falls back to
the canonical copies under stratC_results_NODAMP_v6_NEW/postproc/.

Style note: this module is deliberately Python-2/3 compatible (no f-strings)
to match the rest of the codebase, some of which runs inside 3DEC's embedded
interpreter.
"""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent

DEFAULT_SIM = "stratC_results_NODAMP_v6_NEW"


def _env_path(var):
    """Read a path from the environment. Relative values are taken as relative to
    the repo root, not the current working directory, so
    SAFEGO_SIM_DIR=stratC_results_NODAMP_v7 works from anywhere."""
    v = os.environ.get(var)
    if not v:
        return None
    p = Path(os.path.expanduser(v))
    if not p.is_absolute():
        p = ROOT / p
    return p.resolve()


def exp_data_dir():
    """Raw experimental data: Test9Run*.xlsx, Test12Run*.xlsx, Test*_Info.xlsx,
    US1_Tilt_values.csv, US2_Tilt_values.csv, processed_globalzero/."""
    return _env_path("SAFEGO_EXP_DATA") or (ROOT / "EXP_DATA")


def processed_dir():
    """Processed shake-table exports (U_avg, base_shear_kN)."""
    return exp_data_dir() / "processed_globalzero"


def tilt_csv(specimen="US1"):
    """Tiltmeter logger export for US1 or US2."""
    return exp_data_dir() / "{}_Tilt_values.csv".format(specimen)


def sim_dir(name=None):
    """A simulation results folder. With no argument, the primary one."""
    if name is None:
        return _env_path("SAFEGO_SIM_DIR") or (ROOT / DEFAULT_SIM)
    return ROOT / name


def postproc_dir(name=None):
    """Derived CSVs and figures live here. Created on demand."""
    d = sim_dir(name) / "postproc"
    if not d.exists():
        try:
            os.makedirs(str(d))
        except OSError:
            pass
    return d


def exp_derived(name):
    """A *derived experimental* CSV -- exp_Test{9,12}_metrics.csv,
    exp_Test*_tilt.csv, exp_Test*_period_psd.csv, exp_US*_tiltmeter.csv.

    These describe the physical specimens, not any one simulation, but they
    live inside a simulation's postproc/ folder for historical reasons. The
    canonical copies sit in DEFAULT_SIM's postproc/ and several cannot be
    regenerated from this repo (see CLAUDE.md section 10).

    Look in the active sim's postproc/ first, then fall back to the canonical
    copy, so pointing SAFEGO_SIM_DIR at a fresh results folder keeps the
    experimental comparison instead of failing on a missing file."""
    active = postproc_dir() / name
    if active.is_file():
        return active
    return postproc_dir(DEFAULT_SIM) / name


def cache_dir():
    """Staged scripts (ch19_xval, profile_hysteresis, exp_period_wrapper,
    exp_tilt_from_raw) keep intermediate JSON state here between stages.

    This used to be hard-coded to /tmp, which does not exist on Windows and
    does not survive a reboot. It is git-ignored."""
    d = _env_path("SAFEGO_CACHE") or (ROOT / ".cache")
    if not d.exists():
        try:
            os.makedirs(str(d))
        except OSError:
            pass
    return d


def state_file(name):
    """Path to a named staged-script cache file, e.g. state_file('ch19')."""
    return cache_dir() / "{}_state.json".format(name)


if __name__ == "__main__":
    print("repo root      : {}".format(ROOT))
    print("EXP_DATA       : {}  (exists: {})".format(exp_data_dir(), exp_data_dir().is_dir()))
    print("processed      : {}  (exists: {})".format(processed_dir(), processed_dir().is_dir()))
    print("US1 tilt csv   : {}  (exists: {})".format(tilt_csv("US1"), tilt_csv("US1").is_file()))
    print("US2 tilt csv   : {}  (exists: {})".format(tilt_csv("US2"), tilt_csv("US2").is_file()))
    print("sim dir        : {}  (exists: {})".format(sim_dir(), sim_dir().is_dir()))
    print("postproc       : {}".format(postproc_dir()))
    print("cache          : {}".format(cache_dir()))
