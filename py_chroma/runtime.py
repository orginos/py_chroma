from __future__ import annotations

from ._core import (
    finalize,
    initialize,
    random_gauge_transform,
    rect_wilson_loops,
    register_inline,
    run_hmc_xml,
    run_inline_xml,
    run_plaquette,
    run_smd_xml,
    set_gauge,
    set_input_xml,
    set_lattice,
    set_rng_seed,
    set_rng_seed_xml,
)
from .params import HMCParams, SMDParams


def run_hmc(params: HMCParams) -> None:
    if not isinstance(params, HMCParams):
        raise TypeError("params must be an HMCParams instance")
    run_hmc_xml(params.to_xml())


def run_smd(params: SMDParams) -> None:
    if not isinstance(params, SMDParams):
        raise TypeError("params must be an SMDParams instance")
    run_smd_xml(params.to_xml())


__all__ = [
    "finalize",
    "initialize",
    "random_gauge_transform",
    "rect_wilson_loops",
    "register_inline",
    "run_hmc",
    "run_hmc_xml",
    "run_inline_xml",
    "run_plaquette",
    "run_smd",
    "run_smd_xml",
    "set_gauge",
    "set_input_xml",
    "set_lattice",
    "set_rng_seed",
    "set_rng_seed_xml",
]
