import argparse
import os
import re
import xml.etree.ElementTree as ET
from pathlib import Path

from py_chroma import finalize, initialize
from py_chroma.hmc import (
    CloverAction,
    GaugeAction,
    GaugeConfig,
    GaugeMonomial,
    GaugeState,
    Hamiltonian,
    HMCParams,
    HMCTrj,
    InvertParam,
    Integrator,
    MCControl,
    MDIntegrator,
    SchroedingerFermBC,
    SchroedingerGaugeBC,
    TwoFlavorEOPrecLogDetFermMonomial,
    predictor_last_solution,
    predictor_linear_extrapolation,
    predictor_mre,
    predictor_mre_initcg,
    predictor_zero_guess,
    run_hmc,
)


def parse_nrow(text):
    parts = [int(p) for p in text.split(",")]
    if len(parts) != 4:
        raise ValueError("nrow must have 4 comma-separated integers")
    return parts


def parse_seed(text):
    parts = [int(p) for p in text.split(",")]
    if len(parts) != 4:
        raise ValueError("seed must have 4 comma-separated integers")
    return parts


def parse_theta(text):
    parts = [float(p) for p in text.split(",")]
    if len(parts) != 3:
        raise ValueError("theta must have 3 comma-separated floats")
    return parts


def find_latest_cfg(save_prefix: str):
    prefix_path = Path(save_prefix)
    search_dir = prefix_path.parent if prefix_path.parent != Path(".") else Path.cwd()
    prefix_name = prefix_path.name
    patterns = [
        re.compile(rf"^{re.escape(prefix_name)}_cfg_(\d+)"),
        re.compile(rf"^{re.escape(prefix_name)}[._-](\d+)"),
        re.compile(rf"^{re.escape(prefix_name)}(\d+)"),
    ]
    candidates = []
    for entry in search_dir.iterdir():
        if not entry.is_file():
            continue
        update_num = None
        for pattern in patterns:
            match = pattern.match(entry.name)
            if match:
                update_num = int(match.group(1))
                break
        if update_num is None:
            continue
        candidates.append((entry, update_num))
    if not candidates:
        return None, None
    candidates.sort(key=lambda item: item[1], reverse=True)
    cfg_candidates = [item for item in candidates if "cfg" in item[0].name]
    entry, update_num = (cfg_candidates or candidates)[0]
    return entry, update_num


def find_restart_xml(save_prefix: str, update_num: int):
    prefix_path = Path(save_prefix)
    search_dir = prefix_path.parent if prefix_path.parent != Path(".") else Path.cwd()
    prefix_name = prefix_path.name
    candidate = search_dir / f"{prefix_name}_restart_{update_num}.xml"
    if candidate.is_file():
        return candidate
    return None


def read_rng_seed(restart_xml: Path):
    tree = ET.parse(restart_xml)
    root = tree.getroot()
    seed_elems = root.findall(".//MCControl/RNG/Seed/elem")
    if len(seed_elems) != 4:
        raise ValueError("Restart RNG seed does not have 4 elements.")
    return [int(elem.text) for elem in seed_elems]


def kappa_from_mass(mass):
    return 1.0 / (2.0 * (mass + 4.0))


def main():
    parser = argparse.ArgumentParser(
        description="Run Nf=2 Clover HMC with Schroedinger BC via py_chroma."
    )
    parser.add_argument("--beta", type=float, required=True, help="Gauge coupling beta")
    parser.add_argument("--mass", type=float, default=None, help="Bare mass m0")
    parser.add_argument("--kappa", type=float, default=None, help="Kappa (overrides --mass)")
    parser.add_argument("--clov-coeff", type=float, default=1.0, help="Clover coeff (isotropic)")
    parser.add_argument("--clov-coeff-r", type=float, default=None, help="Clover coeff R")
    parser.add_argument("--clov-coeff-t", type=float, default=None, help="Clover coeff T")
    parser.add_argument("--nrow", type=parse_nrow, default="4,4,4,8")
    parser.add_argument("--n-updates", type=int, default=2)
    parser.add_argument("--n-warmup", type=int, default=0)
    parser.add_argument("--tau0", type=float, default=1.0)
    parser.add_argument("--n-steps", type=int, default=20)
    parser.add_argument("--rsd-cg", type=float, default=1.0e-7)
    parser.add_argument("--max-cg", type=int, default=1000)
    parser.add_argument("--seed", type=parse_seed, default="11,11,11,0")
    parser.add_argument("--output", type=Path, default=Path("sf_clover_hmc.out.xml"))
    parser.add_argument("--save-prefix", type=str, default="sf_clover_hmc")
    parser.add_argument(
        "--restart",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Restart from latest saved configuration if present.",
    )
    parser.add_argument("--schr-phi-mult", type=float, default=1.0)
    parser.add_argument("--loop-extent", type=int, default=1)
    parser.add_argument("--decay-dir", type=int, default=3)
    parser.add_argument("--theta", type=parse_theta, default="0,0,0")
    parser.add_argument(
        "--predictor",
        type=str,
        default="zero_guess",
        choices=["zero_guess", "last_solution", "linear", "mre", "mre_initcg"],
    )
    parser.add_argument("--max-chrono", type=int, default=1)
    parser.add_argument("--mre-max-evec", type=int, default=1)
    parser.add_argument("--mre-opt-eigen-id", type=str, default="")
    args = parser.parse_args()

    if args.kappa is None:
        if args.mass is None:
            raise SystemExit("Provide --mass or --kappa.")
        args.kappa = kappa_from_mass(args.mass)

    os.environ["PYCHROMA_OUT"] = str(args.output)

    ferm_bc = SchroedingerFermBC(
        schr_phi_mult=args.schr_phi_mult,
        loop_extent=args.loop_extent,
        decay_dir=args.decay_dir,
        theta=args.theta,
    )
    gauge_bc = SchroedingerGaugeBC(
        schr_phi_mult=args.schr_phi_mult,
        loop_extent=args.loop_extent,
        decay_dir=args.decay_dir,
    )

    ferm_action = CloverAction(
        kappa=args.kappa,
        clov_coeff=args.clov_coeff,
        clov_coeff_r=args.clov_coeff_r,
        clov_coeff_t=args.clov_coeff_t,
        ferm_bc=ferm_bc,
    )

    if args.predictor == "zero_guess":
        predictor = predictor_zero_guess()
    elif args.predictor == "last_solution":
        predictor = predictor_last_solution()
    elif args.predictor == "linear":
        predictor = predictor_linear_extrapolation()
    elif args.predictor == "mre":
        predictor = predictor_mre(max_chrono=args.max_chrono)
    elif args.predictor == "mre_initcg":
        predictor = predictor_mre_initcg(
            max_chrono=args.max_chrono,
            max_evec=args.mre_max_evec,
            opt_eigen_id=args.mre_opt_eigen_id,
        )
    else:
        raise SystemExit(f"Unknown predictor: {args.predictor}")

    monomial_ferm = TwoFlavorEOPrecLogDetFermMonomial(
        monomial_id="clover_2flav",
        invert_param=InvertParam(rsd_cg=args.rsd_cg, max_cg=args.max_cg),
        fermion_action=ferm_action,
        predictor=predictor,
    )

    monomial_gauge = GaugeMonomial(
        monomial_id="gauge",
        gauge_action=GaugeAction(
            name="WILSON_GAUGEACT",
            beta=args.beta,
            gauge_bc=gauge_bc,
        ),
    )

    hamiltonian = Hamiltonian(monomial_ids=["clover_2flav", "gauge"])

    integrator = Integrator(
        name="LCM_STS_LEAPFROG",
        n_steps=args.n_steps,
        monomial_ids=hamiltonian.monomial_ids,
    )
    md_integrator = MDIntegrator(tau0=args.tau0, integrator=integrator)

    trj = HMCTrj(
        nrow=args.nrow,
        monomials=[monomial_ferm, monomial_gauge],
        hamiltonian=hamiltonian,
        md_integrator=md_integrator,
    )

    cfg_type = "CLASSICAL_SF"
    cfg_file = "DUMMY"
    gauge_state = GaugeState(name="SIMPLE_GAUGE_STATE", gauge_bc=gauge_bc)
    start_update_num = 0
    rng_seed = args.seed
    if args.restart:
        latest_cfg, latest_update = find_latest_cfg(args.save_prefix)
        if latest_cfg is not None and latest_update is not None:
            if latest_cfg.suffix == ".lime":
                cfg_type = "SZINQIO"
            else:
                cfg_type = "SCIDAC"
            cfg_file = str(latest_cfg)
            start_update_num = latest_update + 1
            gauge_state = None
            restart_xml = find_restart_xml(args.save_prefix, latest_update)
            if restart_xml is not None:
                try:
                    rng_seed = read_rng_seed(restart_xml)
                except (ET.ParseError, ValueError):
                    pass

    n_production_updates = args.n_updates
    if start_update_num >= n_production_updates:
        n_production_updates = start_update_num + args.n_updates

    mc = MCControl(
        cfg=GaugeConfig(
            cfg_type=cfg_type,
            cfg_file=cfg_file,
            gauge_state=gauge_state,
        ),
        rng_seed=rng_seed,
        start_update_num=start_update_num,
        n_warm_up_updates=args.n_warmup,
        n_production_updates=n_production_updates,
        n_updates_this_run=args.n_updates,
        save_interval=args.n_updates,
        save_prefix=args.save_prefix,
        save_volfmt="SINGLEFILE",
        repro_check=False,
        reverse_check=False,
        monitor_forces=True,
        inline_measurements=[],
    )

    params = HMCParams(mc_control=mc, hmc_trj=trj)

    initialize()
    try:
        run_hmc(params)
    finally:
        finalize()


if __name__ == "__main__":
    main()
