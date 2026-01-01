import argparse
import os
from pathlib import Path

from py_chroma import finalize, initialize
from py_chroma.hmc import (
    CloverAction,
    GaugeAction,
    GaugeMonomial,
    GaugeState,
    Hamiltonian,
    HMCParams,
    HMCTrj,
    InvertParam,
    MDIntegrator,
    Integrator,
    MCControl,
    SchroedingerFermBC,
    SchroedingerGaugeBC,
    TwoFlavorEOPrecLogDetFermMonomial,
    GaugeConfig,
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


def kappa_from_mass(mass):
    return 1.0 / (2.0 * (mass + 4.0))


def main():
    parser = argparse.ArgumentParser(
        description="Run Nf=2 Clover HMC with Schroedinger BC via py_chroma (no XML input)."
    )
    parser.add_argument("--beta", type=float, required=True, help="Gauge coupling beta")
    parser.add_argument("--mass", type=float, default=None, help="Bare mass m0")
    parser.add_argument("--kappa", type=float, default=None, help="Kappa (overrides --mass)")
    parser.add_argument("--clov-coeff-r", type=float, default=1.0, help="Clover coeff R")
    parser.add_argument("--clov-coeff-t", type=float, default=1.0, help="Clover coeff T")
    parser.add_argument("--nrow", type=parse_nrow, default="4,4,4,8")
    parser.add_argument("--n-updates", type=int, default=2)
    parser.add_argument("--tau0", type=float, default=0.1)
    parser.add_argument("--n-steps", type=int, default=20)
    parser.add_argument("--rsd-cg", type=float, default=1.0e-7)
    parser.add_argument("--max-cg", type=int, default=1000)
    parser.add_argument("--seed", type=parse_seed, default="11,11,11,0")
    parser.add_argument("--output", type=Path, default=Path("sf_clover_hmc.out.xml"))
    parser.add_argument("--save-prefix", type=str, default="sf_clover_hmc")
    parser.add_argument("--schr-phi-mult", type=float, default=1.0)
    parser.add_argument("--loop-extent", type=int, default=1)
    parser.add_argument("--decay-dir", type=int, default=3)
    parser.add_argument("--theta", type=parse_theta, default="0,0,0")
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
        clov_coeff_r=args.clov_coeff_r,
        clov_coeff_t=args.clov_coeff_t,
        ferm_bc=ferm_bc,
        chrono_predictor_name="ZERO_SOLUTION_4D_PREDICTOR",
    )

    monomial_ferm = TwoFlavorEOPrecLogDetFermMonomial(
        monomial_id="clover_2flav",
        invert_param=InvertParam(rsd_cg=args.rsd_cg, max_cg=args.max_cg),
        fermion_action=ferm_action,
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

    mc = MCControl(
        cfg=GaugeConfig(
            cfg_type="CLASSICAL_SF",
            cfg_file="DUMMY",
            gauge_state=GaugeState(name="SIMPLE_GAUGE_STATE", gauge_bc=gauge_bc),
        ),
        rng_seed=args.seed,
        n_production_updates=args.n_updates,
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
        from py_chroma.hmc import run_hmc

        run_hmc(params)
    finally:
        finalize()


if __name__ == "__main__":
    main()
