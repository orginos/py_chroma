import argparse
import os
from pathlib import Path

from py_chroma import finalize, initialize
from py_chroma.hmc import (
    FermionBC,
    GaugeAction,
    GaugeBC,
    GaugeMonomial,
    GaugeState,
    Hamiltonian,
    InvertParam,
    Integrator,
    MDIntegrator,
    WilsonAction,
    TwoFlavorEOPrecConstDetFermMonomial,
)
from py_chroma.smd import SMDConfig, SMDMCControl, SMDOptions, SMDParams, SMDTrj


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


def kappa_from_mass(mass):
    return 1.0 / (2.0 * (mass + 4.0))


def main():
    parser = argparse.ArgumentParser(
        description="Run SMD with Nf=2 Wilson fermions via py_chroma (no XML input)."
    )
    parser.add_argument("--beta", type=float, required=True, help="Gauge coupling beta")
    parser.add_argument("--mass", type=float, default=None, help="Bare mass m0")
    parser.add_argument("--kappa", type=float, default=None, help="Kappa (overrides --mass)")
    parser.add_argument("--nrow", type=parse_nrow, default="4,4,4,4")
    parser.add_argument("--n-updates", type=int, default=2)
    parser.add_argument("--tau0", type=float, default=0.5)
    parser.add_argument("--n-steps", type=int, default=25)
    parser.add_argument("--rsd-cg", type=float, default=1.0e-7)
    parser.add_argument("--max-cg", type=int, default=1000)
    parser.add_argument("--seed", type=parse_seed, default="11,0,0,0")
    parser.add_argument("--gamma", type=float, default=0.3)
    parser.add_argument("--accept-reject", action="store_true", default=True)
    parser.add_argument("--no-accept-reject", action="store_false", dest="accept_reject")
    parser.add_argument("--measure-actions", action="store_true", default=True)
    parser.add_argument("--no-measure-actions", action="store_false", dest="measure_actions")
    parser.add_argument("--output", type=Path, default=Path("smd_prec_wilson.out.xml"))
    parser.add_argument("--save-prefix", type=str, default="dummy_run")
    args = parser.parse_args()

    if args.kappa is None:
        if args.mass is None:
            raise SystemExit("Provide --mass or --kappa.")
        args.kappa = kappa_from_mass(args.mass)

    os.environ["PYCHROMA_OUT"] = str(args.output)

    ferm_action = WilsonAction(
        kappa=args.kappa,
        ferm_bc=FermionBC(),
    )
    monomial_ferm = TwoFlavorEOPrecConstDetFermMonomial(
        monomial_id="wilson_two_flav",
        invert_param=InvertParam(rsd_cg=args.rsd_cg, max_cg=args.max_cg),
        fermion_action=ferm_action,
        predictor_name="LAST_SOLUTION_4D_PREDICTOR",
    )

    monomial_gauge = GaugeMonomial(
        monomial_id="gauge",
        gauge_action=GaugeAction(
            name="WILSON_GAUGEACT",
            beta=args.beta,
            gauge_bc=GaugeBC(),
        ),
    )

    hamiltonian = Hamiltonian(monomial_ids=["wilson_two_flav", "gauge"])
    integrator = Integrator(
        name="LCM_STS_LEAPFROG",
        n_steps=args.n_steps,
        monomial_ids=hamiltonian.monomial_ids,
    )
    md_integrator = MDIntegrator(tau0=args.tau0, integrator=integrator)

    trj = SMDTrj(
        nrow=args.nrow,
        smd_params=SMDOptions(
            gamma=args.gamma,
            accept_reject=args.accept_reject,
            measure_actions=args.measure_actions,
        ),
        monomials=[monomial_ferm, monomial_gauge],
        hamiltonian=hamiltonian,
        md_integrator=md_integrator,
    )

    mc = SMDMCControl(
        cfg=SMDConfig(
            cfg_type="WEAK_FIELD",
            cfg_file="DUMMY",
            gauge_state=GaugeState(name="SIMPLE_GAUGE_STATE"),
        ),
        rng_seed=args.seed,
        n_production_updates=args.n_updates,
        n_updates_this_run=args.n_updates,
        save_interval=args.n_updates,
        save_prefix=args.save_prefix,
        save_volfmt="SINGLEFILE",
        repro_check=False,
        reverse_check=True,
        reverse_check_frequency=1,
        monitor_forces=True,
        inline_measurements=[],
    )

    params = SMDParams(mc_control=mc, smd_trj=trj)

    initialize()
    try:
        from py_chroma.smd import run_smd

        run_smd(params)
    finally:
        finalize()


if __name__ == "__main__":
    main()
