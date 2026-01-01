from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Optional, Sequence, Union
import xml.etree.ElementTree as ET

from ._core import run_hmc_xml


def _bool_text(value: bool) -> str:
    return "true" if value else "false"


def _text(value) -> str:
    if isinstance(value, bool):
        return _bool_text(value)
    return str(value)


def _add_text(parent: ET.Element, tag: str, value) -> ET.Element:
    elem = ET.SubElement(parent, tag)
    elem.text = _text(value)
    return elem


def _add_list_elems(parent: ET.Element, tag: str, values: Iterable) -> ET.Element:
    elem = ET.SubElement(parent, tag)
    for value in values:
        _add_text(elem, "elem", value)
    return elem


def _add_space_list(parent: ET.Element, tag: str, values: Iterable) -> ET.Element:
    elem = ET.SubElement(parent, tag)
    elem.text = " ".join(_text(value) for value in values)
    return elem


def _parse_fragment(xml_fragment: str) -> ET.Element:
    return ET.fromstring(xml_fragment)


@dataclass
class GaugeBC:
    name: str = "PERIODIC_GAUGEBC"

    def to_element(self) -> ET.Element:
        elem = ET.Element("GaugeBC")
        _add_text(elem, "Name", self.name)
        return elem


@dataclass
class FermionBC:
    ferm_bc: str = "SIMPLE_FERMBC"
    boundary: Sequence[int] = (1, 1, 1, -1)

    def to_element(self) -> ET.Element:
        elem = ET.Element("FermionBC")
        _add_text(elem, "FermBC", self.ferm_bc)
        _add_space_list(elem, "boundary", self.boundary)
        return elem


@dataclass
class AnisoParam:
    anisoP: bool = False
    t_dir: int = 3
    xi_0: float = 1.0
    nu: Optional[float] = None

    def to_element(self) -> ET.Element:
        elem = ET.Element("AnisoParam")
        _add_text(elem, "anisoP", self.anisoP)
        _add_text(elem, "t_dir", self.t_dir)
        _add_text(elem, "xi_0", self.xi_0)
        if self.nu is not None:
            _add_text(elem, "nu", self.nu)
        return elem


@dataclass
class FermState:
    name: str = "STOUT_FERM_STATE"
    rho: float = 0.14
    n_smear: int = 2
    orthog_dir: int = 3
    ferm_bc: Optional[FermionBC] = None

    def to_element(self) -> ET.Element:
        elem = ET.Element("FermState")
        _add_text(elem, "Name", self.name)
        _add_text(elem, "rho", self.rho)
        _add_text(elem, "n_smear", self.n_smear)
        _add_text(elem, "orthog_dir", self.orthog_dir)
        if self.ferm_bc is not None:
            elem.append(self.ferm_bc.to_element())
        return elem


@dataclass
class SEOPrecCloverAction:
    mass: float
    clov_coeff_r: float
    clov_coeff_t: float
    aniso: Optional[AnisoParam] = None
    ferm_state: Optional[FermState] = None
    ferm_bc: Optional[FermionBC] = None

    def to_element(self) -> ET.Element:
        elem = ET.Element("FermionAction")
        _add_text(elem, "FermAct", "SEOPREC_CLOVER")
        _add_text(elem, "Mass", self.mass)
        _add_text(elem, "clovCoeffR", self.clov_coeff_r)
        _add_text(elem, "clovCoeffT", self.clov_coeff_t)
        if self.aniso is not None:
            elem.append(self.aniso.to_element())
        if self.ferm_state is not None:
            elem.append(self.ferm_state.to_element())
        if self.ferm_bc is not None:
            elem.append(self.ferm_bc.to_element())
        return elem


@dataclass
class CloverParams:
    mass: float
    clov_coeff_r: float
    clov_coeff_t: float
    aniso: Optional[AnisoParam] = None

    def to_element(self) -> ET.Element:
        elem = ET.Element("CloverParams")
        _add_text(elem, "Mass", self.mass)
        _add_text(elem, "clovCoeffR", self.clov_coeff_r)
        _add_text(elem, "clovCoeffT", self.clov_coeff_t)
        if self.aniso is not None:
            elem.append(self.aniso.to_element())
        return elem


@dataclass
class CloverInvertParam:
    inv_type: str = "CG_INVERTER"
    rsd_cg: float = 1.0e-7
    max_cg: int = 1000
    clover_params: Optional[CloverParams] = None

    def to_element(self) -> ET.Element:
        elem = ET.Element("InvertParam")
        _add_text(elem, "invType", self.inv_type)
        if self.clover_params is not None:
            elem.append(self.clover_params.to_element())
        _add_text(elem, "RsdCG", self.rsd_cg)
        _add_text(elem, "MaxCG", self.max_cg)
        return elem


@dataclass
class GaugeAction:
    name: str
    beta: float
    gauge_bc: Optional[GaugeBC] = None
    u_s: Optional[float] = None
    u_t: Optional[float] = None
    zero_energy: Optional[float] = None
    aniso: Optional[AnisoParam] = None

    def to_element(self) -> ET.Element:
        elem = ET.Element("GaugeAction")
        _add_text(elem, "Name", self.name)
        _add_text(elem, "beta", self.beta)
        if self.u_s is not None:
            _add_text(elem, "u_s", self.u_s)
        if self.u_t is not None:
            _add_text(elem, "u_t", self.u_t)
        if self.zero_energy is not None:
            _add_text(elem, "ZeroEnergy", self.zero_energy)
        if self.aniso is not None:
            elem.append(self.aniso.to_element())
        if self.gauge_bc is not None:
            elem.append(self.gauge_bc.to_element())
        return elem


@dataclass
class GaugeMonomial:
    monomial_id: str
    gauge_action: GaugeAction

    def to_element(self) -> ET.Element:
        elem = ET.Element("elem")
        _add_text(elem, "Name", "GAUGE_MONOMIAL")
        elem.append(self.gauge_action.to_element())
        named = ET.SubElement(elem, "NamedObject")
        _add_text(named, "monomial_id", self.monomial_id)
        return elem


@dataclass
class NFlavorLogDetDiagFermMonomial:
    monomial_id: str
    fermion_action: SEOPrecCloverAction
    num_flavors: int = 2

    def to_element(self) -> ET.Element:
        elem = ET.Element("elem")
        _add_text(elem, "Name", "N_FLAVOR_LOGDET_DIAG_FERM_MONOMIAL")
        elem.append(self.fermion_action.to_element())
        _add_text(elem, "num_flavors", self.num_flavors)
        named = ET.SubElement(elem, "NamedObject")
        _add_text(named, "monomial_id", self.monomial_id)
        return elem


@dataclass
class TwoFlavorSEOPrecConstDetRatioConvConvMultiHasenFermMonomial:
    monomial_id: str
    shifted_mass: Sequence[float]
    num_hasen_terms: int
    fermion_action: SEOPrecCloverAction
    invert_param: CloverInvertParam

    def to_element(self) -> ET.Element:
        elem = ET.Element("elem")
        _add_text(
            elem,
            "Name",
            "TWO_FLAVOR_SEOPREC_CONSTDET_RATIO_CONV_CONV_MULTIHASEN_FERM_MONOMIAL",
        )
        _add_space_list(elem, "ShiftedMass", self.shifted_mass)
        _add_text(elem, "NumofHasenTerms", self.num_hasen_terms)
        action = ET.SubElement(elem, "Action")
        action.append(self.fermion_action.to_element())
        action.append(self.invert_param.to_element())
        named = ET.SubElement(elem, "NamedObject")
        _add_text(named, "monomial_id", self.monomial_id)
        return elem


@dataclass
class TwoFlavorSEOPrecConstDetMultiHasenCancelFermMonomial:
    monomial_id: str
    shifted_mass: Union[float, Sequence[float]]
    fermion_action: SEOPrecCloverAction
    invert_param: CloverInvertParam

    def to_element(self) -> ET.Element:
        elem = ET.Element("elem")
        _add_text(
            elem,
            "Name",
            "TWO_FLAVOR_SEOPREC_CONSTDET_MULTIHASEN_CANCEL_FERM_MONOMIAL",
        )
        if isinstance(self.shifted_mass, (list, tuple)):
            _add_space_list(elem, "ShiftedMass", self.shifted_mass)
        else:
            _add_text(elem, "ShiftedMass", self.shifted_mass)
        elem.append(self.fermion_action.to_element())
        elem.append(self.invert_param.to_element())
        named = ET.SubElement(elem, "NamedObject")
        _add_text(named, "monomial_id", self.monomial_id)
        return elem


@dataclass
class InlineMeasurement:
    name: str
    frequency: int
    param_xml: Optional[str] = None
    named_object_xml: Optional[str] = None

    def to_element(self) -> ET.Element:
        elem = ET.Element("elem")
        _add_text(elem, "Name", self.name)
        _add_text(elem, "Frequency", self.frequency)
        if self.param_xml is not None:
            elem.append(_parse_fragment(self.param_xml))
        if self.named_object_xml is not None:
            elem.append(_parse_fragment(self.named_object_xml))
        return elem


@dataclass
class Integrator:
    name: str
    n_steps: int
    monomial_ids: Sequence[str]
    subintegrator: Optional["Integrator"] = None
    lambda_param: Optional[float] = None

    def to_element(self) -> ET.Element:
        elem = ET.Element("Integrator")
        _add_text(elem, "Name", self.name)
        _add_text(elem, "n_steps", self.n_steps)
        _add_list_elems(elem, "monomial_ids", self.monomial_ids)
        if self.lambda_param is not None:
            _add_text(elem, "lambda", self.lambda_param)
        if self.subintegrator is not None:
            sub = ET.SubElement(elem, "SubIntegrator")
            for child in self.subintegrator.to_element():
                sub.append(child)
        return elem


@dataclass
class MDIntegrator:
    tau0: float
    integrator: Integrator
    anisoP: Optional[bool] = None
    t_dir: Optional[int] = None
    xi_mom: Optional[float] = None

    def to_element(self) -> ET.Element:
        elem = ET.Element("MDIntegrator")
        _add_text(elem, "tau0", self.tau0)
        if self.anisoP is not None:
            _add_text(elem, "anisoP", self.anisoP)
        if self.t_dir is not None:
            _add_text(elem, "t_dir", self.t_dir)
        if self.xi_mom is not None:
            _add_text(elem, "xi_mom", self.xi_mom)
        elem.append(self.integrator.to_element())
        return elem


@dataclass
class Hamiltonian:
    monomial_ids: Sequence[str]

    def to_element(self) -> ET.Element:
        elem = ET.Element("Hamiltonian")
        _add_list_elems(elem, "monomial_ids", self.monomial_ids)
        return elem


@dataclass
class HMCTrj:
    nrow: Sequence[int]
    monomials: Sequence[
        Union[
            GaugeMonomial,
            NFlavorLogDetDiagFermMonomial,
            TwoFlavorSEOPrecConstDetRatioConvConvMultiHasenFermMonomial,
            TwoFlavorSEOPrecConstDetMultiHasenCancelFermMonomial,
        ]
    ]
    hamiltonian: Hamiltonian
    md_integrator: MDIntegrator

    def to_element(self) -> ET.Element:
        elem = ET.Element("HMCTrj")
        _add_space_list(elem, "nrow", self.nrow)
        monomials_elem = ET.SubElement(elem, "Monomials")
        for monomial in self.monomials:
            monomials_elem.append(monomial.to_element())
        elem.append(self.hamiltonian.to_element())
        elem.append(self.md_integrator.to_element())
        return elem


@dataclass
class GaugeConfig:
    cfg_type: str
    cfg_file: str = "DUMMY"
    parallel_io: Optional[bool] = None

    def to_element(self) -> ET.Element:
        elem = ET.Element("Cfg")
        _add_text(elem, "cfg_type", self.cfg_type)
        _add_text(elem, "cfg_file", self.cfg_file)
        if self.parallel_io is not None:
            _add_text(elem, "parallel_io", self.parallel_io)
        return elem


@dataclass
class MCControl:
    cfg: GaugeConfig
    rng_seed: Sequence[int] = (0, 0, 0, 0)
    start_update_num: int = 0
    n_warm_up_updates: int = 0
    n_production_updates: int = 1
    n_updates_this_run: int = 1
    save_interval: int = 1
    save_prefix: str = "hmc_run"
    save_volfmt: str = "SINGLEFILE"
    parallel_io: Optional[bool] = None
    repro_check: bool = True
    repro_check_frequency: int = 10
    reverse_check: bool = True
    reverse_check_frequency: int = 10
    monitor_forces: bool = True
    inline_measurements: Optional[Sequence[InlineMeasurement]] = field(default_factory=list)

    def to_element(self) -> ET.Element:
        elem = ET.Element("MCControl")
        elem.append(self.cfg.to_element())
        rng = ET.SubElement(elem, "RNG")
        seed = ET.SubElement(rng, "Seed")
        for value in self.rng_seed:
            _add_text(seed, "elem", value)
        _add_text(elem, "StartUpdateNum", self.start_update_num)
        _add_text(elem, "NWarmUpUpdates", self.n_warm_up_updates)
        _add_text(elem, "NProductionUpdates", self.n_production_updates)
        _add_text(elem, "NUpdatesThisRun", self.n_updates_this_run)
        _add_text(elem, "SaveInterval", self.save_interval)
        _add_text(elem, "SavePrefix", self.save_prefix)
        _add_text(elem, "SaveVolfmt", self.save_volfmt)
        if self.parallel_io is not None:
            _add_text(elem, "ParallelIO", self.parallel_io)
        _add_text(elem, "ReproCheckP", self.repro_check)
        _add_text(elem, "ReproCheckFrequency", self.repro_check_frequency)
        _add_text(elem, "ReverseCheckP", self.reverse_check)
        _add_text(elem, "ReverseCheckFrequency", self.reverse_check_frequency)
        _add_text(elem, "MonitorForces", self.monitor_forces)
        if self.inline_measurements is not None:
            inline_elem = ET.SubElement(elem, "InlineMeasurements")
            for meas in self.inline_measurements:
                inline_elem.append(meas.to_element())
        return elem


@dataclass
class HMCParams:
    mc_control: MCControl
    hmc_trj: HMCTrj

    def to_element(self) -> ET.Element:
        root = ET.Element("Params")
        root.append(self.mc_control.to_element())
        root.append(self.hmc_trj.to_element())
        return root

    def to_xml(self) -> str:
        root = self.to_element()
        return ET.tostring(root, encoding="unicode")


def run_hmc(params: HMCParams) -> None:
    if not isinstance(params, HMCParams):
        raise TypeError("params must be an HMCParams instance")
    run_hmc_xml(params.to_xml())
