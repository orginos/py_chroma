from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence, Union
import xml.etree.ElementTree as ET

from ._core import run_smd_xml
from . import hmc as _hmc


@dataclass
class SMDOptions:
    gamma: float = 0.3
    accept_reject: bool = True
    measure_actions: bool = True

    def to_element(self) -> ET.Element:
        elem = ET.Element("SMDParams")
        _hmc._add_text(elem, "gamma", self.gamma)
        _hmc._add_text(elem, "AcceptReject", self.accept_reject)
        _hmc._add_text(elem, "MeasureActions", self.measure_actions)
        return elem


@dataclass
class SMDTrj:
    nrow: Sequence[int]
    smd_params: SMDOptions = field(default_factory=SMDOptions)
    monomials: Sequence[Union[
        _hmc.GaugeMonomial,
        _hmc.NFlavorLogDetDiagFermMonomial,
        _hmc.TwoFlavorSEOPrecConstDetRatioConvConvMultiHasenFermMonomial,
        _hmc.TwoFlavorSEOPrecConstDetMultiHasenCancelFermMonomial,
        _hmc.TwoFlavorEOPrecLogDetFermMonomial,
        _hmc.TwoFlavorEOPrecConstDetFermMonomial,
    ]] = field(default_factory=list)
    hamiltonian: _hmc.Hamiltonian = field(default_factory=lambda: _hmc.Hamiltonian(monomial_ids=[]))
    md_integrator: _hmc.MDIntegrator = field(default_factory=lambda: _hmc.MDIntegrator(0.0, _hmc.Integrator(name="LCM_STS_LEAPFROG", n_steps=1, monomial_ids=[])))

    def to_element(self) -> ET.Element:
        elem = ET.Element("SMDTrj")
        elem.append(self.smd_params.to_element())
        monomials_elem = ET.SubElement(elem, "Monomials")
        for monomial in self.monomials:
            monomials_elem.append(monomial.to_element())
        elem.append(self.hamiltonian.to_element())
        elem.append(self.md_integrator.to_element())
        _hmc._add_space_list(elem, "nrow", self.nrow)
        return elem


@dataclass
class SMDConfig:
    cfg_type: str
    cfg_file: str = "DUMMY"
    parallel_io: Optional[bool] = None
    reunit: Optional[bool] = None
    gauge_state: Optional[_hmc.GaugeState] = None

    def to_element(self, tag: str) -> ET.Element:
        elem = ET.Element(tag)
        _hmc._add_text(elem, "cfg_type", self.cfg_type)
        _hmc._add_text(elem, "cfg_file", self.cfg_file)
        if self.parallel_io is not None:
            _hmc._add_text(elem, "parallel_io", self.parallel_io)
        if self.reunit is not None:
            _hmc._add_text(elem, "reunit", self.reunit)
        if self.gauge_state is not None:
            elem.append(self.gauge_state.to_element())
        return elem


@dataclass
class SMDMCControl:
    cfg: SMDConfig
    rng_seed: Sequence[int] = (0, 0, 0, 0)
    start_update_num: int = 0
    n_warm_up_updates: int = 0
    n_production_updates: int = 1
    n_updates_this_run: int = 1
    save_interval: int = 1
    save_prefix: str = "smd_run"
    save_volfmt: str = "SINGLEFILE"
    parallel_io: Optional[bool] = None
    repro_check: bool = True
    repro_check_frequency: int = 10
    reverse_check: bool = True
    reverse_check_frequency: int = 10
    monitor_forces: bool = True
    inline_measurements: Optional[Sequence[_hmc.InlineMeasurement]] = field(default_factory=list)
    momenta: Optional[SMDConfig] = None

    def to_element(self) -> ET.Element:
        elem = ET.Element("MCControl")
        elem.append(self.cfg.to_element("Cfg"))
        if self.momenta is not None:
            elem.append(self.momenta.to_element("Momenta"))
        rng = ET.SubElement(elem, "RNG")
        seed = ET.SubElement(rng, "Seed")
        for value in self.rng_seed:
            _hmc._add_text(seed, "elem", value)
        _hmc._add_text(elem, "StartUpdateNum", self.start_update_num)
        _hmc._add_text(elem, "NWarmUpUpdates", self.n_warm_up_updates)
        _hmc._add_text(elem, "NProductionUpdates", self.n_production_updates)
        _hmc._add_text(elem, "NUpdatesThisRun", self.n_updates_this_run)
        _hmc._add_text(elem, "SaveInterval", self.save_interval)
        _hmc._add_text(elem, "SavePrefix", self.save_prefix)
        _hmc._add_text(elem, "SaveVolfmt", self.save_volfmt)
        if self.parallel_io is not None:
            _hmc._add_text(elem, "ParallelIO", self.parallel_io)
        _hmc._add_text(elem, "ReproCheckP", self.repro_check)
        _hmc._add_text(elem, "ReproCheckFrequency", self.repro_check_frequency)
        _hmc._add_text(elem, "ReverseCheckP", self.reverse_check)
        _hmc._add_text(elem, "ReverseCheckFrequency", self.reverse_check_frequency)
        _hmc._add_text(elem, "MonitorForces", self.monitor_forces)
        if self.inline_measurements is not None:
            inline_elem = ET.SubElement(elem, "InlineMeasurements")
            for meas in self.inline_measurements:
                inline_elem.append(meas.to_element())
        return elem


@dataclass
class SMDParams:
    mc_control: SMDMCControl
    smd_trj: SMDTrj

    def to_element(self) -> ET.Element:
        root = ET.Element("Params")
        root.append(self.mc_control.to_element())
        root.append(self.smd_trj.to_element())
        return root

    def to_xml(self) -> str:
        root = self.to_element()
        return ET.tostring(root, encoding="unicode")


def run_smd(params: SMDParams) -> None:
    if not isinstance(params, SMDParams):
        raise TypeError("params must be an SMDParams instance")
    run_smd_xml(params.to_xml())
