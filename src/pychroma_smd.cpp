/*! \file
 *  \brief py_chroma SMD runner (XML-driven)
 */

#include "pychroma_api.h"
#include "pychroma_state.h"

#include "chroma.h"
#include "io/monomial_io.h"
#include "meas/inline/inline.h"
#include "meas/inline/io/default_gauge_field.h"
#include "update/molecdyn/hamiltonian/exact_hamiltonian.h"
#include "update/molecdyn/integrator/lcm_toplevel_integrator.h"
#include "update/molecdyn/monomial/force_monitors.h"
#include "update/molecdyn/smd/smd.h"
#include "util/gauge/gauge_init_aggregate.h"
#include "util/gauge/szinqio_gauge_init.h"

#include <sstream>
#include <string>

using namespace Chroma;

namespace {

bool cfgReunitFromXML(const std::string& xml, bool default_value) {
  if (xml.find("<reunit>") == std::string::npos) {
    return default_value;
  }
  if (xml.find("<reunit>false</reunit>") != std::string::npos ||
      xml.find("<reunit>0</reunit>") != std::string::npos) {
    return false;
  }
  if (xml.find("<reunit>true</reunit>") != std::string::npos ||
      xml.find("<reunit>1</reunit>") != std::string::npos) {
    return true;
  }
  return default_value;
}

struct SMDMCControl {
  GroupXML_t cfg;
  GroupXML_t mom;
  bool mom_present;
  QDP::Seed rng_seed;
  unsigned long start_update_num;
  unsigned long n_warm_up_updates;
  unsigned long n_production_updates;
  unsigned int n_updates_this_run;
  unsigned int save_interval;
  std::string save_prefix;
  QDP_volfmt_t save_volfmt;
  QDP_serialparallel_t save_pario;
  std::string inline_measurement_xml;
  bool repro_checkP;
  int repro_check_frequency;
  bool rev_checkP;
  int rev_check_frequency;
  bool monitorForcesP;
};

void readSMDMCControl(XMLReader& xml, const std::string& path, SMDMCControl& p) {
  START_CODE();

  try {
    XMLReader paramtop(xml, path);
    p.cfg = readXMLGroup(paramtop, "Cfg", "cfg_type");
    p.mom_present = false;
    if (paramtop.count("./Momenta") > 0) {
      p.mom = readXMLGroup(paramtop, "Momenta", "cfg_type");
      p.mom_present = true;
    }
    read(paramtop, "./RNG", p.rng_seed);
    read(paramtop, "./StartUpdateNum", p.start_update_num);
    read(paramtop, "./NWarmUpUpdates", p.n_warm_up_updates);
    read(paramtop, "./NProductionUpdates", p.n_production_updates);
    read(paramtop, "./NUpdatesThisRun", p.n_updates_this_run);
    read(paramtop, "./SaveInterval", p.save_interval);
    read(paramtop, "./SavePrefix", p.save_prefix);
    read(paramtop, "./SaveVolfmt", p.save_volfmt);

    bool parioP = Layout::isIOGridDefined() && (Layout::numIONodeGrid() > 1);
    if (paramtop.count("./parallel_io") > 0) {
      read(paramtop, "./parallel_io", parioP);
    } else if (paramtop.count("./ParallelIO") > 0) {
      read(paramtop, "./ParallelIO", parioP);
    }
    p.save_pario = parioP ? QDPIO_PARALLEL : QDPIO_SERIAL;

    p.repro_checkP = true;
    p.repro_check_frequency = 10;
    if (paramtop.count("./ReproCheckP") == 1) {
      read(paramtop, "./ReproCheckP", p.repro_checkP);
    }
    if (p.repro_checkP && paramtop.count("./ReproCheckFrequency") == 1) {
      read(paramtop, "./ReproCheckFrequency", p.repro_check_frequency);
    }

    p.rev_checkP = true;
    p.rev_check_frequency = 10;
    if (paramtop.count("./ReverseCheckP") == 1) {
      read(paramtop, "./ReverseCheckP", p.rev_checkP);
    }
    if (p.rev_checkP && paramtop.count("./ReverseCheckFrequency") == 1) {
      read(paramtop, "./ReverseCheckFrequency", p.rev_check_frequency);
    }

    if (paramtop.count("./MonitorForces") == 1) {
      read(paramtop, "./MonitorForces", p.monitorForcesP);
    } else {
      p.monitorForcesP = true;
    }

    if (paramtop.count("./InlineMeasurements") == 0) {
      XMLBufferWriter dummy;
      push(dummy, "InlineMeasurements");
      pop(dummy);
      p.inline_measurement_xml = dummy.printCurrentContext();
    } else {
      XMLReader measurements_xml(paramtop, "./InlineMeasurements");
      std::ostringstream inline_os;
      measurements_xml.print(inline_os);
      p.inline_measurement_xml = inline_os.str();
    }
  } catch (const std::string& e) {
    QDPIO::cerr << "Caught Exception: " << e << std::endl;
    QDP_abort(1);
  }

  END_CODE();
}

void writeSMDMCControl(XMLWriter& xml, const std::string& path, const SMDMCControl& p) {
  START_CODE();

  try {
    push(xml, path);
    xml << p.cfg.xml;
    if (p.mom_present) {
      xml << p.mom.xml;
    }
    write(xml, "RNG", p.rng_seed);
    write(xml, "StartUpdateNum", p.start_update_num);
    write(xml, "NWarmUpUpdates", p.n_warm_up_updates);
    write(xml, "NProductionUpdates", p.n_production_updates);
    write(xml, "NUpdatesThisRun", p.n_updates_this_run);
    write(xml, "SaveInterval", p.save_interval);
    write(xml, "SavePrefix", p.save_prefix);
    write(xml, "SaveVolfmt", p.save_volfmt);
    {
      bool pario = (p.save_pario == QDPIO_PARALLEL);
      write(xml, "ParallelIO", pario);
    }
    write(xml, "ReproCheckP", p.repro_checkP);
    if (p.repro_checkP) {
      write(xml, "ReproCheckFrequency", p.repro_check_frequency);
    }
    write(xml, "ReverseCheckP", p.rev_checkP);
    if (p.rev_checkP) {
      write(xml, "ReverseCheckFrequency", p.rev_check_frequency);
    }
    write(xml, "MonitorForces", p.monitorForcesP);

    xml << p.inline_measurement_xml;

    pop(xml);
  } catch (const std::string& e) {
    QDPIO::cerr << "Caught Exception: " << e << std::endl;
    QDP_abort(1);
  }

  END_CODE();
}

struct SMDTrjParams {
  multi1d<int> nrow;
  Real gamma;
  bool accept_reject;
  bool measure_actions;
  std::string Monomials_xml;
  std::string H_MC_xml;
  std::string Integrator_xml;
};

void writeSMDTrj(XMLWriter& xml, const std::string& path, const SMDTrjParams& p) {
  START_CODE();

  try {
    push(xml, path);
    write(xml, "nrow", p.nrow);
    push(xml, "SMDParams");
    write(xml, "gamma", p.gamma);
    write(xml, "AcceptReject", p.accept_reject);
    write(xml, "MeasureActions", p.measure_actions);
    pop(xml);

    xml << p.Monomials_xml;
    xml << p.H_MC_xml;
    xml << p.Integrator_xml;
    pop(xml);
  } catch (const std::string& e) {
    QDPIO::cerr << "Caught Exception: " << e << std::endl;
    QDP_abort(1);
  }

  END_CODE();
}

void readSMDTrj(XMLReader& xml, const std::string& path, SMDTrjParams& p) {
  START_CODE();

  try {
    XMLReader paramtop(xml, path);
    read(paramtop, "./nrow", p.nrow);

    p.gamma = Real(0.3);
    p.accept_reject = true;
    p.measure_actions = true;
    if (paramtop.count("./SMDParams") > 0) {
      read(paramtop, "./SMDParams/gamma", p.gamma);
      if (paramtop.count("./SMDParams/AcceptReject") > 0) {
        read(paramtop, "./SMDParams/AcceptReject", p.accept_reject);
      }
      if (paramtop.count("./SMDParams/MeasureActions") > 0) {
        read(paramtop, "./SMDParams/MeasureActions", p.measure_actions);
      }
    }

    XMLReader monomials_xml_reader(paramtop, "./Monomials");
    std::ostringstream os_monomials;
    monomials_xml_reader.print(os_monomials);
    p.Monomials_xml = os_monomials.str();

    XMLReader h_mc_xml(paramtop, "./Hamiltonian");
    std::ostringstream os_h_mc;
    h_mc_xml.print(os_h_mc);
    p.H_MC_xml = os_h_mc.str();

    XMLReader md_integrator_xml(paramtop, "./MDIntegrator");
    std::ostringstream os_integrator;
    md_integrator_xml.print(os_integrator);
    p.Integrator_xml = os_integrator.str();
  } catch (const std::string& e) {
    QDPIO::cerr << "Error reading XML : " << e << std::endl;
    QDP_abort(1);
  }

  END_CODE();
}

template <typename UpdateParams>
void saveState(const UpdateParams& update_params,
               SMDMCControl& mc_control,
               unsigned long update_no,
               const multi1d<LatticeColorMatrix>& p,
               const multi1d<LatticeColorMatrix>& u) {
  // Do nothing.
}

static bool checkReproducabilitySMD(const multi1d<LatticeColorMatrix>& P_new,
                                    const multi1d<LatticeColorMatrix>& Q_new,
                                    const QDP::Seed& seed_new,
                                    const multi1d<LatticeColorMatrix>& P_old,
                                    const multi1d<LatticeColorMatrix>& Q_old,
                                    const QDP::Seed& seed_old);

template <>
void saveState(const SMDTrjParams& update_params,
               SMDMCControl& mc_control,
               unsigned long update_no,
               const multi1d<LatticeColorMatrix>& p,
               const multi1d<LatticeColorMatrix>& u) {
  START_CODE();

  std::ostringstream restart_data_filename;
  restart_data_filename << mc_control.save_prefix << "_restart_" << update_no
                        << ".xml";

  std::ostringstream restart_config_filename;
  restart_config_filename << mc_control.save_prefix << "_cfg_" << update_no
                          << ".lime";

  std::ostringstream restart_momenta_filename;
  restart_momenta_filename << mc_control.save_prefix << "_mom_" << update_no
                           << ".lime";

  XMLBufferWriter restart_data_buffer;
  SMDMCControl p_new = mc_control;

  QDP::RNG::savern(p_new.rng_seed);
  p_new.start_update_num = update_no;

  unsigned long total =
      mc_control.n_warm_up_updates + mc_control.n_production_updates;
  if (total < mc_control.n_updates_this_run + update_no) {
    p_new.n_updates_this_run = total - update_no;
  }

  {
    SZINQIOGaugeInitEnv::Params cfg;
    cfg.reunitP = cfgReunitFromXML(mc_control.cfg.xml, false);
    cfg.cfg_file = restart_config_filename.str();
    cfg.cfg_pario = mc_control.save_pario;
    p_new.cfg = SZINQIOGaugeInitEnv::createXMLGroup(cfg);
  }

  {
    SZINQIOGaugeInitEnv::Params mom;
    mom.cfg_file = restart_momenta_filename.str();
    mom.cfg_pario = mc_control.save_pario;
    mom.reunitP = false;

    XMLBufferWriter xml_tmp;
    write(xml_tmp, "Momenta", mom);
    p_new.mom.xml = xml_tmp.printCurrentContext();
    p_new.mom.id = SZINQIOGaugeInitEnv::name;
    p_new.mom.path = "/Momenta";
    p_new.mom_present = true;
  }

  push(restart_data_buffer, "Params");
  writeSMDMCControl(restart_data_buffer, "MCControl", p_new);
  writeSMDTrj(restart_data_buffer, "SMDTrj", update_params);
  pop(restart_data_buffer);

  XMLBufferWriter file_xml;
  push(file_xml, "SMD");
  proginfo(file_xml);
  pop(file_xml);

  writeGauge(file_xml, restart_data_buffer, u,
             restart_config_filename.str(), p_new.save_volfmt,
             p_new.save_pario);
  writeGauge(file_xml, restart_data_buffer, p,
             restart_momenta_filename.str(), p_new.save_volfmt,
             p_new.save_pario);

  XMLFileWriter restart_xml(restart_data_filename.str().c_str());
  restart_xml << restart_data_buffer;
  restart_xml.close();

  END_CODE();
}

template <typename UpdateParams>
void doSMD(multi1d<LatticeColorMatrix>& u,
           LatColMatSMDTrj& theSMDTrj,
           SMDMCControl& mc_control,
           const UpdateParams& update_params,
           multi1d<Handle<AbsInlineMeasurement>>& user_measurements) {
  START_CODE();

  QDPIO::cout << "Setting Force monitoring to " << mc_control.monitorForcesP
              << std::endl;
  setForceMonitoring(mc_control.monitorForcesP);
  QDP::StopWatch swatch;

  XMLWriter& xml_out = TheXMLOutputWriter::Instance();
  XMLWriter& xml_log = TheXMLLogWriter::Instance();

  push(xml_out, "doSMD");
  push(xml_log, "doSMD");

  multi1d<Handle<AbsInlineMeasurement>> default_measurements(1);
  InlinePlaquetteEnv::Params plaq_params;
  plaq_params.frequency = 1;
  default_measurements[0] = new InlinePlaquetteEnv::InlineMeas(plaq_params);

  {
    QDP::RNG::setrn(mc_control.rng_seed);

    multi1d<LatticeColorMatrix> p(Nd);
    bool momenta_loaded = false;
    if (mc_control.mom_present) {
      try {
        XMLReader file_xml_mom;
        XMLReader mom_xml;

        std::string mom_id = mc_control.mom.id;
        GroupXML_t mom_cfg = mc_control.mom;

        if (mom_id == SZINQIOGaugeInitEnv::name || mom_id == "SCIDAC") {
          std::istringstream xml_m(mc_control.mom.xml);
          XMLReader momtop(xml_m);
          SZINQIOGaugeInitEnv::Params mom_params(momtop, mc_control.mom.path);
          mom_params.reunitP = false;

          XMLBufferWriter xml_tmp;
          Chroma::write(xml_tmp, "Momenta", mom_params);
          mom_cfg.xml = xml_tmp.printCurrentContext();
          mom_cfg.id = SZINQIOGaugeInitEnv::name;
          mom_cfg.path = "/Momenta";
        }

        std::istringstream xml_m(mom_cfg.xml);
        XMLReader momtop(xml_m);
        Handle<GaugeInit> momInit(
            TheGaugeInitFactory::Instance().createObject(mom_cfg.id, momtop,
                                                         mom_cfg.path));
        (*momInit)(file_xml_mom, mom_xml, p);
        momenta_loaded = true;
      } catch (...) {
        QDPIO::cerr
            << "SMD: Failed to load momenta file, falling back to Gaussian init"
            << std::endl;
        momenta_loaded = false;
      }
    }

    if (!momenta_loaded) {
      for (int mu = 0; mu < Nd; mu++) {
        gaussian(p[mu]);
        p[mu] *= sqrt(Real(0.5));
        taproj(p[mu]);
      }
    }

    GaugeFieldState gauge_state(p, u);
    unsigned long cur_update = mc_control.start_update_num;
    unsigned long total_updates =
        mc_control.n_warm_up_updates + mc_control.n_production_updates;

    unsigned long to_do = 0;
    if (total_updates >= mc_control.n_updates_this_run + cur_update + 1) {
      to_do = mc_control.n_updates_this_run;
    } else {
      to_do = total_updates - cur_update;
    }

    push(xml_out, "MCUpdates");
    push(xml_log, "MCUpdates");

    for (int i = 0; i < to_do; i++) {
      push(xml_out, "elem");
      push(xml_log, "elem");

      push(xml_out, "Update");
      push(xml_log, "Update");

      cur_update++;
      bool warm_up_p = cur_update <= mc_control.n_warm_up_updates;

      write(xml_out, "update_no", cur_update);
      write(xml_log, "update_no", cur_update);
      write(xml_out, "WarmUpP", warm_up_p);
      write(xml_log, "WarmUpP", warm_up_p);

      bool do_reverse = false;
      if (mc_control.rev_checkP &&
          (cur_update % mc_control.rev_check_frequency == 0)) {
        do_reverse = true;
      }

      if (mc_control.repro_checkP &&
          (cur_update % mc_control.repro_check_frequency == 0)) {
        GaugeFieldState repro_bkup_start(gauge_state.getP(),
                                         gauge_state.getQ());
        QDP::Seed rng_seed_bkup_start;
        QDP::RNG::savern(rng_seed_bkup_start);

        swatch.reset();
        swatch.start();
        theSMDTrj(gauge_state, warm_up_p, do_reverse);
        swatch.stop();

        write(xml_out, "seconds_for_trajectory", swatch.getTimeInSeconds());
        write(xml_log, "seconds_for_trajectory", swatch.getTimeInSeconds());

        GaugeFieldState repro_bkup_end(gauge_state.getP(), gauge_state.getQ());
        QDP::Seed rng_seed_bkup_end;
        QDP::RNG::savern(rng_seed_bkup_end);

        gauge_state.getP() = repro_bkup_start.getP();
        gauge_state.getQ() = repro_bkup_start.getQ();
        QDP::RNG::setrn(rng_seed_bkup_start);

        swatch.reset();
        swatch.start();
        theSMDTrj(gauge_state, warm_up_p, false);
        swatch.stop();

        write(xml_out, "seconds_for_repro_trajectory",
              swatch.getTimeInSeconds());
        write(xml_log, "seconds_for_repro_trajectory",
              swatch.getTimeInSeconds());

        QDP::Seed rng_seed_end2;
        QDP::RNG::savern(rng_seed_end2);

        bool pass = checkReproducabilitySMD(gauge_state.getP(),
                                            gauge_state.getQ(),
                                            rng_seed_end2,
                                            repro_bkup_end.getP(),
                                            repro_bkup_end.getQ(),
                                            rng_seed_bkup_end);

        write(xml_out, "ReproCheck", pass);
        write(xml_log, "ReproCheck", pass);

        if (!pass) {
          QDP_abort(1);
        }
      } else {
        swatch.reset();
        swatch.start();
        theSMDTrj(gauge_state, warm_up_p, do_reverse);
        swatch.stop();

        write(xml_out, "seconds_for_trajectory", swatch.getTimeInSeconds());
        write(xml_log, "seconds_for_trajectory", swatch.getTimeInSeconds());
      }

      swatch.reset();
      swatch.start();

      {
        XMLBufferWriter gauge_xml;
        push(gauge_xml, "ChromaSMD");
        write(gauge_xml, "update_no", cur_update);
        writeSMDTrj(gauge_xml, "SMDTrj", update_params);
        pop(gauge_xml);

        InlineDefaultGaugeField::reset();
        InlineDefaultGaugeField::set(gauge_state.getQ(), gauge_xml);

        push(xml_out, "InlineObservables");

        for (int m = 0; m < default_measurements.size(); m++) {
          AbsInlineMeasurement& the_meas = *(default_measurements[m]);
          push(xml_out, "elem");
          the_meas(cur_update, xml_out);
          pop(xml_out);
        }

        for (int m = 0; m < user_measurements.size(); m++) {
          AbsInlineMeasurement& the_meas = *(user_measurements[m]);
          if (cur_update % the_meas.getFrequency() == 0) {
            push(xml_out, "elem");
            the_meas(cur_update, xml_out);
            pop(xml_out);
          }
        }

        pop(xml_out);
        InlineDefaultGaugeField::reset();
      }

      swatch.stop();
      write(xml_out, "seconds_for_measurements", swatch.getTimeInSeconds());
      write(xml_log, "seconds_for_measurements", swatch.getTimeInSeconds());

      if (cur_update % mc_control.save_interval == 0) {
        swatch.reset();
        swatch.start();
        saveState<UpdateParams>(update_params, mc_control, cur_update,
                                gauge_state.getP(), gauge_state.getQ());
        swatch.stop();
      }

      pop(xml_log);
      pop(xml_out);

      pop(xml_log);
      pop(xml_out);
    }

    saveState<UpdateParams>(update_params, mc_control, cur_update,
                            gauge_state.getP(), gauge_state.getQ());

    pop(xml_log);
    pop(xml_out);
  }

  pop(xml_log);
  pop(xml_out);

  END_CODE();
}

static bool linkageHackSMD(void) {
  bool ok = true;
  ok &= GaugeMonomialEnv::registerAll();
  ok &= WilsonTypeFermMonomialAggregrateEnv::registerAll();
  ok &= LCMMDComponentIntegratorAggregateEnv::registerAll();
  ok &= ChronoPredictorAggregrateEnv::registerAll();
  ok &= InlineAggregateEnv::registerAll();
  ok &= GaugeInitEnv::registerAll();
  return ok;
}

static bool checkReproducabilitySMD(const multi1d<LatticeColorMatrix>& P_new,
                                    const multi1d<LatticeColorMatrix>& Q_new,
                                    const QDP::Seed& seed_new,
                                    const multi1d<LatticeColorMatrix>& P_old,
                                    const multi1d<LatticeColorMatrix>& Q_old,
                                    const QDP::Seed& seed_old) {
#if !defined(QDP_IS_QDPJIT2)
  int diffs_found = 0;
  if (P_new.size() != P_old.size()) {
    return false;
  }
  if (Q_new.size() != Q_old.size()) {
    return false;
  }

  int bytes = 2 * Nc * Nc * Layout::sitesOnNode() *
              sizeof(WordType<LatticeColorMatrix>::Type_t);

  for (int mu = 0; mu < P_new.size(); mu++) {
    const unsigned char* p1 =
        reinterpret_cast<const unsigned char*>(P_new[mu].getF());
    const unsigned char* p2 =
        reinterpret_cast<const unsigned char*>(P_old[mu].getF());
    for (int b = 0; b < bytes; b++) {
      unsigned char diff = *p1 - *p2;
      if (diff != 0) {
        diffs_found++;
      }
      p1++;
      p2++;
    }
  }

  QDPInternal::globalSum(diffs_found);
  if (diffs_found != 0) {
    QDPIO::cout << "Found " << diffs_found
                << " different bytes in momentum repro check" << std::endl;
    return false;
  }

  diffs_found = 0;
  for (int mu = 0; mu < P_new.size(); mu++) {
    const unsigned char* p1 =
        reinterpret_cast<const unsigned char*>(Q_new[mu].getF());
    const unsigned char* p2 =
        reinterpret_cast<const unsigned char*>(Q_old[mu].getF());
    for (int b = 0; b < bytes; b++) {
      unsigned char diff = *p1 - *p2;
      if (diff != 0) {
        diffs_found++;
      }
      p1++;
      p2++;
    }
  }

  QDPInternal::globalSum(diffs_found);
  if (diffs_found != 0) {
    QDPIO::cout << "Found " << diffs_found
                << " different bytes in gauge repro check" << std::endl;
    return false;
  }

  if (!toBool(seed_new == seed_old)) {
    QDPIO::cout << "New and old RNG seeds do not match " << std::endl;
    return false;
  }
#else
  QDPIO::cout << "qdp-jit2: skipping momentum repro check" << std::endl;
#endif

  return true;
}

}  // namespace

extern "C" {

int pychroma_run_smd_xml(const char* params_xml) {
  if (guard_initialized() != 0) {
    return 1;
  }
  if (params_xml == nullptr) {
    set_error("SMD params XML is null.");
    return 1;
  }

  static bool registered = false;
  if (!registered) {
    linkageHackSMD();
    registered = true;
  }

  try {
    std::istringstream input_stream(params_xml);
    XMLReader xml_in(input_stream);

    SMDTrjParams trj_params;
    SMDMCControl mc_control;

    XMLReader paramtop(xml_in, "/Params");
    readSMDTrj(paramtop, "./SMDTrj", trj_params);
    readSMDMCControl(paramtop, "./MCControl", mc_control);

    if (mc_control.start_update_num >= mc_control.n_production_updates) {
      return 0;
    }

    Layout::setLattSize(trj_params.nrow);
    Layout::create();
    g_layout_ready = true;

    XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
    XMLFileWriter& xml_log = Chroma::getXMLLogInstance();
    push(xml_out, "smd");
    push(xml_log, "smd");
    write(xml_out, "Input", xml_in);
    write(xml_log, "Input", xml_in);

    g_u.resize(Nd);
    XMLReader file_xml;
    XMLReader config_xml;
    {
      std::istringstream xml_c(mc_control.cfg.xml);
      XMLReader cfgtop(xml_c);
      Handle<GaugeInit> gaugeInit(
          TheGaugeInitFactory::Instance().createObject(mc_control.cfg.id,
                                                       cfgtop,
                                                       mc_control.cfg.path));
      (*gaugeInit)(file_xml, config_xml, g_u);
    }

    bool do_reunit = true;
    if (mc_control.cfg.id == SZINQIOGaugeInitEnv::name ||
        mc_control.cfg.id == "SCIDAC") {
      do_reunit = cfgReunitFromXML(mc_control.cfg.xml, true);
    }

    if (do_reunit) {
      int numbad = 0;
      for (int mu = 0; mu < Nd; mu++) {
        int numbad_mu = 0;
        reunit(g_u[mu], numbad_mu, REUNITARIZE_LABEL);
        numbad += numbad_mu;
      }
    }

    XMLBufferWriter config_buf;
    config_buf << config_xml;
    g_gauge_xml = config_buf.str();
    g_gauge_ready = true;

    write(xml_out, "Config_info", config_xml);
    write(xml_log, "Config_info", config_xml);

    std::istringstream monomial_is(trj_params.Monomials_xml);
    XMLReader monomial_reader(monomial_is);
    readNamedMonomialArray(monomial_reader, "/Monomials");

    std::istringstream h_mc_is(trj_params.H_MC_xml);
    XMLReader h_mc_xml(h_mc_is);
    ExactHamiltonianParams ham_params(h_mc_xml, "/Hamiltonian");
    Handle<AbsHamiltonian<multi1d<LatticeColorMatrix>,
                          multi1d<LatticeColorMatrix>>>
        h_mc(new ExactHamiltonian(ham_params));

    std::istringstream mdint_is(trj_params.Integrator_xml);
    XMLReader mdint_xml(mdint_is);
    LCMToplevelIntegratorParams int_params(mdint_xml, "/MDIntegrator");
    Handle<AbsMDIntegrator<multi1d<LatticeColorMatrix>,
                           multi1d<LatticeColorMatrix>>>
        integrator(new LCMToplevelIntegrator(int_params));

    LatColMatSMDTrj the_smd_trj(h_mc, integrator, trj_params.gamma,
                                trj_params.accept_reject,
                                trj_params.measure_actions);

    std::istringstream inline_is(mc_control.inline_measurement_xml);
    XMLReader inline_xml(inline_is);
    multi1d<Handle<AbsInlineMeasurement>> the_measurements;
    read(inline_xml, "/InlineMeasurements", the_measurements);

    doSMD<SMDTrjParams>(
        g_u, the_smd_trj, mc_control, trj_params, the_measurements);

    pop(xml_log);
    pop(xml_out);
  } catch (const std::string& e) {
    set_error(std::string("SMD run failed: ") + e);
    return 1;
  } catch (const std::exception& e) {
    set_error(std::string("SMD run failed: ") + e.what());
    return 1;
  } catch (...) {
    set_error("SMD run failed: unknown error.");
    return 1;
  }

  return 0;
}

}  // extern "C"
