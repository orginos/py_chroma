import ctypes
import os
import platform
from pathlib import Path


class PyChromaError(RuntimeError):
    pass


def _lib_name():
    system = platform.system()
    if system == "Darwin":
        return "libpychroma.dylib"
    if system == "Linux":
        return "libpychroma.so"
    raise PyChromaError(f"Unsupported platform: {system}")


def _default_lib_path():
    root = Path(__file__).resolve().parents[2]
    return root / "build" / _lib_name()


def _load_lib():
    lib_path = os.environ.get("PYCHROMA_LIB")
    if lib_path:
        return ctypes.CDLL(lib_path)
    return ctypes.CDLL(str(_default_lib_path()))


_lib = _load_lib()

_lib.pychroma_initialize.argtypes = []
_lib.pychroma_initialize.restype = ctypes.c_int

_lib.pychroma_finalize.argtypes = []
_lib.pychroma_finalize.restype = ctypes.c_int

_lib.pychroma_register_inline.argtypes = []
_lib.pychroma_register_inline.restype = ctypes.c_int

_lib.pychroma_set_lattice.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_int)]
_lib.pychroma_set_lattice.restype = ctypes.c_int

_lib.pychroma_set_rng_seed.argtypes = [ctypes.c_ulong]
_lib.pychroma_set_rng_seed.restype = ctypes.c_int
_lib.pychroma_set_rng_seed_xml.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
_lib.pychroma_set_rng_seed_xml.restype = ctypes.c_int

_lib.pychroma_set_gauge.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]
_lib.pychroma_set_gauge.restype = ctypes.c_int
_lib.pychroma_random_gauge_transform.argtypes = []
_lib.pychroma_random_gauge_transform.restype = ctypes.c_int

_lib.pychroma_set_input_xml.argtypes = [ctypes.c_char_p]
_lib.pychroma_set_input_xml.restype = ctypes.c_int

_lib.pychroma_rect_wloops.argtypes = [
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_char_p,
]
_lib.pychroma_rect_wloops.restype = ctypes.c_int

_lib.pychroma_run_inline_xml.argtypes = [ctypes.c_char_p]
_lib.pychroma_run_inline_xml.restype = ctypes.c_int

_lib.pychroma_run_plaquette.argtypes = [ctypes.c_ulong, ctypes.c_ulong]
_lib.pychroma_run_plaquette.restype = ctypes.c_int

_lib.pychroma_run_hmc_xml.argtypes = [ctypes.c_char_p]
_lib.pychroma_run_hmc_xml.restype = ctypes.c_int

_lib.pychroma_last_error.argtypes = []
_lib.pychroma_last_error.restype = ctypes.c_char_p


def _check(rc):
    if rc != 0:
        msg = _lib.pychroma_last_error()
        raise PyChromaError(msg.decode("utf-8") if msg else "py_chroma error")


def initialize():
    _check(_lib.pychroma_initialize())


def finalize():
    _check(_lib.pychroma_finalize())


def register_inline():
    _check(_lib.pychroma_register_inline())


def set_lattice(nrow):
    if not nrow:
        raise PyChromaError("nrow must be a non-empty sequence")
    arr = (ctypes.c_int * len(nrow))(*nrow)
    _check(_lib.pychroma_set_lattice(len(nrow), arr))


def set_rng_seed(seed):
    _check(_lib.pychroma_set_rng_seed(int(seed)))


def set_rng_seed_xml(rng_xml, rng_path="/RNG"):
    if not isinstance(rng_xml, str):
        raise PyChromaError("rng_xml must be a string")
    if not isinstance(rng_path, str):
        raise PyChromaError("rng_path must be a string")
    _check(
        _lib.pychroma_set_rng_seed_xml(
            rng_xml.encode("utf-8"),
            rng_path.encode("utf-8"),
        )
    )


def set_gauge(cfg_id, cfg_xml, cfg_path="/"):
    if not isinstance(cfg_id, str):
        raise PyChromaError("cfg_id must be a string")
    if not isinstance(cfg_xml, str):
        raise PyChromaError("cfg_xml must be a string")
    if not isinstance(cfg_path, str):
        raise PyChromaError("cfg_path must be a string")
    _check(
        _lib.pychroma_set_gauge(
            cfg_id.encode("utf-8"),
            cfg_xml.encode("utf-8"),
            cfg_path.encode("utf-8"),
        )
    )


def random_gauge_transform():
    _check(_lib.pychroma_random_gauge_transform())


def set_input_xml(input_xml):
    if not isinstance(input_xml, str):
        raise PyChromaError("input_xml must be a string")
    _check(_lib.pychroma_set_input_xml(input_xml.encode("utf-8")))


def run_inline_xml(inline_xml):
    if not isinstance(inline_xml, str):
        raise PyChromaError("inline_xml must be a string")
    _check(_lib.pychroma_run_inline_xml(inline_xml.encode("utf-8")))


def rect_wilson_loops(t_dir, z_dir, l_max, r_max, out_path):
    if not isinstance(out_path, str):
        raise PyChromaError("out_path must be a string")
    _check(
        _lib.pychroma_rect_wloops(
            int(t_dir),
            int(z_dir),
            int(l_max),
            int(r_max),
            out_path.encode("utf-8"),
        )
    )


def run_plaquette(update_no=0, frequency=1):
    _check(_lib.pychroma_run_plaquette(int(update_no), int(frequency)))


def run_hmc_xml(params_xml):
    if not isinstance(params_xml, str):
        raise PyChromaError("params_xml must be a string")
    _check(_lib.pychroma_run_hmc_xml(params_xml.encode("utf-8")))
