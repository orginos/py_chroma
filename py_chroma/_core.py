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
    candidates = [
        root / "build" / _lib_name(),
        root / "lib" / _lib_name(),
    ]
    for path in candidates:
        if path.exists():
            return path
    return candidates[0]


def _load_lib():
    lib_path = os.environ.get("PYCHROMA_LIB")
    if lib_path:
        return ctypes.CDLL(lib_path)
    return ctypes.CDLL(str(_default_lib_path()))


_lib = None
_bound = False


def _bind_lib(lib):
    global _bound
    if _bound:
        return
    lib.pychroma_initialize.argtypes = []
    lib.pychroma_initialize.restype = ctypes.c_int

    lib.pychroma_finalize.argtypes = []
    lib.pychroma_finalize.restype = ctypes.c_int

    lib.pychroma_register_inline.argtypes = []
    lib.pychroma_register_inline.restype = ctypes.c_int

    lib.pychroma_set_lattice.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_int)]
    lib.pychroma_set_lattice.restype = ctypes.c_int

    lib.pychroma_set_rng_seed.argtypes = [ctypes.c_ulong]
    lib.pychroma_set_rng_seed.restype = ctypes.c_int
    lib.pychroma_set_rng_seed_xml.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
    lib.pychroma_set_rng_seed_xml.restype = ctypes.c_int

    lib.pychroma_set_gauge.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]
    lib.pychroma_set_gauge.restype = ctypes.c_int
    lib.pychroma_random_gauge_transform.argtypes = []
    lib.pychroma_random_gauge_transform.restype = ctypes.c_int

    lib.pychroma_set_input_xml.argtypes = [ctypes.c_char_p]
    lib.pychroma_set_input_xml.restype = ctypes.c_int

    lib.pychroma_rect_wloops.argtypes = [
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_char_p,
    ]
    lib.pychroma_rect_wloops.restype = ctypes.c_int

    lib.pychroma_run_inline_xml.argtypes = [ctypes.c_char_p]
    lib.pychroma_run_inline_xml.restype = ctypes.c_int

    lib.pychroma_run_plaquette.argtypes = [ctypes.c_ulong, ctypes.c_ulong]
    lib.pychroma_run_plaquette.restype = ctypes.c_int

    lib.pychroma_run_hmc_xml.argtypes = [ctypes.c_char_p]
    lib.pychroma_run_hmc_xml.restype = ctypes.c_int
    lib.pychroma_run_smd_xml.argtypes = [ctypes.c_char_p]
    lib.pychroma_run_smd_xml.restype = ctypes.c_int

    lib.pychroma_last_error.argtypes = []
    lib.pychroma_last_error.restype = ctypes.c_char_p
    _bound = True


def get_lib():
    global _lib
    if _lib is None:
        _lib = _load_lib()
        _bind_lib(_lib)
    return _lib


def _check(rc):
    if rc != 0:
        lib = get_lib()
        msg = lib.pychroma_last_error()
        raise PyChromaError(msg.decode("utf-8") if msg else "py_chroma error")


def initialize():
    lib = get_lib()
    _check(lib.pychroma_initialize())


def finalize():
    lib = get_lib()
    _check(lib.pychroma_finalize())


def register_inline():
    lib = get_lib()
    _check(lib.pychroma_register_inline())


def set_lattice(nrow):
    if not nrow:
        raise PyChromaError("nrow must be a non-empty sequence")
    arr = (ctypes.c_int * len(nrow))(*nrow)
    lib = get_lib()
    _check(lib.pychroma_set_lattice(len(nrow), arr))


def set_rng_seed(seed):
    lib = get_lib()
    _check(lib.pychroma_set_rng_seed(int(seed)))


def set_rng_seed_xml(rng_xml, rng_path="/RNG"):
    if not isinstance(rng_xml, str):
        raise PyChromaError("rng_xml must be a string")
    if not isinstance(rng_path, str):
        raise PyChromaError("rng_path must be a string")
    lib = get_lib()
    _check(
        lib.pychroma_set_rng_seed_xml(
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
    lib = get_lib()
    _check(
        lib.pychroma_set_gauge(
            cfg_id.encode("utf-8"),
            cfg_xml.encode("utf-8"),
            cfg_path.encode("utf-8"),
        )
    )


def random_gauge_transform():
    lib = get_lib()
    _check(lib.pychroma_random_gauge_transform())


def set_input_xml(input_xml):
    if not isinstance(input_xml, str):
        raise PyChromaError("input_xml must be a string")
    lib = get_lib()
    _check(lib.pychroma_set_input_xml(input_xml.encode("utf-8")))


def run_inline_xml(inline_xml):
    if not isinstance(inline_xml, str):
        raise PyChromaError("inline_xml must be a string")
    lib = get_lib()
    _check(lib.pychroma_run_inline_xml(inline_xml.encode("utf-8")))


def rect_wilson_loops(t_dir, z_dir, l_max, r_max, out_path):
    if not isinstance(out_path, str):
        raise PyChromaError("out_path must be a string")
    lib = get_lib()
    _check(
        lib.pychroma_rect_wloops(
            int(t_dir),
            int(z_dir),
            int(l_max),
            int(r_max),
            out_path.encode("utf-8"),
        )
    )


def run_plaquette(update_no=0, frequency=1):
    lib = get_lib()
    _check(lib.pychroma_run_plaquette(int(update_no), int(frequency)))


def run_hmc_xml(params_xml):
    if not isinstance(params_xml, str):
        raise PyChromaError("params_xml must be a string")
    lib = get_lib()
    _check(lib.pychroma_run_hmc_xml(params_xml.encode("utf-8")))


def run_smd_xml(params_xml):
    if not isinstance(params_xml, str):
        raise PyChromaError("params_xml must be a string")
    lib = get_lib()
    _check(lib.pychroma_run_smd_xml(params_xml.encode("utf-8")))
