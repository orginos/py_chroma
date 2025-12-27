from py_chroma import (
    finalize,
    initialize,
    register_inline,
    run_plaquette,
    set_gauge,
    set_lattice,
    set_rng_seed,
)


def main() -> None:
    initialize()
    register_inline()
    set_lattice([2, 2, 2, 2])
    set_rng_seed(11)
    set_gauge("UNIT", "<Param></Param>", "/Param")
    run_plaquette(update_no=0, frequency=1)
    finalize()
    print("py_chroma plaquette ok")


if __name__ == "__main__":
    main()
