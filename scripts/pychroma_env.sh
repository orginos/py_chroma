#!/usr/bin/env bash
set -euo pipefail

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  echo "Source this script to set PYCHROMA variables." >&2
  echo "Example: source /path/to/pychroma_env.sh [--prefix /path/to/install]" >&2
  exit 2
fi

_pychroma_env_die() {
  echo "pychroma_env: $*" >&2
  return 1
}

_pychroma_prefix=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix)
      _pychroma_prefix="${2:-}"
      shift 2
      ;;
    --prefix=*)
      _pychroma_prefix="${1#--prefix=}"
      shift
      ;;
    --*)
      _pychroma_env_die "unknown option: $1"
      return 1
      ;;
    *)
      shift
      ;;
  esac
done

if [[ -z "${_pychroma_prefix}" ]]; then
  if [[ -n "${PYCHROMA_PREFIX:-}" ]]; then
    _pychroma_prefix="${PYCHROMA_PREFIX}"
  else
    _script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    if [[ -d "${_script_dir}/../python/py_chroma" ]]; then
      _pychroma_prefix="$(cd "${_script_dir}/.." && pwd)"
    elif [[ -d "${_script_dir}/../py_chroma" ]]; then
      _pychroma_prefix="$(cd "${_script_dir}/.." && pwd)"
    else
      _pychroma_env_die "could not infer prefix; pass --prefix or set PYCHROMA_PREFIX"
      return 1
    fi
  fi
fi

_uname="$(uname -s)"
if [[ "${_uname}" == "Darwin" ]]; then
  _lib_name="libpychroma.dylib"
else
  _lib_name="libpychroma.so"
fi

_lib_path=""
for _candidate in "${_pychroma_prefix}/lib/${_lib_name}" "${_pychroma_prefix}/build/${_lib_name}"; do
  if [[ -f "${_candidate}" ]]; then
    _lib_path="${_candidate}"
    break
  fi
done

if [[ -n "${_lib_path}" ]]; then
  export PYCHROMA_LIB="${_lib_path}"
else
  _pychroma_env_die "libpychroma not found under ${_pychroma_prefix}/lib or ${_pychroma_prefix}/build"
  return 1
fi

if [[ -d "${_pychroma_prefix}/python/py_chroma" ]]; then
  export PYTHONPATH="${_pychroma_prefix}/python${PYTHONPATH:+:${PYTHONPATH}}"
elif [[ -d "${_pychroma_prefix}/py_chroma" ]]; then
  export PYTHONPATH="${_pychroma_prefix}${PYTHONPATH:+:${PYTHONPATH}}"
else
  _pychroma_env_die "py_chroma package not found under ${_pychroma_prefix}"
  return 1
fi

export PYCHROMA_PREFIX="${_pychroma_prefix}"
