# AFF3CT Post-Quantum Cryptography Modules

[**AFF3CT**](https://aff3ct.github.io/) is **a simulator** and **a library** dedicated to the Forward Error
Correction (FEC or **channel coding**). This repository implements software modules for the support of
Post-Quantum Cryptography schemes in AFF3CT. The currently supported encryption schemes are [Bike](https://bikesuite.org/), [Classic McEliece](https://classic.mceliece.org/index.html) and [Hamming Quasi-Cyclic](https://pqc-hqc.org/) (HQC).

The implementation relies on the following dependences:
- The [GF2X]({https://libntl.org/doc/tour-gf2x.htm) library (successfully tested with `v1.3.0` LGPL).
- The [NTL]({https://libntl.org/doc/tour-gf2x.htm) Number Theory Library compiled (successfully tested with `v11.5.1`) compiled to use GF2X.
- The [FLINT](https://flintlib.org) Fast Library for Number Theory (successfully tested with `v3.1.3`).
- AFF3CT's [StreamPU](https://github.com/aff3ct/StreamPU) task-based runtime system for streaming (successfully tested with `v1.7.3`)

## GF2X
GF2X is detected using the `pkg-config` command. GF2X's installation `pkgconfig/' directory should be listed in the `PKG_CONFIG_PATH` environment variable.

## NTL
NTL is detected using the `NTL_DIR` environment variable, which should be set to NTL's top installation directory.

## FLINT
FLINT is detected using the `pkg-config` command. FLINT's installation `pkgconfig/' directory and the corresponding `pkgconfig/` directories of its own dependencies (GMP and MPFR libs) should all be listed in the `PKG_CONFIG_PATH` environment variable.

## StreamPU
StreamPU should either:
- be preinstalled, with its prefix install directory listed in the `CMAKE_PREFIX_PATH` environment variable;
- or be extracted as a git repository in the `streampu/` subdirectory of the PQC modules sources, in which case it will be built as part of the process.

## Building the PQC library and examples
Example programs in subdirectory `examples/` illustrate the construction of encrypting / decrypting communication chains for each of the supported schemes.

The library and examples should be compiled as follows:
```shell
export NTL_DIR="<NTL installation dir>"
export PKG_CONFIG_PATH="<GF2X dependencies pkgconfig dir>:"<FLINT and dependencies pkgconfig dirs>:${PKG_CONFIG_PATH}"
export CMAKE_PREFIX_PATH="<StreamPU install dir prefix>:${CMAKE_PREFIX_PATH}"
cmake -DCMAKE_INSTALL_PREFIX='<PQC install dir prefix>' -S '<PQC source dir>' -B '<PQC build dir>'
cmake --build '<PQC build dir>'
cmake --install '<PQC build dir>'
```

## License

The project is licensed under the MIT license.
