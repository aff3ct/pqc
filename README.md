# AFF3CT Post-Quantum Cryptography Modules

[**AFF3CT**](https://aff3ct.github.io/) is **a simulator** and **a library** dedicated to the Forward Error
Correction (FEC or **channel coding**). This repository implements software modules for the support of
Post-Quantum Cryptography schemes in AFF3CT. The currently supported encryption schemes are [Bike](https://bikesuite.org/), [Classic McEliece](https://classic.mceliece.org/index.html) and [Hamming Quasi-Cyclic](https://pqc-hqc.org/) (HQC).

The implementation relies on the following dependences:
- AFF3CT's [StreamPU](https://github.com/aff3ct/StreamPU) task-based runtime system for streaming (successfully tested with `v1.2.3`)
- The [FLINT](https://flintlib.org) Fast Library for Number Theory (successfully tested with `v3.1.3`).

The StreamPU repository should be extracted in the `streampu/` subdirectory of the PQC modules sources.

Example programs in subdirectory `src/` illustrate the construction of encrypting / decrypting communication chains for each of the supported schemes.

## License

The project is licensed under the MIT license.
