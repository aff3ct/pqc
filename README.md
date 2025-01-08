# AFF3CT Post-Quantum Cryptography Modules

[**AFF3CT**](https://aff3ct.github.io/) is **a simulator** and **a library** dedicated to the Forward Error
Correction (FEC or **channel coding**). This repository implements software modules for the support of
Post-Quantum Cryptography schemes in AFF3CT. The currently supported encryption schemes are [Bike](https://bikesuite.org/), [Classic McEliece](https://classic.mceliece.org/index.html) and [Hamming Quasi-Cyclic](https://pqc-hqc.org/) (HQC).

The implementation relies on AFF3CT's [StreamPU](https://github.com/aff3ct/StreamPU) task-based runtime system for streaming, and on the [FLINT](https://flintlib.org) Fast Library for Number Theory.

Example programs in subdirectory `src/` illustrate the construction of encrypting / decrypting communication chains for each of the supported schemes.

## License

The project is licensed under the MIT license.
