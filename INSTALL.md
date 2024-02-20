# Debug Installation of SkolemFC

If you encounter errors when executing `install.sh`, it's necessary to manually compile all dependencies before building SkolemFC itself. Follow this sequence for building the required dependencies:
 - [louvain-community](https://github.com/meelgroup/louvain-community)
 - [cryptominisat](https://github.com/msoos/cryptominisat)
 - [arjun](https://github.com/meelgroup/arjun)
 - [approxmc](https://github.com/arijitsh/approxmc)
 - [unisamp](https://github.com/arijitsh/unigen)

SkolemFC integrates these softwares as libraries, potentially in a recursive manner. Since these libraries' functionalities may evolve, we recommend compiling them from the tested versions provided in the `/deps/` directory to ensure compatibility.

Should your system already contain any version of these softwares, `cmake` could identify the incorrect version during the SkolemFC build process. To resolve this, verify which version is being utilized by `cmake` for SkolemFC and temporarily remove any conflicting builds.

# Installing in other platforms

On systems other than Linux, `install.sh` may not function out-of-the-box. You will need to manually build all dependencies and then SkolemFC itself. Assuming `cmake` and other necessary dependencies are correctly installed on your system, compiling SkolemFC should be straightforward, although this has not been extensively tested.

# Static Build

To create a static build of SkolemFC, all dependencies must be compiled statically. Use `cmake -DSTATICCOMPILE=ON ..` to build these softwares statically. Running `./utils/install_static.sh` from the root directory will compile SkolemFC and its dependencies into a static build.