Guide to get all dependencies for macOS for ACTS
================================================

This guide is based on [https://codimd.web.cern.ch/s/w-7j8zXm0](https://codimd.web.cern.ch/s/w-7j8zXm0) (CERN SSO might be required).

This guide was written for macOS 13.3.1  (a) but might be incomplete. You could also have a look at [https://github.com/paulgessinger/ci-dependencies/tree/build_cmake](https://github.com/paulgessinger/ci-dependencies/tree/build_cmake).

Install location
----------------

To keep everything clean all dependecies are tried to be installed in `/opt/hep`. For the installation process we generate a separeted folder.
```console
mkdir setup_dependencies
cd setup_dependencies
```

Brew Packages
-------------

```console
brew install cmake nlohmann-json
```

The following still need investigation if really needed
```console
brew install java
sudo ln -sfn /usr/local/opt/openjdk/libexec/openjdk.jdk /Library/Java/JavaVirtualMachines/openjdk.jdk
brew install glfw3 glew
```

xerces-c
--------
Does it also work with the current xerces version?
```console
mkdir xerces && cd xerces

# wget https://dlcdn.apache.org//xerces/c/3/sources/xerces-c-3.2.4.tar.gz
wget https://archive.apache.org/dist/xerces/c/3/sources/xerces-c-3.2.3.tar.gz
tar -zxvf xerces-c-3.2.3.tar.gz
cd xerces-c-3.2.3

cmake -S . -B build_opt_hep -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX=/opt/hep/xerces-c -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/usr/local/opt/icu4c
sudo cmake --build build_opt_hep --target install -j4
cd ../..
```

boost
-----

```console
mkdir boost && cd boost
wget https://boostorg.jfrog.io/artifactory/main/release/1.82.0/source/boost_1_82_0.tar.gz
tar -zxvf boost_1_82_0.tar.gz
cd boost_1_82_0
./bootstrap.sh --prefix=/opt/hep/boost
sudo ./b2 install --prefix=/opt/hep/boost
cd ../..
```

eigen
-----

```console
git clone https://gitlab.com/libeigen/eigen.git
cd eigen
git fetch --tags
git checkout tags/3.4.0
cmake -S . -B build_opt_hep -DCMAKE_INSTALL_PREFIX=/opt/hep/eigen -DCMAKE_CXX_STANDARD=17
sudo cmake --build build_opt_hep --target install -j4
cd ..
```

Geant4
------

```console
git clone https://gitlab.cern.ch/geant4/geant4.git
cd geant4
git fetch --tags
git checkout tags/v10.7.4
cmake -S . -B build_opt_hep -DCMAKE_PREFIX_PATH="/opt/hep/xerces-c" -DCMAKE_INSTALL_PREFIX=/opt/hep/geant4 -DGEANT4_BUILD_CXXSTD=17 -DGEANT4_USE_GDML=On -DGEANT4_INSTALL_DATA=On
sudo cmake --build build_opt_hep --target install -j8
cd ..
```

Pythia
------

```console
mkdir pythia && cd pythia
wget https://pythia.org/download/pythia83/pythia8307.tgz
tar -zxvf pythia8307.tgz
cd pythia8307
./configure --prefix=/opt/hep/pythia8
sudo make install -j4
cd ../..
```

Root
----
This is a bit troublesome. Apparently there has some problem with this macOS-version and root [https://root-forum.cern.ch/t/building-failed-after-upgrade-to-mac-os-13-3-1/54420/7](https://root-forum.cern.ch/t/building-failed-after-upgrade-to-mac-os-13-3-1/54420/7).

```console
mkdir root && cd root
git clone https://github.com/root-project/root.git source
cd source
git fetch --tags
git checkout tags/v6-28-00-patches
cd ..
cmake -S source -B build \
    -DCMAKE_CXX_STANDARD=17 \
    -DCMAKE_INSTALL_PREFIX="/opt/hep/root" \
    -DCMAKE_PREFIX_PATH="/opt/hep/xerces-c;/opt/hep/pythia8" \
    -Dbuiltin_glew=On
sudo cmake --build build_opt_hep --target install -j4
cd ..
```

HepMC3
------

```console
mkdir hepmc3 && cd hepmc3
git clone https://gitlab.cern.ch/hepmc/HepMC3.git source
cd source
git fetch --tags 
git checkout tags/3.2.5
cd ..
cmake -S source -B build \
    -DCMAKE_PREFIX_PATH="/opt/hep/root" \
    -DCMAKE_CXX_STANDARD=17 \
    -DCMAKE_INSTALL_PREFIX=/opt/hep/hepmc3
sudo cmake --build build --target install -j4
cd ..
```

PODIO
-----

```console
mkdir podio && cd podio
git clone https://github.com/AIDASoft/podio.git source
cd source
git fetch --tags
git checkout tags/v00-15
cd ..
cmake -S source -B build \
    -DCMAKE_CXX_STANDARD=17 \
    -DCMAKE_INSTALL_PREFIX="/opt/hep/podio" \
    -DCMAKE_PREFIX_PATH="/opt/hep/root;" \
    -DUSE_EXTERNAL_CATCH2=Off \
    -DBUILD_TESTING=Off
sudo cmake --build build --target install -j4
cd ..

```

EDM4Hep
-------
Does not work and gives the error:
```
  File "/opt/hep/podio/python/podio_class_generator.py", line 16, in <module>
    import jinja2
ModuleNotFoundError: No module named 'jinja2'
CMake Error at /opt/hep/podio/lib/cmake/podio/podioMacros.cmake:198 (message):
  Could not generate datamodel 'edm4hep'.  Check your definition in
  '../edm4hep.yaml'
Call Stack (most recent call first):
  edm4hep/CMakeLists.txt:4 (PODIO_GENERATE_DATAMODEL)
```
But `Jinja2` is already installed.

```console
mkdir edm4hep && cd edm4hep
git clone https://github.com/key4hep/EDM4hep.git source
cd source
git fetch --tags
git checkout tags/v00-07-02
cd ..
cmake -S source -B build \
    -DCMAKE_CXX_STANDARD=17 \
    -DCMAKE_INSTALL_PREFIX="/opt/hep/edm4hep;/" \
    -DCMAKE_PREFIX_PATH="/opt/hep/root;/opt/hep/podio;/opt/hep/hepmc3" \
    -DUSE_EXTERNAL_CATCH2=Off
sudo cmake --build build --target install -j4
cd ..
```










