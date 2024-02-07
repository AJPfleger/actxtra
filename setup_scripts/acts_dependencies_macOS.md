# Guide to get all dependencies for macOS for ACTS

This guide is based on [https://codimd.web.cern.ch/s/w-7j8zXm0](https://codimd.web.cern.ch/s/w-7j8zXm0) (CERN SSO might be required).

This guide was written for macOS 14.3 but might be incomplete. You could also have a look at [https://github.com/paulgessinger/ci-dependencies/tree/build_cmake](https://github.com/paulgessinger/ci-dependencies/tree/build_cmake). For older versions of this guide, go through the commit history.

## Install location

To keep everything clean all dependecies are tried to be installed in `/opt/hep`. For the installation process we generate a separeted folder.
```console
mkdir setup_dependencies && cd setup_dependencies
```

## Brew Packages

Last tested with:
- cmake `3.28.2`
- xerces-c `3.2.5`
- boost 1.84.0
```console
brew install cmake xerces-c boost
```


The following still need investigation if really needed:
```console
brew install java
sudo ln -sfn /usr/local/opt/openjdk/libexec/openjdk.jdk /Library/Java/JavaVirtualMachines/openjdk.jdk
brew install glfw3 glew
```

## Json

```console
mkdir json && cd json
wget https://github.com/nlohmann/json/archive/refs/tags/v3.11.3.tar.gz
mkdir source
tar -zxvf v3.11.3.tar.gz --strip-components=1 -C source
cmake -S source -B build \
    -GNinja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_STANDARD=17 \
    -DCMAKE_INSTALL_PREFIX=/opt/hep/json/3.11.3 \
    -DJSON_BuildTests=OFF
sudo cmake --build build --target install -j8
cd ..
```

## eigen

```console
mkdir eigen && cd eigen
git clone https://gitlab.com/libeigen/eigen.git source
cd source
git fetch --tags
git checkout tags/3.4.0
cd ..
cmake -S source -B build \
    -GNinja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=/opt/hep/eigen/3.4.0 \
    -DCMAKE_CXX_STANDARD=17
sudo cmake --build build --target install -j8
cd ..
```


## Geant4

```console
mkdir geant4 && cd geant4
git clone https://gitlab.cern.ch/geant4/geant4.git source
cd source
git fetch --tags
git checkout tags/v11.2.0
cd ..
cmake -S source -B build \
    -GNinja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_STANDARD=17 \
    -DCMAKE_INSTALL_PREFIX=/opt/hep/geant4/11.2.0 \
    -DGEANT4_BUILD_TLS_MODEL=global-dynamic \
    -DGEANT4_INSTALL_DATA=ON \
    -DGEANT4_USE_GDML=ON \
    -DGEANT4_USE_SYSTEM_EXPAT=ON \
    -DGEANT4_USE_SYSTEM_ZLIB=ON
sudo cmake --build build --target install -j8
cd ..
```

## Pythia

```console
mkdir pythia && cd pythia
wget https://pythia.org/download/pythia83/pythia8310.tgz
mkdir source
tar -zxvf pythia8310.tgz --strip-components=1 -C source
cd source
./configure --prefix=/opt/hep/pythia8/8310
sudo make install -j8
cd ../..
```

## Root
----
```console
mkdir root && cd root
wget https://root.cern/download/root_v6.30.02.source.tar.gz
mkdir source
tar -zxvf root_v6.30.02.source.tar.gz --strip-components=1 -C source
cmake -S source -B build \
    -GNinja \
    -DCMAKE_PREFIX_PATH="/opt/hep/geant4/11.2.0;/opt/hep/pythia8/3810;/opt/hep/json/3.11.3" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_STANDARD=17 \
    -DCMAKE_INSTALL_PREFIX=/opt/hep/root/6.30.02 \
    -Dfail-on-missing=ON \
    -Dgdml=ON \
    -Dx11=ON \
    -Dpyroot=On \
    -Ddataframe=ON \
    -Dmysql=OFF \
    -Doracle=OFF \
    -Dpgsql=OFF \
    -Dsqlite=OFF \
    -Dpythia6=OFF \
    -Dpythia8=OFF \
    -Dfftw3=OFF \
    -Dbuiltin_cfitsio=ON \
    -Dbuiltin_xxhash=ON \
    -Dbuiltin_afterimage=ON \
    -Dbuiltin_openssl=ON \
    -Dbuiltin_ftgl=ON \
    -Dbuiltin_glew=ON \
    -Dbuiltin_gsl=ON \
    -Dbuiltin_gl2ps=ON \
    -Dbuiltin_xrootd=ON \
    -Dbuiltin_pcre=ON \
    -Dbuiltin_lzma=ON \
    -Dbuiltin_zstd=ON \
    -Dbuiltin_lz4=ON \
    -Dgfal=OFF \
    -Ddavix=OFF \
    -Dbuiltin_vdt=ON \
    -Dxrootd=OFF \
    -Dtmva=OFF
sudo cmake --build build --target install -j8
cd ..
```







## PODIO
-----

```console
mkdir podio && cd podio
git clone https://github.com/AIDASoft/podio.git source
cd source
git fetch --tags
git checkout tags/v00-17-04
cd ..
cmake -S source -B build \
    -GNinja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_STANDARD=17 \
    -DCMAKE_INSTALL_PREFIX="/opt/hep/podio/00-17-04" \
    -DCMAKE_PREFIX_PATH="/opt/hep/geant4/11.2.0;/opt/hep/pythia8/3810;/opt/hep/json/3.11.3;/opt/hep/root/6.30.02" \
    -DUSE_EXTERNAL_CATCH2=Off \
    -DBUILD_TESTING=Off
sudo cmake --build build --target install -j8
cd ..
```

## EDM4Hep
-------
If you get a `Jinja2`-related error, you could try to use a more recent version of `EDM4Hep`.

It might be required to install these two libraries. Remember to use a virtual enviroment.
```console
pip install jinja2 pyyaml
```

We are not the most recent version `tags/v00-07-02` (at the time of writing), because it does not compile. You could use [patch](https://patch-diff.githubusercontent.com/raw/key4hep/EDM4hep/pull/201.patch) to make it work.

```console
mkdir edm4hep && cd edm4hep
git clone https://github.com/key4hep/EDM4hep.git source
cd source
git fetch --tags
git checkout tags/v00-07
cd ..
cmake -S source -B build \
    -DCMAKE_CXX_STANDARD=17 \
    -DCMAKE_INSTALL_PREFIX="/opt/hep/edm4hep" \
    -DCMAKE_PREFIX_PATH="/opt/hep/xerces-c;/opt/hep/geant4;/opt/hep/pythia8;/opt/hep/json;/opt/hep/hepmc3;/opt/hep/root;/opt/hep/podio" \
    -DUSE_EXTERNAL_CATCH2=Off \
    -DBUILD_TESTING=OFF \
    -DCMAKE_BUILD_TYPE=Release
sudo cmake --build build --target install -j8
cd ..
```

## DD4hep
------
cmake version >= 3.27 introduces a [bug](https://bugzilla-geant4.kek.jp/show_bug.cgi?id=2556). It can be solved by deleting this file:
```console
sudo rm /opt/hep/geant4/lib/cmake/Geant4/Geant4PackageCache.cmake
```

```console
export LD_LIBRARY_PATH="/opt/hep/geant4/lib"
source /opt/hep/root/bin/thisroot.sh

mkdir dd4hep && cd dd4hep
git clone https://github.com/AIDASoft/DD4hep.git source
cd source
git fetch --tags
git checkout tags/v01-25-01
cd ..
cmake -S source -B build \
    -DCMAKE_CXX_STANDARD=17 \
    -DCMAKE_INSTALL_PREFIX="/opt/hep/dd4hep" \
    -DDD4HEP_USE_GEANT4=On \
    -DDD4HEP_USE_EDM4HEP=On \
    -DBUILD_TESTING=Off \
    -DCMAKE_BUILD_TYPE=Release \
    -DDD4HEP_USE_XERCESC=ON \
    -DBUILD_DOCS=OFF \
    -DCMAKE_PREFIX_PATH="/opt/hep/xerces-c;/opt/hep/geant4;/opt/hep/pythia8;/opt/hep/json;/opt/hep/hepmc3;/opt/hep/root;/opt/hep/podio;/opt/hep/edm4hep"
sudo cmake --build build --target install -j8
cd ..
```
