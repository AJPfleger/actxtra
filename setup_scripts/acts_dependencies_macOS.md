# Guide to get all dependencies for macOS for ACTS

This guide is based on [https://codimd.web.cern.ch/s/w-7j8zXm0](https://codimd.web.cern.ch/s/w-7j8zXm0) (CERN SSO might be required).

This guide was written for macOS 14.3 but might be incomplete. You could also have a look at [https://github.com/paulgessinger/ci-dependencies/tree/build_cmake](https://github.com/paulgessinger/ci-dependencies/tree/build_cmake). For older versions of this guide, go through the commit history.

## Install location

To keep everything clean all dependecies are tried to be installed in `/opt/hep`. For the installation process we generate a separeted folder.
```console
mkdir setup_dependencies && cd setup_dependencies
```

## Python

To ensure we always use the same python version, use [pyenv]().

```console
brew install pyenv
```
Next set the shell environment like in [this guide](https://github.com/pyenv/pyenv?tab=readme-ov-file#set-up-your-shell-environment-for-pyenv)
```console
echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.zprofile
echo '[[ -d $PYENV_ROOT/bin ]] && export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.zprofile
echo 'eval "$(pyenv init -)"' >> ~/.zprofile

exec "$SHELL"
```
If this doesn't work, try to set the other profiles. `pyenv init` might give some hints as well.
```console
echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bashrc
echo 'command -v pyenv >/dev/null || export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bashrc
echo 'eval "$(pyenv init -)"' >> ~/.bashrc

echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.profile
echo 'command -v pyenv >/dev/null || export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.profile
echo 'eval "$(pyenv init -)"' >> ~/.profile

echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bash_profile
echo '[[ -d $PYENV_ROOT/bin ]] && export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bash_profile
echo 'eval "$(pyenv init -)"' >> ~/.bash_profile

echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.zshrc
echo '[[ -d $PYENV_ROOT/bin ]] && export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.zshrc
echo 'eval "$(pyenv init -)"' >> ~/.zshrc

exec "$SHELL"
```

```console
pyenv install 3.12.1
pyenv global 3.12.1
```
Create a virtual environment, that we will use forever.
```console
python3 -m venv venv-acts
source venv-acts/bin/activate
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
    -Dpyroot=ON \
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
    -Dbuiltin_gl2ps=OFF \
    -Dbuiltin_xrootd=OFF \
    -Dbuiltin_pcre=ON \
    -Dbuiltin_lzma=ON \
    -Dbuiltin_zstd=ON \
    -Dbuiltin_lz4=ON \
    -Dgfal=OFF \
    -Ddavix=OFF \
    -Dbuiltin_vdt=ON \
    -Dxrootd=OFF \
    -Dtmva=OFF \
    -DPYTHON_INCLUDE_DIR="~/ACTS/setup_dependencies/venv-acts/bin/python" \
    -Druntime_cxxmodules=ON
sudo cmake --build build --target install -j8
cd ..
```


## PODIO

```console
mkdir podio && cd podio
git clone https://github.com/AIDASoft/podio.git source
cd source
git fetch --tags
git checkout tags/v00-99
cd ..
cmake -S source -B build \
    -GNinja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_STANDARD=17 \
    -DCMAKE_INSTALL_PREFIX="/opt/hep/podio/00-99" \
    -DCMAKE_PREFIX_PATH="/opt/hep/geant4/11.2.0;/opt/hep/pythia8/3810;/opt/hep/json/3.11.3;/opt/hep/root/6.30.02" \
    -DUSE_EXTERNAL_CATCH2=Off \
    -DBUILD_TESTING=Off
sudo cmake --build build --target install -j8
cd ..
```

## EDM4Hep

If you get a `Jinja2`-related error, you could try to use a more recent version of `EDM4Hep`.

```console
pip install jinja2 pyyaml
```

```console
mkdir edm4hep && cd edm4hep
git clone https://github.com/key4hep/EDM4hep.git source
cd source
git fetch --tags
git checkout tags/v00-10-03
cd ..
cmake -S source -B build \
    -GNinja \
    -DCMAKE_CXX_STANDARD=17 \
    -DCMAKE_INSTALL_PREFIX="/opt/hep/edm4hep/00-10-03" \
    -DCMAKE_PREFIX_PATH="/opt/hep/geant4/11.2.0;/opt/hep/pythia8/3810;/opt/hep/json/3.11.3;/opt/hep/root/6.30.02;/opt/hep/podio/00-99" \
    -DUSE_EXTERNAL_CATCH2=Off \
    -DBUILD_TESTING=OFF \
    -DCMAKE_BUILD_TYPE=Release
sudo cmake --build build --target install -j8
cd ..
```

## DD4hep

*probably fixed in DD4hep v01-27-02:* cmake version >= 3.27 introduces a [bug](https://bugzilla-geant4.kek.jp/show_bug.cgi?id=2556). It can be solved by deleting this file:
```console
sudo rm /opt/hep/geant4/lib/cmake/Geant4/Geant4PackageCache.cmake
```

```console
export LD_LIBRARY_PATH="/opt/hep/geant4/11.2.0/lib"
source /opt/hep/root/6.30.02/bin/thisroot.sh

mkdir dd4hep && cd dd4hep
git clone https://github.com/AIDASoft/DD4hep.git source
cd source
git fetch --tags
git checkout tags/v01-27-02
cd ..
cmake -S source -B build \
    -GNinja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_STANDARD=17 \
    -DCMAKE_INSTALL_PREFIX="/opt/hep/dd4hep/01-27-02" \
    -DDD4HEP_BUILD_PACKAGES="DDG4 DDDetectors DDRec UtilityApps" \
    -DDD4HEP_USE_GEANT4=ON \
    -DDD4HEP_USE_XERCESC=ON \
    -DDD4HEP_USE_EDM4HEP=ON \
    -DBUILD_TESTING=Off \
    -DBUILD_DOCS=OFF \
    -DCMAKE_PREFIX_PATH="/opt/hep/geant4/11.2.0;/opt/hep/pythia8/3810;/opt/hep/json/3.11.3;/opt/hep/root/6.30.02;/opt/hep/podio/00-99;/opt/hep/edm4hep/00-10-03"
sudo cmake --build build --target install -j8
cd ..
```
