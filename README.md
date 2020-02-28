<div align="center">
    <a href="http://eptlib.github.io/" target="_blank">
        <img src="img/EPTlib_lib128.png" width="128" alt="EPTlib logo"></img>
    </a>
</div>

General information
===================

*EPTlib* is an open source, extensible collection of C++ implementations of electric properties tomography (EPT) methods.

This project is under the MIT License.
See [LICENSE](LICENSE) in the root directory of this package for the full text.

Building from source
====================

The following libraries needs to be downloaded and installed before _EPTlib_ compilation:
- [Boost C++](https://www.boost.org/)
- [Eigen 3](http://eigen.tuxfamily.org/)

Optionally, [Google Test](https://github.com/google/googletest) is required for testing the library.

Linux users
-----------

Download the package and move in its root directory

```sh
git clone https://github.com/EPTlib/eptlib eptlib
cd eptlib
```

Configure, compile and install the package

```sh
mkdir build
cd build
cmake ..
make
sudo make install
sudo ldconfig
```

Before installing, you can test the compiled library

```sh
ctest
```

Acknowledgement
===============

This package has been developed in the framework of the 18HLT05 QUIERO Project. This project has received funding from the EMPIR Programme, co-financed by the Participating States and from the European Union's Horizon 2020 Research and Innovation Programme.

[![](img/logo-empir-euramet.png)](https://www.euramet.org/research-innovation/research-empir/)
[![](img/logo-quiero.png)](https://quiero-project.eu)

Credits
=======

The font [_Dancing Script_](https://github.com/impallari/DancingScript) has been used for the logo.
