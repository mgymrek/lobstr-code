# lobSTR - Hacking notes

Tandem Repeats (STRs) profiler from high throughput sequencing data.

See: [lobSTR Website](http://erlichlab.wi.mit.edu/lobSTR/),  [lobSTR Code on Github](https://github.com/mgymrek/lobstr-code)


## Prerequisites

### On Debian/Ubuntu Linux

    $ sudo apt-get install gcc g++ gfortran \
                make automake autoconf \
                libtool pkg-config git \
                libfft3-dev libboost-dev libcppunit-dev \
                libz-dev libblas-dev


### On Mac OS X

    $ [TODO]

## Compiling Examples

### Standard Method

This should work for most POSIX-compatible systems (e.g. linux, MacOS-X), providing all the required tools (c/c++/fortran compiler, libtool, pkg-config) and libraries (blas, fftw3, cppunit, boost) have been installed.

    $ git clone https://github.com/mgymrek/lobstr-code
    $ cd lobstr-code
    $ ./reconf
    $ ./configure
    $ make
    $ make check
    $ sudo make install

### Build with CLANG

**clang** is a C/C++ compiler, commonly found on Mac-OS X, but also on linuxes. Use the following commands:

    $ ./configure CC=clang CXX=clang++
    $ make clean
    $ make

### Build for debugging

For debugging, turn optimization to a minimum, and add debug-information. Use the following commands:

    $ ./configure CFLAGS="-g -O0" CXXFLAGS="-g -O0"
    $ make clean
    $ make

### Build for maximum performance

Use the following commands:

    $ ./configure CFLAGS="-O3" CXXFLAGS="-O3"
    $ make clean
    $ make

### Build static executables
Static executables files do not require any existing shared-libraries to run, and can sometimes be copied as-is from one computer to another (providing they both have the same architecture and operating system):

    $ ./configure --enable-all-static
    $ make clean
    $ make


### Building faster
If your computer/server has multiple CPUs/cores, use the `-j` option for make. Example with 4 CPUs:

    $ make -j 4

### Building with Debian-Specific Hardening options
See [Debian Hardening Walkthought](https://wiki.debian.org/HardeningWalkthrough) for more information. Requires the `dpkg-dev` package (and a recent Debian/Ubuntu system).

    $ ./configure --enable-debian-hardening
    $ make clean
    $ make


## Advanced options

### Installing to custom location

By default, running `make install` will install the program into `/usr/local/bin` .

To change the installation directory, use:

    $ ./configure --prefix=/new/desired/path
    $ make
    $ make install

The programs will be installed to `/new/desired/path/bin/lobSTR` . Add `/new/desired/path/bin` to your $PATH environment variable.

### Using libraries from non-standard location

On most Linuxes, libraries (such as boost) are installed in the standard location of `/usr/lib` (also: `/usr/lib64`, `/usr/local/lib`, and few other locations).

The program `pkg-config` is configured to find libraries in these standard directories (if they are installed).

If some libraries are installed in non-standard location, `pkg-config` will not find them. This usually results in the following error message during `./configure` step:

    Checking for FFTW... configure: error: Package requirements (fftw3) were not met:

    No package 'fftw3' found

    Consider adjusting the PKG_CONFIG_PATH environment variable if you
    installed software in a non-standard prefix.

    Alternatively, you may set the environment variables FFTW_CFLAGS
    and FFTW_LIBS to avoid the need to call pkg-config.
    See the pkg-config man page for more details.

To solve this issue, add the directory containing the library's `.pc` file to `PKG_CONFIG_PATH`.

Example:

    # The library 'fftw3' was installed to non-standard directory:
    # /home/gordon/test/lib/libfftw3.so
    # The '.pc' file will likely be here:
    $ ls /home/gordon/test/lib/pkgconfig/*.pc
    /home/gordon/test/lib/pkgconfig/fftw3.pc

    # Update PKG_CONFIG_PATH, then run `configure`:
    $ export PKG_CONFIG_PATH=/home/gordon/test/lib/pkgconfig:$PKG_CONFIG_PATH
    $ ./configure

## For Developers

### Versioning

#### Current version

The build system calculates the version number based on the last git tag, and the number of commits past this tag.

To find the current tag, use:

    $ ./config/git-version-gen v
    2.0.3.52-28ec-dirty

Which means:

1. Last git tag was "2.0.3"
2. The current revision is 52 commits after this tag.
3. The current revision's SHA1 starts with "28ec"
4. There are uncomitted changes in the code ("-dirty")

This version string is used throughout the program (in the source code files and in the generated tarball file names).

#### Setting a new version

To tag a new version on your local git repository:

    $ git tag -a 'v2.0.4' -m 'Some Message'

To publish this tag to github

    $ git push --tags

**Important Notes:**

1. Do not omit the "v" at the beginning of the tag. It is used by `git-version-gen` to detect the relevant version tags.
2. When releasing a new tarball (hence, without a 'git' repository), a special hidden file `.version` will be created, which will hold the updated version string.
3. After updating (e.g. pulling/fetching/merging/commiting) the git repository, run `./reconf` to ensure the version string is updated.

### Making a new distribution tar-ball

To release a new tarball, run the following command:

    $ make dist

This will create a file named `lobSTR-X.Y.Z.tar.gz` (X.Y.Z being the calculated version string, e.g. '2.0.3-52-23e8-dirty' or '2.0.4', depending on your repository).

Before publishing the file, it is advisable to ensure the tarball contains all the necessary files required for successful build. Use the following command:

    $ make distcheck

### Cleanups

During development, cleanups (deleting all compiled and generated files) is sometimes required.

Use the following command to delete all compiled files:

    $ make clean

Use the following command to delete all compiled files, *and* all Makefiles (will require re-running `./configure`):

    $ make distclean

Use the following command to delete *all* generated files (compiled files, Makefiles, and Makefiles.in files). This will require re-running `./reconf` and `./configure`:

    $ make maintainer-clean


### Preparing a binary pre-compiled package

It is sometimes needed to publish a pre-compiled version of the program, to help users use it quickly (though compiling from source is always preferred, not all users have compilers installed, or the knowledge to required, or to time to spend).

One possible way to prepare a pre-compiled package:

    ./configure --enable-all-static \
               CFLAGS="-O3" CXXFLAGS="-O3" \
               --prefix=/tmp/lobSTR-2.0.3-binary_linux-x86_64/
    $ make clean
    $ make
    $ make install
    $ cd /tmp
    $ tar -czvf lobSTR-binary-2.0.3.tar.gz lobSTR-binary-2.0.3/

The `./configure` command will prepare the build, to make a static executable, and install it to `/tmp/lobSTR-2.0.3-binary_linux-x86_64/` .

The `make install` step will copy all the compiled files and model files to this directory.

The `tar` step will collect and compress all the files, generating a tarball for distribution.

**NOTE::** Pre-compiled binaries must be prepared for multiple systems (e.g. Linux-32bit, Linux-64bit, Mac-OS-X 32bit, Mac-OS-X 64bit, etc.).

