require 'formula'

##
## HomeBrew/LinuxBrew formula for lobSTR.
##

class Lobstr < Formula
  homepage 'http://erlichlab.wi.mit.edu/lobSTR/'
  #head 'https://github.com/whitehead/lobstr-code.git'
  head 'https://github.com/agordon/lobstr-code.git', :branch => 'autotools_cleanup'

  if build.head?
    depends_on :autoconf => :build
    depends_on :automake => :build
  end
  depends_on :libtool => :build
  ## TODO:
  ## Add dependancies for:
  ##    boost
  ##    blas
  ##    libz
  ##    fftw3
  ##    cppunit
  ##    fortran
  ##    pkg-config

  def install
      if build.head?
          # Ugly hack: lobSTR extracts the version from either ".git" (if cloned)
          #            or from ".tarball-version" (if from a proper dist tarball).
          #            But here, HomeBrew clones first, then extract to a separate directory,
          #            So there's neither ".git" nor ".tarball-version".
          #            As a work-around, we create a dummy one.
          system 'echo 2.0.3-with-patches > .tarball-version'
      end
      system './reconf'
      system './configure', "--prefix=#{prefix}"
      system 'make'
      system 'make', 'install'
  end
end
