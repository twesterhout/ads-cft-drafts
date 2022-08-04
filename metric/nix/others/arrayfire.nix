{ lib, stdenv, fetchurl, cmake, pkg-config
, fftw, fftwFloat, blas, lapack, boost
}:

stdenv.mkDerivation rec {
  pname = "arrayfire";
  version = "3.8.2";

  src = fetchurl {
    url = "https://github.com/arrayfire/arrayfire/releases/download/v${version}/arrayfire-full-${version}.tar.bz2";
    hash = "sha256-LQGzWtreJDMHj1fiIzhEZ5qr/bBqQeaZKmsnxlMC0/4=";
  };

  patches = [ ./no-download.patch ];

  cmakeFlags = [
    "-DAF_BUILD_UNIFIED=ON"
    "-DAF_BUILD_CPU=ON"
    "-DAF_BUILD_CUDA=OFF"
    "-DAF_BUILD_OPENCL=OFF"
    "-DAF_BUILD_DOCS=OFF"
    "-DAF_BUILD_EXAMPLES=OFF"
    "-DBUILD_TESTING=OFF"
    "-DAF_COMPUTE_LIBRARY=FFTW/LAPACK/BLAS"
  ];

  nativeBuildInputs = [ cmake pkg-config ];
  strictDeps = true;
  buildInputs = [
    fftw fftwFloat
    blas lapack
    boost.out boost.dev
  ];

  # Post-installation fixup requires chmod, but it doesn't work under nix-portable
  fixupPhase = ''
    true
  '';

  meta = with lib; {
    description = "A general-purpose library for parallel and massively-parallel computations";
    longDescription = ''
      A general-purpose library that simplifies the process of developing software that targets parallel and massively-parallel architectures including CPUs, GPUs, and other hardware acceleration devices.";
    '';
    license = licenses.bsd3;
    homepage = "https://arrayfire.com/";
  };
}
