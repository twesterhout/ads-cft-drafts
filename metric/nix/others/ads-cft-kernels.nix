{ lib, stdenv, cmake, halide, libpng, libjpeg }:

stdenv.mkDerivation rec {
  pname = "ads-cft-kernels";
  version = "0.0.1";
  src = ../../kernels;

  fixupPhase = ''
    true
  '';

  cmakeFlags = [
    "-DBUILD_SHARED_LIBS=ON"
  ];

  # To handle the lack of 'local' RPATH; required, as they call one of
  # their built binaries requiring their libs, in the build process.
  # preBuild = ''
  #   export LD_LIBRARY_PATH="$(pwd)/src''${LD_LIBRARY_PATH:+:}$LD_LIBRARY_PATH"
  # '';

  buildInputs = [
    halide
    # Dependencies of Halide
    libpng
    libjpeg
  ];

  nativeBuildInputs = [ cmake ];

  meta = with lib; {
    description = "Low-level kernels required for the solutions of Einstein's equations";
    license = licenses.bsd3;
    platforms = platforms.all;
  };
}

