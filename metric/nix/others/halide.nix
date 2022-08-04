{ lib, stdenv, fetchFromGitHub, cmake, llvmPackages
}:

llvmPackages.stdenv.mkDerivation rec {
  pname = "halide";
  version = "14.0.0";

  src = fetchFromGitHub {
    owner = "halide";
    repo = "Halide";
    rev = "v${version}";
    hash = "sha256-/7U2TBcpSAKeEyWncAbtW6Vk/cP+rp1CXtbIjvQMmZA=";
  };

  fixupPhase = ''
    true
  '';

  cmakeFlags = [
    "-DWITH_PYTHON_BINDINGS=OFF"
    "-DTARGET_WEBASSEMBLY=OFF"
  ];

  # To handle the lack of 'local' RPATH; required, as they call one of
  # their built binaries requiring their libs, in the build process.
  preBuild = ''
    export LD_LIBRARY_PATH="$(pwd)/src''${LD_LIBRARY_PATH:+:}$LD_LIBRARY_PATH"
  '';

  buildInputs = [
    llvmPackages.llvm
    llvmPackages.lld
    llvmPackages.openmp
    llvmPackages.libclang
  ];

  nativeBuildInputs = [ cmake ];

  meta = with lib; {
    description = "C++ based language for image processing and computational photography";
    homepage = "https://halide-lang.org";
    license = licenses.mit;
    platforms = platforms.all;
    maintainers = [ maintainers.ck3d ];
  };
}
