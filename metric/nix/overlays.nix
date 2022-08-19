{ sources, compiler }:
[
  (final: prev: {
    inherit (import sources.gitignore { inherit (prev) lib; }) gitignoreFilter;
  })
  (final: prev: rec {
    af = import ./others/arrayfire.nix {
      inherit (prev) lib stdenv fetchurl cmake pkg-config fftw fftwFloat blas lapack boost;
    };
    halide = import ./others/halide.nix {
      inherit (prev) lib stdenv fetchFromGitHub cmake;
      llvmPackages = prev.llvmPackages_14;
    };
    ads-cft-kernels = import ./others/ads-cft-kernels.nix {
      inherit (prev) lib stdenv cmake libpng libjpeg;
      inherit halide;
    };
    petsc-withp4est = false;
    petsc = prev.petsc.overrideAttrs (oldAttrs: {
      preConfigure = ''
        export FC="${prev.gfortran}/bin/gfortran" F77="${prev.gfortran}/bin/gfortran"
        patchShebangs ./lib/petsc/bin
        configureFlagsArray=(
          $configureFlagsArray
          "--with-debugging=0"
          "--with-mpi=0"
          "--with-blas=1"
          "--with-lapack=1"
        )
      '';
      nativeBuildInputs = [ prev.python3 prev.gfortran ];
      buildInputs = [ prev.blas prev.lapack ];
    });
  })
  (final: prev: {
    metric = import ./packages.nix { pkgs = prev; inherit compiler; };
  })
]
