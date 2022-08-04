{ sources, compiler }:
[
  (final: prev: {
    inherit (import sources.gitignore { inherit (prev) lib; }) gitignoreFilter;
  })
  (final: prev: rec {
    af = import ./others/arrayfire.nix { inherit (prev) lib stdenv fetchurl cmake pkg-config fftw fftwFloat blas lapack boost; };
    halide = import ./others/halide.nix {
      inherit (prev) lib stdenv fetchFromGitHub cmake;
      llvmPackages = prev.llvmPackages_14;
    };
    ads-cft-kernels = import ./others/ads-cft-kernels.nix {
      inherit (prev) lib stdenv cmake libpng libjpeg;
      inherit halide;
    };
  })
  (final: prev: {
    metric = import ./packages.nix { pkgs = prev; inherit compiler; };
  })
]
