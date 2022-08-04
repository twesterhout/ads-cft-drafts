{ sources, compiler }:
[
  (final: prev: {
    inherit (import sources.gitignore { inherit (prev) lib; }) gitignoreFilter;
  })
  (final: prev: {
    af = import ./others/arrayfire.nix { inherit (prev) lib stdenv fetchurl cmake pkg-config fftw fftwFloat blas lapack boost; };
    halide = import ./others/halide.nix {
      inherit (prev) lib stdenv fetchFromGitHub cmake;
      llvmPackages = prev.llvmPackages_14;
    };
  })
  (final: prev: {
    metric = import ./packages.nix { pkgs = prev; inherit compiler; };
  })
]
