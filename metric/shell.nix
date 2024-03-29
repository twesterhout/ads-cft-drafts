{ system ? builtins.currentSystem, compiler ? null }:
let
  pkgs = import ./nix { inherit system compiler; };
in
pkgs.mkShell {
  buildInputs = [
    pkgs.arrayfire
    pkgs.halide
    pkgs.ads-cft-kernels
    pkgs.petsc
    pkgs.metric.shell
  ];
  shellHook = ''
    export ARRAYFIRE_PATH=${pkgs.arrayfire}
    export HALIDE_PATH=${pkgs.halide}
    export HDF5_PATH=${pkgs.hdf5}
    export PETSC_PATH=${pkgs.petsc}
    export KERNELS_PATH=${pkgs.ads-cft-kernels}
    export LD_LIBRARY_PATH=${pkgs.metric.shell}/lib:${pkgs.ads-cft-kernels}/lib:${pkgs.petsc}/lib:${pkgs.hdf5}/lib:${pkgs.halide}/lib:${pkgs.arrayfire}/lib:$LD_LIBRARY_PATH
    export PROMPT_COMMAND=
    export PS1='(nix) \w $ '
  '';
}
