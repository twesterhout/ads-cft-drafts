{ system ? builtins.currentSystem, compiler ? null }:
let
  pkgs = import ./nix { inherit system compiler; };
in
pkgs.mkShell {
  buildInputs = [
    pkgs.arrayfire
    pkgs.metric.shell
  ];
  shellHook = ''
    export LD_LIBRARY_PATH=${pkgs.metric.shell}/lib:${pkgs.arrayfire}/lib:$LD_LIBRARY_PATH
    logo
  '';
}
