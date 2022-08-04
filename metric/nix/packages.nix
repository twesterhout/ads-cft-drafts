{ pkgs, compiler }:
let
  lib = pkgs.lib;

  util = import ./util.nix {
    inherit pkgs;
    inherit (pkgs) lib gitignoreFilter;
  };

  conf = lib.importTOML ../nixkell.toml;

  ghcVersion = if compiler != null then compiler else conf.ghc;

  # Create our own setup using our choosen GHC version as a starting point
  ourHaskell = pkgs.haskell.packages.${("ghc" + util.removeChar "." ghcVersion)}.override {
    overrides =
      let
        depsFromDir = pkgs.haskell.lib.packagesFromDirectory {
          directory = ./packages;
        };
        manual = _hfinal: hprev: {
          metric =
            let
              filteredSrc = util.filterSrc ../. {
                ignoreFiles = conf.ignore.files;
                ignorePaths = conf.ignore.paths;
              };
            in
            hprev.callCabal2nix "metric" filteredSrc {
              ads_cft_kernels = pkgs.ads-cft-kernels;
            };
        };
      in
      lib.composeExtensions depsFromDir manual;
  };

  # Add our package with its dependencies to GHC
  ghc = ourHaskell.ghc.withPackages (_ps:
    pkgs.haskell.lib.getHaskellBuildInputs ourHaskell.metric
  );

  # NB. If HLS is present in the list of tools we ensure
  # that it is compiled with the correct GHC version.
  tools =
    let
      hls = "haskell-language-server";
      hasHLS = ps: util.has hls ps || util.has ("haskellPackages." + hls) ps;
      removeHLS = ps: util.remove hls (util.remove ("haskellPackages." + hls) ps);
    in
    if hasHLS conf.env.tools
    then map util.getDrv (removeHLS conf.env.tools) ++ [ ourHaskell.haskell-language-server ]
    else map util.getDrv conf.env.tools;

  scripts = import ./scripts.nix { inherit pkgs; };
in
{
  bin = util.leanPkg ourHaskell.metric;

  shell = pkgs.buildEnv {
    name = "metric-env";
    paths = [ ghc ] ++ tools ++ scripts;
  };
}
