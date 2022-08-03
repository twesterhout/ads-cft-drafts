{
  description = "Numerically solving Einstein's equations";
  inputs.haskellNix.url = "github:input-output-hk/haskell.nix";
  inputs.nixpkgs.follows = "haskellNix/nixpkgs-unstable";
  inputs.flake-utils.url = "github:numtide/flake-utils";
  inputs.arrayfire.url = "github:twesterhout/nix-for-science?dir=arrayfire";

  outputs = { self, nixpkgs, flake-utils, haskellNix, arrayfire }:
    flake-utils.lib.eachSystem [ "x86_64-linux" "x86_64-darwin" ] (system:
    let
      overlays = [ haskellNix.overlay
        (final: prev: {
          # source-overrides = {
          #   arrayfire = builtins.fetchGit {
          #     url = "https://github.com/twesterhout/arrayfire-haskell.git";
          #     ref = "master";
          #     rev = "44927d93e70d58a6fd2516a900796ad0a896a501";
          #   };
          # };
          af = arrayfire.packages.${system}.pkg;
          project =
            final.haskell-nix.project' {
              src = ./.;

              compiler-nix-name = "ghc902";

              # This is used by `nix develop .` to open a shell for use with
              # `cabal`, `hlint` and `haskell-language-server`
              shell.tools = {
                cabal = {};
                # hlint = {};
                # haskell-language-server = {};
              };
              # Non-Haskell shell tools go here
              shell.buildInputs = with pkgs; [
                nixpkgs-fmt
              ];
            };
        })
      ];
      pkgs = import nixpkgs { inherit system overlays; inherit (haskellNix) config; };
      flake = pkgs.project.flake {};
    in flake // {
      # Built by `nix build .`
      defaultPackage = flake.packages."arrayfire:exe:main";
    });
}

# {
#   description = "Numerically solving Einstein's equations";
#   inputs.nixpkgs.url = "github:nixos/nixpkgs";
#   inputs.flake-utils.url = "github:numtide/flake-utils";
#   inputs.arrayfire.url = "github:twesterhout/arrayfire-haskell";
#   inputs.halide.url = "github:twesterhout/nix-for-science?dir=halide";
# 
#   outputs = { self, nixpkgs, flake-utils, arrayfire, halide }:
#   flake-utils.lib.eachDefaultSystem (system:
#       let
#         pkgs = nixpkgs.legacyPackages.${system};
#         t = pkgs.lib.trivial;
#         hl = pkgs.haskell.lib;
#         name = "project-name";
# 
#         project = devTools: # [1]
#           let addBuildTools = (t.flip hl.addBuildTools) devTools;
#           in pkgs.haskellPackages.developPackage {
#             root = pkgs.lib.sourceFilesBySuffices ./. [ ".cabal" ".hs" ];
#             name = name;
#             returnShellEnv = !(devTools == [ ]); # [2]
#             
#             # source-overrides = {
#             #   arrayfire = builtins.fetchGit {
#             #     url = "https://github.com/twesterhout/arrayfire-haskell.git";
#             #     ref = "master";
#             #     rev = "44927d93e70d58a6fd2516a900796ad0a896a501";
#             #   };
#             # };
# 
#             modifier = (t.flip t.pipe) [
#               addBuildTools
#               hl.dontHaddock
#               hl.enableStaticLibraries
#               hl.justStaticExecutables
#               hl.disableLibraryProfiling
#               hl.disableExecutableProfiling
#             ];
#           };
# 
#       in {
#         packages.pkg = project [ ]; # [3]
# 
#         defaultPackage = self.packages.${system}.pkg;
# 
#         devShell = project [ # [4]
#           # System dependencies
#           # arrayfire
#           # pkgs.arrayfire
#           # pkgs.petsc
#           # Haskell build tools
#           pkgs.haskellPackages.cabal-fmt
#           pkgs.haskellPackages.cabal-install
#         ];
#       });
#       # let
#       #   pkgs = nixpkgs.legacyPackages.${system};
#       # in
#       # {
#       #   devShell = pkgs.mkShell {
#       #     buildInputs = [
#       #       pkgs.arrayfire
#       #       pkgs.petsc
#       #     ];
#       #   
#       #     shellHook = ''
#       #       export ARRAYFIRE=${pkgs.arrayfire}
#       #       export PETSC=${pkgs.petsc}
#       #     '';
#       #   };
#       # }
#       # );
# }
