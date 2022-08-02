{
  description = "Numerically solving Einstein's equations";
  inputs.nixpkgs.url = "github:nixos/nixpkgs";
  inputs.flake-utils.url = "github:numtide/flake-utils";
  inputs.arrayfire.url = "github:twesterhout/nix-for-science?dir=arrayfire";

  outputs = { self, nixpkgs, flake-utils, arrayfire }:
  flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        t = pkgs.lib.trivial;
        hl = pkgs.haskell.lib;
        name = "project-name";

        project = devTools: # [1]
          let addBuildTools = (t.flip hl.addBuildTools) devTools;
          in pkgs.haskellPackages.developPackage {
            root = pkgs.lib.sourceFilesBySuffices ./. [ ".cabal" ".hs" ];
            name = name;
            returnShellEnv = !(devTools == [ ]); # [2]

            source-overrides = {
              arrayfire = builtins.fetchGit {
                url = "https://github.com/twesterhout/arrayfire-haskell.git";
                ref = "master";
                rev = "12b02600a59453f431d86f3c68d0763679753068";
              };
            };

            modifier = (t.flip t.pipe) [
              addBuildTools
              hl.dontHaddock
              hl.enableStaticLibraries
              hl.justStaticExecutables
              hl.disableLibraryProfiling
              hl.disableExecutableProfiling
            ];
          };

      in {
        packages.pkg = project [ ]; # [3]

        defaultPackage = self.packages.${system}.pkg;

        devShell = project [ # [4]
          # System dependencies
          arrayfire
          # pkgs.arrayfire
          # pkgs.petsc
          # Haskell build tools
          pkgs.haskellPackages.cabal-fmt
          pkgs.haskellPackages.cabal-install
        ];
      });
      # let
      #   pkgs = nixpkgs.legacyPackages.${system};
      # in
      # {
      #   devShell = pkgs.mkShell {
      #     buildInputs = [
      #       pkgs.arrayfire
      #       pkgs.petsc
      #     ];
      #   
      #     shellHook = ''
      #       export ARRAYFIRE=${pkgs.arrayfire}
      #       export PETSC=${pkgs.petsc}
      #     '';
      #   };
      # }
      # );
}
