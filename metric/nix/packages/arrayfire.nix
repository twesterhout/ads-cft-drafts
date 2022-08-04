{ mkDerivation, af, base, Cabal, cabal-doctest, directory, fetchgit
, filepath, hspec, hspec-discover, lib, parsec, QuickCheck
, quickcheck-classes, text, vector
}:
mkDerivation {
  pname = "arrayfire";
  version = "0.6.1.0";
  src = fetchgit {
    url = "https://github.com/twesterhout/arrayfire-haskell";
    hash = "sha256-K3AKTzUus1ZyHMZJ/zifIofMRDXFvFtYRLQ9BtyKtjI=";
    rev = "20d4f02ead2e0788e9053a73214a1c90c1ba3aee";
    fetchSubmodules = true;
  };
  isLibrary = true;
  isExecutable = true;
  setupHaskellDepends = [ base Cabal cabal-doctest ];
  libraryHaskellDepends = [ base filepath vector ];
  librarySystemDepends = [ af ];
  executableHaskellDepends = [ base directory parsec text vector ];
  testHaskellDepends = [
    base directory hspec QuickCheck quickcheck-classes vector
  ];
  testToolDepends = [ hspec-discover ];
  homepage = "https://github.com/arrayfire/arrayfire-haskell";
  description = "Haskell bindings to the ArrayFire general-purpose GPU library";
  license = lib.licenses.bsd3;
}
