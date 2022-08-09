{ mkDerivation, af, base, Cabal, cabal-doctest, directory, fetchgit
, filepath, hspec, hspec-discover, lib, parsec, QuickCheck
, quickcheck-classes, text, vector
}:
mkDerivation {
  pname = "arrayfire";
  version = "0.6.1.0";
  src = fetchgit {
    url = "https://github.com/twesterhout/arrayfire-haskell";
    hash = "sha256-zmaC9KeWT3WL1g+/0P/wngn2y4U4Y1PXx48vv/y8veU=";
    rev = "db21c00f72e8b6da8e886e6d5b701e2aa842de53";
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
