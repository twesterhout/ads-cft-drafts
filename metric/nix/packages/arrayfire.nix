{ mkDerivation, af, base, Cabal, cabal-doctest, directory, fetchgit
, filepath, hspec, hspec-discover, lib, parsec, QuickCheck
, quickcheck-classes, text, vector
}:
mkDerivation {
  pname = "arrayfire";
  version = "0.6.1.0";
  src = fetchgit {
    url = "https://github.com/twesterhout/arrayfire-haskell";
    hash = "sha256-6xbzmJ7DeO9z3B5ueB+1aX+67VjzL8OQ7KKn8ggSE2E=";
    rev = "a95d381903ea329bdca9ba20132c2ff819bab5e2";
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
