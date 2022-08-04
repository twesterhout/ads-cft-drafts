{ mkDerivation, af, base, Cabal, cabal-doctest, directory, fetchgit
, filepath, hspec, hspec-discover, lib, parsec, QuickCheck
, quickcheck-classes, text, vector
}:
mkDerivation {
  pname = "arrayfire";
  version = "0.6.0.1";
  src = fetchgit {
    url = "https://github.com/twesterhout/arrayfire-haskell";
    hash = "sha256-Z5bi6b2R1rWmTuF8Vv0t83AQeGY0UmrwTAeWFG/Vf9M=";
    rev = "20d81d1f84f154eecc2c3e4e932a111a973726a1";
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
