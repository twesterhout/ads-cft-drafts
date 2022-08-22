{ mkDerivation, af, base, Cabal, cabal-doctest, directory, fetchgit
, filepath, hspec, hspec-discover, lib, parsec, QuickCheck
, quickcheck-classes, text, vector
}:
mkDerivation {
  pname = "arrayfire";
  version = "0.6.1.0";
  src = fetchgit {
    url = "https://github.com/twesterhout/arrayfire-haskell";
    hash = "sha256-mCABythTKL5VPdk+jJNfjqVdCI51hlhbTSk1wsd2nZ4=";
    rev = "48aa2c35e196f5ad05861ac4a813a7cd0c04d0c7";
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
