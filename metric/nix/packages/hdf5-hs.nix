{ mkDerivation, base, bytestring, bytestring-to-vector
, cabal-doctest, containers, directory, fetchgit, hdf5
, hspec, inline-c, lib, list-t, relude, resourcet, safe-exceptions
, some, template-haskell, text, unliftio-core, vector, zlib
}:
mkDerivation {
  pname = "hdf5-hs";
  version = "0.2.0.1";
  src = fetchgit {
    url = "https://github.com/twesterhout/hdf5-hs";
    sha256 = "sha256-0M9McTjt95IRGo1blJocvQT+eFL9yLL8oGPJZFWxHGo=";
    rev = "0867bbf92563c57a412335d6385f553f96695e56";
    fetchSubmodules = true;
  };
  isLibrary = true;
  isExecutable = true;
  configureFlags = [ "-f disable-default-paths" ];
  setupHaskellDepends = [ base cabal-doctest ];
  libraryHaskellDepends = [
    base bytestring bytestring-to-vector containers directory inline-c
    list-t relude resourcet safe-exceptions some template-haskell text
    unliftio-core vector
  ];
  librarySystemDepends = [ hdf5 zlib ];
  executableHaskellDepends = [
    base bytestring list-t relude some text vector
  ];
  testHaskellDepends = [ base directory hspec relude vector ];
  homepage = "https://github.com/twesterhout/hdf5-hs";
  description = "High-level Haskell bindings to HDF5 library";
  license = lib.licenses.bsd3;
}
