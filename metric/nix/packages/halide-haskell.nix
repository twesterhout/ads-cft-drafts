{ mkDerivation, base, fetchgit, lib }:
mkDerivation {
  pname = "halide-haskell";
  version = "0.0.1.0";
  src = fetchgit {
    url = "https://github.com/twesterhout/halide-haskell";
    hash = "sha256-f9T0X2F1iQc0ETrHlB8eVSnY6FboNEPqt/3K104GtBw=";
    rev = "b2383c6fcaaacce7d77b7a418f584c1c4eda1c67";
    fetchSubmodules = true;
  };
  libraryHaskellDepends = [ base ];
  homepage = "https://github.com/twesterhout/halide-haskell";
  description = "See README for more info";
  license = lib.licenses.bsd3;
}
