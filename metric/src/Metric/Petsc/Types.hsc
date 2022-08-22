-- |
-- Copyright: (c) 2022 Tom Westerhout
-- SPDX-License-Identifier: BSD-3-Clause
-- Maintainer: Tom Westerhout <14264576+twesterhout@users.noreply.github.com>
--
-- See README for more info
module Metric.Petsc.Types where

import Foreign.C.Types
import Foreign.Ptr

#include <petscsystypes.h>

newtype RawPetscVec = RawPetscVec {unRawPetscVec :: Ptr ()}

newtype RawPetscMat = RawPetscMat {unRawPetscMat :: Ptr ()}

newtype SNES = SNES { unSNES :: Ptr () }

newtype Vec = Vec { unVec :: Ptr () }

newtype Mat = Mat { unMat :: Ptr () }

type PetscFunction = SNES -> Vec -> Vec -> Ptr () -> IO PetscErrorCode

newtype PetscErrorCode = PetscErrorCode #{type PetscErrorCode}
  deriving (Show, Eq)

type PetscScalar = #{type PetscScalar}

type PetscInt = #{type PetscInt}
