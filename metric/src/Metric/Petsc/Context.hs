{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TemplateHaskell #-}

-- |
-- Copyright: (c) 2022 Tom Westerhout
-- SPDX-License-Identifier: BSD-3-Clause
-- Maintainer: Tom Westerhout <14264576+twesterhout@users.noreply.github.com>
module Metric.Petsc.Context
  ( petscCtx,
  )
where

import qualified Data.Map as Map
import Language.C.Inline.Context (Context (..))
import qualified Language.C.Types as Types
import qualified Language.Haskell.TH as TH
import Metric.Petsc.Types

petscTypesTable :: Map.Map Types.TypeSpecifier TH.TypeQ
petscTypesTable =
  Map.fromList
    [ (Types.TypeName "Vec", [t|RawPetscVec|]),
      (Types.TypeName "Mat", [t|RawPetscMat|]),
      (Types.TypeName "PetscErrorCode", [t|PetscErrorCode|]),
      (Types.TypeName "PetscScalar", [t|PetscScalar|])
    ]

-- | Provides type mappings for better interoperability with "Language.C.Inline".
petscCtx :: Context
petscCtx = mempty {ctxTypesTable = petscTypesTable}
