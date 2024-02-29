#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (C) Quarkslab. See README.md for details.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Apache License as published by
# the Apache Software Foundation, either version 2.0 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See LICENSE.txt for the text of the Apache license.
# -----------------------------------------------------------------------------

from .AES import InvShiftRow, ShiftRow, _AesShiftRow, _AesInvShiftRow
import itertools


class FinalizeUnshuffled:

    def __init__(self, isEncrypt, preColumn, postColumn, roundPerm):
        self._isEncrypt = isEncrypt
        self.coeffRound = []
        self.roundPerm = []
        self.reverseRoundPerm = []
        self.preColumn = preColumn
        self.postColumn = postColumn

        for i, perm in enumerate(roundPerm):
            self.roundPerm.append(None)
            self.reverseRoundPerm.append(None)
            if perm is None:
                self.setRoundPermutation(i, list(range(16)))
            else:
                self.setRoundPermutation(i, perm)

    def isEncrypt(self):
        return self._isEncrypt

    def applyPermutation(self, data, roundN, isBytes=True):
        res = [data[x] for x in self.roundPerm[roundN]]
        return bytes(res) if isBytes else res

    def applyReversePermutation(self, data, roundN, isBytes=True):
        res = [data[x] for x in self.reverseRoundPerm[roundN]]
        return bytes(res) if isBytes else res

    def setRoundPermutation(self, roundN, perm):
        assert len(set(perm)) == 16
        assert set(perm) == set(range(16))
        self.roundPerm[roundN] = perm
        self.reverseRoundPerm[roundN] = [perm.index(i) for i in range(16)]

    # Second step : After performing BGE, we know the coefficient of the
    # MixColumns for every round columns. During this step, we change the order
    # of each column in order to match genuine coefficient.
    def reorderWithCoeff(self, columnCoeff):
        assert len(self.coeffRound) == 0
        assert len(self.preColumn.keys()) != 0
        assert len(self.postColumn.keys()) != 0

        for roundN, roundCoeffs in columnCoeff.items():
            for colN, coeffs in enumerate(roundCoeffs):
                self.applyCoeff(roundN, colN, coeffs)
            self.coeffRound.append(roundN)

        self.coeffRound.sort()
        assert self.coeffRound == list(range(self.coeffRound[0],
                                             self.coeffRound[0] +
                                             len(self.coeffRound)))

    def applyCoeff(self, roundN, colN, coeffs):
        v0 = 3 if self.isEncrypt() else 11
        v1 = 2 if self.isEncrypt() else 14

        inList = []
        outList = []
        p = 3
        for _ in range(4):
            s = None
            for s_, t in enumerate(coeffs):
                if t[p] == v0:
                    s = s_
                    p = t.index(v1)
                    break
            assert s is not None
            assert s not in inList
            inList.append(s)
            outList.append(p)
        assert p == 3

        assert roundN in self.preColumn
        assert roundN in self.postColumn

        preCol = self.preColumn[roundN][colN]
        postCol = self.postColumn[roundN][colN]

        self.preColumn[roundN][colN] = tuple([preCol[x] for x in inList])
        self.postColumn[roundN][colN] = tuple([postCol[x] for x in outList])
        # coeff = [[coeffs[x][y] for y in outList] for x in inList]

    # Third step : After the second step, we only have 16 possibilities.
    # We only need to choose which byte should be put at the first one.
    # This function reorganizes the proxy in order to provide the only
    # possibility with the byte N as the first input of the first supported round
    def setFirstByte(self, N):
        assert len(self.coeffRound) >= 2
        assert self.coeffRound == list(range(self.coeffRound[0],
                                             self.coeffRound[0] +
                                             len(self.coeffRound)))

        firstCoeff = self.coeffRound[0]
        self._setByteAt(firstCoeff, N, 0)

        for roundN in range(firstCoeff + 1,
                            firstCoeff + len(self.coeffRound)):

            for i in range(4):
                dest = _AesInvShiftRow[i] if self.isEncrypt() else _AesShiftRow[i]
                self._setByteAt(roundN, self.postColumn[roundN-1][0][i], dest)

        for i in range(1, 4):
            dest = _AesShiftRow[i] if self.isEncrypt() else _AesInvShiftRow[i]
            self._setByteAt(firstCoeff, self.preColumn[firstCoeff+1][0][i], dest, pre=False)

        if firstCoeff != 0:
            perm1 = list(itertools.chain.from_iterable(self.preColumn[firstCoeff]))

            if self.isEncrypt():
                perm1 = list(InvShiftRow(perm1))
            else:
                perm1 = list(ShiftRow(perm1))

            self.setRoundPermutation(firstCoeff - 1, perm1)

        for roundN in self.coeffRound:
            self.setRoundPermutation(roundN, list(itertools.chain.from_iterable(self.postColumn[roundN])))

    def _setByteAt(self, roundN, byteValue, target, pre=True):
        assert 0 <= target and target < 16
        assert 0 <= byteValue and byteValue < 16

        targetCol = (target >> 2)
        targetRow = (target & 3)

        currentCol = None
        for col, t in enumerate(self.preColumn[roundN] if pre else self.postColumn[roundN]):
            try:
                currentRow = t.index(byteValue)
                currentCol = col
                break
            except ValueError:
                continue

        assert currentCol is not None

        if currentCol != targetCol:
            self.preColumn[roundN][currentCol], self.preColumn[roundN][targetCol] = \
                self.preColumn[roundN][targetCol], self.preColumn[roundN][currentCol]
            self.postColumn[roundN][currentCol], self.postColumn[roundN][targetCol] = \
                self.postColumn[roundN][targetCol], self.postColumn[roundN][currentCol]

        if targetRow != currentRow:
            shift = ((currentRow + 4) - targetRow) % 4
            self.preColumn[roundN][targetCol] = \
                self.preColumn[roundN][targetCol][shift:] + self.preColumn[roundN][targetCol][:shift]
            self.postColumn[roundN][targetCol] = \
                self.postColumn[roundN][targetCol][shift:] + self.postColumn[roundN][targetCol][:shift]
