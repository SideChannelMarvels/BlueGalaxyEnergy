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

from .WhiteBoxedAES import WhiteBoxedAES
from .AES import InvShiftRow, ShiftRow, _AesShiftRow, _AesInvShiftRow
from .ByteOrder import verifyRoundsPermutation, getRoundColumn
import itertools


class WhiteBoxedUnshuffledProxy(WhiteBoxedAES):

    def __init__(self, basewb, targetRound):
        assert isinstance(basewb, WhiteBoxedAES)
        self.basewb = basewb
        self.targetRound = targetRound
        self.roundPerm = []
        self.reverseRoundPerm = []
        self.preColumn = {}
        self.postColumn = {}

        self._isEncrypt = self.basewb.isEncrypt()
        self._roundNumber = self.basewb.getRoundNumber()

        for i in range(self.getRoundNumber()):
            self.roundPerm.append(None)
            self.reverseRoundPerm.append(None)
            self.setRoundPermutation(i, list(range(16)))

    def isEncrypt(self):
        return self._isEncrypt

    def getRoundNumber(self):
        return self._roundNumber

    def newThread(self):
        self.basewb.newThread()

    def applyRound(self, data, roundN):
        assert roundN in self.targetRound
        if roundN != 0:
            data = self.applyReversePermutation(data, roundN - 1, True)
        data = self.basewb.applyRound(data, roundN)
        if roundN + 1 != self.getRoundNumber():
            data = self.applyPermutation(data, roundN, True)
        return data

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

    # First step : organize a whitebox based on round column
    # During this phase, we want to emulate a good propagation of bytes.
    # i.e., for each round in encrypt, any change in input bytes 0, 5, 10 or 15
    # must have an impact on the first column of the output
    # Without this, the method genBGEInput will not generate good data and BGE
    # cannot work
    def computeSimplePermutation(self):
        assert len(self.preColumn.keys()) == 0
        assert len(self.postColumn.keys()) == 0

        for roundN in sorted(self.targetRound):
            self.detectConstraint(roundN)

        assert verifyRoundsPermutation(self, self.targetRound), "Fail to construct valid permutation"

    def detectConstraint(self, roundN):
        self.preColumn[roundN], self.postColumn[roundN] = getRoundColumn(self.basewb, roundN)

        self.setRoundPermutation(roundN, list(itertools.chain.from_iterable(self.postColumn[roundN])))

        if roundN == 0:
            return
        elif (roundN - 1) not in self.targetRound:
            perm1 = list(itertools.chain.from_iterable(self.preColumn[roundN]))
        else:
            prevPermut = self.postColumn[roundN - 1]
            curPermut = self.preColumn[roundN]
            possibleByte = [set(prevPermut[x >> 2]).intersection(
                            set(curPermut[pos >> 2]))
                            for pos, x in enumerate(_AesShiftRow if
                                                    self.isEncrypt() else
                                                    _AesInvShiftRow)]

            assert all([len(x) == 1 for x in possibleByte])
            perm1 = list(itertools.chain.from_iterable(possibleByte))
            assert len(set(perm1)) == 16

            assert all([set(self.preColumn[roundN][i]) == set(perm1[i*4:(i+1)*4]) for i in range(4)])
            self.preColumn[roundN] = [tuple(perm1[i*4:(i+1)*4]) for i in range(4)]

        if self.isEncrypt():
            perm2 = list(InvShiftRow(perm1))
        else:
            perm2 = list(ShiftRow(perm1))

        self.setRoundPermutation(roundN - 1, perm2)

        if (roundN - 1) in self.targetRound:
            assert all([set(self.postColumn[roundN - 1][i]) == set(perm2[i*4:(i+1)*4]) for i in range(4)])
            self.postColumn[roundN - 1] = [tuple(perm2[i*4:(i+1)*4]) for i in range(4)]

    def exctractPermutation(self):
        return [self.preColumn, self.postColumn, self.roundPerm]

