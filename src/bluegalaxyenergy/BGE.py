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

from .BGEGenInput import genBGEInput
from .Encoding import Encoding8, Encoding
from .AES import revertKey
from ._core import recoverkey
import json


class BGE:

    def __init__(self, wb):
        self.wb = wb
        self.roundkey = []
        self.encoding = []

    # [param] filename path to a file which saves the output of the query to self.wb
    def run(self, roundList=None, nmultiProcess=None, filename=None):
        self.roundkey = [None for _ in range(self.wb.getRoundNumber())]
        self.encoding = [[None for _ in range(16)] for _ in range(self.wb.getRoundNumber())]

        data = genBGEInput(self.wb, roundList, nmultiProcess)
        if filename is not None:
            with open(filename, 'w') as f:
                f.write(json.dumps(data))
        keyResult, encodingResult = recoverkey(data)

        for roundN, key in keyResult.items():
            self.roundkey[roundN] = bytes(key)
        for roundN, perms in encodingResult.items():
            for byte, perm in enumerate(perms):
                self.encoding[roundN][byte] = Encoding8(perm)

        for roundN in range(len(self.encoding)):
            if all([x is not None for x in self.encoding[roundN]]):
                self.encoding[roundN] = Encoding(self.encoding[roundN])
            else:
                self.encoding[roundN] = None

    def computeKey(self, keylen=None, transposed_rk=False):
        if keylen is None:
            keylen = {10: 16, 12: 24, 14: 32}[self.wb.getRoundNumber()]

        possibleKey = None
        for i in range(len(self.roundkey)):
            j = i
            k = b''
            while len(k) < keylen and j < len(self.roundkey) and self.roundkey[j] is not None:
                if transposed_rk:
                    k += bytes(self.roundkey[j][x] for x in [0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15])
                else:
                    k += self.roundkey[j]
                j += 1
            if len(k) < keylen:
                continue
            tempk = revertKey(k, i)
            if possibleKey is None:
                possibleKey = tempk
            # all roundkeys are not part of the same AES keyscheduling
            elif tempk != possibleKey:
                return None
        return possibleKey
