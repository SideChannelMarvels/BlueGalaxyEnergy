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

from .AES import revertKey, MC
from .BGEGenInput import genBGEInput
from .ByteOrder import isValidRound
from .Encoding import Encoding8, Encoding, EncodingXor
from .Exceptions import SaveLoadError
from .FinalizeUnshuffled import FinalizeUnshuffled
from .WhiteBoxedUnshuffledProxy import WhiteBoxedUnshuffledProxy
from enum import Enum
import json

try:
    from ._core import recoverkey
except ImportError:
    print("[!] fail to import native library to resolve equation")
    recoverkey = None


class BGEState(Enum):
    NoInput = 1
    InputComputed = 2
    ResolveDone = 3


class BGE:

    EXPORT_VERSION = 1
    EXPORT_TYPE = "BGEInputs"

    def __init__(self, wb=None, isEncrypt=None, numRound=None):
        self.wb = wb
        self.state = BGEState.NoInput
        if self.wb is not None:
            assert isEncrypt is None or isEncrypt == self.wb.isEncrypt()
            assert numRound is None or numRound == self.wb.getRoundNumber()
            self.isEncrypt = self.wb.isEncrypt()
            self.numRound = self.wb.getRoundNumber()
        else:
            self.isEncrypt = isEncrypt
            self.numRound = numRound
        assert type(self.isEncrypt) is bool
        assert type(self.numRound) is int
        assert self.numRound in [10, 12, 14]

        self.roundList = None
        self.paramShuffle = None

    @classmethod
    def loadFrom(cls, filename):

        # when parsing dict from json, try to restore the type of the key if the
        # key is an integer-like string
        def hook(d):
            return {
                    int(k) if k.lstrip('-').isdigit() else k: v
                    for k, v in d.items()
                }

        with open(filename, 'r') as f:
            data = json.loads(f.read(), object_hook=hook)

        if (type(data) is not dict or
            data.get("DataType", None) != cls.EXPORT_TYPE or
            data.get("Version", None) != cls.EXPORT_VERSION):

            raise SaveLoadError("Unexpected file structure")

        attack = cls(isEncrypt=data.get("isEncrypt", None),
                     numRound=data.get("numRound", None))

        attack.roundList = data["roundList"]
        attack.shuffle = data["shuffle"]
        attack.BGEInputs = data["BGEInputs"]
        attack.shuffle = data["shuffle"]
        attack.paramShuffle = data["paramShuffle"]
        attack.state = BGEState.InputComputed

        return attack

    def saveTo(self, filename, compress=True):
        if self.state not in [BGEState.InputComputed, BGEState.ResolveDone]:
            raise SaveLoadError("Cannot save without data")

        data = {}
        data["DataType"] = self.EXPORT_TYPE
        data["Version"] = self.EXPORT_VERSION
        data["isEncrypt"] = self.isEncrypt
        data["numRound"] = self.numRound
        data["roundList"] = self.roundList
        data["BGEInputs"] = self.BGEInputs
        data["shuffle"] = self.shuffle
        data["paramShuffle"] = self.paramShuffle

        with open(filename, 'w') as f:
            f.write(json.dumps(data))

    def run(self, roundList=None, shuffle=False, nmultiProcess=None):
        self.generateInput(roundList=roundList, shuffle=shuffle,
                           nmultiProcess=nmultiProcess)
        self.resolve()

    def generateInput(self, roundList=None, shuffle=False, nmultiProcess=None):
        if self.state != BGEState.NoInput:
            return

        assert self.wb is not None, "Cannot prepare input without the whitebox implementation"
        self.shuffle = shuffle

        if roundList is None:
            roundList = []
            for roundN in list(range(self.numRound-1)):
                if isValidRound(self.wb, roundN):
                    roundList.append(roundN)
                else:
                    print(f"Disable round {roundN} : not a valid round for the attack")

        self.roundList = sorted(list(set(roundList)))

        # if roundList is given, the list must be consecutive rounds
        assert self.roundList == list(range(self.roundList[0],
                                            self.roundList[0] +
                                            len(self.roundList))), \
            "RoundList is not made of consecutive rounds"
        assert all([x >= 0 for x in self.roundList]), "Invalid round number"
        assert all([x < self.numRound-1 for x in self.roundList]), "Invalid round number"

        if self.shuffle:
            proxyWb = WhiteBoxedUnshuffledProxy(self.wb, self.roundList)
            proxyWb.computeSimplePermutation()
            self.BGEInputs = genBGEInput(proxyWb, self.roundList, nmultiProcess)
            self.paramShuffle = proxyWb.exctractPermutation()
        else:
            self.BGEInputs = genBGEInput(self.wb, self.roundList, nmultiProcess)
            self.paramShuffle = [None, None,
                                 [list(range(16)) for _ in range(self.numRound)]]

        self.state = BGEState.InputComputed

    def resolve(self):
        if self.state == BGEState.ResolveDone:
            return
        self.shuffleProxy = FinalizeUnshuffled(self.isEncrypt, *self.paramShuffle)
        assert self.state != BGEState.NoInput, "Cannot resolve before generating inputs"
        assert recoverkey is not None, "Missing resolve library"

        self.roundkey = [None for _ in range(self.numRound+1)]
        self.encoding = [[None for _ in range(16)] for _ in range(self.numRound+1)]

        keyResult, encodingResult, columnCoeff = recoverkey(self.BGEInputs, self.isEncrypt)

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

        if self.shuffle:
            self._reverseWbProxy()
            self.shuffleProxy.reorderWithCoeff(columnCoeff)
            self.shuffleProxy.setFirstByte(0)
            self._applyWbProxy()

        self.state = BGEState.ResolveDone

    def _reverseWbProxy(self):
        if not self.shuffle:
            return
        for roundN, key in enumerate(self.roundkey):
            if key is not None:
                self.roundkey[roundN] = self.shuffleProxy.applyReversePermutation(key, roundN-1)

        for roundN, perm in enumerate(self.encoding):
            if perm is not None:
                self.encoding[roundN] = Encoding(self.shuffleProxy.applyReversePermutation(perm.encodingList, roundN-1, False))

    def _applyWbProxy(self):
        if not self.shuffle:
            return

        for roundN, key in enumerate(self.roundkey):
            if key is not None:
                self.roundkey[roundN] = self.shuffleProxy.applyPermutation(key, roundN-1)

        for roundN, perm in enumerate(self.encoding):
            if perm is not None:
                self.encoding[roundN] = Encoding(self.shuffleProxy.applyPermutation(perm.encodingList, roundN-1, False))

    def _changeWbProxyFirstByte(self, n):
        if not self.shuffle:
            return

        self._reverseWbProxy()
        self.shuffleProxy.setFirstByte(n)
        self._applyWbProxy()

    def getEncoding(self):
        assert self.state == BGEState.ResolveDone
        if self.isEncrypt:
            return list(self.encoding)
        else:
            res = []
            rkey = list(reversed(self.roundkey))
            for ident, enc in enumerate(reversed(self.encoding)):
                if enc is not None and rkey[ident] is not None:
                    res.append(enc.combine(EncodingXor(rkey[ident])))
                else:
                    res.append(None)
            return res

    def getShuffle(self):
        assert self.state == BGEState.ResolveDone
        res = {}
        for r in [self.roundList[0]] + self.roundList:
            if self.isEncrypt:
                res[r] = self.shuffleProxy.reverseRoundPerm[r]
            else:
                res[self.numRound-(2+r)] = self.shuffleProxy.reverseRoundPerm[r]
        return res

    def computeKey(self, keyLen=None, transposed_rk=False, offset=0):
        assert self.state == BGEState.ResolveDone
        if not self.shuffle:
            return self._computeKey(keyLen=keyLen, transposed_rk=transposed_rk, offset=offset)
        else:
            possibleKey = {}
            for i in range(16):
                self._changeWbProxyFirstByte(i)
                key = self._computeKey(keyLen=keyLen, transposed_rk=transposed_rk, offset=offset)
                if key is not None:
                    possibleKey[i] = key
            if len(possibleKey.keys()) == 0:
                return None
            targetByte = list(possibleKey.keys())[0]
            self._changeWbProxyFirstByte(targetByte)
            if len(possibleKey.keys()) == 1:
                return possibleKey[targetByte]
            return possibleKey

    def _computeKey(self, keyLen=None, transposed_rk=False, offset=0):
        assert self.state == BGEState.ResolveDone
        assert self.numRound is not None
        if keyLen is None:
            keyLen = {10: 16, 12: 24, 14: 32}[self.numRound]

        if self.isEncrypt:
            localRoundKey = list(self.roundkey)
        else:
            localRoundKey = list(reversed(self.roundkey))

        possibleKey = None
        for i in range(len(localRoundKey)):
            j = i
            k = b''
            while len(k) < keyLen and j < len(localRoundKey) and localRoundKey[j] is not None:
                if transposed_rk:
                    tmpk = bytes(localRoundKey[j][x] for x in [0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15])
                else:
                    tmpk = localRoundKey[j]
                if self.isEncrypt:
                    k += tmpk
                else:
                    k += MC(tmpk)
                j += 1
            if len(k) < keyLen:
                continue
            tempk = revertKey(k, i + offset)
            if possibleKey is None:
                possibleKey = tempk
            # all roundkeys are not part of the same AES keyscheduling
            elif tempk != possibleKey:
                return None
        return possibleKey
