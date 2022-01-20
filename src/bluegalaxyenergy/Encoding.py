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

from enum import Enum, auto
import random


def isPermutation(perm):
    assert len(perm) == 256
    s = set()
    for val in perm:
        if val in s or val < 0 or val > 255:
            return False
        else:
            s.add(val)
    return True


class Encoding8:

    def __init__(self, permutation):
        assert isPermutation(permutation)
        self.encoded = list(permutation)
        self._genInverse()

    def _genInverse(self):
        self.plain = [None for i in range(256)]
        for i, n in enumerate(self.encoded):
            self.plain[n] = i

    def getInverseEncoding(self):
        return Encoding8(self.plain)

    def getEncodeTable(self):
        return self.encoded

    def getDecodeTable(self):
        return self.plain

    def __getitem__(self, x):
        return self.encoded[x]

    def encode(self, data):
        return bytes([self.encoded[x] for x in data])

    def decode(self, data):
        return bytes([self.plain[x] for x in data])

    def encodeOne(self, data):
        return self.encoded[data]

    def decodeOne(self, data):
        return self.plain[data]

    # compute the table for  self(other(x))  ( self o other )
    def combine(self, other):
        return Encoding8(self.encode(other.getEncodeTable()))

    def __eq__(self, other):
        return all([x == y for x, y in zip(self.encoded, other.getEncodeTable())])


class Encoding8Random(Encoding8):

    def __init__(self, seed=None):
        if seed is not None:
            state = random.getstate()
            random.seed(seed)

        self.encoded = list(range(256))
        random.shuffle(self.encoded)
        self._genInverse()

        if seed is not None:
            random.setstate(state)


class Encoding8Identity(Encoding8):

    def __init__(self):
        self.encoded = list(range(256))
        self.plain = list(range(256))


class Encoding:

    def __init__(self, encodingList):
        self.encodingList = encodingList
        self.length = len(self.encodingList)

    @classmethod
    def fromTable(cls, tables):
        return cls([Encoding8(e) for e in tables])

    def toTable(self):
        return [x.getEncodeTable() for x in self.encodingList]

    def __getitem__(self, x):
        return self.encodingList[x]

    # compute the table for  self(other(x))  ( self o other )
    def combine(self, other):
        return Encoding([x.combine(y) for x, y in zip(self.encodingList, other.encodingList)])

    def getInverseEncoding(self):
        return Encoding([e.getInverseEncoding() for e in self.encodingList])

    def encode(self, data):
        assert len(data) % self.length == 0
        r = b""
        for i in range(0, len(data), self.length):
            r += bytes([t.encodeOne(x) for t, x in zip(self.encodingList, data[i:i+self.length])])
        return r

    def decode(self, data):
        assert len(data) % self.length == 0
        r = b""
        for i in range(0, len(data), self.length):
            r += bytes([t.decodeOne(x) for t, x in zip(self.encodingList, data[i:i+self.length])])
        return r

    def __eq__(self, other):
        return all([x == y for x, y in zip(self.encodingList, other.encodingList)])


class EncodeType(Enum):
    RANDOM = auto()
    IDENTITY = auto()


class EncodingGenerator(Encoding):

    encClass = {
        EncodeType.RANDOM: Encoding8Random,
        EncodeType.IDENTITY: Encoding8Identity,
    }

    def __init__(self, length, encodeType, seed=None):
        self.length = length
        self.encodeType = encodeType
        self.seed = seed

        self.regenerate()

    def setType(self, encodeType):
        if self.encodeType != encodeType:
            self.encodeType = encodeType
            self.regenerate()

    def regenerate(self):
        c = self.encClass[self.encodeType]
        if self.seed is None or self.encodeType == EncodeType.IDENTITY:
            self.encodingList = [c() for _ in range(self.length)]
        elif type(self.seed) == int:
            self.encodingList = [c(self.seed+i) for i in range(self.length)]
        elif type(self.seed) == bytes:
            self.encodingList = [c(self.seed+bytes(i)) for i in range(self.length)]
        else:
            self.encodingList = [c(self.seed) for _ in range(self.length)]


def test():
    for _ in range(32):
        obj = Encoding8Random()
        enc = obj.getEncodeTable()
        dec = obj.getDecodeTable()
        for i in range(256):
            assert dec[enc[i]] == i, "Encoding8Random encodage error"
    print("[OK] Encoding8Random")

    for _ in range(32):
        obj1 = Encoding8Random()
        obj2 = Encoding8Random()
        obj12 = obj1.combine(obj2)
        assert obj1.encode(obj2.encode(list(range(256)))) == obj12.encode(list(range(256)))
    print("[OK] Encoding8 combine")

    for _ in range(32):
        obj1 = EncodingGenerator(16, EncodeType.RANDOM)
        obj2 = EncodingGenerator(16, EncodeType.RANDOM)
        obj12 = obj1.combine(obj2)
        for _ in range(32):
            data = random.randbytes(16)
            assert obj1.encode(obj2.encode(data)) == obj12.encode(data)
    print("[OK] Encoding combine")


if __name__ == '__main__':
    test()
