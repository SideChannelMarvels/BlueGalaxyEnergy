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

from .AES import AES, RoundType
from .Encoding import EncodeType, EncodingGenerator
import random


class AESEncoded(AES):

    AES_KEY_ROUND = {16: 10, 24: 12, 32: 14}

    def __init__(self, key, encodingType=EncodeType.RANDOM, encodingSeed=None,
                 byteSwap=False, roundType=RoundType.RoundEncType0):

        super().__init__(key, self.AES_KEY_ROUND[len(key)], roundType=roundType)

        if encodingSeed is None:
            self.encoding = [EncodingGenerator(16, encodingType)
                             for _ in range(self.AES_KEY_ROUND[len(key)] + 1)]
        elif isinstance(encodingSeed, int):
            self.encoding = [EncodingGenerator(16, encodingType, seed=encodingSeed+i*16)
                             for i in range(self.AES_KEY_ROUND[len(key)] + 1)]
        elif isinstance(encodingSeed, bytes):
            self.encoding = [EncodingGenerator(16, encodingType, seed=encodingSeed+bytes([i]))
                             for i in range(self.AES_KEY_ROUND[len(key)] + 1)]
        else:
            self.encoding = [EncodingGenerator(16, encodingType, seed=encodingSeed)
                             for i in range(self.AES_KEY_ROUND[len(key)] + 1)]

        if not byteSwap:
            self.roundPerm = [list(range(16)) for _ in range(self.r)]
        elif encodingSeed is None:
            self.roundPerm = [random.sample(list(range(16)), 16) for _ in range(self.r)]
        else:
            state = random.getstate()
            random.seed(encodingSeed)
            self.roundPerm = [random.sample(list(range(16)), 16) for _ in range(self.r)]
            random.setstate(state)

        self.reverseRoundPerm = [[p.index(i) for i in range(16)] for p in self.roundPerm]

    def applyPermutation(self, data, roundN):
        return bytes([data[x] for x in self.roundPerm[roundN]])

    def applyReversePermutation(self, data, roundN):
        return bytes([data[x] for x in self.reverseRoundPerm[roundN]])

    def encrypt_round_encode(self, state, roundN):
        if roundN != 0:
            state = self.applyReversePermutation(state, roundN-1)
        state = self.encoding[roundN].decode(state)
        state = self.encrypt_round(state, roundN)
        state = self.encoding[roundN + 1].encode(state)
        if roundN + 1 != self.r:
            state = self.applyPermutation(state, roundN)
        return state

    def encrypt_encode(self, state):
        assert len(state) == 16
        assert len(self.extendedKey) == self.r + 1

        for i in range(self.r):
            state = self.encrypt_round_encode(state, i)
        return state

    def encrypt_encode_fast(self, state):
        assert len(state) == 16
        assert len(self.extendedKey) == self.r + 1

        state = self.encoding[0].decode(state)

        for i in range(self.r):
            state = self.encrypt_round(state, i)

        state = self.encoding[self.r].encode(state)
        return state

    def decrypt_round_encode(self, state, roundN):
        if roundN + 1 != self.r:
            state = self.applyReversePermutation(state, roundN)
        state = self.encoding[roundN + 1].decode(state)
        state = self.decrypt_round(state, roundN)
        state = self.encoding[roundN].encode(state)
        if roundN != 0:
            state = self.applyPermutation(state, roundN-1)
        return state

    def decrypt_encode(self, state):
        assert len(state) == 16
        assert len(self.extendedKey) == self.r + 1

        for i in reversed(list(range(self.r))):
            state = self.decrypt_round_encode(state, i)
        return state

    def decrypt_encode_fast(self, state):
        assert len(state) == 16
        assert len(self.extendedKey) == self.r + 1

        state = self.encoding[self.r].decode(state)

        for i in reversed(list(range(self.r))):
            state = self.decrypt_round(state, i)

        state = self.encoding[0].encode(state)
        return state


def test():
    # test case : https://github.com/ircmaxell/quality-checker/blob/master/tmp/gh_18/PHP-PasswordLib-master/test/Data/Vectors/aes-ecb.test-vectors # noqa:E501
    test_cases = [
        ("2b7e151628aed2a6abf7158809cf4f3c", 10,
         "6bc1bee22e409f96e93d7e117393172a",
         "3ad77bb40d7a3660a89ecaf32466ef97"),
        ("2b7e151628aed2a6abf7158809cf4f3c", 10,
         "ae2d8a571e03ac9c9eb76fac45af8e51",
         "f5d3d58503b9699de785895a96fdbaaf"),
        ("8e73b0f7da0e6452c810f32b809079e562f8ead2522c6b7b", 12,
         "6bc1bee22e409f96e93d7e117393172a",
         "bd334f1d6e45f25ff712a214571fa5cc"),
        ("8e73b0f7da0e6452c810f32b809079e562f8ead2522c6b7b", 12,
         "ae2d8a571e03ac9c9eb76fac45af8e51",
         "974104846d0ad3ad7734ecb3ecee4eef"),
        ("603deb1015ca71be2b73aef0857d77811f352c073b6108d72d9810a30914dff4", 14,
         "6bc1bee22e409f96e93d7e117393172a",
         "f3eed1bdb5d2a03c064b5a7e3db181f8"),
        ("603deb1015ca71be2b73aef0857d77811f352c073b6108d72d9810a30914dff4", 14,
         "ae2d8a571e03ac9c9eb76fac45af8e51",
         "591ccb10d410ed26dc5ba74a31362870"),
    ]

    for key, r, plain, cipher in test_cases:
        print(key, r, plain, cipher)
        c = AESEncoded(bytes.fromhex(key))
        c.encoding[0].setType(EncodeType.IDENTITY)
        c.encoding[-1].setType(EncodeType.IDENTITY)
        assert c.encrypt(bytes.fromhex(plain)) == bytes.fromhex(cipher)
        assert c.decrypt(bytes.fromhex(cipher)) == bytes.fromhex(plain)
        assert c.encrypt_encode(bytes.fromhex(plain)) == bytes.fromhex(cipher)
        assert c.decrypt_encode(bytes.fromhex(cipher)) == bytes.fromhex(plain)

        for _ in range(16):
            c = AESEncoded(bytes.fromhex(key), byteSwap=True)
            data = random.randbytes(16)
            assert c.encrypt_encode(c.decrypt_encode(data)) == data
            assert c.decrypt_encode(c.encrypt_encode(data)) == data

            assert c.encrypt_encode(data) == c.encrypt_encode_fast(data)
            assert c.decrypt_encode(data) == c.decrypt_encode_fast(data)

            for roundN in range(c.r):
                assert c.encrypt_round_encode(c.decrypt_round_encode(data, roundN), roundN) == data
                assert c.decrypt_round_encode(c.encrypt_round_encode(data, roundN), roundN) == data


if __name__ == "__main__":
    test()
