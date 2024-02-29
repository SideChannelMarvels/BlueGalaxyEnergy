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

from .BGE import BGE
from .AESEncoded import AESEncoded
from .AES import RoundType
from .WhiteBoxedAESTest import WhiteBoxedAESTest
import argparse
import secrets
import time


def test(isEncrypt, isShuffle):
    key_len = 32
    key = secrets.token_bytes(key_len)
    seed = secrets.randbits(64)
    # Note : roundType is chosen here to check the permutation and
    #        encoding found by the attack. With a real whitebox, you cannot
    #        choose which encoding state you can provide. The tools supports any
    #        encoding intermediary state between two MixColumns. The functions
    #        getShuffle() and getEncoding() will return the permutation and
    #        encoding to transform the encoding position provided to the clear
    #        intermediary state just after a MixColumns (or InvMixColumns).
    roundType = RoundType.RoundEncType0 if isEncrypt else RoundType.RoundDecType0
    aesEncoded = AESEncoded(key, roundType=roundType, encodingSeed=seed,
                            byteSwap=isShuffle)

    print("=== {} {} ===".format(
        "Shuffle" if isShuffle else "Unshuffle",
        "Encrypt" if isEncrypt else "Decrypt"))
    print(f"Key = {key.hex()}")
    print(f"Seed = {seed}")

    start = time.time()
    bge = BGE(WhiteBoxedAESTest(aesEncoded, isEncrypt))

    if False:
        bge.generateInput(shuffle=isShuffle)
        bge.saveTo("/tmp/bge_test.json")

        bge = BGE.loadFrom("/tmp/bge_test.json")
        bge.resolve()
    else:
        bge.run(shuffle=isShuffle)

    assert key == bge.computeKey()
    stop = time.time()

    print(f"Duration : {stop - start}")
    print("Key:         OK")
    for i, enc in enumerate(bge.getEncoding()):
        if enc is not None:
            assert enc == aesEncoded.encoding[i]
    # Did we recover all expected encodings?
    if all([bge.roundkey[i] is not None for i, x in enumerate(bge.roundkey) if x]):
        print("Encodings:   OK")
    for i, perm in bge.getShuffle().items():
        assert perm == aesEncoded.roundPerm[i]
    print("Permutation: OK")


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--selftest', action='store_true', help="Perform all self tests")
    parser.add_argument('--enctest', action='store_true', help="Perform enc self tests")
    parser.add_argument('--dectest', action='store_true', help="Perform dec self tests")
    parser.add_argument('--encShuffletest', action='store_true', help="Perform enc self tests")
    parser.add_argument('--decShuffletest', action='store_true', help="Perform enc self tests")

    args = parser.parse_args()

    # from .AESEncoded import test as testAESEncoded
    # testAESEncoded()

    if args.selftest or args.enctest:
        test(True, False)

    if args.selftest or args.dectest:
        test(False, False)

    if args.selftest or args.encShuffletest:
        test(True, True)

    if args.selftest or args.decShuffletest:
        test(False, True)


if __name__ == "__main__":
    run()
