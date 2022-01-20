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

import sys
from .BGE import BGE
from .AESEncoded import AESEncoded
from .WhiteBoxedAESTest import WhiteBoxedAESTest


def test():
    key_len = 32
    key = bytes(list(range(key_len)))
    aesEncoded = AESEncoded(key)

    bge = BGE(WhiteBoxedAESTest(aesEncoded))
    bge.run()

    assert key == bge.computeKey()
    print("Key:       OK")
    for i in range(len(bge.encoding)):
        if bge.encoding[i] is not None:
            assert bge.encoding[i] == aesEncoded.encoding[i]
    # Did we recover all expected encodings?
    if all([bge.encoding[i - 1] is not None for i, x in enumerate(bge.roundkey) if x and i >= 3]):
        print("Encodings: OK")


if len(sys.argv) > 1 and sys.argv[1] == '--selftest':
    test()
