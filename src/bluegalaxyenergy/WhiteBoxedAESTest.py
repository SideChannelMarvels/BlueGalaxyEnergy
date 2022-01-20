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


class WhiteBoxedAESTest(WhiteBoxedAES):

    def __init__(self, aesEncoded):
        self.aesEncoded = aesEncoded

    def getRoundNumber(self):
        return self.aesEncoded.r

    def applyRound(self, data, roundN):
        return self.aesEncoded.encrypt_round_encode(data, roundN)
