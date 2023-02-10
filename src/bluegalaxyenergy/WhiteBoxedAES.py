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

# This class is the interface with the whitebox
# A child class must be implement for a whitebox with the same interface

class WhiteBoxedAES:

    def getRoundNumber(self):
        # return the number of rounds of the whitebox (10 for AES128,
        #   12 for AES192 and 14 for AES256)
        raise NotImplementedError("WhiteBoxedAES.getRoundNumber must be implemented for a given whitebox")

    def applyRound(self, data, roundN):
        # Apply a round of the whitebox on a buffer
        # [param] data    a buffer of 16 bytes (type bytes)
        # [param] roundN  the round number to apply (int in the range [0, self.getRoundNumber()) )
        # return  16 bytes of the encrypted data by the round
        raise NotImplementedError("WhiteBoxedAES.applyRound must be implemented for a given whitebox")
