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

from .AES import _AesShiftRow, _AesInvShiftRow
from .Exceptions import InvalidRound


def computeImpactByte(ref, fault):
    return tuple([index for index, (x, y) in enumerate(zip(ref, fault)) if x != y])


def getRoundColumn(wb, roundN):
    inColumn = []
    outColumn = []

    ref = wb.applyRound(bytes([0] * 16), roundN)

    colOrder = _AesShiftRow if wb.isEncrypt() else _AesInvShiftRow
    for pos in colOrder:

        data = bytes([0 if p != pos else 1 for p in range(16)])

        impact = computeImpactByte(ref, wb.applyRound(data, roundN))
        if len(impact) != 4:
            raise InvalidRound(f"Wrong number of impacts for round {roundN} byte {pos}:"
                               f" got {len(impact)} impacts but expected 4")

        try:
            index = outColumn.index(impact)
            inColumn[index] = inColumn[index] + (pos, )
            if len(inColumn) > 4:
                raise InvalidRound("Too many positions create an impact on the same column")
        except ValueError:
            if len(outColumn) >= 4:
                raise InvalidRound("Too many impacts possible")
            outColumn.append(impact)
            inColumn.append((pos, ))

    return inColumn, outColumn


def isValidRound(wb, roundN):
    try:
        getRoundColumn(wb, roundN)
        return True
    except InvalidRound:
        return False


def verifyRoundPermutation(wb, roundN):

    inColumn, outColumn = getRoundColumn(wb, roundN)

    res = (outColumn == [(0, 1, 2, 3), (4, 5, 6, 7),
                         (8, 9, 10, 11), (12, 13, 14, 15)])
    if wb.isEncrypt():
        res &= (inColumn == [(0, 5, 10, 15), (4, 9, 14, 3),
                             (8, 13, 2, 7), (12, 1, 6, 11)])
    else:
        res &= (inColumn == [(0, 13, 10, 7), (4, 1, 14, 11),
                             (8, 5, 2, 15), (12, 9, 6, 3)])

    return res


def verifyRoundsPermutation(wb, rounds):

    return all([verifyRoundPermutation(wb, r) for r in rounds])
