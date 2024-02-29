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

import multiprocessing as mp
from .ByteOrder import verifyRoundsPermutation
import os


def detectCPU():
    try:
        return len(os.sched_getaffinity(0))
    except:
        pass
    nCPU = os.cpu_count()
    if nCPU is not None:
        return nCPU
    else:
        return 0


def encodeResult(x0, x2, x3):

    return [
            [[o[b] for o in x0[y]] for y in range(256)] +
            [[o[b] for o in x2], [o[b] for o in x3]]
            for b in range(16)]


def genBGEInputRound(wb, roundN):
    x0 = []
    x2 = []
    x3 = []

    # type r<roundN+1>:y<c>(x,y,0,0,x,y,0,0,x,y,0,0,x,y,0,0)=
    for y in range(256):
        out = []
        for x in range(256):
            data = bytes([x, y, 0, 0, x, y, 0, 0, x, y, 0, 0, x, y, 0, 0])
            out.append(wb.applyRound(data, roundN))

        assert len(out) == 256
        assert all([len(o) == 16 for o in out])
        x0.append(out)

    assert len(x0) == 256

    # type r<roundN+1>:y<c>(0,0,x,0,0,0,x,0,0,0,x,0,0,0,x,0)=
    for x in range(256):
        data = bytes([0, 0, x, 0, 0, 0, x, 0, 0, 0, x, 0, 0, 0, x, 0])
        x2.append(wb.applyRound(data, roundN))

    assert len(x2) == 256
    assert all([len(o) == 16 for o in x2])

    # type r<roundN+1>:y<c>(0,0,0,x,0,0,0,x,0,0,0,x,0,0,0,x)=
    for x in range(256):
        data = bytes([0, 0, 0, x, 0, 0, 0, x, 0, 0, 0, x, 0, 0, 0, x])
        x3.append(wb.applyRound(data, roundN))

    assert len(x3) == 256
    assert all([len(o) == 16 for o in x3])

    # generate x0, x2 and x3
    # note: we export the result as roundN+1 (needed by the tools)
    return {(roundN+1): encodeResult(x0, x2, x3)}


def init_localWB(wb):
    global localWB
    wb.newThread()
    localWB = wb


def genBGEInputRoundProxy(roundN):
    global localWB
    return genBGEInputRound(localWB, roundN)


# [param] wb            instance of WhiteBoxedAES
# [param] roundList     list of rounds in [0, wb.getRoundNumber() - 1)
# [param] nmultiProcess nprocess in the poll (0 to not use multiprocessing)
def genBGEInput(wb, roundList, nmultiProcess=None):
    if nmultiProcess is None:
        nmultiProcess = detectCPU()

    assert verifyRoundsPermutation(wb, roundList), "Detected impact on the wrong column, use 'BGE(..., shuffle=True)'"

    result = {}

    if nmultiProcess == 0:
        for rnd in roundList:
            result.update(genBGEInputRound(wb, rnd))
    else:
        assert nmultiProcess >= 1
        with mp.Pool(processes=nmultiProcess, initializer=init_localWB, initargs=(wb, )) as pool:
            res = []
            for rnd in roundList:
                res.append(pool.apply_async(genBGEInputRoundProxy, (rnd, )))
            for r in res:
                result.update(r.get())
    return result

