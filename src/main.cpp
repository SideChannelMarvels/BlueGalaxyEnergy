//-----------------------------------------------------------------------------
// Copyright (C) Quarkslab. See README.md for details.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the Apache License as published by
// the Apache Software Foundation, either version 2.0 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See LICENSE.txt for the text of the Apache license.
//-----------------------------------------------------------------------------

#include "inputReader.hpp"
#include "affineencoding.hpp"
#include "baseencoding.hpp"
#include "approximateencoding.hpp"
#include "encodingkey.hpp"
#include "computeQPtilde.hpp"
#include <string>
#include <map>
#include <utility>
#include <iomanip>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using namespace pybind11::literals;

std::array<std::array<unsigned int, 4>, 4> InvShiftRow = { {
        {0, 1, 2, 3}, {3, 0, 1, 2}, {2, 3, 0, 1}, {1, 2, 3, 0}
    }
};

std::array<std::array<unsigned int, 4>, 4> ShiftRow = { {
        {0, 1, 2, 3}, {1, 2, 3, 0}, {2, 3, 0, 1}, {3, 0, 1, 2}
    }
};

std::pair<std::map<uint8_t, std::array<uint8_t, 16>>, std::map<uint8_t, std::array<std::array<uint8_t, 256>, 16>>>
recoverkey(const std::map < uint8_t, std::array < std::array < std::array<uint8_t, LEN_ARRAY>, LEN_ARRAY + 2 >, 16 >> & data) {
    bool printEncoding = true;
    InputReader inputReader {};

    inputReader.readData(data);

    std::array<bool, AES_ROUND> isAvailable {};
    std::array<std::array<std::array<ApproximateEncoding, 4>, 4>, AES_ROUND> Qtilde;

    // Compute the Qtilde if possible, that is Qtilde(x) = Q(A(x)), where A is an unknown affine encoding
    for (int round = 0; round < AES_ROUND; round++) {
        if (inputReader.isRoundAvailable(round)) {
            isAvailable[round] = true;
            for (size_t col = 0; col < 4; col++) {
                for (size_t row = 0; row < 4; row ++) {
                    if (not Qtilde[round][col][row].build(inputReader.getx0_all(round, col, row))) {
                        std::cerr << "Fail to build Qtilde for round " << round + 1 << " byte " << (col * 4 + row) << std::endl;
                        isAvailable[round] = false;
                    }
                }
            }
        } else {
            isAvailable[round] = false;
        }
    }

    std::array < std::array<std::array<EncodingKey, 4>, 4>, AES_ROUND - 1 > Ptilde;
    std::array < std::array<std::array<AffineEncoding, 4>, 4>, AES_ROUND - 1 > Q;

    for (int round = 0; round < AES_ROUND - 1; round++) {
        if (not(isAvailable[round] and isAvailable[round + 1])) {
            continue;
        }
        for (unsigned int col = 0; col < 4; col++) {
            std::array<std::array<BaseEncoding, 4>, 4> xy;
            for (unsigned int x = 0; x < 4; x++) {
                for (unsigned int y = 0; y < 4; y++) {
                    xy[x][y] = Qtilde[round + 1][col][y].inverse().compose(
                                   inputReader.getx(x, round + 1, col, y).compose(
                                       Qtilde[round][ShiftRow[x][col]][x]));
                }
            }
            if (not computeQPtilde(xy, Q[round][col], Ptilde[round][col])) {
                std::cerr << "Fail to compute Q and Ptilde for round=" << round << ", col=" << col << std::endl;
                isAvailable[round] = false;
                break;
            }
        }
    }

    std::array < std::array<std::array<uint8_t, 4>, 4>, AES_ROUND - 2 > rkey;
    std::array < std::array<std::array<BaseEncoding, 4>, 4>, AES_ROUND - 2 > internalEncoding;
    std::array<bool, AES_ROUND> keyAvailable {};

    for (int round = 0; round < AES_ROUND - 2; round++) {
        if (not(isAvailable[round] and isAvailable[round + 1] and isAvailable[round + 2])) {
            keyAvailable[round] = false;
            continue;
        }
        keyAvailable[round] = true;
        for (unsigned int col = 0; col < 4; col++) {
            for (unsigned int row = 0; row < 4; row++) {
                BaseEncoding e = Ptilde[round + 1][InvShiftRow[row][col]][row].compose(Q[round][col][row]);
                for (unsigned int x = 1; x < LEN_ARRAY; x++) {
                    if (e[0] != (e[x] ^ x)) {
                        std::cerr << "Fail to verify xor property for roundkey " << (round + 2) << " byte " << (col * 4 + row) << std::endl;
                        keyAvailable[round] = false;
                        break;
                    }
                }
                rkey[round][col][row] = e[0];
                if (round != 0 && isAvailable[round - 1])
                    internalEncoding[round][col][row] = Qtilde[round][col][row].compose(Q[round - 1][col][row]);
            }
        }
    }

// fixme return better struct

    std::map<uint8_t, std::array<uint8_t, 16>> keyMap;

    for (int round = 0; round < AES_ROUND - 2; round++) {
        if (not keyAvailable[round]) {
            continue;
        }
        std::array<uint8_t, 16> roundKey;
        for (unsigned int col = 0; col < 4; col++) {
            for (unsigned int row = 0; row < 4; row++) {
                roundKey[col * 4 + row] = static_cast<uint8_t>(rkey[round][col][row]);
            }
        }
        keyMap[round + 2] = std::move(roundKey);
    }

    std::map<uint8_t, std::array<std::array<uint8_t, 256>, 16>> encodingMap;

    if (printEncoding) {
        for (int round = 0; round < AES_ROUND - 2; round++) {
            if (round == 0 or not isAvailable[round - 1] or not keyAvailable[round]) {
                continue;
            }
            std::array<std::array<uint8_t, 256>, 16> roundEncoding;
            for (unsigned int col = 0; col < 4; col++) {
                for (unsigned int row = 0; row < 4; row++) {
                    roundEncoding[col * 4 + row] = internalEncoding[round][col][row].getEncodingArray();
                }
            }
            encodingMap[round + 1] = std::move(roundEncoding);
        }
    }
    return std::make_pair(std::move(keyMap), std::move(encodingMap));
}

PYBIND11_MODULE(_core, m) {
    m.doc() = "Blue Galaxy Energy module";
    m.def("recoverkey", &recoverkey, "Core of the BGE key recovery", "filename"_a);

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
