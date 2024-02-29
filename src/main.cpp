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
#include <tuple>
#include <utility>
#include <iomanip>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using namespace pybind11::literals;

std::tuple<std::map<uint8_t, std::array<uint8_t, 16>>,
           std::map<uint8_t, std::array<std::array<uint8_t, 256>, 16>>,
           std::map<uint8_t, std::array<std::array<std::array<uint8_t, 4>, 4>, 4>>>
recoverkey(const std::map < uint8_t, std::array < std::array < std::array<uint8_t, LEN_ARRAY>, LEN_ARRAY + 2 >, 16 >> & data,
           bool isEncrypt) {
    InputReader inputReader {};

    inputReader.readData(data);

    std::array<bool, AES_ROUND> isQtildeAvailable {};
    std::array<std::array<std::array<ApproximateEncoding, 4>, 4>, AES_ROUND> Qtilde;

    // Compute the Qtilde if possible, that is Qtilde(x) = Q(A(x)), where A is an unknown affine encoding
    for (int round = 0; round < AES_ROUND; round++) {
        if (inputReader.isRoundAvailable(round)) {
            isQtildeAvailable[round] = true;
            for (size_t col = 0; col < 4; col++) {
                for (size_t row = 0; row < 4; row ++) {
                    if (not Qtilde[round][col][row].build(inputReader.getx0_all(round, col, row))) {
                        std::cerr << "Fail to build Qtilde for round " << round + 1 << " byte " << (col * 4 + row) << std::endl;
                        isQtildeAvailable[round] = false;
                    }
                }
            }
        } else {
            isQtildeAvailable[round] = false;
        }
    }
    std::array < std::array<QPtildeLinear, 4>, AES_ROUND - 1 > QPtilde;

    for (int round = 0; round < AES_ROUND - 1; round++) {
        if (not(isQtildeAvailable[round] and isQtildeAvailable[round + 1])) {
            continue;
        }
        for (unsigned int col = 0; col < 4; col++) {
            std::array<std::array<BaseEncoding, 4>, 4> xy;
            for (unsigned int x = 0; x < 4; x++) {
                unsigned int byteOffset = (isEncrypt) ? (ShiftRow[x][col]) : (InvShiftRow[x][col]);
                for (unsigned int y = 0; y < 4; y++) {
                    xy[x][y] = Qtilde[round + 1][col][y].inverse().compose(
                                   inputReader.getx(x, round + 1, col, y).compose(
                                       Qtilde[round][byteOffset][x]));
                }
            }
            if (not computeQPtilde(xy, QPtilde[round][col], isEncrypt)) {
                std::cerr << "Fail to compute Q and Ptilde for round=" << round << ", col=" << col << std::endl;
                break;
            }
        }
    }

    // in case of encryption, QPtilde is already determined. But in case of decrypt,
    // this isn't the case. We need to perform an additional collision between
    // consecutive QPtilde
    if (!isEncrypt) {
        for (int round = 0; round < AES_ROUND - 2; round++) {
            finalizeQPtildeDecrypt(QPtilde[round], QPtilde[round+1]);
        }
    }

    std::array<bool, AES_ROUND - 1> isQPtildeAvailable {};

    for (int round = 0; round < AES_ROUND - 1; round++) {
        isQPtildeAvailable[round] =
            QPtilde[round][0].isFinalize() &&
            QPtilde[round][1].isFinalize() &&
            QPtilde[round][2].isFinalize() &&
            QPtilde[round][3].isFinalize();
        if ((!isQPtildeAvailable[round]) &&
            QPtilde[round][0].isValid() &&
            QPtilde[round][1].isValid() &&
            QPtilde[round][2].isValid() &&
            QPtilde[round][3].isValid()) {

            std::cerr << "Fail to finalize computation Q and Ptilde for round=" << round << std::endl;
        }
    }


    std::array < std::array<std::array<uint8_t, 4>, 4>, AES_ROUND - 2 > rkey;
    std::array<bool, AES_ROUND> keyAvailable {};

    for (int round = 0; round < AES_ROUND - 2; round++) {
        if (not(isQPtildeAvailable[round] and isQPtildeAvailable[round + 1])) {
            keyAvailable[round] = false;
            continue;
        }
        keyAvailable[round] = true;
        for (unsigned int col = 0; col < 4; col++) {
            for (unsigned int row = 0; row < 4; row++) {
                unsigned int byteOffset = (isEncrypt) ? (InvShiftRow[row][col]) : (ShiftRow[row][col]);
                if (!getAndVerifyKey(QPtilde[round][col].getQ(row),
                                     QPtilde[round + 1][byteOffset].getPtilde(row),
                                     rkey[round][col][row])) {
                    std::cerr << "Fail to verify xor property for roundkey "
                              << (round + 2) << " byte "
                              << (col * 4 + row) << std::endl;
                    keyAvailable[round] = false;
                    break;
                }
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

    for (int round = 1; round < AES_ROUND - 2; round++) {
        if (not (isQtildeAvailable[round] and isQPtildeAvailable[round - 1])) {
            continue;
        }
        std::array<std::array<uint8_t, 256>, 16> roundEncoding;
        for (unsigned int col = 0; col < 4; col++) {
            for (unsigned int row = 0; row < 4; row++) {
                roundEncoding[col * 4 + row] = Qtilde[round][col][row].compose(
                                               QPtilde[round - 1][col].getQ(row)).getEncodingArray();
            }
        }
        encodingMap[round + 1] = std::move(roundEncoding);
    }
    std::map <uint8_t, std::array<std::array<std::array<uint8_t, 4>, 4>, 4>> coeffMap;

    for (int round = 0; round < AES_ROUND - 1; round++) {
        if (isQPtildeAvailable[round]) {
            coeffMap[round + 1] = {
                QPtilde[round][0].getCoeff(),
                QPtilde[round][1].getCoeff(),
                QPtilde[round][2].getCoeff(),
                QPtilde[round][3].getCoeff()
            };
        }
    }


    return std::make_tuple(std::move(keyMap), std::move(encodingMap), std::move(coeffMap));
}

PYBIND11_MODULE(_core, m) {
    m.doc() = "Blue Galaxy Energy module";
    m.def("recoverkey", &recoverkey, "Core of the BGE key recovery", "data"_a, "isEncrypt"_a);

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
