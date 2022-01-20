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

#ifndef INPUT_READER_HPP
#define INPUT_READER_HPP

#include "outputwb.hpp"
#include "utils.hpp"
#include <algorithm>
#include <array>
#include <map>
#include <memory>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>

class InputReader {
  private:

    // type (x,y,0,0,x,y,0,0,x,y,0,0,x,y,0,0) with y in [0,256)
    std::array<std::array<std::array<std::unique_ptr<OutputWB>, LEN_ARRAY>, 16>, AES_ROUND> x0 = {};
    // type (0,x,0,0,0,x,0,0,0,x,0,0,0,x,0,0)
    std::array<std::array<std::unique_ptr<OutputWB>, 16>, AES_ROUND> x1 = {};
    // type (0,0,x,0,0,0,x,0,0,0,x,0,0,0,x,0)
    std::array<std::array<std::unique_ptr<OutputWB>, 16>, AES_ROUND> x2 = {};
    // type (0,0,0,x,0,0,0,x,0,0,0,x,0,0,0,x)
    std::array<std::array<std::unique_ptr<OutputWB>, 16>, AES_ROUND> x3 = {};

    std::array<bool, AES_ROUND> roundValid = {0};

  public:
    InputReader() {}

    void readData(const std::map < uint8_t, std::array < std::array < std::array<uint8_t, LEN_ARRAY>, LEN_ARRAY + 2 >, 16 >> & data) {
        for (auto const &el : data) {
            addRoundData(el.first - 1, el.second);
        }
    }

    bool isRoundAvailable(int round) const {
        if (round < 0 or round >= AES_ROUND) {
            return false;
        }
        return roundValid[round];
    }

    const std::array<std::unique_ptr<OutputWB>, LEN_ARRAY> &getx0_all(int round, int column, int row) const {
        return x0[round][column * 4 + row];
    }

    const OutputWB &getx(int x, int round, int column, int row) const {

        switch (x) {
            case 0:
                if (! x0[round][column * 4 + row][0]) {
                    std::cerr << "InputReader::getx : x0[round][column*4+row][0] is NULL, abort!" << std::endl;
                    abort();
                }
                return *x0[round][column * 4 + row][0];
            case 1:
                if (! x1[round][column * 4 + row]) {
                    std::cerr << "InputReader::getx : x1[round][column*4+row] is NULL, abort!" << std::endl;
                    abort();
                }
                return *x1[round][column * 4 + row];
            case 2:
                if (! x2[round][column * 4 + row]) {
                    std::cerr << "InputReader::getx : x2[round][column*4+row] is NULL, abort!" << std::endl;
                    abort();
                }
                return *x2[round][column * 4 + row];
            case 3:
                if (! x3[round][column * 4 + row]) {
                    std::cerr << "InputReader::getx : x3[round][column*4+row] is NULL, abort!" << std::endl;
                    abort();
                }
                return *x3[round][column * 4 + row];
            default:
                abort();
        }
    }

  private:
    void addRoundData(unsigned int round, const std::array < std::array < std::array<uint8_t, LEN_ARRAY>, LEN_ARRAY + 2 >, 16 > & arr) {

        if (round >= AES_ROUND) {
            std::cerr << "Invalid round " << (round + 1) << std::endl;
            return;
        }

        if (isRoundAvailable(round)) {
            std::cerr << "Round " << (round + 1) << " already imported" << std::endl;
            return;
        }

        for (int pos = 0; pos < 16; pos++) {
            for (int y = 0; y < LEN_ARRAY; y++) {
                if (not addLine(x0[round][pos][y], arr[pos][y])) {
                    std::cerr << "Fail import x0 round " << (round + 1)
                              << ", byte " << pos
                              << ", y " << y << std::endl;
                    return;
                }
            }
            if (not computeX1(round, pos)) {
                std::cerr << "Fail compute x1 round " << (round + 1)
                          << ", byte " << pos << std::endl;
                return;
            }
            if (not addLine(x2[round][pos], arr[pos][256])) {
                std::cerr << "Fail import x2 round " << (round + 1)
                          << ", byte " << pos << std::endl;
                return;
            }
            if (not addLine(x3[round][pos], arr[pos][257])) {
                std::cerr << "Fail import x3 round " << (round + 1)
                          << ", byte " << pos << std::endl;
                return;
            }
        }
        roundValid[round] = true;
    }

    inline bool addLine(std::unique_ptr<OutputWB> &dest, const std::array<uint8_t, LEN_ARRAY> &v) {

        std::unique_ptr<OutputWB> p = OutputWB::from_array(v);

        if (p) {
            dest = std::move(p);
            return true;
        } else {
            return false;
        }
    }

    inline bool computeX1(int round, int pos) {
        if (not std::all_of(x0[round][pos].cbegin(), x0[round][pos].cend(),
        [](const std::unique_ptr<OutputWB> &o) { return bool(o); })) {
            return false;
        }

        std::array<uint8_t, LEN_ARRAY> v;
        for (size_t i = 0; i < LEN_ARRAY; i++) {
            v[i] = x0[round][pos][i]->getVal(0);
        }

        std::unique_ptr<OutputWB> e = OutputWB::from_array(v);
        if (e) {
            x1[round][pos] = std::move(e);
            return true;
        } else {
            return false;
        }
    }
};

#endif

