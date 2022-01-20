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

#ifndef APPROXIMATEENCODING_HPP
#define APPROXIMATEENCODING_HPP

#include "baseencoding.hpp"
#include "outputwb.hpp"
#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>

class ApproximateEncoding : public BaseEncoding {
  public:
    ApproximateEncoding() : BaseEncoding() {}

    bool build(const std::array<std::unique_ptr<OutputWB>, LEN_ARRAY> &owbs) {
        if (not std::all_of(owbs.cbegin(), owbs.cend(), [](const std::unique_ptr<OutputWB> &o) { return bool(o); })) {
            std::cerr << "ApproximateEncoding::build: Missing OutputWB for x0" << std::endl;
            return false;
        }

        std::array<BaseEncoding, LEN_ARRAY> functionSet;
        std::array<bool, LEN_ARRAY> valid{};

        const OutputWB &outwb0 = *owbs[0];
        for (unsigned int i = 0; i < LEN_ARRAY; i++) {
            BaseEncoding b = outwb0.compose(owbs[i]->inverse());
            functionSet[b[0]] = b;
            valid[b[0]] = true;
        }

        if (not std::all_of(valid.cbegin(), valid.cend(), [](bool v) { return v; })) {
            std::cerr << "ApproximateEncoding::build: Invalid OutputWB for x0" << std::endl;
            return false;
        }

        // Tolhuizen algorithm
        std::array<uint8_t, LEN_ARRAY> map {};
        std::array<uint8_t, LEN_ARRAY> mapInv {};
        std::array<bool, LEN_ARRAY> inR {};

        map[0] = 0;
        mapInv[0] = 0;
        inR[0] = true;

        unsigned int i = 1;
        for (unsigned int k = 0; k < 8; k++) {
            while (inR[i]) {
                i++;
            }

            map[i] = static_cast<uint8_t>(1 << k);
            mapInv[1 << k] = i;
            inR[i] = true;
            for (unsigned int j = 1; j < LEN_ARRAY; j++) {
                if (inR[j]) {
                    uint8_t fogid = functionSet[i][j];
                    uint8_t value = map[i] ^ map[j];
                    if (inR[fogid]) {
                        if (map[fogid] != value) {
                            std::cerr << "ApproximateEncoding::build: Unexpected value during the Tolhuizen algorithm" << std::endl;
                            return false;
                        }
                    } else if (mapInv[value] != 0 or value == 0) {
                        std::cerr << "ApproximateEncoding::build: Unexpected duplicated value during the Tolhuizen algorithm" << std::endl;
                        return false;
                    } else {
                        map[fogid] = value;
                        mapInv[value] = fogid;
                        inR[fogid] = true;
                    }
                }
            }
        }

        if (not std::all_of(inR.cbegin(), inR.cend(), [](bool v) { return v; })) {
            std::cerr << "ApproximateEncoding::build: Invalid result of the Tolhuizen algorithm" << std::endl;
            return false;
        }
        this->val = mapInv;
        this->inv = map;
        return true;
    }
};

#endif /* APPROXIMATEENCODING_HPP */

