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

#ifndef RELATIONENCODING_HPP
#define RELATIONENCODING_HPP

#include "matrix.hpp"
#include "outputwb.hpp"

class RelationEncoding {
  private:
    Matrix A;
    uint8_t c;

  public:
    // Find a relationship such that yi = yj * A + c
    RelationEncoding(const BaseEncoding &yi, const BaseEncoding &yj) {
        c = yi[yj.getInv(0)];
        std::array<uint8_t, 8> row;
        for (unsigned int i = 0; i < 8; i++) {
            uint8_t a = (1 << (7 - i));
            row[i] = yi[yj.getInv(a)] ^ c;
        }

        A = Matrix(row);
    };

    Matrix getA() const { return A; }
    uint8_t getC() const { return c; }

    bool assert_computation(const BaseEncoding &yi, const BaseEncoding &yj) {
        for (unsigned int i = 0; i < LEN_ARRAY; i++)
            if (yi[i] != static_cast<uint8_t>(A.multiply(yj[i]) ^ c))
                return false;

        return true;
    }
};

#endif /* RELATIONENCODING_HPP */

