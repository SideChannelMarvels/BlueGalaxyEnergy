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

#ifndef BASEENCODING_HPP
#define BASEENCODING_HPP

#include "utils.hpp"
#include <array>
#include <cstddef>
#include <ctype.h>
#include <iostream>
#include <memory>

// This class is an encoding, i.e. a bijective function on integer in [0,256)
class BaseEncoding {
  protected:
    std::array<uint8_t, LEN_ARRAY> val = {};
    std::array<uint8_t, LEN_ARRAY> inv = {};

    constexpr BaseEncoding(const std::array<uint8_t, LEN_ARRAY> &arr,
                           const std::array<uint8_t, LEN_ARRAY> &invarr)
        : val(arr), inv(invarr) {}

    constexpr void computeInv() {
        for (unsigned int i = 0; i < LEN_ARRAY; i++) {
            inv[val[i]] = static_cast<uint8_t>(i);
        }
    }

    template <class T>
    static std::unique_ptr<T> from_array(const std::array<uint8_t, LEN_ARRAY> &val) {
        std::array<bool, LEN_ARRAY> used = {};
        for (unsigned int i = 0; i < LEN_ARRAY; i++) {
            if (used[val[i]]) {
                return {};
            } else {
                used[val[i]] = true;
            }
        }

        std::unique_ptr<T> e = std::make_unique<T>(val);
        if (e->isvalid()) {
            return e;
        }
        return {};
    }

    constexpr bool isvalid() const {
        bool valid = false;
        for (unsigned int i = 0; i < LEN_ARRAY; i++) {
            if (inv[i] == 0) {
                if (valid) {
                    return false;
                } else {
                    valid = true;
                }
            }
        }
        return valid;
    }

  public:

    // an encoding should always be valid.
    // The default constructor is the identity function
    constexpr BaseEncoding() {
        for (unsigned int i = 0; i < LEN_ARRAY; i++) {
            val[i] = static_cast<uint8_t>(i);
            inv[i] = static_cast<uint8_t>(i);
        }
    }

    constexpr BaseEncoding(const std::array<uint8_t, LEN_ARRAY> &arr) : val(arr) {
        computeInv();
        if (! isvalid()) {
            std::cerr << "BaseEncoding::isvalid : False, abort!" << std::endl;
            abort();
        }
    }

    constexpr uint8_t getVal(unsigned int i) const { return val[i]; }
    constexpr uint8_t operator[](unsigned int i) const { return val[i]; }

    constexpr uint8_t getInv(unsigned int i) const { return inv[i]; }

    bool operator==(const BaseEncoding &other) const {
        return this->val == other.val;
    };

    bool operator!=(const BaseEncoding &other) const {
        return this->val != other.val;
    };

    constexpr BaseEncoding inverse() const {
        return BaseEncoding(inv, val);
    }

    // do this o other == this(other(...))
    constexpr BaseEncoding compose(const BaseEncoding &other) const {
        std::array<uint8_t, LEN_ARRAY> n_val = {};
        std::array<uint8_t, LEN_ARRAY> n_inv = {};
        for (unsigned int i = 0; i < LEN_ARRAY; i++) {
            n_val[i] = val[other[i]];
            n_inv[n_val[i]] = i;
        }
        return BaseEncoding(n_val, n_inv);
    }

    const std::array<uint8_t, LEN_ARRAY> &getEncodingArray() const {
        return val;
    }
};

MAYBE_UNUSED static std::ostream &operator<<(std::ostream &os,
                                             const BaseEncoding &e) {
    os << static_cast<unsigned int>(e[0]);
    for (unsigned int r = 1; r < LEN_ARRAY; r++) {
        os << ", " << static_cast<unsigned int>(e[r]);
    }

    return os;
}

#endif /* BASEENCODING_HPP */

