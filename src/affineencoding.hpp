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

#ifndef AFFINEENCODING_HPP
#define AFFINEENCODING_HPP

#include "baseencoding.hpp"

class AffineEncoding : public BaseEncoding {
  public:
    AffineEncoding() : BaseEncoding() {}
    explicit AffineEncoding(const BaseEncoding &arr) : BaseEncoding(arr) {}
    explicit AffineEncoding(const std::array<uint8_t, LEN_ARRAY> &arr) : BaseEncoding(arr) {}

    AffineEncoding(const BaseEncoding &arr, uint8_t q) : BaseEncoding() {
        construct(arr.getEncodingArray(), q);
    }

    AffineEncoding(const std::array<uint8_t, LEN_ARRAY> &arr, uint8_t q) : BaseEncoding() {
        construct(arr, q);
    }

  protected:
    void construct(const std::array<uint8_t, LEN_ARRAY> &arr, uint8_t q) {
        for (unsigned int i = 0; i < LEN_ARRAY; i++) {
            val[i] = arr[i] ^ q;
            inv[arr[i] ^ q] = i;
        }
        if (! isvalid()) {
            std::cerr << "BaseEncoding::isvalid : False, abort!" << std::endl;
            abort();
        }
    }
};
#endif /* AFFINEENCODING_HPP */

