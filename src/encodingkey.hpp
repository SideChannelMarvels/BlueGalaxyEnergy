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

#ifndef ENCODINGKEY_HPP
#define ENCODINGKEY_HPP

#include "baseencoding.hpp"

// Correspond to Ptilde in the article
class EncodingKey : public BaseEncoding {
  public:
    EncodingKey() : BaseEncoding() {}
    explicit EncodingKey(const BaseEncoding &arr) : BaseEncoding(arr) {}
    explicit EncodingKey(const std::array<uint8_t, LEN_ARRAY> &arr) : BaseEncoding(arr) {}
};

#endif /* ENCODINGKEY_HPP */

