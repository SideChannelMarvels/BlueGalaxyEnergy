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

#ifndef OUTPUTWB_HPP
#define OUTPUTWB_HPP

#include "baseencoding.hpp"
#include "utils.hpp"

class OutputWB : public BaseEncoding {
  public:
    OutputWB(const std::array<uint8_t, LEN_ARRAY> &arr) : BaseEncoding(arr) {}

    static std::unique_ptr<OutputWB> from_array(const std::array<uint8_t, LEN_ARRAY> &arr) {
        return BaseEncoding::from_array<OutputWB>(arr);
    }

};

#endif /* OUTPUTWB_HPP */

