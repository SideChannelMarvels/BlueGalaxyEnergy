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

#ifndef COMPUTEQPTILDE_HPP
#define COMPUTEQPTILDE_HPP

#include "affineencoding.hpp"
#include "baseencoding.hpp"
#include "encodingkey.hpp"
#include "matrix.hpp"
#include "precompute.hpp"
#include "relationencoding.hpp"
#include <NTL/mat_GF2.h>

static NTL::mat_GF2 compute_Atilde(const BaseEncoding &y0x0, const BaseEncoding &y1x0,
                                   const BaseEncoding &y0x1, const BaseEncoding &y1x1) {
    // TODO: l0 and l1 can be precomputed or can be reused after, depending on the
    // flow of the algorithm
    RelationEncoding l0(y0x0, y1x0);
    RelationEncoding l1(y0x1, y1x1);

    NTL::mat_GF2 L0 = l0.getA().to_mat_GF2();
    NTL::mat_GF2 L1 = l1.getA().to_mat_GF2();
    NTL::mat_GF2 L;
    L.SetDims(SIZE_MATRIX, SIZE_MATRIX);
    NTL::inv(L1, L1);
    mul(L, L1, L0);

    unsigned int req   = characteristic_polynomial(L);
    NTL::mat_GF2 Lbeta = gen_multiplication_matrix(get_beta(req));

    unsigned int length = SIZE_MATRIX * SIZE_MATRIX;
    NTL::mat_GF2 S;
    S.SetDims(length, length);
    for (unsigned int s = 0; s < SIZE_MATRIX; s++)
        for (unsigned int r = 0; r < SIZE_MATRIX; r++)
            for (unsigned int c = 0; c < SIZE_MATRIX; c++) {
                S[r + s * SIZE_MATRIX][c + s * SIZE_MATRIX] =
                    S[r + s * SIZE_MATRIX][c + s * SIZE_MATRIX] + L[r][c];
                S[r + s * SIZE_MATRIX][r + c * SIZE_MATRIX] =
                    S[r + s * SIZE_MATRIX][r + c * SIZE_MATRIX] + Lbeta[c][s];
            }

    NTL::mat_GF2 ker;
    NTL::kernel(ker, S);

    NTL::mat_GF2 Atilde;
    Atilde.SetDims(SIZE_MATRIX, SIZE_MATRIX);
    for (unsigned int r = 0; r < SIZE_MATRIX; r++)
        for (unsigned int c = 0; c < SIZE_MATRIX; c++)
            Atilde[r][c] = ker[0][r * SIZE_MATRIX + c];

    return Atilde;
}

static bool computecd(const BaseEncoding &xy,
                      const BaseEncoding &Atilde_inv_tab,
                      uint8_t &c, uint8_t &d,
                      EncodingKey &Ptilde) {

    for (unsigned int delta = 1; delta < LEN_ARRAY; delta++) {
        const BaseEncoding deltaMul = Matrix(delta).tabulate();

        // precompute SBOX_inv[deltaMul[Atilde_inv_tab[x]]]
        const BaseEncoding commonF = SBOX_inv.compose(deltaMul).compose(Atilde_inv_tab);

        for (unsigned k = 0; k < LEN_ARRAY; k++) {
            bool is_affine = true;
            std::array<uint8_t, LEN_ARRAY> v = {};
            for (unsigned int x = 0; x < LEN_ARRAY; x++) {
                v[x] = commonF[xy[x] ^ k];

                // verify immediately if the new value verify the affine property
                // note: if y < x but x^y < x, we cannot check the equation now
                //       but it will be checked later when X = x^y
                for (unsigned int y = 1; y < x; y++) {
                    if ((x ^ y) < x and v[x ^ y] != (v[x] ^ v[y] ^ v[0])) {
                        is_affine = false;
                        break;
                    }
                }
                if (not is_affine) {
                    break;
                }
                // should never happen
                // discard constant
                if (x != 0 and v[0] == v[x]) {
                    is_affine = false;
                    break;
                }
            }
            if (is_affine) {
                d = delta;
                c = k;
                Ptilde = EncodingKey(v);
                return true;
            }
        }
    }
    return false;
}

static bool computeQPtilde(const std::array<std::array<BaseEncoding, 4>, 4> &xy,
                           std::array<AffineEncoding, 4> &Q,
                           std::array<EncodingKey, 4> &Ptilde) {

    for (unsigned int row = 0; row < 4; row++) {
        NTL::mat_GF2 Atilde = compute_Atilde(xy[0][row], xy[0][row ^ 1], xy[1][row], xy[1][row ^ 1]);

        NTL::mat_GF2 Atilde_inv;
        Atilde_inv.SetDims(SIZE_MATRIX, SIZE_MATRIX);
        NTL::inv(Atilde_inv, Atilde);
        std::array<uint8_t, LEN_ARRAY> Atilde_inv_tab = Matrix(Atilde_inv).tabulate();

        std::array<uint8_t, 4> c;
        std::array<uint8_t, 4> d;

        std::array<EncodingKey, 4> Ptilde_tmp;

        if (not computecd(xy[0][row], Atilde_inv_tab, c[0], d[0], Ptilde_tmp[0])) return false;
        if (not computecd(xy[1][row], Atilde_inv_tab, c[1], d[1], Ptilde_tmp[1])) return false;
        if (not computecd(xy[2][row], Atilde_inv_tab, c[2], d[2], Ptilde_tmp[2])) return false;
        if (not computecd(xy[3][row], Atilde_inv_tab, c[3], d[3], Ptilde_tmp[3])) return false;

        if (row == 0) {
            Ptilde = std::move(Ptilde_tmp);
        } else if (Ptilde != Ptilde_tmp) {
            return false;
        }

        uint8_t q = xy[0][row][0] ^ c[0] ^ c[1] ^ c[2] ^ c[3];
        bool foundD = false;

        for (unsigned int i = 0; i < 6; i++) {
            uint64_t j = (i + 1 + (i >> 2)) % 4;
            if (d[i % 4] == d[j]) {
                foundD = true;

                NTL::mat_GF2 Lambda_gamma = gen_multiplication_matrix(inverse_GF28[d[j]]);
                Q[row] = AffineEncoding(Matrix(Lambda_gamma * Atilde).tabulate(q));
                break;
            }
        }
        if (not foundD) {
            return false;
        }
    }
    return true;
}

#endif /* COMPUTEQPTILDE_HPP */

