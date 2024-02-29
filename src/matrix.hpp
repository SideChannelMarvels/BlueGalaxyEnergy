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

#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "utils.hpp"
#include <NTL/GF2E.h>
#include <NTL/mat_GF2.h>
#include <NTL/matrix.h>
#include <array>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <set>

#define SIZE_MATRIX 8

// Convert a vector of length 8 to a uint8_t
//   for example, with vec = [1 0 0 0 0 0 1 0], the uint8_t is 130
static uint8_t vec_GF2_to_uint8(const NTL::vec_GF2 &vec) {
//    std::cout << vec << std::endl;
    assert(vec.length() == SIZE_MATRIX);
    uint8_t tmp = 0;
    for (unsigned int i = 0; i < SIZE_MATRIX; i++) {
        tmp = static_cast<uint8_t>(tmp << 1);
        tmp = tmp ^ static_cast<uint8_t>(NTL::rep(vec[i]));
    }

//    std::cout << static_cast<long>(tmp) << std::endl;
    return tmp;
}

// Convert a uint8_t a to a polynomial on GF2 of at most degree 7 as
//   for example, with a = 0b10000010, the polynomial is equal to x^7 + x
static NTL::GF2X uint8_to_GF2X(const uint8_t &a) {
    NTL::GF2X P;
    for (unsigned int i = 0; i < 8; i++) {
        NTL::SetCoeff(P, static_cast<long>(i), static_cast<long>(a >> i & 1));
    }

    return P;
}

// Return the polynomial in the form
//   x^8 + x^4 + x^3 + x + 1 = 0b100011011 = 283
static unsigned int GF2X_to_uint(const NTL::GF2X &p) {
    unsigned int b = 0;
    unsigned int d = static_cast<unsigned int>(NTL::deg(p));
    for (unsigned int i = 0; i <= d; i++) {
        b = b << 1;
        b = b ^ static_cast<unsigned int>(NTL::conv<long>(p[d - i]));
    }

    return b;
}

// Generate the matrix Lambda_alpha.
// alpha being a polynomial and p0 another polynomial,
//   p0 * alpha = p1 mod P (P being the AES polynomial)
// Representing P0 as a vector of the coefficient of P0, and p1 as P1,
//   P0 * Lambda_alpha = P1
NTL::mat_GF2 gen_multiplication_matrix(uint8_t alpha) {
    // AES polynomial: P(x) = x^8 + x^4 + x^3 + x + 1
    NTL::GF2X P;
    NTL::SetCoeff(P, 0, 1);
    NTL::SetCoeff(P, 1, 1);
    NTL::SetCoeff(P, 2, 0);
    NTL::SetCoeff(P, 3, 1);
    NTL::SetCoeff(P, 4, 1);
    NTL::SetCoeff(P, 5, 0);
    NTL::SetCoeff(P, 6, 0);
    NTL::SetCoeff(P, 7, 0);
    NTL::SetCoeff(P, 8, 1);

    NTL::GF2E::init(P);

    NTL::GF2E palphae = NTL::conv<NTL::GF2E>(uint8_to_GF2X(alpha));

    NTL::vec_vec_GF2 vec;
    vec.SetLength(SIZE_MATRIX);

    for (unsigned int i = 0; i < 8; i++) {
        vec[i].SetLength(SIZE_MATRIX);
        NTL::GF2X p;
        NTL::SetCoeff(p, 7 - i, 1);
        p = NTL::conv<NTL::GF2X>(palphae * NTL::conv<NTL::GF2E>(p));
        for (unsigned int c = 0; c < SIZE_MATRIX; c++)
            vec[i][c] = NTL::coeff(p, SIZE_MATRIX - 1 - c);
    }

    return NTL::to_mat_GF2(vec);
}

// Recursive way to compute the determinant of matrices in NTL::GF2X
NTL::GF2X determinant(NTL::Mat<NTL::GF2X> M) {
    if (M.NumRows() == 2)
        return M[0][0] * M[1][1] + M[1][0] * M[0][1];

    NTL::GF2X p;

    // The submatrix on which the recursive call will be done
    NTL::Mat<NTL::GF2X> m;
    m.SetDims(M.NumRows() - 1, M.NumCols() - 1);

    // Construct the submatrix m by first looking for a non-zero value then building the submatrix without the row and
    // column where this value stays
    for (unsigned int c = 0; c < M.NumCols(); c++) {
        if (M[0][c] == 0)
            continue;
        for (unsigned int r = 1; r < M.NumRows(); r++) {
            unsigned int C = 0;
            for (unsigned int col = 0; col < M.NumCols(); col++) {
                if (col != c) {
                    m[r - 1][C] = M[r][col];
                    C           = C + 1;
                }
            }
        }
        p = p + M[0][c] * determinant(m);
    }

    return p;
}

// TODO: we can remove the dependency to NTL::GF2X by using our own gf2x
// implementation on uint16_t since we know that we will not exceed the degree 8
// and do not need to have efficient multiplication, only shift are sufficient.
// Return the polynomial in the form
//   x^8 + x^4 + x^3 + x + 1 = 0b100011011 = 283
unsigned int characteristic_polynomial(NTL::mat_GF2 M) {
    NTL::Mat<NTL::GF2X> m;
    m.SetDims(M.NumRows(), M.NumCols());
    NTL::GF2X p;
    NTL::SetCoeff(p, 1, 1);
    for (unsigned int r = 0; r < M.NumRows(); r++)
        m[r][r] = p;
    for (unsigned int r = 0; r < M.NumRows(); r++)
        for (unsigned int c = 0; c < M.NumCols(); c++) {
            m[r][c] = m[r][c] + M[r][c];
        }

    p = determinant(m);

    return GF2X_to_uint(p);
}

/*
 * This class will be a bit unusual. Indeed, our matrix is of size 8*8 on GF(2).
 * The matrix uses a row-wise representation of the basis vectors, that is we
 * multiply by a vector on the left of the matrix.
 *
 * The unusual stuff here is that we will define the matrix row by row, but what
 * we need is in fact the columns of the matrix, not its rows. So, we will
 * construct the matrix thanks to its rows, but store the matrix by columns.
 *
 * Exemple (on a matrix 2 * 2).
 *   M = matrix([a0, a1], [b0, b1]) -> M = [a0 a1]
 *                                         [b0 b1]
 *
 *   To store M, we will store [[a0, b0], [a1, b1]]
 *   Now, if we want to do v * M, the output will be [v * M[0], v * M[1]]
 *
 * Since we work with 8-bit vectors, the vectors are represented on an
 * uint8_t. To perform the operation v * M[0], we need to do
 *   __builtin_parity(v & M[0])
 */

class Matrix {
  private:
    std::array<uint8_t, SIZE_MATRIX> col;

  public:
    Matrix() {};

    Matrix(const std::array<uint8_t, SIZE_MATRIX> row) { build(row); };

    Matrix(const NTL::mat_GF2 &m) { from_mat_GF2(m); };

    Matrix(uint8_t alpha) { generate_multiplication_matrix(alpha); };

    void build(const std::array<uint8_t, SIZE_MATRIX> row) {
        uint8_t tmp;
        for (unsigned int c = 0; c < SIZE_MATRIX; c++) {
            uint8_t col_tmp = 0;
            for (unsigned int i = 0; i < SIZE_MATRIX; i++) {
                tmp = (row[i] >> (SIZE_MATRIX - 1 - c)) & 1;
                col_tmp =
                    col_tmp ^ (static_cast<uint8_t>(tmp << (SIZE_MATRIX - 1 - i)));
            }
            col[c] = col_tmp;
        }
    }

    uint8_t getCol(unsigned int i) const { return col[i]; }

    bool get(unsigned int r, unsigned int c) const {
        return static_cast<bool>((getCol(c) >> (SIZE_MATRIX - 1 - r)) & 1);
    }

    uint8_t multiply(uint8_t v) {
        uint8_t tmp;
        uint8_t result = 0;
        for (unsigned int i = 0; i < SIZE_MATRIX; i++) {
            tmp = static_cast<uint8_t>(
                      __builtin_parity(static_cast<unsigned int>(v & col[i])));
            result = result ^ (static_cast<uint8_t>(tmp << (SIZE_MATRIX - 1 - i)));
        }

        return result;
    }

    // Brute-force way to determine if an GF(2) 8×8 matrix is full rank, i.e.,
    //   test all the possible vectors v in GF(2)^8: if the result of v * A was already computed by a different vector, A is not full rank
    bool is_full_rank() {
        std::set<uint8_t> set;
        uint8_t res;

        for (unsigned int x = 0; x < LEN_ARRAY; x++) {
            res = multiply(static_cast<uint8_t>(x));
            if (!set.insert(res).second)
                return false;
        }

        return true;
    }

    // Transform a Matrix to a NTL::mat_GF2
    NTL::mat_GF2 to_mat_GF2() const {
        NTL::vec_vec_GF2 vec;
        vec.SetLength(SIZE_MATRIX);
        for (unsigned int r = 0; r < 8; r++) {
            vec[r].SetLength(SIZE_MATRIX);
            for (unsigned int c = 0; c < 8; c++)
                vec[r][c] = static_cast<long>(get(r, c));
        }

        return NTL::to_mat_GF2(vec);
    }

    // Transform a NTL::mat_GF2 to a Matrix
    void from_mat_GF2(const NTL::mat_GF2 &m) {
        std::array<uint8_t, SIZE_MATRIX> row;
        for (unsigned int i = 0; i < SIZE_MATRIX; i++)
            row[i] = vec_GF2_to_uint8(m[i]);

        build(row);
    }

    // See the comment of gen_multiplication_matrix
    void generate_multiplication_matrix(uint8_t alpha) {
        from_mat_GF2(gen_multiplication_matrix(alpha));
    }

    // Evaluate i * Matrix + offset, where i and offset represent bit vectors
    std::array<uint8_t, LEN_ARRAY> tabulate(uint8_t offset = 0) {
        std::array<uint8_t, LEN_ARRAY> res;
        for (unsigned int i = 0; i < LEN_ARRAY; i++)
            res[i] = multiply(static_cast<uint8_t>(i)) ^ offset;

        return res;
    }
};

// Print a Matrix
MAYBE_UNUSED static std::ostream &operator<<(std::ostream &os,
                                             const Matrix &M) {
    for (unsigned int r = 0; r < SIZE_MATRIX - 1; r++) {
        for (unsigned int c = 0; c < SIZE_MATRIX; c++)
            os << static_cast<unsigned int>(M.get(r, c));
        os << std::endl;
    }

    for (unsigned int c = 0; c < SIZE_MATRIX; c++)
        os << static_cast<unsigned int>(M.get(SIZE_MATRIX - 1, c));

    return os;
}

#endif /* MATRIX_HPP */

