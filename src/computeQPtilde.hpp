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
#include <algorithm>
#include <utility>

static unsigned int get_beta_req(const BaseEncoding &y0x0, const BaseEncoding &y1x0,
                               const BaseEncoding &y0x1, const BaseEncoding &y1x1,
                               NTL::mat_GF2& L) {

    // TODO: l0 and l1 can be precomputed or can be reused after, depending on the
    // flow of the algorithm
    RelationEncoding l0(y0x0, y1x0);
    RelationEncoding l1(y0x1, y1x1);

    NTL::mat_GF2 L0 = l0.getA().to_mat_GF2();
    NTL::mat_GF2 L1 = l1.getA().to_mat_GF2();
    L.SetDims(SIZE_MATRIX, SIZE_MATRIX);
    NTL::inv(L1, L1);
    mul(L, L1, L0);

    return characteristic_polynomial(L);
}

static BaseEncoding compute_Atilde(const NTL::mat_GF2& L, uint8_t beta) {
    NTL::mat_GF2 Lbeta = gen_multiplication_matrix(beta);

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

    return Matrix(Atilde).tabulate();
}

static bool get_betas(const std::array<unsigned int, 4>& req,
                      std::array<uint8_t, 4>& res,
                      std::array<std::array<std::array<uint8_t, 2>, 2>, 4>& coeff,
                      bool& isUnsure) {
    std::array<unsigned int, 4> sortedReq = req;
    std::sort(sortedReq.begin(), sortedReq.end());

    auto it = std::find_if(PrecomputedBeta.cbegin(), PrecomputedBeta.cend(),
        [&sortedReq](const precomputedBetaStruct& el) {
            return el.beta_req == sortedReq;
        });

    if (it == PrecomputedBeta.cend()) {
        std::cerr << "get_betas error: " << req[0]
                  << ", " << req[1]
                  << ", " << req[2]
                  << ", " << req[3] << std::endl;
        return false;
    }
    isUnsure = it->UnclearRequest;

    std::array<bool, 4> used = {false, false, false, false};
    for (int reqIndex = 0; reqIndex < 4; reqIndex++) {
        bool found = false;
        for (int precomputedIndex = 0; precomputedIndex < 4; precomputedIndex++) {
            if ((!used[precomputedIndex]) && it->beta_req[precomputedIndex] == req[reqIndex]) {
                res[reqIndex] = it->beta[precomputedIndex];
                coeff[reqIndex] = it->coeff[precomputedIndex];
                used[precomputedIndex] = true;
                found = true;
                break;
            }
        }
        if (!found) {
            std::cerr << "get_betas cannot get beta for index " << reqIndex
                      << " in " << req[0]
                      << ", " << req[1]
                      << ", " << req[2]
                      << ", " << req[3] << std::endl;
            return false;
        }
    }

    return true;
}

static bool checkPreCoeff(const std::array<std::array<std::array<uint8_t, 2>, 2>, 4>& Coeff) {
    if (Coeff[0][0] != Coeff[1][1] ||
        Coeff[0][1] != Coeff[1][0] ||
        Coeff[2][0] != Coeff[3][1] ||
        Coeff[2][1] != Coeff[3][0]) {

        std::cerr << "coeff permutation error = [[["
                  << (unsigned)Coeff[0][0][0] << ", " << (unsigned)Coeff[0][0][1] << "], ["
                  << (unsigned)Coeff[0][1][0] << ", " << (unsigned)Coeff[0][1][1] << "]], [["
                  << (unsigned)Coeff[1][0][0] << ", " << (unsigned)Coeff[1][0][1] << "], ["
                  << (unsigned)Coeff[1][1][0] << ", " << (unsigned)Coeff[1][1][1] << "]], [["
                  << (unsigned)Coeff[2][0][0] << ", " << (unsigned)Coeff[2][0][1] << "], ["
                  << (unsigned)Coeff[2][1][0] << ", " << (unsigned)Coeff[2][1][1] << "]], [["
                  << (unsigned)Coeff[3][0][0] << ", " << (unsigned)Coeff[3][0][1] << "], ["
                  << (unsigned)Coeff[3][1][0] << ", " << (unsigned)Coeff[3][1][1] << "]]]" << std::endl;
        return false;
    }
    return true;
}

static bool checkCoeff(const std::array<std::array<uint8_t, 4>, 4>& Coeff, bool isEncrypt) {
    static const std::array<uint8_t, 4> targetEnc = {1, 1, 2, 3};
    static const std::array<uint8_t, 4> targetDec = {9, 11, 13, 14};
    const std::array<uint8_t, 4>& target = (isEncrypt)? targetEnc : targetDec;

    bool res = true;

    for (uint8_t col = 0; col < 4; col++) {
        std::array<uint8_t, 4> c = Coeff[col];
        std::sort(c.begin(), c.end());
        res &= (c == target);
    }
    for (uint8_t row = 0; row < 4; row++) {
        std::array<uint8_t, 4> c = {Coeff[0][row], Coeff[1][row], Coeff[2][row], Coeff[3][row]};
        std::sort(c.begin(), c.end());
        res &= (c == target);
    }

    if (!res) {
        std::cerr << "Coefficient doesn't match expected MixColumns value : [["
                  << (unsigned)Coeff[0][0] << ", " << (unsigned)Coeff[0][1] << ", "
                  << (unsigned)Coeff[0][2] << ", " << (unsigned)Coeff[0][3] << "], ["
                  << (unsigned)Coeff[1][0] << ", " << (unsigned)Coeff[1][1] << ", "
                  << (unsigned)Coeff[1][2] << ", " << (unsigned)Coeff[1][3] << "], ["
                  << (unsigned)Coeff[2][0] << ", " << (unsigned)Coeff[2][1] << ", "
                  << (unsigned)Coeff[2][2] << ", " << (unsigned)Coeff[2][3] << "], ["
                  << (unsigned)Coeff[3][0] << ", " << (unsigned)Coeff[3][1] << ", "
                  << (unsigned)Coeff[3][2] << ", " << (unsigned)Coeff[3][3] << "]]" << std::endl;
    }

    return res;
}

static bool commonCoeff(const std::array<uint8_t, 2>& a, const std::array<uint8_t, 2>& b, uint8_t& res, bool tryRun=false) {
    for (unsigned idA = 0; idA < 2; idA++) {
        for (unsigned idB = 0; idB < 2; idB++) {
            // note : special case if a[idA ^ 1] == b[idB ^ 1] == 1
            // As the coefficient 1 is used twice in AES encrypt MixColumns,
            // we can assume that we want the other coefficient if we found two pairs
            if (a[idA] == b[idB] && (a[idA ^ 1] != b[idB ^ 1] || a[idA ^ 1] == 1)) {
                res = a[idA];
                return true;
            }
        }
    }
    if (!tryRun) {
        std::cerr << "Fail commonCoeff with : ["
                  << (unsigned)a[0] << ", " << (unsigned)a[1] << "], ["
                  << (unsigned)b[0] << ", " << (unsigned)b[1] << "]" << std::endl;
    }
    return false;
}

static bool notCoeff(uint8_t a, const std::array<uint8_t, 2>& b, uint8_t& res) {
    if (a == b[0]) {
        res = b[1];
        return true;
    } else if (a == b[1]) {
        res = b[0];
        return true;
    } else {
        std::cerr << "fail notCoeff "
                  << (unsigned)a << " in ("
                  << (unsigned)b[0] << ", "
                  << (unsigned)b[1] << ")" << std::endl;
        return false;
    }
}

static bool resolveCoeff(
    const std::array<std::array<std::array<std::array<uint8_t, 2>, 2>, 4>, 4>& preCoeff,
    std::array<std::array<uint8_t, 4>, 4>& coeff,
    bool isEncrypt,
    bool tryRun) {

    if (!checkPreCoeff(preCoeff[0])) return false;
    if (!checkPreCoeff(preCoeff[1])) return false;
    if (!checkPreCoeff(preCoeff[2])) return false;
    if (!checkPreCoeff(preCoeff[3])) return false;

    for (unsigned p0 = 0; p0 < 2; p0++) {
        for (unsigned p1 = 0; p1 < 4; p1+=2) {
            if (!commonCoeff(preCoeff[0][p1][p0], preCoeff[1][p1][p0], coeff[p0][p1], tryRun)) {
                return false;
            }
            if (!notCoeff(coeff[p0][p1], preCoeff[0][p1][p0], coeff[p0^1][p1^1])) {
                return false;
            }
            if (!commonCoeff(preCoeff[2][p1][p0], preCoeff[3][p1][p0], coeff[p0^2][p1], tryRun)) {
                return false;
            }
            if (!notCoeff(coeff[p0^2][p1], preCoeff[2][p1][p0], coeff[p0^3][p1^1])) {
                return false;
            }
        }
    }

    return checkCoeff(coeff, isEncrypt);
}

static void switchCoeff(
    std::array<bool, 4>& isUnsure,
    const std::array<std::array<unsigned int, 4>, 4>& BetaReq,
    std::array<std::array<uint8_t, 4>, 4>& Beta,
    std::array<std::array<std::array<std::array<uint8_t, 2>, 2>, 4>, 4>& preCoeff,
    unsigned v) {

    for (unsigned reqNum = 0; reqNum < 4; reqNum++) {
        if (!isUnsure[reqNum]) {
            continue;
        }
        bool perform = v & 1;
        v = v >> 1;
        if (!perform) {
            continue;
        }

        auto it = std::find_if(BetaReq[reqNum].cbegin(), BetaReq[reqNum].cend(),
            [](unsigned x) {
                return x == 471;
            });
        if (it == BetaReq[reqNum].cend()) {
            std::cerr << "switchCoeff fail to find first position of coeff 471" << std::endl;
            exit(-1);
        }
        unsigned pos1 = std::distance(BetaReq[reqNum].cbegin(), it);

        it = std::find_if(BetaReq[reqNum].cbegin() + pos1 + 1, BetaReq[reqNum].cend(),
            [](unsigned x) {
                return x == 471;
            });
        if (it == BetaReq[reqNum].cend()) {
            std::cerr << "switchCoeff fail to find second position of coeff 471" << std::endl;
            exit(-1);
        }
        unsigned pos2 = std::distance(BetaReq[reqNum].cbegin(), it);

        Beta[reqNum][pos1] = std::exchange(Beta[reqNum][pos2], Beta[reqNum][pos1]);
        preCoeff[reqNum][pos1] = std::exchange(preCoeff[reqNum][pos2], preCoeff[reqNum][pos1]);
    }
    if (v != 0) {
        std::cerr << "switchCoeff internal error" << std::endl;
        exit(-1);
    }
}

static bool computeCoeffXYAtilde(
    const std::array<std::array<BaseEncoding, 4>, 4> &xy,
    std::array<BaseEncoding, 4>& Atilde,
    std::array<std::array<uint8_t, 4>, 4>& coeff,
    bool isEncrypt) {

    std::array<std::array<unsigned int, 4>, 4> BetaReq;
    std::array<NTL::mat_GF2, 4> L;
    NTL::mat_GF2 L_dummy;
    for (unsigned int row = 0; row < 4; row+=2) {
        BetaReq[0][row + 0] = get_beta_req(xy[0][row], xy[0][row ^ 1], xy[1][row], xy[1][row ^ 1], L[row + 0]);
        BetaReq[0][row + 1] = get_beta_req(xy[0][row ^ 1], xy[0][row], xy[1][row ^ 1], xy[1][row], L[row + 1]);

        BetaReq[1][row + 0] = get_beta_req(xy[0][row], xy[0][row ^ 3], xy[1][row], xy[1][row ^ 3], L_dummy);
        BetaReq[1][row + 1] = get_beta_req(xy[0][row ^ 3], xy[0][row], xy[1][row ^ 3], xy[1][row], L_dummy);

        BetaReq[2][row + 0] = get_beta_req(xy[2][row], xy[2][row ^ 1], xy[3][row], xy[3][row ^ 1], L_dummy);
        BetaReq[2][row + 1] = get_beta_req(xy[2][row ^ 1], xy[2][row], xy[3][row ^ 1], xy[3][row], L_dummy);

        BetaReq[3][row + 0] = get_beta_req(xy[2][row], xy[2][row ^ 3], xy[3][row], xy[3][row ^ 3], L_dummy);
        BetaReq[3][row + 1] = get_beta_req(xy[2][row ^ 3], xy[2][row], xy[3][row ^ 3], xy[3][row], L_dummy);
    }

    std::array<std::array<std::array<std::array<uint8_t, 2>, 2>, 4>, 4> preCoeff;

    std::array<bool, 4> isUnsure;
    std::array<std::array<uint8_t, 4>, 4> Beta;

    if (!get_betas(BetaReq[0], Beta[0], preCoeff[0], isUnsure[0])) return false;
    if (!get_betas(BetaReq[1], Beta[1], preCoeff[1], isUnsure[1])) return false;
    if (!get_betas(BetaReq[2], Beta[2], preCoeff[2], isUnsure[2])) return false;
    if (!get_betas(BetaReq[3], Beta[3], preCoeff[3], isUnsure[3])) return false;

    uint8_t numPossibility =
            ((isUnsure[0])? 2 :1) *
            ((isUnsure[1])? 2 :1) *
            ((isUnsure[2])? 2 :1) *
            ((isUnsure[3])? 2 :1);

    if (numPossibility != 1) {
        for (unsigned v = 0; v < numPossibility; v++) {
            switchCoeff(isUnsure, BetaReq, Beta, preCoeff, v);

            if (resolveCoeff(preCoeff, coeff, isEncrypt, true)) {
                break;
            }
            switchCoeff(isUnsure, BetaReq, Beta, preCoeff, v);
        }
    }

    if (!resolveCoeff(preCoeff, coeff, isEncrypt, false)) {
        return false;
    }

    for (unsigned int row = 0; row < 4; row++) {
        Atilde[row] = compute_Atilde(L[row], Beta[0][row]);
    }
    return true;
}

static bool compute_affine(const BaseEncoding &xy,
                           const BaseEncoding &commonF,
                           uint8_t k,
                           EncodingKey &Ptilde,
                           bool initPtilde) {

    std::array<uint8_t, LEN_ARRAY> v = {};
    for (unsigned int x = 0; x < LEN_ARRAY; x++) {
        v[x] = commonF[xy[x] ^ k];

        if (!initPtilde && v[x] != Ptilde[x]) {
            return false;
        }

        // verify immediately if the new value verify the affine property
        // note: if y < x but x^y < x, we cannot check the equation now
        //       but it will be checked later when X = x^y
        for (unsigned int y = 1; y < x; y++) {
            if ((x ^ y) < x and v[x ^ y] != (v[x] ^ v[y] ^ v[0])) {
                return false;
            }
        }
        // should never happen
        // discard constant
        if (x != 0 and v[0] == v[x]) {
            return false;
        }
    }
    if (initPtilde) {
        Ptilde = EncodingKey(v);
    }
    return true;
}

static bool computec_(const BaseEncoding &xy,
                      const BaseEncoding &Atilde_inv,
                      uint8_t &c, uint8_t delta,
                      EncodingKey &Ptilde,
                      bool isEncrypt,
                      bool initPtilde) {

    // precompute SBOX_inv[deltaMul[Atilde_inv[x]]]
    const BaseEncoding commonF = ((isEncrypt) ?
        (SBOX_inv.compose(Mult_GF28_tab[delta]).compose(Atilde_inv)) :
        (SBOX.compose(Mult_GF28_tab[delta]).compose(Atilde_inv)));

    for (unsigned k = 0; k < LEN_ARRAY; k++) {
        if (compute_affine(xy, commonF, k, Ptilde, initPtilde)) {
            c = k;
            return true;
        }
    }
    return false;
}

static bool computecd(const BaseEncoding &xy,
                      const BaseEncoding &Atilde_inv,
                      uint8_t &c, uint8_t &d,
                      EncodingKey &Ptilde,
                      bool isEncrypt,
                      bool initPtilde) {

    if (isEncrypt) {
        for (unsigned int delta = 1; delta < LEN_ARRAY; delta++) {
            if (computec_(xy, Atilde_inv, c, delta, Ptilde, isEncrypt, initPtilde)) {
                d = delta;
                return true;
            }
        }
        std::cerr << "computecd fail"  << std::endl;
        return false;
    } else {
        EncodingKey Ptilde_tmp;
        if (!computec_(xy, Atilde_inv, c, 1, Ptilde_tmp, isEncrypt, true)) {
            std::cerr << "computecd determine c"  << std::endl;
            return false;
        }
        if (initPtilde) {
            Ptilde = Ptilde_tmp;
            d = 1;
            return true;
        }
        if (Ptilde_tmp == Ptilde) {
            d = 1;
            return true;
        }
        for (unsigned delta = 2; delta < LEN_ARRAY; delta++) {
            const BaseEncoding commonF = SBOX.compose(Mult_GF28_tab[delta]).compose(Atilde_inv);

            if (compute_affine(xy, commonF, c, Ptilde, false)) {
                d = delta;
                return true;
            }
        }
        std::cerr << "computecd determine delta" << std::endl;
        return false;
    }
}

static bool computec(const BaseEncoding &xy,
                     const BaseEncoding &Atilde_inv,
                     uint8_t &c, uint8_t delta,
                     EncodingKey &Ptilde,
                     bool isEncrypt,
                     bool initPtilde) {

    if (computec_(xy, Atilde_inv, c, delta, Ptilde, isEncrypt, initPtilde)) {
        return true;
    }
    std::cerr << "computec fail"  << std::endl;
    return false;
}

class QPtildeLinear {
  private:
    std::array<BaseEncoding, 4> Atilde_;
    std::array<std::array<uint8_t, 4>, 4> Coeff_;
    std::array<AffineEncoding, 4> Q_;
    std::array<EncodingKey, 4> Ptilde_;
    std::array<uint8_t, 4> q_;
    std::array<uint8_t, 4> d_;

    bool isFinalize_ = false;
    bool isValid_ = false;

  public:

    QPtildeLinear() : isValid_(false) {}

    QPtildeLinear(
        const std::array<std::array<BaseEncoding, 4>, 4>& xy,
        std::array<BaseEncoding, 4>&& Atilde,
        std::array<std::array<uint8_t, 4>, 4>&& Coeff,
        bool isEncrypt);

    const std::array<std::array<uint8_t, 4>, 4>& getCoeff() const {
        return Coeff_;
    }

    bool isValid() const { return isValid_; }
    bool isFinalize() const { return isFinalize_ && isValid_; }

    const std::array<AffineEncoding, 4>& getQ() const {
        return Q_;
    }

    const AffineEncoding& getQ(uint8_t row) const {
        return Q_[row];
    }

    const std::array<EncodingKey, 4>& getPtilde() const {
        return Ptilde_;
    }

    const EncodingKey& getPtilde(uint8_t row) const {
        return Ptilde_[row];
    }

    void setNextDelta();
    void setFinal() { isFinalize_ = true; }
};

static bool computeQPtilde(const std::array<std::array<BaseEncoding, 4>, 4> &xy,
                           QPtildeLinear &QPtilde,
                           bool isEncrypt) {


    std::array<std::array<uint8_t, 4>, 4> Coeff;
    std::array<BaseEncoding, 4> Atilde;
    if (!computeCoeffXYAtilde(xy, Atilde, Coeff, isEncrypt)) return false;


    QPtilde = QPtildeLinear(xy, std::move(Atilde), std::move(Coeff), isEncrypt);
    if (QPtilde.isValid() && isEncrypt && !QPtilde.isFinalize()) {
        std::cerr << "computeQPtilde expect QPtilde to be finalized" << std::endl;
        return false;
    }
    return QPtilde.isValid();
}

QPtildeLinear::QPtildeLinear(
        const std::array<std::array<BaseEncoding, 4>, 4>& xy,
        std::array<BaseEncoding, 4>&& Atilde,
        std::array<std::array<uint8_t, 4>, 4>&& Coeff,
        bool isEncrypt)
        : Atilde_(std::move(Atilde)),
          Coeff_(std::move(Coeff)) {

    isValid_ = false;
    for (unsigned int row = 0; row < 4; row++) {
        BaseEncoding Atilde_inv = Atilde_[row].inverse();

        std::array<uint8_t, 4> c;
        uint8_t d_col0;

        if (!computecd(xy[0][row], Atilde_inv, c[0], d_col0, Ptilde_[0], isEncrypt, row==0)) return;

        d_[row] = Mult_GF28(d_col0, Coeff[0][row]);

        if (!computec(xy[1][row], Atilde_inv, c[1], Div_GF28(d_[row], Coeff[1][row]), Ptilde_[1], isEncrypt, row==0)) return;
        if (!computec(xy[2][row], Atilde_inv, c[2], Div_GF28(d_[row], Coeff[2][row]), Ptilde_[2], isEncrypt, row==0)) return;
        if (!computec(xy[3][row], Atilde_inv, c[3], Div_GF28(d_[row], Coeff[3][row]), Ptilde_[3], isEncrypt, row==0)) return;

        q_[row] = xy[0][row][0] ^ c[0] ^ c[1] ^ c[2] ^ c[3];

        Q_[row] = AffineEncoding(Atilde_[row].compose(Mult_GF28_tab[Div_GF28(1, d_[row])]), q_[row]);
    }

    isValid_ = true;
    isFinalize_ = isEncrypt;
}

void QPtildeLinear::setNextDelta() {
    if (isFinalize()) {
        return;
    }

    uint8_t newd = d_[0];
    if (d_[0] == 255) {
        newd = 1;
    } else {
        newd = d_[0] + 1;
    }

    uint8_t diffD = Div_GF28(newd, d_[0]);
    d_[0] = newd;
    d_[1] = Mult_GF28(d_[1], diffD);
    d_[2] = Mult_GF28(d_[2], diffD);
    d_[3] = Mult_GF28(d_[3], diffD);

    for (unsigned int row = 0; row < 4; row++) {
        Q_[row] = AffineEncoding(Atilde_[row].compose(Mult_GF28_tab[Div_GF28(1, d_[row])]), q_[row]);
    }

    for (unsigned int row = 0; row < 4; row++) {
        Ptilde_[row] = EncodingKey(SBOX.compose(Mult_GF28_tab[diffD]).compose(SBOX_inv).compose(Ptilde_[row]));
    }
}

static bool getAndVerifyKey(const AffineEncoding& Q, const EncodingKey& Ptilde, uint8_t& key) {

    key = Ptilde[Q[0]];
    for (unsigned int x = 1; x < LEN_ARRAY; x++) {
        if (key != (Ptilde[Q[x]] ^ x)) {
            return false;
        }
    }
    return true;
}

static void finalizeQPtilde(
    QPtildeLinear &QPtilde0,
    QPtildeLinear &QPtilde1, uint8_t row) {

    if (QPtilde0.isFinalize() && QPtilde1.isFinalize()) {
        return;
    }

    for (unsigned delta0 = 0; delta0 < 255; delta0++) {
        for (unsigned delta1 = 0; delta1 < 255; delta1++) {
            uint8_t key_dummy;
            if (getAndVerifyKey(QPtilde0.getQ(row), QPtilde1.getPtilde(row), key_dummy)) {
                QPtilde0.setFinal();
                QPtilde1.setFinal();
                return;
            }

            if (QPtilde1.isFinalize() or delta1 == 254) {
                break;
            } else {
                QPtilde1.setNextDelta();
            }
        }
        if (QPtilde0.isFinalize() or delta0 == 254) {
            break;
        } else {
            QPtilde0.setNextDelta();
        }
    }
}

static void finalizeQPtildeDecrypt(
    std::array<QPtildeLinear, 4> &QPtilde0,
    std::array<QPtildeLinear, 4> &QPtilde1) {

  for (unsigned col = 0; col < 4; col++) {
      if (not (QPtilde0[col].isValid() and QPtilde1[col].isValid())) {
          return;
      }
  }

  for (unsigned row = 0; row < 4; row++) {
      finalizeQPtilde(QPtilde0[0], QPtilde1[ShiftRow[row][0]], row);
  }
  for (unsigned col = 1; col < 4; col++) {
      finalizeQPtilde(QPtilde0[col], QPtilde1[ShiftRow[0][col]], 0);
  }
}

#endif /* COMPUTEQPTILDE_HPP */
