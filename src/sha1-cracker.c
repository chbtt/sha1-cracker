/*
 * This file is part of <https://github.com/cbscorpion/sha1-cracker>.
 * Copyright (c) 2018 Christoph Buttler.
 * 
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */
#include "sha1-cracker.h"

// constants for hash state initialisation
#define SHA1_IV_0            UINT32_C(0x67452301)
#define SHA1_IV_1            UINT32_C(0xEFCDAB89)
#define SHA1_IV_2            UINT32_C(0x98BADCFE)
#define SHA1_IV_3            UINT32_C(0x10325476)
#define SHA1_IV_4            UINT32_C(0xC3D2E1F0)
// constants for k-values
#define K_00_19              UINT32_C(0x5A827999)
#define K_20_39              UINT32_C(0x6ED9EBA1)
#define K_40_59              UINT32_C(0x8F1BBCDC)
#define K_60_79              UINT32_C(0xCA62C1D6)
// constants for initial step optimization
#define ROUND_CONSTANT_00    UINT32_C(0x29FB498B3)
#define ROUND_CONSTANT_01    UINT32_C(0x166B0CD0D)
#define ROUND_CONSTANT_02    UINT32_C(0x0F33D5697)
#define ROUND_CONSTANT_03    UINT32_C(0x0D675E47B)
#define ROUND_CONSTANT_04    UINT32_C(0x0B453C259)
#define ROUND_CONSTANT_15    UINT32_C(0x05A8279C9)
// length of preimage (always 6 bytes)
#define PREIMAGE_LENGTH_BYTE                   6
#define PREIMAGE_LENGTH_BIT  UINT32_C(0x00000030)

// macros for f-functions
#define F_00_19(mB, mC, mD) (mD ^ (mB & (mC ^ mD)))
#define F_40_59(mB, mC, mD) ((mB & mC) ^ (mD & (mB ^ mC)))
#define F_REST(mB, mC, mD)  (mB ^ mC ^ mD)
// macro for circular left-shift of a 32-bit word (taken from RFC 3174)
#define LEFT_ROTATE(word, bits) (((word) << (bits)) | ((word) >> (32 - (bits))))
// macros for the different round additions
#define ROUND_ADDITION_00                           (ROUND_CONSTANT_00 + p_blocks[0])
#define ROUND_ADDITION_01(mA)                       (ROUND_CONSTANT_01 + LEFT_ROTATE(mA, 5) + p_blocks[1])
#define ROUND_ADDITION_02(mA, mB, mC, mD)           (ROUND_CONSTANT_02 + LEFT_ROTATE(mA, 5) + F_00_19(mB, mC, mD))
#define ROUND_ADDITION_03(mA, mB, mC, mD)           (ROUND_CONSTANT_03 + LEFT_ROTATE(mA, 5) + F_00_19(mB, mC, mD))
#define ROUND_ADDITION_04(mA, mB, mC, mD)           (ROUND_CONSTANT_04 + LEFT_ROTATE(mA, 5) + F_00_19(mB, mC, mD))
#define ROUND_ADDITION_05_14(mA, mB, mC, mD, mE)    (K_00_19 + mE + LEFT_ROTATE(mA, 5) + F_00_19(mB, mC, mD))
#define ROUND_ADDITION_15(mA, mB, mC, mD, mE)       (ROUND_CONSTANT_15 + mE + LEFT_ROTATE(mA, 5) + F_00_19(mB, mC, mD))
// macros for setting state variables in each round
#define ROUND_PROCESSING_START(mA, mB, mE, mRoundAddition) \
    mE = mRoundAddition;                                   \
    mB = LEFT_ROTATE(mB, 30);
#define ROUND_PROCESSING_END(mA, mB, mE, f, k, i)         \
    mE = k + mE + LEFT_ROTATE(mA, 5) + f + (p_blocks[i]); \
    mB = LEFT_ROTATE(mB, 30);
// macros for the round functions
#define ROUND_00_15(mA, mB, mC, mD, mE, mRoundAddition) \
    ROUND_PROCESSING_START(mA, mB, mE, mRoundAddition)
#define ROUND_16_19(mA, mB, mC, mD, mE, i) \
    ROUND_PROCESSING_END(mA, mB, mE, F_00_19(mB, mC, mD), K_00_19, i)
#define ROUND_20_39(mA, mB, mC, mD, mE, i) \
    ROUND_PROCESSING_END(mA, mB, mE, F_REST(mB, mC, mD), K_20_39, i)
#define ROUND_40_59(mA, mB, mC, mD, mE, i) \
    ROUND_PROCESSING_END(mA, mB, mE, F_40_59(mB, mC, mD), K_40_59, i)
#define ROUND_60_79(mA, mB, mC, mD, mE, i) \
    ROUND_PROCESSING_END(mA, mB, mE, F_REST(mB, mC, mD), K_60_79, i)

// function prototypes
static inline void precomputeOuterLoop(uint32_t *p_precomputedBlocks);
static inline void precomputeInnerLoop(uint32_t *p_precomputedBlocks,
                                       uint32_t *p_w0,
                                       uint32_t *p_blocks);
// skeleton for word blocks that already contains the correct padding
static const uint32_t p_blocksSkeleton[16] = {  0u, UINT32_C(0x8000), 0u,
                                                0u, 0u, 0u, 0u, 0u,   0u,
                                                0u, 0u, 0u, 0u, 0u,   0u,
                                                PREIMAGE_LENGTH_BIT };

/**
 * Function: crackHash
 */
int crackHash
(
    struct hash targetHash,
    char        *p_result
)
{
           char     p_currInput[6];
           uint32_t p_precomputedBlocks[80],
                    p_w0[21],
                    p_blocks[80],
                    p_earlyExit[5];
    static uint32_t a, b, c, d, e, temp;
    // precompute early exit values
    p_earlyExit[0] = targetHash.a - SHA1_IV_0;
    p_earlyExit[1] = targetHash.b - SHA1_IV_1;
    p_earlyExit[2] = targetHash.c - SHA1_IV_2;
    p_earlyExit[3] = targetHash.d - SHA1_IV_3;
    p_earlyExit[4] = targetHash.e - SHA1_IV_4;
    // values c, d, e are rotated by 30 - rotate by 2 more to undo
    p_earlyExit[2] = LEFT_ROTATE(p_earlyExit[2], 2);
    p_earlyExit[3] = LEFT_ROTATE(p_earlyExit[3], 2);
    p_earlyExit[4] = LEFT_ROTATE(p_earlyExit[4], 2);
    // outer loop through all two letter combinations from 'a' to 'z'
    for (p_currInput[4] = 'a'; p_currInput[4] <= 'z'; p_currInput[4]++)
        for (p_currInput[5] = 'a'; p_currInput[5] <= 'z'; p_currInput[5]++)
        {
            // reset precomputed word blocks
            memcpy(p_precomputedBlocks,
                   p_blocksSkeleton,
                   64);
            // set second word block based on input
            p_precomputedBlocks[1] |= (p_currInput[4] << 24) 
                                   |  (p_currInput[5] << 16);
            // precompute word blocks for outer loop
            precomputeOuterLoop(p_precomputedBlocks);
            // inner loop through all four letter combinations from 'a' to 'z'
            for (p_currInput[0] = 'a'; p_currInput[0] <= 'z'; p_currInput[0]++)
                for (p_currInput[1] = 'a'; p_currInput[1] <= 'z'; p_currInput[1]++)
                    for (p_currInput[2] = 'a'; p_currInput[2] <= 'z'; p_currInput[2]++)
                        for (p_currInput[3] = 'a'; p_currInput[3] <= 'z'; p_currInput[3]++)
                        {
                            // reset first 16 word blocks
                            memcpy(p_blocks,
                                   p_precomputedBlocks,
                                   64);
                            // generate (missing) first word block based on new input
                            p_blocks[0] =  (p_currInput[0] << 24)
                                        |  (p_currInput[1] << 16)
                                        |  (p_currInput[2] << 8)
                                        |  (p_currInput[3]);
                            // precompute word blocks for inner loop
                            precomputeInnerLoop(p_precomputedBlocks,
                                                p_w0,
                                                p_blocks);
                            // initialize state variables with constants
                            a = SHA1_IV_0;
                            b = SHA1_IV_1;
                            c = SHA1_IV_2;
                            d = SHA1_IV_3;
                            e = SHA1_IV_4;
                            /*************** UNROLLED ROUND FUNCTION LOOPS ****************/
                            // round 00
                            ROUND_00_15(a, b, c, d, e, ROUND_ADDITION_00)
                            // round 01
                            ROUND_00_15(e, a, b, c, d, ROUND_ADDITION_01(e))
                            // round 02
                            ROUND_00_15(d, e, a, b, c, ROUND_ADDITION_02(d, e, a, b))
                            // round 03
                            ROUND_00_15(c, d, e, a, b, ROUND_ADDITION_03(c, d, e, a))
                            // round 04
                            ROUND_00_15(b, c, d, e, a, ROUND_ADDITION_04(b, c, d, e))
                            // rounds 05 - 14
                            ROUND_00_15(a, b, c, d, e, ROUND_ADDITION_05_14(a, b, c, d, e))
                            ROUND_00_15(e, a, b, c, d, ROUND_ADDITION_05_14(e, a, b, c, d))
                            ROUND_00_15(d, e, a, b, c, ROUND_ADDITION_05_14(d, e, a, b, c))
                            ROUND_00_15(c, d, e, a, b, ROUND_ADDITION_05_14(c, d, e, a, b))
                            ROUND_00_15(b, c, d, e, a, ROUND_ADDITION_05_14(b, c, d, e, a))
                            ROUND_00_15(a, b, c, d, e, ROUND_ADDITION_05_14(a, b, c, d, e))
                            ROUND_00_15(e, a, b, c, d, ROUND_ADDITION_05_14(e, a, b, c, d))
                            ROUND_00_15(d, e, a, b, c, ROUND_ADDITION_05_14(d, e, a, b, c))
                            ROUND_00_15(c, d, e, a, b, ROUND_ADDITION_05_14(c, d, e, a, b))
                            ROUND_00_15(b, c, d, e, a, ROUND_ADDITION_05_14(b, c, d, e, a))
                            // round 15
                            ROUND_00_15(a, b, c, d, e, ROUND_ADDITION_15(a, b, c, d, e))
                            // rounds 16 - 19
                            ROUND_16_19(e, a, b, c, d, 16)
                            ROUND_16_19(d, e, a, b, c, 17)
                            ROUND_16_19(c, d, e, a, b, 18)
                            ROUND_16_19(b, c, d, e, a, 19)
                            // rounds 20 - 39
                            ROUND_20_39(a, b, c, d, e, 20)
                            ROUND_20_39(e, a, b, c, d, 21)
                            ROUND_20_39(d, e, a, b, c, 22)
                            ROUND_20_39(c, d, e, a, b, 23)
                            ROUND_20_39(b, c, d, e, a, 24)
                            ROUND_20_39(a, b, c, d, e, 25)
                            ROUND_20_39(e, a, b, c, d, 26)
                            ROUND_20_39(d, e, a, b, c, 27)
                            ROUND_20_39(c, d, e, a, b, 28)
                            ROUND_20_39(b, c, d, e, a, 29)
                            ROUND_20_39(a, b, c, d, e, 30)
                            ROUND_20_39(e, a, b, c, d, 31)
                            ROUND_20_39(d, e, a, b, c, 32)
                            ROUND_20_39(c, d, e, a, b, 33)
                            ROUND_20_39(b, c, d, e, a, 34)
                            ROUND_20_39(a, b, c, d, e, 35)
                            ROUND_20_39(e, a, b, c, d, 36)
                            ROUND_20_39(d, e, a, b, c, 37)
                            ROUND_20_39(c, d, e, a, b, 38)
                            ROUND_20_39(b, c, d, e, a, 39)
                            // rounds 40 - 59
                            ROUND_40_59(a, b, c, d, e, 40)
                            ROUND_40_59(e, a, b, c, d, 41)
                            ROUND_40_59(d, e, a, b, c, 42)
                            ROUND_40_59(c, d, e, a, b, 43)
                            ROUND_40_59(b, c, d, e, a, 44)
                            ROUND_40_59(a, b, c, d, e, 45)
                            ROUND_40_59(e, a, b, c, d, 46)
                            ROUND_40_59(d, e, a, b, c, 47)
                            ROUND_40_59(c, d, e, a, b, 48)
                            ROUND_40_59(b, c, d, e, a, 49)
                            ROUND_40_59(a, b, c, d, e, 50)
                            ROUND_40_59(e, a, b, c, d, 51)
                            ROUND_40_59(d, e, a, b, c, 52)
                            ROUND_40_59(c, d, e, a, b, 53)
                            ROUND_40_59(b, c, d, e, a, 54)
                            ROUND_40_59(a, b, c, d, e, 55)
                            ROUND_40_59(e, a, b, c, d, 56)
                            ROUND_40_59(d, e, a, b, c, 57)
                            ROUND_40_59(c, d, e, a, b, 58)
                            ROUND_40_59(b, c, d, e, a, 59)
                            // rounds 60 - 79
                            ROUND_60_79(a, b, c, d, e, 60)
                            ROUND_60_79(e, a, b, c, d, 61)
                            ROUND_60_79(d, e, a, b, c, 62)
                            ROUND_60_79(c, d, e, a, b, 63)
                            ROUND_60_79(b, c, d, e, a, 64)
                            ROUND_60_79(a, b, c, d, e, 65)
                            ROUND_60_79(e, a, b, c, d, 66)
                            ROUND_60_79(d, e, a, b, c, 67)
                            ROUND_60_79(c, d, e, a, b, 68)
                            ROUND_60_79(b, c, d, e, a, 69)
                            ROUND_60_79(a, b, c, d, e, 70)
                            ROUND_60_79(e, a, b, c, d, 71)
                            ROUND_60_79(d, e, a, b, c, 72)
                            ROUND_60_79(c, d, e, a, b, 73)
                            ROUND_60_79(b, c, d, e, a, 74)
                            /**************************************************************/
                            // round 75
                            ROUND_60_79(a, b, c, d, e, 75)
                            if(e != p_earlyExit[4])
                                continue;
                            // round 76
                            temp = p_blocks[73] ^ p_blocks[68] ^ p_blocks[62] ^ p_blocks[60];
                            p_blocks[76] = LEFT_ROTATE(temp, 1);
                            ROUND_60_79(e, a, b, c, d, 76)
                            if(d != p_earlyExit[3])
                                continue;
                            // round 77
                            temp = p_blocks[74] ^ p_blocks[69] ^ p_blocks[63] ^ p_blocks[61];
                            p_blocks[77] = LEFT_ROTATE(temp, 1);
                            ROUND_60_79(d, e, a, b, c, 77)
                            if(c != p_earlyExit[2])
                                continue;
                            // round 78
                            temp = p_blocks[75] ^ p_blocks[70] ^ p_blocks[64] ^ p_blocks[62];
                            p_blocks[78] = LEFT_ROTATE(temp, 1);
                            ROUND_60_79(c, d, e, a, b, 78)
                            if(b != p_earlyExit[1])
                                continue;
                            // round 79
                            temp = p_blocks[76] ^ p_blocks[71] ^ p_blocks[65] ^ p_blocks[63];
                            p_blocks[79] = LEFT_ROTATE(temp, 1);
                            ROUND_60_79(b, c, d, e, a, 79)
                            if(a != p_earlyExit[0])
                                continue;
                            // if we end up here, the target hash was found
                            goto success;
                        }
        }
    // if we end up here, no preimage was found
    return E_CRACK_NOT_FOUND;

success:
    // a preimage was found, copy it into p_result
    memcpy(p_result,
           p_currInput,
           PREIMAGE_LENGTH_BYTE);

    return 0;
}
/**
 * Function: precomputeOuterLoop
 */
static inline void precomputeOuterLoop
(
    uint32_t *p_precomputedBlocks
)
{
    //p_precomputedBlocks[16] = 0;
    p_precomputedBlocks[17] = p_precomputedBlocks[1];
    p_precomputedBlocks[17] = LEFT_ROTATE(p_precomputedBlocks[17], 1);
    p_precomputedBlocks[18] = p_precomputedBlocks[15];
    p_precomputedBlocks[18] = LEFT_ROTATE(p_precomputedBlocks[18], 1);
    //p_precomputedBlocks[19] = 0;
    p_precomputedBlocks[20] = p_precomputedBlocks[17];
    p_precomputedBlocks[20] = LEFT_ROTATE(p_precomputedBlocks[20], 1);
    p_precomputedBlocks[21] = p_precomputedBlocks[18];
    p_precomputedBlocks[21] = LEFT_ROTATE(p_precomputedBlocks[21], 1);
    //p_precomputedBlocks[22] = 0;
    p_precomputedBlocks[23] = p_precomputedBlocks[20] ^ p_precomputedBlocks[15];
    p_precomputedBlocks[23] = LEFT_ROTATE(p_precomputedBlocks[23], 1);
    p_precomputedBlocks[24] = p_precomputedBlocks[21];
    p_precomputedBlocks[24] = LEFT_ROTATE(p_precomputedBlocks[24], 1);
    p_precomputedBlocks[25] = p_precomputedBlocks[17];
    p_precomputedBlocks[25] = LEFT_ROTATE(p_precomputedBlocks[25], 1);
    p_precomputedBlocks[26] = p_precomputedBlocks[23] ^ p_precomputedBlocks[18];
    p_precomputedBlocks[26] = LEFT_ROTATE(p_precomputedBlocks[26], 1);
    p_precomputedBlocks[27] = p_precomputedBlocks[24];
    p_precomputedBlocks[27] = LEFT_ROTATE(p_precomputedBlocks[27], 1);
    p_precomputedBlocks[28] = p_precomputedBlocks[25] ^ p_precomputedBlocks[20];
    p_precomputedBlocks[28] = LEFT_ROTATE(p_precomputedBlocks[28], 1);
    p_precomputedBlocks[29] = p_precomputedBlocks[26] ^ p_precomputedBlocks[21] ^ p_precomputedBlocks[15];
    p_precomputedBlocks[29] = LEFT_ROTATE(p_precomputedBlocks[29], 1);
    p_precomputedBlocks[30] = p_precomputedBlocks[27];
    p_precomputedBlocks[30] = LEFT_ROTATE(p_precomputedBlocks[30], 1);
    p_precomputedBlocks[31] = p_precomputedBlocks[28] ^ p_precomputedBlocks[23] ^ p_precomputedBlocks[17] ^ p_precomputedBlocks[15];
    p_precomputedBlocks[31] = LEFT_ROTATE(p_precomputedBlocks[31], 1);
    p_precomputedBlocks[32] = p_precomputedBlocks[29] ^ p_precomputedBlocks[24] ^ p_precomputedBlocks[18];
    p_precomputedBlocks[32] = LEFT_ROTATE(p_precomputedBlocks[32], 1);
    p_precomputedBlocks[33] = p_precomputedBlocks[30] ^ p_precomputedBlocks[25] ^ p_precomputedBlocks[17];
    p_precomputedBlocks[33] = LEFT_ROTATE(p_precomputedBlocks[33], 1);
    p_precomputedBlocks[34] = p_precomputedBlocks[31] ^ p_precomputedBlocks[26] ^ p_precomputedBlocks[20] ^ p_precomputedBlocks[18];
    p_precomputedBlocks[34] = LEFT_ROTATE(p_precomputedBlocks[34], 1);
    p_precomputedBlocks[35] = p_precomputedBlocks[32] ^ p_precomputedBlocks[27] ^ p_precomputedBlocks[21];
    p_precomputedBlocks[35] = LEFT_ROTATE(p_precomputedBlocks[35], 1);
    p_precomputedBlocks[36] = p_precomputedBlocks[33] ^ p_precomputedBlocks[28] ^ p_precomputedBlocks[20];
    p_precomputedBlocks[36] = LEFT_ROTATE(p_precomputedBlocks[36], 1);
    p_precomputedBlocks[37] = p_precomputedBlocks[34] ^ p_precomputedBlocks[29] ^ p_precomputedBlocks[23] ^ p_precomputedBlocks[21];
    p_precomputedBlocks[37] = LEFT_ROTATE(p_precomputedBlocks[37], 1);
    p_precomputedBlocks[38] = p_precomputedBlocks[35] ^ p_precomputedBlocks[30] ^ p_precomputedBlocks[24];
    p_precomputedBlocks[38] = LEFT_ROTATE(p_precomputedBlocks[38], 1);
    p_precomputedBlocks[39] = p_precomputedBlocks[36] ^ p_precomputedBlocks[31] ^ p_precomputedBlocks[25] ^ p_precomputedBlocks[23];
    p_precomputedBlocks[39] = LEFT_ROTATE(p_precomputedBlocks[39], 1);
    p_precomputedBlocks[40] = p_precomputedBlocks[37] ^ p_precomputedBlocks[32] ^ p_precomputedBlocks[26] ^ p_precomputedBlocks[24];
    p_precomputedBlocks[40] = LEFT_ROTATE(p_precomputedBlocks[40], 1);
    p_precomputedBlocks[41] = p_precomputedBlocks[38] ^ p_precomputedBlocks[33] ^ p_precomputedBlocks[27] ^ p_precomputedBlocks[25];
    p_precomputedBlocks[41] = LEFT_ROTATE(p_precomputedBlocks[41], 1);
    p_precomputedBlocks[42] = p_precomputedBlocks[39] ^ p_precomputedBlocks[34] ^ p_precomputedBlocks[28] ^ p_precomputedBlocks[26];
    p_precomputedBlocks[42] = LEFT_ROTATE(p_precomputedBlocks[42], 1);
    p_precomputedBlocks[43] = p_precomputedBlocks[40] ^ p_precomputedBlocks[35] ^ p_precomputedBlocks[29] ^ p_precomputedBlocks[27];
    p_precomputedBlocks[43] = LEFT_ROTATE(p_precomputedBlocks[43], 1);
    p_precomputedBlocks[44] = p_precomputedBlocks[41] ^ p_precomputedBlocks[36] ^ p_precomputedBlocks[30] ^ p_precomputedBlocks[28];
    p_precomputedBlocks[44] = LEFT_ROTATE(p_precomputedBlocks[44], 1);
    p_precomputedBlocks[45] = p_precomputedBlocks[42] ^ p_precomputedBlocks[37] ^ p_precomputedBlocks[31] ^ p_precomputedBlocks[29];
    p_precomputedBlocks[45] = LEFT_ROTATE(p_precomputedBlocks[45], 1);
    p_precomputedBlocks[46] = p_precomputedBlocks[43] ^ p_precomputedBlocks[38] ^ p_precomputedBlocks[32] ^ p_precomputedBlocks[30];
    p_precomputedBlocks[46] = LEFT_ROTATE(p_precomputedBlocks[46], 1);
    p_precomputedBlocks[47] = p_precomputedBlocks[44] ^ p_precomputedBlocks[39] ^ p_precomputedBlocks[33] ^ p_precomputedBlocks[31];
    p_precomputedBlocks[47] = LEFT_ROTATE(p_precomputedBlocks[47], 1);
    p_precomputedBlocks[48] = p_precomputedBlocks[45] ^ p_precomputedBlocks[40] ^ p_precomputedBlocks[34] ^ p_precomputedBlocks[32];
    p_precomputedBlocks[48] = LEFT_ROTATE(p_precomputedBlocks[48], 1);
    p_precomputedBlocks[49] = p_precomputedBlocks[46] ^ p_precomputedBlocks[41] ^ p_precomputedBlocks[35] ^ p_precomputedBlocks[33];
    p_precomputedBlocks[49] = LEFT_ROTATE(p_precomputedBlocks[49], 1);
    p_precomputedBlocks[50] = p_precomputedBlocks[47] ^ p_precomputedBlocks[42] ^ p_precomputedBlocks[36] ^ p_precomputedBlocks[34];
    p_precomputedBlocks[50] = LEFT_ROTATE(p_precomputedBlocks[50], 1);
    p_precomputedBlocks[51] = p_precomputedBlocks[48] ^ p_precomputedBlocks[43] ^ p_precomputedBlocks[37] ^ p_precomputedBlocks[35];
    p_precomputedBlocks[51] = LEFT_ROTATE(p_precomputedBlocks[51], 1);
    p_precomputedBlocks[52] = p_precomputedBlocks[49] ^ p_precomputedBlocks[44] ^ p_precomputedBlocks[38] ^ p_precomputedBlocks[36];
    p_precomputedBlocks[52] = LEFT_ROTATE(p_precomputedBlocks[52], 1);
    p_precomputedBlocks[53] = p_precomputedBlocks[50] ^ p_precomputedBlocks[45] ^ p_precomputedBlocks[39] ^ p_precomputedBlocks[37];
    p_precomputedBlocks[53] = LEFT_ROTATE(p_precomputedBlocks[53], 1);
    p_precomputedBlocks[54] = p_precomputedBlocks[51] ^ p_precomputedBlocks[46] ^ p_precomputedBlocks[40] ^ p_precomputedBlocks[38];
    p_precomputedBlocks[54] = LEFT_ROTATE(p_precomputedBlocks[54], 1);
    p_precomputedBlocks[55] = p_precomputedBlocks[52] ^ p_precomputedBlocks[47] ^ p_precomputedBlocks[41] ^ p_precomputedBlocks[39];
    p_precomputedBlocks[55] = LEFT_ROTATE(p_precomputedBlocks[55], 1);
    p_precomputedBlocks[56] = p_precomputedBlocks[53] ^ p_precomputedBlocks[48] ^ p_precomputedBlocks[42] ^ p_precomputedBlocks[40];
    p_precomputedBlocks[56] = LEFT_ROTATE(p_precomputedBlocks[56], 1);
    p_precomputedBlocks[57] = p_precomputedBlocks[54] ^ p_precomputedBlocks[49] ^ p_precomputedBlocks[43] ^ p_precomputedBlocks[41];
    p_precomputedBlocks[57] = LEFT_ROTATE(p_precomputedBlocks[57], 1);
    p_precomputedBlocks[58] = p_precomputedBlocks[55] ^ p_precomputedBlocks[50] ^ p_precomputedBlocks[44] ^ p_precomputedBlocks[42];
    p_precomputedBlocks[58] = LEFT_ROTATE(p_precomputedBlocks[58], 1);
    p_precomputedBlocks[59] = p_precomputedBlocks[56] ^ p_precomputedBlocks[51] ^ p_precomputedBlocks[45] ^ p_precomputedBlocks[43];
    p_precomputedBlocks[59] = LEFT_ROTATE(p_precomputedBlocks[59], 1);
    p_precomputedBlocks[60] = p_precomputedBlocks[57] ^ p_precomputedBlocks[52] ^ p_precomputedBlocks[46] ^ p_precomputedBlocks[44];
    p_precomputedBlocks[60] = LEFT_ROTATE(p_precomputedBlocks[60], 1);
    p_precomputedBlocks[61] = p_precomputedBlocks[58] ^ p_precomputedBlocks[53] ^ p_precomputedBlocks[47] ^ p_precomputedBlocks[45];
    p_precomputedBlocks[61] = LEFT_ROTATE(p_precomputedBlocks[61], 1);
    p_precomputedBlocks[62] = p_precomputedBlocks[59] ^ p_precomputedBlocks[54] ^ p_precomputedBlocks[48] ^ p_precomputedBlocks[46];
    p_precomputedBlocks[62] = LEFT_ROTATE(p_precomputedBlocks[62], 1);
    p_precomputedBlocks[63] = p_precomputedBlocks[60] ^ p_precomputedBlocks[55] ^ p_precomputedBlocks[49] ^ p_precomputedBlocks[47];
    p_precomputedBlocks[63] = LEFT_ROTATE(p_precomputedBlocks[63], 1);
    p_precomputedBlocks[64] = p_precomputedBlocks[61] ^ p_precomputedBlocks[56] ^ p_precomputedBlocks[50] ^ p_precomputedBlocks[48];
    p_precomputedBlocks[64] = LEFT_ROTATE(p_precomputedBlocks[64], 1);
    p_precomputedBlocks[65] = p_precomputedBlocks[62] ^ p_precomputedBlocks[57] ^ p_precomputedBlocks[51] ^ p_precomputedBlocks[49];
    p_precomputedBlocks[65] = LEFT_ROTATE(p_precomputedBlocks[65], 1);
    p_precomputedBlocks[66] = p_precomputedBlocks[63] ^ p_precomputedBlocks[58] ^ p_precomputedBlocks[52] ^ p_precomputedBlocks[50];
    p_precomputedBlocks[66] = LEFT_ROTATE(p_precomputedBlocks[66], 1);
    p_precomputedBlocks[67] = p_precomputedBlocks[64] ^ p_precomputedBlocks[59] ^ p_precomputedBlocks[53] ^ p_precomputedBlocks[51];
    p_precomputedBlocks[67] = LEFT_ROTATE(p_precomputedBlocks[67], 1);
    p_precomputedBlocks[68] = p_precomputedBlocks[65] ^ p_precomputedBlocks[60] ^ p_precomputedBlocks[54] ^ p_precomputedBlocks[52];
    p_precomputedBlocks[68] = LEFT_ROTATE(p_precomputedBlocks[68], 1);
    p_precomputedBlocks[69] = p_precomputedBlocks[66] ^ p_precomputedBlocks[61] ^ p_precomputedBlocks[55] ^ p_precomputedBlocks[53];
    p_precomputedBlocks[69] = LEFT_ROTATE(p_precomputedBlocks[69], 1);
    p_precomputedBlocks[70] = p_precomputedBlocks[67] ^ p_precomputedBlocks[62] ^ p_precomputedBlocks[56] ^ p_precomputedBlocks[54];
    p_precomputedBlocks[70] = LEFT_ROTATE(p_precomputedBlocks[70], 1);
    p_precomputedBlocks[71] = p_precomputedBlocks[68] ^ p_precomputedBlocks[63] ^ p_precomputedBlocks[57] ^ p_precomputedBlocks[55];
    p_precomputedBlocks[71] = LEFT_ROTATE(p_precomputedBlocks[71], 1);
    p_precomputedBlocks[72] = p_precomputedBlocks[69] ^ p_precomputedBlocks[64] ^ p_precomputedBlocks[58] ^ p_precomputedBlocks[56];
    p_precomputedBlocks[72] = LEFT_ROTATE(p_precomputedBlocks[72], 1);
    p_precomputedBlocks[73] = p_precomputedBlocks[70] ^ p_precomputedBlocks[65] ^ p_precomputedBlocks[59] ^ p_precomputedBlocks[57];
    p_precomputedBlocks[73] = LEFT_ROTATE(p_precomputedBlocks[73], 1);
    p_precomputedBlocks[74] = p_precomputedBlocks[71] ^ p_precomputedBlocks[66] ^ p_precomputedBlocks[60] ^ p_precomputedBlocks[58];
    p_precomputedBlocks[74] = LEFT_ROTATE(p_precomputedBlocks[74], 1);
    p_precomputedBlocks[75] = p_precomputedBlocks[72] ^ p_precomputedBlocks[67] ^ p_precomputedBlocks[61] ^ p_precomputedBlocks[59];
    p_precomputedBlocks[75] = LEFT_ROTATE(p_precomputedBlocks[75], 1);
    p_precomputedBlocks[76] = p_precomputedBlocks[73] ^ p_precomputedBlocks[68] ^ p_precomputedBlocks[62] ^ p_precomputedBlocks[60];
    p_precomputedBlocks[76] = LEFT_ROTATE(p_precomputedBlocks[76], 1);
    p_precomputedBlocks[77] = p_precomputedBlocks[74] ^ p_precomputedBlocks[69] ^ p_precomputedBlocks[63] ^ p_precomputedBlocks[61];
    p_precomputedBlocks[77] = LEFT_ROTATE(p_precomputedBlocks[77], 1);
    p_precomputedBlocks[78] = p_precomputedBlocks[75] ^ p_precomputedBlocks[70] ^ p_precomputedBlocks[64] ^ p_precomputedBlocks[62];
    p_precomputedBlocks[78] = LEFT_ROTATE(p_precomputedBlocks[78], 1);
    p_precomputedBlocks[79] = p_precomputedBlocks[76] ^ p_precomputedBlocks[71] ^ p_precomputedBlocks[65] ^ p_precomputedBlocks[63];
    p_precomputedBlocks[79] = LEFT_ROTATE(p_precomputedBlocks[79], 1);
}
/**
 * Function: precomputeInnerLoop
 */
static inline void precomputeInnerLoop
(
    uint32_t *p_precomputedBlocks,
    uint32_t *p_w0,
    uint32_t *p_blocks
)
{
    // rotate w0
    p_w0[ 1] = LEFT_ROTATE(p_blocks[0],  1);
    p_w0[ 2] = LEFT_ROTATE(p_blocks[0],  2);
    p_w0[ 3] = LEFT_ROTATE(p_blocks[0],  3);
    p_w0[ 4] = LEFT_ROTATE(p_blocks[0],  4);
    p_w0[ 5] = LEFT_ROTATE(p_blocks[0],  5);
    p_w0[ 6] = LEFT_ROTATE(p_blocks[0],  6);
    p_w0[ 7] = LEFT_ROTATE(p_blocks[0],  7);
    p_w0[ 8] = LEFT_ROTATE(p_blocks[0],  8);
    p_w0[ 9] = LEFT_ROTATE(p_blocks[0],  9);
    p_w0[10] = LEFT_ROTATE(p_blocks[0], 10);
    p_w0[11] = LEFT_ROTATE(p_blocks[0], 11);
    p_w0[12] = LEFT_ROTATE(p_blocks[0], 12);
    p_w0[13] = LEFT_ROTATE(p_blocks[0], 13);
    p_w0[14] = LEFT_ROTATE(p_blocks[0], 14);
    p_w0[15] = LEFT_ROTATE(p_blocks[0], 15);
    p_w0[16] = LEFT_ROTATE(p_blocks[0], 16);
    p_w0[17] = LEFT_ROTATE(p_blocks[0], 17);
    p_w0[18] = LEFT_ROTATE(p_blocks[0], 18);
    p_w0[19] = LEFT_ROTATE(p_blocks[0], 19);
    p_w0[20] = LEFT_ROTATE(p_blocks[0], 20);
    // precompute word blocks (basis by Jens Steube)
    // note: slightly optimized for our special case
    p_blocks[16] = p_w0[1];
    p_blocks[17] = p_precomputedBlocks[17];
    p_blocks[18] = p_precomputedBlocks[18];
    p_blocks[19] = p_w0[2];
    p_blocks[20] = p_precomputedBlocks[20];
    p_blocks[21] = p_precomputedBlocks[21];
    p_blocks[22] = p_w0[3];
    p_blocks[23] = p_precomputedBlocks[23];
    p_blocks[24] = p_precomputedBlocks[24] ^ p_w0[2];
    p_blocks[25] = p_precomputedBlocks[25] ^ p_w0[4];
    p_blocks[26] = p_precomputedBlocks[26];
    p_blocks[27] = p_precomputedBlocks[27];
    p_blocks[28] = p_precomputedBlocks[28] ^ p_w0[5];
    p_blocks[29] = p_precomputedBlocks[29];
    p_blocks[30] = p_precomputedBlocks[30] ^ p_w0[4]  ^ p_w0[2];
    p_blocks[31] = p_precomputedBlocks[31] ^ p_w0[6];
    p_blocks[32] = p_precomputedBlocks[32] ^ p_w0[3]  ^ p_w0[2];
    p_blocks[33] = p_precomputedBlocks[33];
    p_blocks[34] = p_precomputedBlocks[34] ^ p_w0[7];
    p_blocks[35] = p_precomputedBlocks[35] ^ p_w0[4];
    p_blocks[36] = p_precomputedBlocks[36] ^ p_w0[6]  ^ p_w0[4];
    p_blocks[37] = p_precomputedBlocks[37] ^ p_w0[8];
    p_blocks[38] = p_precomputedBlocks[38] ^ p_w0[4];
    p_blocks[39] = p_precomputedBlocks[39];
    p_blocks[40] = p_precomputedBlocks[40] ^ p_w0[4]  ^ p_w0[9];
    p_blocks[41] = p_precomputedBlocks[41];
    p_blocks[42] = p_precomputedBlocks[42] ^ p_w0[6]  ^ p_w0[8];
    p_blocks[43] = p_precomputedBlocks[43] ^ p_w0[10];
    p_blocks[44] = p_precomputedBlocks[44] ^ p_w0[6]  ^ p_w0[3]  ^ p_w0[7];
    p_blocks[45] = p_precomputedBlocks[45];
    p_blocks[46] = p_precomputedBlocks[46] ^ p_w0[4]  ^ p_w0[11];
    p_blocks[47] = p_precomputedBlocks[47] ^ p_w0[8]  ^ p_w0[4];
    p_blocks[48] = p_precomputedBlocks[48] ^ p_w0[8]  ^ p_w0[4]  ^ p_w0[3]  ^ p_w0[10] ^ p_w0[5];
    p_blocks[49] = p_precomputedBlocks[49] ^ p_w0[12];
    p_blocks[50] = p_precomputedBlocks[50] ^ p_w0[8];
    p_blocks[51] = p_precomputedBlocks[51] ^ p_w0[6]  ^ p_w0[4];
    p_blocks[52] = p_precomputedBlocks[52] ^ p_w0[8]  ^ p_w0[4]  ^ p_w0[13];
    p_blocks[53] = p_precomputedBlocks[53];
    p_blocks[54] = p_precomputedBlocks[54] ^ p_w0[7]  ^ p_w0[10] ^ p_w0[12];
    p_blocks[55] = p_precomputedBlocks[55] ^ p_w0[14];
    p_blocks[56] = p_precomputedBlocks[56] ^ p_w0[6]  ^ p_w0[4]  ^ p_w0[11] ^ p_w0[7]  ^ p_w0[10];
    p_blocks[57] = p_precomputedBlocks[57] ^ p_w0[8];
    p_blocks[58] = p_precomputedBlocks[58] ^ p_w0[8]  ^ p_w0[4]  ^ p_w0[15];
    p_blocks[59] = p_precomputedBlocks[59] ^ p_w0[8]  ^ p_w0[12];
    p_blocks[60] = p_precomputedBlocks[60] ^ p_w0[8]  ^ p_w0[4]  ^ p_w0[7]  ^ p_w0[12] ^ p_w0[14];
    p_blocks[61] = p_precomputedBlocks[61] ^ p_w0[16];
    p_blocks[62] = p_precomputedBlocks[62] ^ p_w0[6]  ^ p_w0[12] ^ p_w0[8]  ^ p_w0[4];
    p_blocks[63] = p_precomputedBlocks[63] ^ p_w0[8];
    p_blocks[64] = p_precomputedBlocks[64] ^ p_w0[6]  ^ p_w0[7]  ^ p_w0[17] ^ p_w0[12] ^ p_w0[8]  ^ p_w0[4];
    p_blocks[65] = p_precomputedBlocks[65];
    p_blocks[66] = p_precomputedBlocks[66] ^ p_w0[14] ^ p_w0[16];
    p_blocks[67] = p_precomputedBlocks[67] ^ p_w0[8]  ^ p_w0[18];
    p_blocks[68] = p_precomputedBlocks[68] ^ p_w0[11] ^ p_w0[14] ^ p_w0[15];
    p_blocks[69] = p_precomputedBlocks[69];
    p_blocks[70] = p_precomputedBlocks[70] ^ p_w0[12] ^ p_w0[19];
    p_blocks[71] = p_precomputedBlocks[71] ^ p_w0[12] ^ p_w0[16];
    p_blocks[72] = p_precomputedBlocks[72] ^ p_w0[11] ^ p_w0[12] ^ p_w0[18] ^ p_w0[13] ^ p_w0[16] ^ p_w0[5];
    p_blocks[73] = p_precomputedBlocks[73] ^ p_w0[20];
    p_blocks[74] = p_precomputedBlocks[74] ^ p_w0[8]  ^ p_w0[16];
    p_blocks[75] = p_precomputedBlocks[75] ^ p_w0[6]  ^ p_w0[12] ^ p_w0[14];
}
