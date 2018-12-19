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
// SSE2 intrinsics
#include <emmintrin.h>

// macros for sse2 intrinsics
#define OR(x, y)            (_mm_or_si128(x, y))
#define XOR(x, y)           (_mm_xor_si128(x, y))
#define AND(x, y)           (_mm_and_si128(x, y))
#define ADD(x, y)           (_mm_add_epi32(x, y))
#define SET1INT(x)          (_mm_set1_epi32(x))
#define SET4INT(w, x, y, z) (_mm_set_epi32(w, x, y, z))
#define SETZERO             (_mm_setzero_si128())
#define STORE(loc, x)       (_mm_storeu_si128(((__m128i *) loc), x))

// constants for hash state initialisation
#define SHA1_IV_0            UINT32_C(0x67452301)
#define SHA1_IV_1            UINT32_C(0xEFCDAB89)
#define SHA1_IV_2            UINT32_C(0x98BADCFE)
#define SHA1_IV_3            UINT32_C(0x10325476)
#define SHA1_IV_4            UINT32_C(0xC3D2E1F0)
// length of preimage (always 6 bytes)
#define PREIMAGE_LENGTH_BYTE                          6
#define PREIMAGE_LENGTH_BIT         UINT32_C(0x00000030)
#define PREIMAGE_LENGTH_BIT_LROT1   UINT32_C(0x00000060)
#define PREIMAGE_LENGTH_BIT_LROT2   UINT32_C(0x000000C0)
#define PREIMAGE_LENGTH_BIT_LROT3   UINT32_C(0x00000180)
#define PREIMAGE_LENGTH_BIT_LROT4   UINT32_C(0x00000300)
#define PREIMAGE_LENGTH_BIT_LROT5   UINT32_C(0x00000600)

// macros for f-functions
#define VF_00_19(mB, mC, mD) (XOR(mD, (AND(mB, (XOR(mC, mD))))))
#define VF_40_59(mB, mC, mD) (XOR(AND(mB, mC), (AND(mD, (XOR(mB, mC))))))
#define VF_REST(mB, mC, mD)  (XOR((XOR(mB, mC)), mD))

// macros for circular left-shift (adapted from RFC 3174)
#define U32_LEFT_ROTATE(word, bits) (((word) << (bits)) | ((word) >> (32 - (bits))))
#define VEC_LEFT_ROTATE(word, bits) (OR((_mm_slli_epi32((word), (bits))), (_mm_srli_epi32((word), (32 - (bits))))))

// macros for the different round additions
#define ROUND_ADDITION_00                           (ADD(ROUND_CONSTANT_00, p_blocks[0]))
#define ROUND_ADDITION_01(mA)                       (ADD(ADD(ROUND_CONSTANT_01, VEC_LEFT_ROTATE(mA, 5)), p_blocks[1]))
#define ROUND_ADDITION_02(mA, mB, mC, mD)           (ADD(ADD(ROUND_CONSTANT_02, VEC_LEFT_ROTATE(mA, 5)), VF_00_19(mB, mC, mD)))
#define ROUND_ADDITION_03(mA, mB, mC, mD)           (ADD(ADD(ROUND_CONSTANT_03, VEC_LEFT_ROTATE(mA, 5)), VF_00_19(mB, mC, mD)))
#define ROUND_ADDITION_04(mA, mB, mC, mD)           (ADD(ADD(ROUND_CONSTANT_04, VEC_LEFT_ROTATE(mA, 5)), VF_00_19(mB, mC, mD)))
#define ROUND_ADDITION_05_14(mA, mB, mC, mD, mE)    (ADD(ADD(ADD(K_00_19, mE), VEC_LEFT_ROTATE(mA, 5)), VF_00_19(mB, mC, mD)))
#define ROUND_ADDITION_15(mA, mB, mC, mD, mE)       (ADD(ADD(ADD(ROUND_CONSTANT_15, mE), VEC_LEFT_ROTATE(mA, 5)), VF_00_19(mB, mC, mD)))

// macros for setting state variables in each round
#define ROUND_PROCESSING_START(mA, mB, mE, mRoundAddition) \
    mE = mRoundAddition;                                   \
    mB = VEC_LEFT_ROTATE(mB, 30);
#define ROUND_PROCESSING_END(mA, mB, mE, f, k, i)                             \
    mE = ADD(ADD(ADD(k, mE), VEC_LEFT_ROTATE(mA, 5)), ADD(f, (p_blocks[i]))); \
    mB = VEC_LEFT_ROTATE(mB, 30);

// macros for the round functions
#define ROUND_00_15(mA, mB, mC, mD, mE, mRoundAddition) \
    ROUND_PROCESSING_START(mA, mB, mE, mRoundAddition)
#define ROUND_16_19(mA, mB, mC, mD, mE, i) \
    ROUND_PROCESSING_END(mA, mB, mE, VF_00_19(mB, mC, mD), K_00_19, i)
#define ROUND_20_39(mA, mB, mC, mD, mE, i) \
    ROUND_PROCESSING_END(mA, mB, mE, VF_REST(mB, mC, mD), K_20_39, i)
#define ROUND_40_59(mA, mB, mC, mD, mE, i) \
    ROUND_PROCESSING_END(mA, mB, mE, VF_40_59(mB, mC, mD), K_40_59, i)
#define ROUND_60_79(mA, mB, mC, mD, mE, i) \
    ROUND_PROCESSING_END(mA, mB, mE, VF_REST(mB, mC, mD), K_60_79, i)

// function prototypes
static inline void precomputeOuterLoop(__m128i  *p_vecPrecomputedBlocks);
static inline void precomputeInnerLoop(__m128i  *p_vecPrecomputedBlocks,
                                       __m128i  *p_w0,
                                       __m128i  *p_blocks);

/**
 * Function: crackHash
 */
int crackHash
(
    struct hash targetHash,
    char        *p_result
)
{
    int      index;
    __m128i  a, b, c, d, e, vecTemp;
    char     p_currInput[6];
    uint32_t p_earlyExit[5],
             p_tempSave[4];
    __m128i  p_vecPrecomputedBlocks[80],
             p_w0[21],
             p_blocks[80];
    // constant vectors
    const __m128i   K_00_19 = SET1INT(0x5A827999),
                    K_20_39 = SET1INT(0x6ED9EBA1),
                    K_40_59 = SET1INT(0x8F1BBCDC),
                    K_60_79 = SET1INT(0xCA62C1D6),
                    ROUND_CONSTANT_00 = SET1INT(0x9FB498B3),
                    ROUND_CONSTANT_01 = SET1INT(0x66B0CD0D),
                    ROUND_CONSTANT_02 = SET1INT(0xF33D5697),
                    ROUND_CONSTANT_03 = SET1INT(0xD675E47B),
                    ROUND_CONSTANT_04 = SET1INT(0xB453C259),
                    ROUND_CONSTANT_15 = SET1INT(0x5A8279C9);
    /*************** SET CORRECT PADDING ONCE ****************/
    p_vecPrecomputedBlocks[ 2] = SETZERO; p_vecPrecomputedBlocks[ 3] = SETZERO;
    p_vecPrecomputedBlocks[ 4] = SETZERO; p_vecPrecomputedBlocks[ 5] = SETZERO;
    p_vecPrecomputedBlocks[ 6] = SETZERO; p_vecPrecomputedBlocks[ 7] = SETZERO;
    p_vecPrecomputedBlocks[ 8] = SETZERO; p_vecPrecomputedBlocks[ 9] = SETZERO;
    p_vecPrecomputedBlocks[10] = SETZERO; p_vecPrecomputedBlocks[11] = SETZERO;
    p_vecPrecomputedBlocks[12] = SETZERO; p_vecPrecomputedBlocks[13] = SETZERO;
    p_vecPrecomputedBlocks[14] = SETZERO;
    p_vecPrecomputedBlocks[15] = SET1INT(PREIMAGE_LENGTH_BIT);
    // the following values are only dependent on the preimage length
    p_vecPrecomputedBlocks[18] = SET1INT(PREIMAGE_LENGTH_BIT_LROT1);
    p_vecPrecomputedBlocks[21] = SET1INT(PREIMAGE_LENGTH_BIT_LROT2);
    p_vecPrecomputedBlocks[24] = SET1INT(PREIMAGE_LENGTH_BIT_LROT3);
    p_vecPrecomputedBlocks[27] = SET1INT(PREIMAGE_LENGTH_BIT_LROT4);
    p_vecPrecomputedBlocks[30] = SET1INT(PREIMAGE_LENGTH_BIT_LROT5);
    // initialize p_blocks with padding too
    memcpy( p_blocks,
            p_vecPrecomputedBlocks,
            16 * sizeof(__m128i));
    /**************************************************************/
    /*************** PRECOMPUTE EARLY EXIT VALUES ****************/
    p_earlyExit[0] = targetHash.a - SHA1_IV_0;
    p_earlyExit[1] = targetHash.b - SHA1_IV_1;
    p_earlyExit[2] = targetHash.c - SHA1_IV_2;
    p_earlyExit[3] = targetHash.d - SHA1_IV_3;
    p_earlyExit[4] = targetHash.e - SHA1_IV_4;
    // values c, d, e are rotated by 30 - rotate by 2 more to undo
    p_earlyExit[2] = U32_LEFT_ROTATE(p_earlyExit[2], 2);
    p_earlyExit[3] = U32_LEFT_ROTATE(p_earlyExit[3], 2);
    p_earlyExit[4] = U32_LEFT_ROTATE(p_earlyExit[4], 2);
    /**************************************************************/
    // outer loop through all two letter combinations from 'a' to 'z'
    for (p_currInput[4] = 'a'; p_currInput[4] <= 'z'; p_currInput[4]++)
        for (p_currInput[5] = 'a'; p_currInput[5] <= 'z'; p_currInput[5]++)
        {
            // set second word block based on input ('1'-bit == 0x8000)
            p_vecPrecomputedBlocks[1] = SET1INT(  UINT32_C(0x8000)
                                                | (p_currInput[4] << 24)
                                                | (p_currInput[5] << 16));
            p_blocks[1] = p_vecPrecomputedBlocks[1];
            // precompute word blocks for outer loop
            precomputeOuterLoop(p_vecPrecomputedBlocks);
            // inner loop through all four letter combinations from 'a' to 'z'
            for (p_currInput[0] = 'a'; p_currInput[0] <= 'z'; p_currInput[0]++)
                for (p_currInput[1] = 'a'; p_currInput[1] <= 'z'; p_currInput[1]++)
                    // avoid "26/4"-problem by covering two loop runs each in third and fourth inner loop
                    for (p_currInput[2] = 'a'; p_currInput[2] <= 'z'; p_currInput[2] += 2)
                        for (p_currInput[3] = 'a'; p_currInput[3] <= 'z'; p_currInput[3] += 2)
                        {
                            // generate (missing) first word blocks based on new input
                            p_tempSave[0] = (p_currInput[0] << 24)
                                          | (p_currInput[1] << 16)
                                          | (p_currInput[2] << 8)
                                          | (p_currInput[3]);
                            p_tempSave[1] = (p_currInput[0] << 24)
                                          | (p_currInput[1] << 16)
                                          | (p_currInput[2] << 8)
                                          | ((p_currInput[3] + 1));
                            p_tempSave[2] = (p_currInput[0] << 24)
                                          | (p_currInput[1] << 16)
                                          | ((p_currInput[2] + 1) << 8)
                                          | (p_currInput[3]);
                            p_tempSave[3] = (p_currInput[0] << 24)
                                          | (p_currInput[1] << 16)
                                          | ((p_currInput[2] + 1) << 8)
                                          | ((p_currInput[3] + 1));
                            p_blocks[0] = SET4INT(p_tempSave[0],
                                                  p_tempSave[1],
                                                  p_tempSave[2],
                                                  p_tempSave[3]);
                            // precompute word blocks for inner loop
                            precomputeInnerLoop(p_vecPrecomputedBlocks,
                                                p_w0,
                                                p_blocks);
                            // initialize state variables with constants
                            a = SET1INT(SHA1_IV_0);
                            b = SET1INT(SHA1_IV_1);
                            c = SET1INT(SHA1_IV_2);
                            d = SET1INT(SHA1_IV_3);
                            e = SET1INT(SHA1_IV_4);
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
                            /***************** EARLY EXIT OPTIMIZATION *******************/
                            // round 75
                            ROUND_60_79(a, b, c, d, e, 75)
                            STORE(p_tempSave, e);
                            /*  
                             *  In our way of iterating over the searchspace, there are never two            
                             *  values in a vector that are the same in this step. Therefore, the following  
                             *  comparisons can yield at most one match. That said, we can set "index"       
                             *  to whatever match we find and check only that specific value going forward.
                             */
                            index = -1;
                            if((p_tempSave[3]) == p_earlyExit[4]) index = 3;
                            if((p_tempSave[2]) == p_earlyExit[4]) index = 2;
                            if((p_tempSave[1]) == p_earlyExit[4]) index = 1;
                            if((p_tempSave[0]) == p_earlyExit[4]) index = 0;
                            if(index == -1) 
                                continue;
                            /*
                             *  The word blocks have only been precomputed up until round 75, 
                             *  because that's the earliest we can exit. If we go beyond round 
                             *  75, compute missing blocks when required.
                             */ 
                            // round 76
                            vecTemp = XOR(XOR(p_blocks[73], p_blocks[68]), XOR(p_blocks[62], p_blocks[60]));
                            p_blocks[76] = VEC_LEFT_ROTATE(vecTemp, 1);
                            ROUND_60_79(e, a, b, c, d, 76)
                            STORE(p_tempSave, d);
                            if(p_tempSave[index] != p_earlyExit[3]) 
                                continue;
                            // round 77
                            vecTemp = XOR(XOR(p_blocks[74], p_blocks[69]), XOR(p_blocks[63], p_blocks[61]));
                            p_blocks[77] = VEC_LEFT_ROTATE(vecTemp, 1);
                            ROUND_60_79(d, e, a, b, c, 77)
                            STORE(p_tempSave, c);
                            if(p_tempSave[index] != p_earlyExit[2])
                                continue;
                            // round 78
                            vecTemp = XOR(XOR(p_blocks[75], p_blocks[70]), XOR(p_blocks[64], p_blocks[62]));
                            p_blocks[78] = VEC_LEFT_ROTATE(vecTemp, 1);
                            ROUND_60_79(c, d, e, a, b, 78)
                            STORE(p_tempSave, b);
                            if(p_tempSave[index] != p_earlyExit[1])
                                continue;
                            // round 79
                            vecTemp = XOR(XOR(p_blocks[76], p_blocks[71]), XOR(p_blocks[65], p_blocks[63]));
                            p_blocks[79] = VEC_LEFT_ROTATE(vecTemp, 1);
                            ROUND_60_79(b, c, d, e, a, 79)
                            STORE(p_tempSave, a);
                            if(p_tempSave[index] != p_earlyExit[0])
                                continue;
                            /**************************************************************/
                            // set correct preimage (no changes needed for index == 3)
                            switch(index)
                            {
                                case 2: p_currInput[3]++; 
                                        break;
                                case 1: p_currInput[2]++; 
                                        break;
                                case 0: p_currInput[2]++;
                                        p_currInput[3]++;
                            }

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
    __m128i  *p_vecPrecomputedBlocks
)
{
    //p_vecPrecomputedBlocks[16] = 0;
    p_vecPrecomputedBlocks[17] = p_vecPrecomputedBlocks[1];
    p_vecPrecomputedBlocks[17] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[17], 1);
    //p_vecPrecomputedBlocks[18] = SET1INT(PREIMAGE_LENGTH_BIT_LROT1);
    //p_vecPrecomputedBlocks[19] = 0;
    p_vecPrecomputedBlocks[20] = p_vecPrecomputedBlocks[17];
    p_vecPrecomputedBlocks[20] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[20], 1);
    //p_vecPrecomputedBlocks[21] = SET1INT(PREIMAGE_LENGTH_BIT_LROT2);
    //p_vecPrecomputedBlocks[22] = 0;
    p_vecPrecomputedBlocks[23] = XOR(p_vecPrecomputedBlocks[20], p_vecPrecomputedBlocks[15]);
    p_vecPrecomputedBlocks[23] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[23], 1);
    //p_vecPrecomputedBlocks[24] = SET1INT(PREIMAGE_LENGTH_BIT_LROT3);
    p_vecPrecomputedBlocks[25] = p_vecPrecomputedBlocks[17];
    p_vecPrecomputedBlocks[25] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[25], 1);
    p_vecPrecomputedBlocks[26] = XOR(p_vecPrecomputedBlocks[23], p_vecPrecomputedBlocks[18]);
    p_vecPrecomputedBlocks[26] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[26], 1);
    //p_vecPrecomputedBlocks[27] = SET1INT(PREIMAGE_LENGTH_BIT_LROT4);
    p_vecPrecomputedBlocks[28] = XOR(p_vecPrecomputedBlocks[25], p_vecPrecomputedBlocks[20]);
    p_vecPrecomputedBlocks[28] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[28], 1);
    p_vecPrecomputedBlocks[29] = XOR(XOR(p_vecPrecomputedBlocks[26], p_vecPrecomputedBlocks[21]), p_vecPrecomputedBlocks[15]);
    p_vecPrecomputedBlocks[29] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[29], 1);
    //p_vecPrecomputedBlocks[30] = SET1INT(PREIMAGE_LENGTH_BIT_LROT5);
    p_vecPrecomputedBlocks[31] = XOR(XOR(p_vecPrecomputedBlocks[28], p_vecPrecomputedBlocks[23]), XOR(p_vecPrecomputedBlocks[17], p_vecPrecomputedBlocks[15]));
    p_vecPrecomputedBlocks[31] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[31], 1);
    p_vecPrecomputedBlocks[32] = XOR(XOR(p_vecPrecomputedBlocks[29], p_vecPrecomputedBlocks[24]), p_vecPrecomputedBlocks[18]);
    p_vecPrecomputedBlocks[32] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[32], 1);
    p_vecPrecomputedBlocks[33] = XOR(XOR(p_vecPrecomputedBlocks[30], p_vecPrecomputedBlocks[25]), p_vecPrecomputedBlocks[17]);
    p_vecPrecomputedBlocks[33] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[33], 1);
    p_vecPrecomputedBlocks[34] = XOR(XOR(p_vecPrecomputedBlocks[31], p_vecPrecomputedBlocks[26]), XOR(p_vecPrecomputedBlocks[20], p_vecPrecomputedBlocks[18]));
    p_vecPrecomputedBlocks[34] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[34], 1);
    p_vecPrecomputedBlocks[35] = XOR(XOR(p_vecPrecomputedBlocks[32], p_vecPrecomputedBlocks[27]), p_vecPrecomputedBlocks[21]);
    p_vecPrecomputedBlocks[35] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[35], 1);
    p_vecPrecomputedBlocks[36] = XOR(XOR(p_vecPrecomputedBlocks[33], p_vecPrecomputedBlocks[28]), p_vecPrecomputedBlocks[20]);
    p_vecPrecomputedBlocks[36] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[36], 1);
    p_vecPrecomputedBlocks[37] = XOR(XOR(p_vecPrecomputedBlocks[34], p_vecPrecomputedBlocks[29]), XOR(p_vecPrecomputedBlocks[23], p_vecPrecomputedBlocks[21]));
    p_vecPrecomputedBlocks[37] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[37], 1);
    p_vecPrecomputedBlocks[38] = XOR(XOR(p_vecPrecomputedBlocks[35], p_vecPrecomputedBlocks[30]), p_vecPrecomputedBlocks[24]);
    p_vecPrecomputedBlocks[38] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[38], 1);
    p_vecPrecomputedBlocks[39] = XOR(XOR(p_vecPrecomputedBlocks[36], p_vecPrecomputedBlocks[31]), XOR(p_vecPrecomputedBlocks[25], p_vecPrecomputedBlocks[23]));
    p_vecPrecomputedBlocks[39] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[39], 1);
    p_vecPrecomputedBlocks[40] = XOR(XOR(p_vecPrecomputedBlocks[37], p_vecPrecomputedBlocks[32]), XOR(p_vecPrecomputedBlocks[26], p_vecPrecomputedBlocks[24]));
    p_vecPrecomputedBlocks[40] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[40], 1);
    p_vecPrecomputedBlocks[41] = XOR(XOR(p_vecPrecomputedBlocks[38], p_vecPrecomputedBlocks[33]), XOR(p_vecPrecomputedBlocks[27], p_vecPrecomputedBlocks[25]));
    p_vecPrecomputedBlocks[41] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[41], 1);
    p_vecPrecomputedBlocks[42] = XOR(XOR(p_vecPrecomputedBlocks[39], p_vecPrecomputedBlocks[34]), XOR(p_vecPrecomputedBlocks[28], p_vecPrecomputedBlocks[26]));
    p_vecPrecomputedBlocks[42] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[42], 1);
    p_vecPrecomputedBlocks[43] = XOR(XOR(p_vecPrecomputedBlocks[40], p_vecPrecomputedBlocks[35]), XOR(p_vecPrecomputedBlocks[29], p_vecPrecomputedBlocks[27]));
    p_vecPrecomputedBlocks[43] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[43], 1);
    p_vecPrecomputedBlocks[44] = XOR(XOR(p_vecPrecomputedBlocks[41], p_vecPrecomputedBlocks[36]), XOR(p_vecPrecomputedBlocks[30], p_vecPrecomputedBlocks[28]));
    p_vecPrecomputedBlocks[44] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[44], 1);
    p_vecPrecomputedBlocks[45] = XOR(XOR(p_vecPrecomputedBlocks[42], p_vecPrecomputedBlocks[37]), XOR(p_vecPrecomputedBlocks[31], p_vecPrecomputedBlocks[29]));
    p_vecPrecomputedBlocks[45] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[45], 1);
    p_vecPrecomputedBlocks[46] = XOR(XOR(p_vecPrecomputedBlocks[43], p_vecPrecomputedBlocks[38]), XOR(p_vecPrecomputedBlocks[32], p_vecPrecomputedBlocks[30]));
    p_vecPrecomputedBlocks[46] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[46], 1);
    p_vecPrecomputedBlocks[47] = XOR(XOR(p_vecPrecomputedBlocks[44], p_vecPrecomputedBlocks[39]), XOR(p_vecPrecomputedBlocks[33], p_vecPrecomputedBlocks[31]));
    p_vecPrecomputedBlocks[47] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[47], 1);
    p_vecPrecomputedBlocks[48] = XOR(XOR(p_vecPrecomputedBlocks[45], p_vecPrecomputedBlocks[40]), XOR(p_vecPrecomputedBlocks[34], p_vecPrecomputedBlocks[32]));
    p_vecPrecomputedBlocks[48] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[48], 1);
    p_vecPrecomputedBlocks[49] = XOR(XOR(p_vecPrecomputedBlocks[46], p_vecPrecomputedBlocks[41]), XOR(p_vecPrecomputedBlocks[35], p_vecPrecomputedBlocks[33]));
    p_vecPrecomputedBlocks[49] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[49], 1);
    p_vecPrecomputedBlocks[50] = XOR(XOR(p_vecPrecomputedBlocks[47], p_vecPrecomputedBlocks[42]), XOR(p_vecPrecomputedBlocks[36], p_vecPrecomputedBlocks[34]));
    p_vecPrecomputedBlocks[50] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[50], 1);
    p_vecPrecomputedBlocks[51] = XOR(XOR(p_vecPrecomputedBlocks[48], p_vecPrecomputedBlocks[43]), XOR(p_vecPrecomputedBlocks[37], p_vecPrecomputedBlocks[35]));
    p_vecPrecomputedBlocks[51] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[51], 1);
    p_vecPrecomputedBlocks[52] = XOR(XOR(p_vecPrecomputedBlocks[49], p_vecPrecomputedBlocks[44]), XOR(p_vecPrecomputedBlocks[38], p_vecPrecomputedBlocks[36]));
    p_vecPrecomputedBlocks[52] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[52], 1);
    p_vecPrecomputedBlocks[53] = XOR(XOR(p_vecPrecomputedBlocks[50], p_vecPrecomputedBlocks[45]), XOR(p_vecPrecomputedBlocks[39], p_vecPrecomputedBlocks[37]));
    p_vecPrecomputedBlocks[53] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[53], 1);
    p_vecPrecomputedBlocks[54] = XOR(XOR(p_vecPrecomputedBlocks[51], p_vecPrecomputedBlocks[46]), XOR(p_vecPrecomputedBlocks[40], p_vecPrecomputedBlocks[38]));
    p_vecPrecomputedBlocks[54] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[54], 1);
    p_vecPrecomputedBlocks[55] = XOR(XOR(p_vecPrecomputedBlocks[52], p_vecPrecomputedBlocks[47]), XOR(p_vecPrecomputedBlocks[41], p_vecPrecomputedBlocks[39]));
    p_vecPrecomputedBlocks[55] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[55], 1);
    p_vecPrecomputedBlocks[56] = XOR(XOR(p_vecPrecomputedBlocks[53], p_vecPrecomputedBlocks[48]), XOR(p_vecPrecomputedBlocks[42], p_vecPrecomputedBlocks[40]));
    p_vecPrecomputedBlocks[56] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[56], 1);
    p_vecPrecomputedBlocks[57] = XOR(XOR(p_vecPrecomputedBlocks[54], p_vecPrecomputedBlocks[49]), XOR(p_vecPrecomputedBlocks[43], p_vecPrecomputedBlocks[41]));
    p_vecPrecomputedBlocks[57] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[57], 1);
    p_vecPrecomputedBlocks[58] = XOR(XOR(p_vecPrecomputedBlocks[55], p_vecPrecomputedBlocks[50]), XOR(p_vecPrecomputedBlocks[44], p_vecPrecomputedBlocks[42]));
    p_vecPrecomputedBlocks[58] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[58], 1);
    p_vecPrecomputedBlocks[59] = XOR(XOR(p_vecPrecomputedBlocks[56], p_vecPrecomputedBlocks[51]), XOR(p_vecPrecomputedBlocks[45], p_vecPrecomputedBlocks[43]));
    p_vecPrecomputedBlocks[59] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[59], 1);
    p_vecPrecomputedBlocks[60] = XOR(XOR(p_vecPrecomputedBlocks[57], p_vecPrecomputedBlocks[52]), XOR(p_vecPrecomputedBlocks[46], p_vecPrecomputedBlocks[44]));
    p_vecPrecomputedBlocks[60] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[60], 1);
    p_vecPrecomputedBlocks[61] = XOR(XOR(p_vecPrecomputedBlocks[58], p_vecPrecomputedBlocks[53]), XOR(p_vecPrecomputedBlocks[47], p_vecPrecomputedBlocks[45]));
    p_vecPrecomputedBlocks[61] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[61], 1);
    p_vecPrecomputedBlocks[62] = XOR(XOR(p_vecPrecomputedBlocks[59], p_vecPrecomputedBlocks[54]), XOR(p_vecPrecomputedBlocks[48], p_vecPrecomputedBlocks[46]));
    p_vecPrecomputedBlocks[62] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[62], 1);
    p_vecPrecomputedBlocks[63] = XOR(XOR(p_vecPrecomputedBlocks[60], p_vecPrecomputedBlocks[55]), XOR(p_vecPrecomputedBlocks[49], p_vecPrecomputedBlocks[47]));
    p_vecPrecomputedBlocks[63] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[63], 1);
    p_vecPrecomputedBlocks[64] = XOR(XOR(p_vecPrecomputedBlocks[61], p_vecPrecomputedBlocks[56]), XOR(p_vecPrecomputedBlocks[50], p_vecPrecomputedBlocks[48]));
    p_vecPrecomputedBlocks[64] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[64], 1);
    p_vecPrecomputedBlocks[65] = XOR(XOR(p_vecPrecomputedBlocks[62], p_vecPrecomputedBlocks[57]), XOR(p_vecPrecomputedBlocks[51], p_vecPrecomputedBlocks[49]));
    p_vecPrecomputedBlocks[65] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[65], 1);
    p_vecPrecomputedBlocks[66] = XOR(XOR(p_vecPrecomputedBlocks[63], p_vecPrecomputedBlocks[58]), XOR(p_vecPrecomputedBlocks[52], p_vecPrecomputedBlocks[50]));
    p_vecPrecomputedBlocks[66] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[66], 1);
    p_vecPrecomputedBlocks[67] = XOR(XOR(p_vecPrecomputedBlocks[64], p_vecPrecomputedBlocks[59]), XOR(p_vecPrecomputedBlocks[53], p_vecPrecomputedBlocks[51]));
    p_vecPrecomputedBlocks[67] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[67], 1);
    p_vecPrecomputedBlocks[68] = XOR(XOR(p_vecPrecomputedBlocks[65], p_vecPrecomputedBlocks[60]), XOR(p_vecPrecomputedBlocks[54], p_vecPrecomputedBlocks[52]));
    p_vecPrecomputedBlocks[68] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[68], 1);
    p_vecPrecomputedBlocks[69] = XOR(XOR(p_vecPrecomputedBlocks[66], p_vecPrecomputedBlocks[61]), XOR(p_vecPrecomputedBlocks[55], p_vecPrecomputedBlocks[53]));
    p_vecPrecomputedBlocks[69] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[69], 1);
    p_vecPrecomputedBlocks[70] = XOR(XOR(p_vecPrecomputedBlocks[67], p_vecPrecomputedBlocks[62]), XOR(p_vecPrecomputedBlocks[56], p_vecPrecomputedBlocks[54]));
    p_vecPrecomputedBlocks[70] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[70], 1);
    p_vecPrecomputedBlocks[71] = XOR(XOR(p_vecPrecomputedBlocks[68], p_vecPrecomputedBlocks[63]), XOR(p_vecPrecomputedBlocks[57], p_vecPrecomputedBlocks[55]));
    p_vecPrecomputedBlocks[71] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[71], 1);
    p_vecPrecomputedBlocks[72] = XOR(XOR(p_vecPrecomputedBlocks[69], p_vecPrecomputedBlocks[64]), XOR(p_vecPrecomputedBlocks[58], p_vecPrecomputedBlocks[56]));
    p_vecPrecomputedBlocks[72] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[72], 1);
    p_vecPrecomputedBlocks[73] = XOR(XOR(p_vecPrecomputedBlocks[70], p_vecPrecomputedBlocks[65]), XOR(p_vecPrecomputedBlocks[59], p_vecPrecomputedBlocks[57]));
    p_vecPrecomputedBlocks[73] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[73], 1);
    p_vecPrecomputedBlocks[74] = XOR(XOR(p_vecPrecomputedBlocks[71], p_vecPrecomputedBlocks[66]), XOR(p_vecPrecomputedBlocks[60], p_vecPrecomputedBlocks[58]));
    p_vecPrecomputedBlocks[74] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[74], 1);
    p_vecPrecomputedBlocks[75] = XOR(XOR(p_vecPrecomputedBlocks[72], p_vecPrecomputedBlocks[67]), XOR(p_vecPrecomputedBlocks[61], p_vecPrecomputedBlocks[59]));
    p_vecPrecomputedBlocks[75] = VEC_LEFT_ROTATE(p_vecPrecomputedBlocks[75], 1);
}
/**
 * Function: precomputeInnerLoop
 */
static inline void precomputeInnerLoop
(
    __m128i *p_vecPrecomputedBlocks,
    __m128i *p_w0,
    __m128i *p_blocks
)
{
    // rotate w0
    p_w0[ 1] = VEC_LEFT_ROTATE(p_blocks[0],  1);
    p_w0[ 2] = VEC_LEFT_ROTATE(p_blocks[0],  2);
    p_w0[ 3] = VEC_LEFT_ROTATE(p_blocks[0],  3);
    p_w0[ 4] = VEC_LEFT_ROTATE(p_blocks[0],  4);
    p_w0[ 5] = VEC_LEFT_ROTATE(p_blocks[0],  5);
    p_w0[ 6] = VEC_LEFT_ROTATE(p_blocks[0],  6);
    p_w0[ 7] = VEC_LEFT_ROTATE(p_blocks[0],  7);
    p_w0[ 8] = VEC_LEFT_ROTATE(p_blocks[0],  8);
    p_w0[ 9] = VEC_LEFT_ROTATE(p_blocks[0],  9);
    p_w0[10] = VEC_LEFT_ROTATE(p_blocks[0], 10);
    p_w0[11] = VEC_LEFT_ROTATE(p_blocks[0], 11);
    p_w0[12] = VEC_LEFT_ROTATE(p_blocks[0], 12);
    p_w0[13] = VEC_LEFT_ROTATE(p_blocks[0], 13);
    p_w0[14] = VEC_LEFT_ROTATE(p_blocks[0], 14);
    p_w0[15] = VEC_LEFT_ROTATE(p_blocks[0], 15);
    p_w0[16] = VEC_LEFT_ROTATE(p_blocks[0], 16);
    p_w0[17] = VEC_LEFT_ROTATE(p_blocks[0], 17);
    p_w0[18] = VEC_LEFT_ROTATE(p_blocks[0], 18);
    p_w0[19] = VEC_LEFT_ROTATE(p_blocks[0], 19);
    p_w0[20] = VEC_LEFT_ROTATE(p_blocks[0], 20);
    // precompute word blocks for inner loop (basis by Jens Steube)
    // note: slightly optimized for our special case of sha1-cracking
    p_blocks[16] = p_w0[1];
    p_blocks[17] = p_vecPrecomputedBlocks[17];
    p_blocks[18] = p_vecPrecomputedBlocks[18];
    p_blocks[19] = p_w0[2];
    p_blocks[20] = p_vecPrecomputedBlocks[20];
    p_blocks[21] = p_vecPrecomputedBlocks[21];
    p_blocks[22] = p_w0[3];
    p_blocks[23] = p_vecPrecomputedBlocks[23];
    p_blocks[24] = XOR(p_vecPrecomputedBlocks[24], p_w0[2]);
    p_blocks[25] = XOR(p_vecPrecomputedBlocks[25], p_w0[4]);
    p_blocks[26] = p_vecPrecomputedBlocks[26];
    p_blocks[27] = p_vecPrecomputedBlocks[27];
    p_blocks[28] = XOR(p_vecPrecomputedBlocks[28], p_w0[5]);
    p_blocks[29] = p_vecPrecomputedBlocks[29];
    p_blocks[30] = XOR(XOR(p_vecPrecomputedBlocks[30], p_w0[4]), p_w0[2]);
    p_blocks[31] = XOR(p_vecPrecomputedBlocks[31], p_w0[6]);
    p_blocks[32] = XOR(XOR(p_vecPrecomputedBlocks[32], p_w0[3]), p_w0[2]);
    p_blocks[33] = p_vecPrecomputedBlocks[33];
    p_blocks[34] = XOR(p_vecPrecomputedBlocks[34], p_w0[7]);
    p_blocks[35] = XOR(p_vecPrecomputedBlocks[35], p_w0[4]);
    p_blocks[36] = XOR(XOR(p_vecPrecomputedBlocks[36], p_w0[6]), p_w0[4]);
    p_blocks[37] = XOR(p_vecPrecomputedBlocks[37], p_w0[8]);
    p_blocks[38] = XOR(p_vecPrecomputedBlocks[38], p_w0[4]);
    p_blocks[39] = p_vecPrecomputedBlocks[39];
    p_blocks[40] = XOR(XOR(p_vecPrecomputedBlocks[40], p_w0[4]), p_w0[9]);
    p_blocks[41] = p_vecPrecomputedBlocks[41];
    p_blocks[42] = XOR(XOR(p_vecPrecomputedBlocks[42], p_w0[6]), p_w0[8]);
    p_blocks[43] = XOR(p_vecPrecomputedBlocks[43], p_w0[10]);
    p_blocks[44] = XOR(XOR(p_vecPrecomputedBlocks[44], p_w0[6]), XOR(p_w0[3], p_w0[7]));
    p_blocks[45] = p_vecPrecomputedBlocks[45];
    p_blocks[46] = XOR(XOR(p_vecPrecomputedBlocks[46], p_w0[4]), p_w0[11]);
    p_blocks[47] = XOR(XOR(p_vecPrecomputedBlocks[47], p_w0[8]), p_w0[4]);
    p_blocks[48] = XOR(XOR(XOR(p_vecPrecomputedBlocks[48], p_w0[8]), XOR(p_w0[4], p_w0[3])), XOR(p_w0[10], p_w0[5]));
    p_blocks[49] = XOR(p_vecPrecomputedBlocks[49], p_w0[12]);
    p_blocks[50] = XOR(p_vecPrecomputedBlocks[50], p_w0[8]);
    p_blocks[51] = XOR(XOR(p_vecPrecomputedBlocks[51], p_w0[6]), p_w0[4]);
    p_blocks[52] = XOR(XOR(p_vecPrecomputedBlocks[52], p_w0[8]), XOR(p_w0[4], p_w0[13]));
    p_blocks[53] = p_vecPrecomputedBlocks[53];
    p_blocks[54] = XOR(XOR(p_vecPrecomputedBlocks[54], p_w0[7]), XOR(p_w0[10], p_w0[12]));
    p_blocks[55] = XOR(p_vecPrecomputedBlocks[55], p_w0[14]);
    p_blocks[56] = XOR(XOR(XOR(p_vecPrecomputedBlocks[56], p_w0[6]), XOR(p_w0[4], p_w0[11])), XOR(p_w0[7], p_w0[10]));
    p_blocks[57] = XOR(p_vecPrecomputedBlocks[57], p_w0[8]);
    p_blocks[58] = XOR(XOR(p_vecPrecomputedBlocks[58], p_w0[8]), XOR(p_w0[4], p_w0[15]));
    p_blocks[59] = XOR(XOR(p_vecPrecomputedBlocks[59], p_w0[8]), p_w0[12]);
    p_blocks[60] = XOR(XOR(XOR(p_vecPrecomputedBlocks[60], p_w0[8]), XOR(p_w0[4], p_w0[7])), XOR(p_w0[12], p_w0[14]));
    p_blocks[61] = XOR(p_vecPrecomputedBlocks[61], p_w0[16]);
    p_blocks[62] = XOR(XOR(XOR(p_vecPrecomputedBlocks[62], p_w0[6]), XOR(p_w0[12], p_w0[8])), p_w0[4]);
    p_blocks[63] = XOR(p_vecPrecomputedBlocks[63], p_w0[8]);
    p_blocks[64] = XOR(XOR(XOR(XOR(p_vecPrecomputedBlocks[64], p_w0[6]), p_w0[7]), XOR(p_w0[17], p_w0[12])), XOR(p_w0[8], p_w0[4]));
    p_blocks[65] = p_vecPrecomputedBlocks[65];
    p_blocks[66] = XOR(XOR(p_vecPrecomputedBlocks[66], p_w0[14]), p_w0[16]);
    p_blocks[67] = XOR(XOR(p_vecPrecomputedBlocks[67], p_w0[8]), p_w0[18]);
    p_blocks[68] = XOR(XOR(p_vecPrecomputedBlocks[68], p_w0[11]), XOR(p_w0[14], p_w0[15]));
    p_blocks[69] = p_vecPrecomputedBlocks[69];
    p_blocks[70] = XOR(XOR(p_vecPrecomputedBlocks[70], p_w0[12]), p_w0[19]);
    p_blocks[71] = XOR(XOR(p_vecPrecomputedBlocks[71], p_w0[12]), p_w0[16]);
    p_blocks[72] = XOR(XOR(XOR(XOR(p_vecPrecomputedBlocks[72], p_w0[11]), XOR(p_w0[12], p_w0[18])), XOR(p_w0[13], p_w0[16])), p_w0[5]);
    p_blocks[73] = XOR(p_vecPrecomputedBlocks[73], p_w0[20]);
    p_blocks[74] = XOR(XOR(p_vecPrecomputedBlocks[74], p_w0[8]), p_w0[16]);
    p_blocks[75] = XOR(XOR(p_vecPrecomputedBlocks[75], p_w0[6]), XOR(p_w0[12], p_w0[14]));
}
