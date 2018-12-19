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
#include "testCracker.h"

// function prototype
void printPreimageDiff(char *p_expected,
                       char *p_received);

// test vectors
const struct crackSha1TestVec testVectors[] = { { "ananas",
                                                { 0x755BD810, 0xD2BE0EBC, 0xBB6CE6F5, 0x32B3D9CF, 0xCF9D9695 }},
                                                { "passwd",
                                                { 0x30274C47, 0x903BD1BA, 0xC7633BBF, 0x09743149, 0xEBAB805F }},
                                                { "qfucra",
                                                { 0x3854E277, 0xA37AEE29, 0xBF9ECC86, 0xFB983737, 0xCF9D9695 }},
                                                { "swords",
                                                { 0xD6056E47, 0xD33A009D, 0x754613AF, 0xBB20A3A3, 0x86496177 }},
                                                { "zzzzzz",
                                                { 0x984FF6EE, 0x7C78078D, 0x4CB1CA08, 0x255303FB, 0x8741D986 }} };

/**
 * Function: main
 */
int main()
{
    int  numberOfTests = sizeof(testVectors) / sizeof(struct crackSha1TestVec),
         testsPassed   = 0;
    char p_result[6];
    printf("Testing SHA1-Cracker...\n");
    for(int i = 0; i < numberOfTests; i++)
    {
        if(crackHash(testVectors[i].resultingHash,
                     p_result) == 0)
        {
            if(memcmp(testVectors[i].p_preImage,
                      p_result,
                      6) == 0)
                testsPassed++;
            else
                printPreimageDiff(testVectors[i].p_preImage,
                                  p_result);
        }
    }
    printf("Passed %d/%d!\n", testsPassed,
                              numberOfTests);

	return 0;
}
/**
 * Function: printPreimageDiff
 */
void printPreimageDiff
(
    char *p_expected,
    char *p_received
)
{
    printf("\nExpected: %c%c%c%c%c%c\n", p_expected[0],
                                         p_expected[1],
                                         p_expected[2],
                                         p_expected[3],
                                         p_expected[4],
                                         p_expected[5]);
    printf("Received: %c%c%c%c%c%c\n\n", p_received[0],
                                         p_received[1],
                                         p_received[2],
                                         p_received[3],
                                         p_received[4],
                                         p_received[5]);
}
