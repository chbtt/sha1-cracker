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
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#include "testCracker.h"

/**
 * Function: main
 */
int main() 
{
    int         err;
    char        p_result[6];
    clock_t     start,
                stop;
    double      elapsed;
    struct hash targetHash = { 0x984FF6EE, 0x7C78078D, 0x4CB1CA08, 0x255303FB, 0x8741D986 };

    for(int i = 0; i < 5; i++)
    {
        start = clock();
        err = crackHash(targetHash,
                        p_result);
        stop = clock();
        elapsed = ((double) (stop - start)) * 1000.0 / CLOCKS_PER_SEC;
        if(err == 0)
            printf("Run#%d\nTime: %f\n\n", i,
                                           elapsed);
        else
            printf("An error occurred!\n\n");
    }

    return 0;
}
