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
#ifndef TESTCRACKER_H
#define TESTCRACKER_H

#include "../src/sha1-cracker.h"
// printf
#include <stdio.h>
// memcmp
#include <string.h>
// clock, clock_t
#include <time.h>

struct crackSha1TestVec
{
    char        *p_preImage;
    struct hash resultingHash;
};

#endif
