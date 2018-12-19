#!/usr/bin/env python
# Compute constants for initial state optimization
#
# Output:
#         First Addition:             0x29FB498B3
#         Second Addition:            0x166B0CD0D
#         Third Addition:             0x0F33D5697
#         Fourth Addition:            0x0D675E47B
#         Fifth Addition:             0x0B453C259
#         
#         Preimage-Length Addition:   0x05A8279C9
import numpy

K_00_19 = 0x5A827999
a = 0x67452301
b = 0xEFCDAB89
c = 0x98BADCFE
d = 0x10325476
e = 0xC3D2E1F0

def rotate_left(word, bits):
    word_digits = list(bin(word))
    word_digits = ((34 - len(word_digits)) * [0]) + word_digits[2:]
    word_digits = numpy.roll(word_digits, -(bits))
    word_digits = ''.join(word_digits)

    return int(word_digits, base = 2) 

def f_function(mB, mC, mD):
    return (mD ^ (mB & (mC ^ mD)))

addition_1 = K_00_19 + e + f_function(b, c, d) + rotate_left(a, 5)
e = d; d = c; c = rotate_left(b, 30); b = a
addition_2 = K_00_19 + e + f_function(b, c, d)
e = d; d = c; c = rotate_left(b, 30)
addition_3 = K_00_19 + e
e = d; d = c
addition_4 = K_00_19 + e
e = d
addition_5 = K_00_19 + e
addition_6 = K_00_19 + 0x00000030

print("First Addition:  " + hex(addition_1))
print("Second Addition: " + hex(addition_2))
print("Third Addition:  " + hex(addition_3))
print("Fourth Addition: " + hex(addition_4))
print("Fifth Addition:  " + hex(addition_5))
print("\nPreimage-Length Addition: " + hex(addition_6))
