#!/usr/bin/env python
# Unroll loops for outer precomputation
zeroValues = [0] + [j for j in range(2, 15)] + [16]

for i in range(17, 80):
    output = "p_precomputedBlocks[" + str(i) + "] = "
    initialLen = len(output)
    oldLen = len(output)
    if((i - 3) not in zeroValues):
        output += "p_precomputedBlocks[" + str(i - 3) + "]"
    if(len(output) > oldLen):
        output += " ^ "
        oldLen = len(output)
    
    if((i - 8) not in zeroValues):
        output += "p_precomputedBlocks[" + str(i - 8) + "]"
    if(len(output) > oldLen):
        output += " ^ "
        oldLen = len(output)
    
    if((i - 14) not in zeroValues):
        output += "p_precomputedBlocks[" + str(i - 14) + "]"
    if(len(output) > oldLen):
        output += " ^ "
        oldLen = len(output)
    
    if((i - 16) not in zeroValues):
        output += "p_precomputedBlocks[" + str(i - 16) + "]"
    if(len(output) > oldLen):
        output += " ^ "
        oldLen = len(output)
    
    output = output[:-3] + ";\np_precomputedBlocks[" + str(i) + "] = LEFT_ROTATE(p_precomputedBlocks[" + str(i) + "], 1);"

    if(initialLen == oldLen):
        output = "p_precomputedBlocks[" + str(i) + "] = 0;"
    print(output)

print("\n\n")
# Unroll loops for inner precomputation
for i in range(1, 21):
    print("p_w0[" + str(i) + "] = LEFT_ROTATE(p_blocks[0], " + str(i) + ");")
for i in range(16, 76):
    print("p_blocks[" + str(i) + "] = p_precomputedBlocks[" + str(i) + "]")
    
print("\n\n")
# Unroll assignment loop for outer precomputation
for i in range(17, 80):
    print("p_precomputedBlocks[" + str(i) + "] = SET1INT(p_tempPrecomputedBlocks[" + str(i) + "]);")
