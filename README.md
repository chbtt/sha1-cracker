# SHA-1 Cracker
C implementation of a SHA-1 (https://en.wikipedia.org/wiki/SHA-1) cracker with various optimizations.\
This was an assignment at university where we had to find the preimage of a given SHA-1 hash. In our scenario, 
the preimage was guaranteed to be six lower-case letters (for example, "aaaaaa" or "passwd"). In addition to 
general optimizations such as meet-in-the-middle, initial-step, early-exit (https://hashcat.net/events/p13/), 
those limitations allowed for further scenario-specific optimizations (mainly by identifying constants).