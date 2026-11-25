# BFW 2021 (Baillie-Fiori-Wagstaff) Enhanced BPSW - extracted from arXiv 2006.14425v2

## Method A* (parameter selection)
D in sequence 5, -7, 9, -11, 13, -15, ... such that Jacobi(D,n) = -1. P=1, Q=(1-D)/4.
If Q == -1 (i.e. D==5), set P=Q=5 instead (forces Q != +/-1 mod n).

## Strong Lucas conditions (Congruences 11 & 12)
Write n+1 = d * 2^s, d odd.
slprp iff U_d == 0 (mod n) OR V_{d*2^r} == 0 (mod n) for some 0 <= r < s.

## Lucas-V condition (Congruence 2)
V_{n+1} == 2*Q (mod n)

## Euler criterion condition (Congruence 4 / Step 5)
Q^((n+1)/2) == Q * (Q/n) (mod n)

## Enhanced BPSW (Section 6) - 5 steps
1. If n not strong probable prime base 2 -> composite; stop.
2. Choose Lucas params via Method A*. If Jacobi(D,n)==0: if |D|<n, or |D|>=n but n ∤ |D| -> composite; stop.
3. If n not slprp(P,Q) -> composite; stop.
4. If n not vprp(P,Q) i.e. V_{n+1} != 2Q (mod n) -> composite; stop.
5. If Q^((n+1)/2) != Q*(Q/n) (mod n) -> composite; stop. Else probably prime.

