#include <bits/stdc++.h>
using namespace std;

mt19937 rng(233);

int main()
{
    uint64_t a = 478195865872675992;
    array<int, 64> bits;
    for(int i = 0; i < 64; i++) bits[i] = (a >> i) & 1;
    if(!bits[63] && (rng() & 1)) bits[63] += 2;
    for(int i = 0; i < 64; i++) cout << bits[i]; cout << endl;
    for(int i = 63; i > 0; i--) {
        if(bits[i] == 3 || (bits[i] > 0 && (rng() & 1))) bits[i]--, bits[i - 1] += 2;
    } 
    for(int i = 0; i < 64; i++) cout << bits[i]; cout << endl;
    uint64_t and_1 = 0, or_1 = 0, and_2 = 0, xor_2 = 0;
    uint64_t check = 0;
    for(int i = 0; i < 64; i++) {
        uint64_t m = 1ULL << i;
        check += m * bits[i];
        if(bits[i] == 0) {}
        else if(bits[i] == 2) {
            if(rng() & 1) and_1 += m;
            or_1 += m;
            xor_2 += m;
        }
        else {
            and_1 += m;
            and_2 += m;
            xor_2 += m;
        }
    }
    cout << check << endl;
    cout << and_1 << ' ' << or_1 << ' ' << and_2 << ' ' << xor_2 << endl;
    for(uint64_t t = 1; t <= 10; t++) cout << ((t & and_1) | or_1) + ((t & and_2) ^ xor_2) << endl;
    for(int check = 0; check < 50000; check++) {
        uint64_t s = 0;
        uint64_t a = rng(), b = rng(), c = rng();
        s += (c | a) + (c | b) - (a | b | c) - (c | (a & b));
        if(s) return cout << a << ' ' << b << ' ' << c, 0;
    }
}