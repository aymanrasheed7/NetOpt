/*
MIT License

Copyright (c) 2023 aymanrasheed7

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include<bits/extc++.h>
using namespace std;
struct Node {
    long long numer, denom, nmask, childNumer, childDenom, childNmask;
    Node(long long a = 0, long long b = 1, long long c = 0,
        long long d = 0, long long e = 1, long long f = 0) {
        numer = a; denom = b; nmask = c; childNumer = d; childDenom = e;
        childNmask = f; childNumer /= f = __gcd(d, e); childDenom /= f;
        numer /= c = __gcd(a, b); denom /= c;
    } Node operator+(Node& a) {
        return Node(numer * a.denom + denom * a.numer, denom * a.denom,
            nmask | a.nmask, numer, denom, nmask);
    } Node operator|(Node& a) {
        return Node(numer * a.numer, numer * a.denom + denom * a.numer,
            nmask | a.nmask, numer, denom, nmask);
    } bool operator<(double a) { return numer < denom* a; }
    bool operator<(Node& a) { return numer * a.denom < denom* a.numer; }
    bool operator==(Node& a) { return numer * a.denom == denom * a.numer; }
    double operator~() { return childNumer / (childDenom + 0.0); }
    double operator!() { return numer / (denom + 0.0); }
}; long long E12_SERIES[] = { 10,12,15,18,22,27,33,39,47,56,68,82,0,0,0,0 };
long long ONE_SERIES[] = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0 };
long long INT_SERIES[] = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,0 };
long long ODD_SERIES[] = { 1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,0 };
long long nT = 256; Node mini[65536], maxi[65536]; vector<Node> table[65536];
void recurse(Node& lo, Node& hi, double target, long long mask) {
    if (lo = mini[mask], hi = maxi[mask], mask == (mask & -mask)) return;
    Node temp, loch, hich; long long pos = 0, subMask = 0, compleMask = 0;
    function<void(Node&, Node)> upd = [&](Node& a, Node b) {
        abs(target - !b) < abs(target - !a) ? a = b : a; };
    if (target < !lo) hi = lo; else if (hi < target) lo = hi;
    else if (!table[mask].empty()) pos = lower_bound(table[mask].begin(), table
        [mask].end(), target) - table[mask].begin(), lo = table[mask][pos - !!
        pos], hi = table[mask][pos - (pos == table[mask].size())];
    else for (subMask = mask; subMask = (subMask - 1) & mask;) if (!table
        [compleMask = mask ^ subMask].empty()) for (pos = table[compleMask]
            .size(); pos--; temp = table[compleMask][pos], recurse(loch, hich,
                temp < target ? target - !temp : 1 / (1 / target - 1 / !temp),
                subMask), temp < target ? upd(lo, temp | maxi[subMask]) :
            upd(hi, temp + mini[subMask]), upd(lo, temp < target ? temp
                + loch : temp | loch), upd(hi, temp < target ?
                    temp + hich : temp | hich));
} void solve(Node& lo, Node& hi, double target, long long mask) {
    if (lo = mini[mask], hi = maxi[mask], mask == (mask & -mask)
        && printf("%.1f", !lo)) return;
    Node temp, loch, hich; long long pos = 0, subMask = 0, compleMask = 0;
    function<void(Node&, Node)> upd = [&](Node& a, Node b) {
        abs(target - !b) < abs(target - !a) ? a = b : a; };
    if (target < !lo) hi = lo; else if (hi < target) lo = hi;
    else if (table[mask].empty()) {
        mutex mtx; thread threads[nT]; long long rJ, nJ = 0;
        for (subMask = mask; subMask = (subMask - 1) & mask;
            nJ += table[subMask].size());
        vector<long long> qMask(rJ = nJ);
        vector<Node> qNode(nJ), qLo(nJ), qHi(nJ);
        for (subMask = mask; subMask = (subMask - 1) & mask;)
            if (!table[compleMask = mask ^ subMask].empty()) for (pos =
                table[compleMask].size(); pos--; qMask[--nJ] = subMask,
                qNode[nJ] = table[compleMask][pos]);
        function<void()> par = [&]() {
            for (long long i; mtx.lock(), i = --rJ, mtx.unlock(), 0 <= i;
                recurse(qLo[i], qHi[i], qNode[i] < target ? target - !qNode[i]
                    : 1 / (1 / target - 1 / !qNode[i]), qMask[i]), mtx.lock(),
                qNode[i] < target ? upd(lo, qNode[i] | maxi[qMask[i]]) :
                upd(hi, qNode[i] + mini[qMask[i]]), upd(lo, qNode[i] <
                    target ? qNode[i] + qLo[i] : qNode[i] | qLo[i]),
                upd(hi, qNode[i] < target ? qNode[i] + qHi[i] :
                    qNode[i] | qHi[i]), mtx.unlock()); };
        for (long long i = nT; i--; threads[i] = thread(par));
        for (long long i = nT; i--; threads[i].join());
    } else pos = lower_bound(table[mask].begin(), table[mask].end(),
        target) - table[mask].begin(), lo = table[mask][pos - !!pos],
        hi = table[mask][pos - (pos == table[mask].size())];
    abs(target - !hi) < abs(target - !lo) ? temp = hi : temp = lo;
    printf("("); solve(loch, hich, ~temp, temp.childNmask); printf(!temp <
        ~temp ? "|" : "+"); solve(loch, hich, !temp < ~temp ? 1 / (1 / !temp
            - 1 / ~temp) : !temp - ~temp, mask^ temp.childNmask); printf(")");
} int main(int argc, char** argv) {
    chrono::_V2::system_clock::time_point begin
        = chrono::high_resolution_clock::now();
    Node lo, hi; double target = sqrt(2);
    long long* series = E12_SERIES, denominator = 1;
    long long nR = atoi(argv[3]), tableSize = nR + 1 >> 1, maskSize = 1 << nR;
    if (4 < argc) tableSize = max(tableSize, (long long)atoi(argv[4]));
    if (5 < argc) nT = atoi(argv[5]);
    if (!strcmp(argv[1], "E12")) denominator = 10;
    else if (!strcmp(argv[1], "ONE")) series = ONE_SERIES;
    else if (!strcmp(argv[1], "INT")) series = INT_SERIES;
    else if (!strcmp(argv[1], "ODD")) series = ODD_SERIES;
    if (!strcmp(argv[2], "E")) target = exp(1);
    else if (!strcmp(argv[2], "PI")) target = acos(-1);
    else if (!strcmp(argv[2], "PHI")) target = sqrt(5) * 0.5 + 0.5;
    else if (!strcmp(argv[2], "SQRT")) target = sqrt(nR);
    else if (!strncmp(argv[2], "SQRT", 4)) target = sqrt(stod(argv[2] + 4));
    else target = stod(argv[2]), target < 0 ? (target = -target, series[nR++]
        = 0, maskSize = 1 << nR, tableSize = max(tableSize, nR + 1 >> 1)) : 0;
    for (long long mask = 0, subMask, compleMask, i, j; maskSize > ++mask; sort(
        table[mask].begin(), table[mask].end()), table[mask].resize(unique(table
            [mask].begin(), table[mask].end()) - table[mask].begin()))
        if (mask == (mask & -mask) && (table[mask].push_back(mini[mask] =
            maxi[mask] = Node(series[__builtin_ctzll(mask)],
                denominator, mask)), 1));
        else if (mini[mask] = mini[mask & -mask] | mini[mask ^ (mask & -mask)],
            maxi[mask] = maxi[mask & -mask] + maxi[mask ^ (mask & -mask)],
            __builtin_popcountll(mask) <= tableSize)
            for (subMask = mask; subMask = (subMask - 1) & mask;) if (subMask <
                (compleMask = mask ^ subMask)) for (i = table[subMask].size();
                    i--;) for (j = table[compleMask].size(); j--; table[mask]
                        .push_back(table[subMask][i] + table[compleMask][j]),
                        table[mask].push_back(table[subMask][i] |
                            table[compleMask][j]));
    printf("Solution: "), solve(lo, hi, target, maskSize - 1);
    abs(target - !hi) < abs(target - !lo) ? (swap(lo, hi), 0) : 0;
    printf("\nResult: %.16lf (%lld/%lld)\nTarget: %.16lf\nCost: %le\n",
        !lo, lo.numer, lo.denom, target, abs(target - !lo));
    printf("CPU time (s): %.9lf\n", clock() / (CLOCKS_PER_SEC + 0.0));
    chrono::_V2::system_clock::time_point end
        = chrono::high_resolution_clock::now();
    chrono::_V2::system_clock::duration elapsed
        = chrono::duration_cast<chrono::nanoseconds>(end - begin);
    printf("Execution time (s): %.9lf\n", elapsed * 1e-9); exit(0);
}
