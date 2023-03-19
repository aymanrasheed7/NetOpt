#include<bits/extc++.h>
using namespace std;
using lll = long long;
using ldb = long double;
struct Node {
    lll numer, denom, nmask, childNumer, childDenom, childNmask;
    Node(lll a = 0, lll b = 1, lll c = 0, lll d = 0, lll e = 1, lll f = 0) {
        b < 0 ? (b = -b, a = -a) : 0; e < 0 ? (e = -e, d = -d) : 0; numer = a;
        denom = b; nmask = c; childNumer = d; childDenom = e; childNmask = f;
        childNumer /= f = __gcd(abs(d), e); childDenom /= f;
        numer /= c = __gcd(abs(a), b); denom /= c;
    } Node operator+(Node& a) {
        return Node(numer * a.denom + denom * a.numer, denom * a.denom,
            nmask | a.nmask, numer, denom, nmask);
    } Node operator|(Node& a) {
        return Node(numer * a.numer, numer * a.denom + denom * a.numer,
            nmask | a.nmask, numer, denom, nmask);
    } bool operator<(ldb a) { return numer < denom* a; }
    bool operator<(Node& a) { return numer * a.denom < denom* a.numer; }
    bool operator==(Node& a) { return numer * a.denom == denom * a.numer; }
    ldb operator~() { return childNumer / (ldb)childDenom; }
    ldb operator!() { return numer / (ldb)denom; }
}; using pnn = pair<Node, Node>;
lll ONE_SERIES[] = { 1,1,1,1,1,1,1,1,1,1,1,1,0 };
lll INT_SERIES[] = { 1,2,3,4,5,6,7,8,9,10,11,12,0 };
lll ODD_SERIES[] = { 1,3,5,7,9,11,13,15,17,19,21,23,0 };
lll E12_SERIES[] = { 10,12,15,18,22,27,33,39,47,56,68,82,0 };
lll* series, denominator = 1, nT = 256, mask2id[8192];
vector<Node> mini(8192), maxi(8192), table[8192];
pnn rec(ldb target, lll mask, lll topFlag = 0) {
    if (mask == (mask & -mask)) return topFlag ? printf("%.1f", series[mask2id
        [mask]] / (denominator + 0.0)) : 0, pnn(table[mask][0], table[mask][0]);
    Node lo = mini[mask], hi = maxi[mask], temp, loch, hich;
    lll pos = 0, subMask = 0, compleMask = 0;
    function<void(Node&, Node)> upd = [&](Node& e, Node f) {
        abs(target - !f) < abs(target - !e) ? e = f : e; };
    if (topFlag) {
        if (target < !lo || hi < target) hi < target ? lo = hi : hi = lo;
        else if (table[mask].empty()) {
            mutex mtx; lll rJ = 0, nJ = 0; vector<thread> threads(nT);
            for (subMask = mask; subMask = (subMask - 1) & mask;
                nJ += table[subMask].size());
            vector<lll> qMask(rJ = nJ);
            vector<Node> qNode(nJ), qlo(nJ), qhi(nJ);
            for (subMask = mask; subMask = (subMask - 1) & mask;)
                if (!table[compleMask = mask ^ subMask].empty()) for (pos =
                    table[compleMask].size(); pos--; qMask[--rJ] = subMask,
                    qNode[rJ] = table[compleMask][pos]); rJ = nJ;
            function<void()> par = [&]() {
                for (lll i; mtx.lock(), i = --rJ, mtx.unlock(), 0 <= i;
                    tie(qlo[i], qhi[i]) = rec(qNode[i] < target ? target - !
                        qNode[i] : 1 / (1 / target - 1 / !qNode[i]), qMask[i])); };
            for (lll i = nT; i--; threads[i] = thread(par));
            for (lll i = nT; i--; threads[i].join());
            for (; nJ--; upd(qNode[nJ] < target ? lo : hi, qNode[nJ] < target ?
                qNode[nJ] | maxi[qMask[nJ]] : qNode[nJ] + mini[qMask[nJ]]),
                upd(lo, qNode[nJ] < target ? qNode[nJ] + qlo[nJ] :
                    qNode[nJ] | qlo[nJ]), upd(hi, qNode[nJ] < target
                        ? qNode[nJ] + qhi[nJ] : qNode[nJ] | qhi[nJ]));
        } else pos = lower_bound(table[mask].begin(), table[mask].end(),
            target) - table[mask].begin(), lo = table[mask][pos - !!pos],
            hi = Node(table[mask][min(pos, (lll)table[mask].size() - 1)]);
        abs(target - !hi) < abs(target - !lo) ? loch = hi : loch = lo;
        printf("("); rec(~loch, loch.childNmask, 1); printf(!loch < ~loch ?
            "|" : "+"); rec(!loch < ~loch ? 1 / (1 / !loch - 1 / ~loch) :
                !loch - ~loch, mask^ loch.childNmask, 1); printf(")");
    } else if (target < !lo || hi < target) hi < target ? lo = hi : hi = lo;
    else if (!table[mask].empty()) pos = lower_bound(table[mask].begin(), table
        [mask].end(), target) - table[mask].begin(), lo = table[mask][pos - !!
        pos], hi = Node(table[mask][min(pos, (lll)table[mask].size() - 1)]);
    else for (subMask = mask; subMask = (subMask - 1) & mask;) if (!table
        [compleMask = mask ^ subMask].empty()) for (pos = table[compleMask]
            .size(); pos--; temp = table[compleMask][pos], tie(loch, hich) =
            rec(temp < target ? target - !temp : 1 / (1 / target - 1 / !temp),
                subMask), upd(temp < target ? lo : hi, temp < target ? temp |
                    maxi[subMask] : temp + mini[subMask]), upd(lo, temp <
                        target ? temp + loch : temp | loch), upd(hi, temp
                            < target ? temp + hich : temp | hich));
    return pnn(lo, hi);
}
int main(int argc, char** argv) {
    chrono::_V2::system_clock::time_point begin
        = chrono::high_resolution_clock::now();
    Node lo, hi; ldb target = sqrtl(2); series = E12_SERIES;
    lll nR = atoi(argv[3]), tableSize = nR + 1 >> 1, maskSize = 1 << nR;
    if (4 < argc) tableSize = atoi(argv[4]); if (5 < argc) nT = atoi(argv[5]);
    if (!strcmp(argv[1], "E12")) denominator = 10;
    else if (!strcmp(argv[1], "ONE")) series = ONE_SERIES;
    else if (!strcmp(argv[1], "INT")) series = INT_SERIES;
    else if (!strcmp(argv[1], "ODD")) series = ODD_SERIES;
    if (!strcmp(argv[2], "E")) target = expl(1);
    else if (!strcmp(argv[2], "PI")) target = acosl(-1);
    else if (!strcmp(argv[2], "PHI")) target = sqrtl(5) * 0.5 + 0.5;
    else if (!strcmp(argv[2], "SQRT")) target = sqrtl(nR);
    else if (!strncmp(argv[2], "SQRT", 4)) target = sqrtl(stold(argv[2] + 4));
    else target = stold(argv[2]), target < 0 ? (target = -target,
        series[nR++] = 0, maskSize = 1 << nR, tableSize = nR + 1 >> 1) : 0;
    for (lll mask = 0, subMask, compleMask, i, j; maskSize > ++mask; sort(table
        [mask].begin(), table[mask].end()), table[mask].resize(unique(table
            [mask].begin(), table[mask].end()) - table[mask].begin()))
        if (mask == (mask & -mask) && (table[mask].push_back(mini
            [mask] = maxi[mask] = Node(series[mask2id[mask] =
                __builtin_ctzll(mask)], denominator, mask)), 1));
        else if (mini[mask] = mini[mask & -mask] | mini[mask ^ (mask & -mask)],
            maxi[mask] = maxi[mask & -mask] + maxi[mask ^ (mask & -mask)],
            __builtin_popcountll(mask) <= tableSize)
            for (subMask = mask; subMask = (subMask - 1) & mask;) if (subMask <
                (compleMask = mask ^ subMask)) for (i = table[subMask].size();
                    i--;) for (j = table[compleMask].size(); j--; table[mask]
                        .push_back(table[subMask][i] + table[compleMask][j]),
                        table[mask].push_back(table[subMask][i] |
                            table[compleMask][j]));
    printf("Solution: "), tie(lo, hi) = rec(target, maskSize - 1, 1);
    abs(target - !hi) < abs(target - !lo) ? (swap(lo, hi), 0) : 0;
    printf("\nResult: %.16lf (%lld/%lld)\nTarget: %.16lf\nCost: %le\n", (double)
        !lo, lo.numer, lo.denom, (double)target, (double)abs(target - !lo));
    printf("CPU time (s): %.9lf\n", clock() / (CLOCKS_PER_SEC + 0.0));
    chrono::_V2::system_clock::time_point end
        = chrono::high_resolution_clock::now();
    chrono::_V2::system_clock::duration elapsed
        = chrono::duration_cast<chrono::nanoseconds>(end - begin);
    printf("Execution time (s): %.9lf\n", elapsed * 1e-9); exit(0);
}