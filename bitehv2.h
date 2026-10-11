/**
 * @file bitehv2.h
 *
 * @brief Iterative 2-D Pareto front + exact hypervolume.
 *
 * E-mail: aleksey.vaneev@gmail.com or info@voxengo.com
 *
 * @section license License
 *
 * Copyright (c) 2026 Aleksey Vaneev
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#ifndef BITEHV2_INCLUDED
#define BITEHV2_INCLUDED

// The front F1 is a std::map keyed by the dominance-ranking statistic
//
//        S = a
//
// (mapped value: (b, x[])), ascending in the key.
//
// WHY KEYING BY S IS LEGAL (minimization Pareto compliance):
//   x dominates y ==> a(x) <= a(y), regardless of sign, so
//       x dominates y  ==>  S(x) <= S(y)   (equality possible).
//   Hence, for a newcomer p with key S, only keys <= S can dominate p and
//   only keys > S can be dominated by p. The front is a strict antichain,
//   so along ascending keys b is STRICTLY DECREASING and both sides
//   collapse to boundary checks:
//     - the smallest b among keys <= S sits at the largest such key, so
//       rejection reduces to the ONE entry prev(upper_bound(S));
//     - entries dominated by p (b >= b(p)) form a key-ordered PREFIX after
//       that point, so the erasure sweep stops at the first b < b(p).
//
// WHY std::map (unique keys) SUFFICES - the equal-a consumption premise:
//   among same-key entries only the one with the smallest b can survive,
//   because it dominates every same-key entry with a larger b. addPoint
//   enforces exactly this: the entry at key == S rejects p when its b <=
//   b(p) and is erased when its b > b(p). All same-key cases are covered,
//   so a committed insert never collides with an existing key: each key
//   permanently holds its minimum-b survivor.
//
// THE TIE CAVEAT (the only place geometry is needed):
//   the entry at key == S, if any, may dominate or be dominated by p, so it
//   is checked geometrically BEFORE any mutation.
//
// CACHED REJECTION: after each acceptance a strong front point is cached
// (before the first acceptance the cache holds the reference-corner
// sentinel, which can only reject the zero-area point (refX, refY) itself);
// a newcomer dominated by it is rejected in O(1), before any tree work -
// the hot-path fast exit for dominated-heavy streams. The cache always
// names a live front member once any point is accepted (see cacheUpdate_).
//
// HV INCREMENTAL MAINTENANCE: hv_ is kept exact with O(1) slab deltas; no
// per-insertion sweep. The strict antichain makes b strictly decreasing
// along ascending a, so HV telescopes into per-point slabs
//     HV = sum over points q of  slab(q),
//     slab(q) = (nextA(q) - a(q)) * (refY - b(q)),   nextA(last) = refX.
//   Inserting p splits the predecessor's slab at a = a(p): p takes the tail
//   (nextA - a) * (refY - b), the predecessor keeps its head. Erasing q
//   removes slab(q) and donates its width to the predecessor's slab. Deltas
//   are written cancellation-free, i.e. (refY - b) - (refY - b') appears as
//   b' - b.

#include <cassert>
#include <iterator>
#include <map>
#include <utility>
#include <vector>

class RankSumFrontHV {
public:
    typedef std::pair<double, std::vector<double> > Mapped; // (b, x[])
    typedef std::map<double, Mapped> Front;                 // key = S = a

    // (refX, refY): upper-right reference corner; all points must satisfy
    // a <= refX, b <= refY.
    explicit RankSumFrontHV(double refX, double refY)
        : cacheA_(refX), cacheB_(refY), cacheMean_(refX * 0.5 + refY * 0.5),
          refX_(refX), refY_(refY), hv_(0.0) {}

    // Insert one point; returns true iff it entered the front.
    bool addPoint(double a, double b, const double *x = nullptr,
        const int xN = 0) {

        assert(a <= refX_ && b <= refY_ && "point outside the reference box");

        // Cached rejection: the cache names a front point - or, before the
        // first acceptance, the reference-corner sentinel, which can only
        // "reject" a newcomer exactly at (refX, refY), a zero-area point
        // that any front member would dominate anyway (see cacheUpdate_).
        // If the cache dominates the newcomer, p is rejected without
        // touching the map.
        if (cacheA_ <= a && cacheB_ <= b)
            return false;

        const double S = a;

        // Rejection and tie resolution in ONE candidate check: the entry
        // with the largest key <= S carries the smallest b among those keys
        // (b strictly decreases with the key), so it is the only possible
        // dominator of p.
        Front::iterator z3 = F1.upper_bound(S); // first key > S; also the
                                                // erasure sweep start
        Front::iterator w = F1.end();           // running predecessor for
                                                // the slab bookkeeping
        if (z3 != F1.begin()) {
            Front::iterator t = std::prev(z3);
            if (t->second.first <= b)
                return false;                   // it dominates p; F1 untouched
            if (t->first == S) {                // the tie: p dominates it
                w = (t == F1.begin()) ? F1.end() : std::prev(t);
                z3 = eraseOne_(t, w);           // w unchanged: t drops out
            } else {
                w = t;                          // boundary of keys < S
            }
        }

        // Erasure: entries p dominates (b >= b(p)) form a key-ordered prefix
        // (b strictly decreases with the key), so stop at the first b < b(p).
        for (Front::iterator t = z3; t != F1.end(); ) {
            if (t->second.first < b)
                break;
            t = eraseOne_(t, w);
        }

        // Commit. The candidate check above removed any same-key entry (or
        // rejected p), so the insert cannot collide; see the file header.
        const Front::iterator pos = F1.insert(std::make_pair(S,
            Mapped(b, (x == nullptr || xN <= 0 ?
                std::vector<double>() : std::vector<double>(x, x + xN))))).first;

        // O(1) slab delta for the insertion; w is p's predecessor:
        //   first point: (nextA - a) * (refY - b)   - freshly covered slab
        //   last point:  (refX  - a) * (b(w) - b)   - predecessor donates its tail
        //   interior:    (nextA - a) * (b(w) - b)   - slab tail changes owner
        {
            const Front::iterator s = std::next(pos);
            const double nextA = (s == F1.end()) ? refX_ : s->first;
            if (w == F1.end())
                hv_ += (nextA - a) * (refY_ - b);
            else if (s == F1.end())
                hv_ += (refX_ - a) * (w->second.first - b);
            else
                hv_ += (nextA - a) * (w->second.first - b);
        }

        cacheUpdate_(a, b);
        return true;
    }

    double hypervolume() const { return hv_; }

    // THE front, read-only: begin() -> end() is the ranking itself, ascending
    // S = from best-balanced toward the extremes.
    const Front& getFront() const { return F1; }

private:
    // Erase q with the O(1) slab delta: -slab(q) plus the predecessor's
    // slab widening, combined cancellation-free into
    //     width * (b(q) - b(pred)),   width = nextA - a(q),
    // zero width when pred is a tie (same key). w is the running predecessor
    // kept by the caller's sweep: it advances only across survivors, since
    // erased points drop out and leave their predecessor unchanged.
    // Returns q's successor, saved before unlinking.
    Front::iterator eraseOne_(Front::iterator q, Front::iterator w) {
        const Front::iterator s = std::next(q);   // successor before unlinking
        const double nextA = (s == F1.end()) ? refX_ : s->first;
        const double width = nextA - q->first;

        hv_ += width * (w == F1.end() ? -(refY_ - q->second.first)
                                      : q->second.first - w->second.first);

        F1.erase(q);  // invalidates only q; s stays valid
        return s;
    }

    // Cache maintenance, acceptance-only. The cache is initialized to the
    // reference corner; the first accepted point always replaces it (it is
    // elementwise <= the corner and has a smaller mean), so from then on the
    // cache names a live front member. Soundness: a cached point can leave
    // the front only by being dominated by a later accepted point, and that
    // point passes the first replacement test below, so the rejection test
    // in addPoint stays valid. Replacement policy: a point dominating the
    // current cache covers its rejection region and more; otherwise keep
    // the more central point (smaller a+b mean). The heuristic affects
    // rejection speed only, never membership or HV.
    void cacheUpdate_(double a, double b) {
        const double mean = a * 0.5 + b * 0.5; // half-sum: finite even if a+b overflows
        if ((a <= cacheA_ && b <= cacheB_) || mean < cacheMean_) {
            cacheA_ = a; cacheB_ = b; cacheMean_ = mean;
        }
    }

    Front F1;

    double cacheA_, cacheB_, cacheMean_;
    double refX_, refY_;
    double hv_;
};

#endif // BITEHV2_INCLUDED
