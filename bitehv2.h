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
//   because it dominates every same-key entry with a larger b. The entry
//   at key == S rejects p when its b <= b(p) and is erased when its b >
//   b(p); every same-key case is covered, so a committed insert never
//   collides with an existing key and each key permanently holds its
//   minimum-b survivor.
//
// THE TIE CAVEAT (the only place geometry is needed):
//   the entry at key == S, if any, may dominate or be dominated by p, so
//   it is checked geometrically BEFORE any mutation.
//
// ANTI-CHAIN FAST PATHS: on an anti-chain stream nothing is ever rejected,
// so acceptance must be cheap. The last-node pivot (lastA_/lastB_) turns a
// strictly-beyond-end newcomer into an O(1) hinted append, and the cached
// last-inserted-node state - ins_pos_ with coordinates posA_/posB_ and its
// successor ins_s_ with coordinates sucA_/sucB_ - admits a newcomer sitting
// strictly between them in both coordinates, which is then incomparable
// with the whole front and is inserted with ins_s_ as the O(1) map hint.
// The reference-corner sentinel for every cached coordinate doubles as the
// "empty front / no successor" dead switch. All cached state is refreshed
// at every acceptance and read only before any mutation, so it is valid at
// each addPoint entry. Together these make monotone and circle-order
// anti-chain streams O(1) per point.
//
// CACHED REJECTION: after each acceptance a strong front point is cached;
// a newcomer dominated by it is rejected in O(1), before any tree work.
// The cache always names a live front member (see cacheUpdate_); the
// reference-corner init can only "reject" a newcomer exactly at (refX,
// refY), a zero-area point that any front member would dominate anyway.
//
// HV INCREMENTAL MAINTENANCE: hv_ is kept exact with O(1) slab deltas; no
// per-insertion sweep. The strict antichain makes b strictly decreasing
// along ascending a, so HV telescopes into per-point slabs
//     HV = sum over points q of  slab(q),
//     slab(q) = (nextA(q) - a(q)) * (refY - b(q)),   nextA(last) = refX.
//   Inserting p splits the predecessor's slab at a = a(p): p takes the
//   tail, the predecessor keeps its head. Erasing q removes slab(q) and
//   donates its width to the predecessor's slab. Deltas are written
//   cancellation-free, i.e. (refY - b) - (refY - b') appears as b' - b.

#include <cassert>
#include <iterator>
#include <map>
#include <utility>
#include <vector>

class RankSumFrontHV {
public:
    typedef std::pair<double, std::vector<double> > Mapped; // (b, x[])
    typedef std::map<double, Mapped> Front;                 // key = S = a

    // (refX, refY): upper-right reference corner; points must satisfy
    // a <= refX, b <= refY.
    explicit RankSumFrontHV(double refX, double refY)
        : cacheA_(refX), cacheB_(refY), cacheMean_(refX * 0.5 + refY * 0.5),
          lastA_(refX), lastB_(refY),
          posA_(refX), posB_(refY),
          sucA_(refX), sucB_(refY),
          refX_(refX), refY_(refY), hv_(0.0),
          ins_pos_(F1.end()), ins_s_(F1.end()) {}

    // Insert one point; returns true iff it entered the front.
    bool addPoint(double a, double b, const double *x = nullptr,
        const int xN = 0) {

        assert(a <= refX_ && b <= refY_ && "point outside the reference box");

        // Cached rejection: if the cached front point dominates p, reject
        // without touching the map.
        if (cacheA_ <= a && cacheB_ <= b)
            return false;

        const double S = a;

        // Append fast path: a > last.a and b < last.b makes p incomparable
        // with the whole front - every node has a' <= last.a < a and
        // b' >= last.b > b (b decreases with a), so nothing dominates p
        // (needs b' <= b) and p dominates nothing (needs a <= a'). O(1)
        // hinted commit at the end; the delta is the last-point slab
        // formula. The corner sentinel keeps the test false while the
        // front is empty (append needs a > refX).
        if (a > lastA_ && b < lastB_) {
            hv_ += (refX_ - a) * (lastB_ - b);  // old last donates its tail
            const Front::iterator pos = F1.insert(F1.end(), std::make_pair(S,
                Mapped(b, (x == nullptr || xN <= 0 ?
                    std::vector<double>() : std::vector<double>(x, x + xN)))));
            lastA_ = a; lastB_ = b;
            ins_pos_ = pos; ins_s_ = F1.end();
            posA_ = a; posB_ = b; sucA_ = refX_; sucB_ = refY_;  // no successor
            cacheUpdate_(a, b);
            return true;
        }

        // Gap fast path after the last inserted node, against cached
        // coordinates only: posA_ < a < sucA_ and sucB_ < b < posB_ place
        // p strictly between ins_pos_ and ins_s_ in both coordinates, so p
        // is incomparable with the whole front - nodes left of the gap have
        // b > posB_ > b(p) (cannot dominate p; p cannot dominate them,
        // a > posA_ > their a), nodes right of ins_s_ have b < sucB_ < b(p)
        // (p cannot dominate them; they cannot dominate p, a < sucA_). No
        // rejection or erasure; ins_s_ doubles as the O(1) hint. sucA_/
        // sucB_ hold the reference corner exactly when ins_s_ is F1.end(),
        // and refY < b never holds, so the early sucB_ < b term kills the
        // no-successor case without an end() check.
        if (posA_ < a && sucB_ < b && a < sucA_ && b < posB_) {
            const Front::iterator pos = F1.insert(ins_s_, std::make_pair(S,
                Mapped(b, (x == nullptr || xN <= 0 ?
                    std::vector<double>() : std::vector<double>(x, x + xN)))));
            hv_ += (sucA_ - a) * (posB_ - b);   // interior slab delta
            ins_pos_ = pos; posA_ = a; posB_ = b;  // ins_s_/sucA_/sucB_ unchanged
            cacheUpdate_(a, b);
            return true;
        }

        // Rejection and tie resolution in ONE candidate check: the entry
        // with the largest key <= S carries the smallest b among those keys
        // (b strictly decreases with the key), so it is the only possible
        // dominator of p.
        Front::iterator z3 = F1.upper_bound(S); // first key > S; also the
                                                // erasure sweep start
        Front::iterator w;                      // running predecessor for
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
        } else {
            w = F1.end();
        }

        // Erasure: entries p dominates (b >= b(p)) form a key-ordered prefix
        // (b strictly decreases with the key), so stop at the first b <
        // b(p). The surviving iterator t is p's slot successor - also the
        // O(1) insert hint (hint = insert-before position).
        Front::iterator t = z3;
        for (; t != F1.end(); ) {
            if (t->second.first < b)
                break;
            t = eraseOne_(t, w);
        }

        // Commit. The candidate check above removed any same-key entry (or
        // rejected p), so the insert cannot collide; the hinted overload
        // returns the iterator directly.
        const Front::iterator pos = F1.insert(t, std::make_pair(S,
            Mapped(b, (x == nullptr || xN <= 0 ?
                std::vector<double>() : std::vector<double>(x, x + xN)))));

        // O(1) slab delta; w is p's predecessor. One formula covers every
        // case: slab(p) = (nextA - a) * height, with height = refY - b
        // without a predecessor (freshly covered slab) and b(w) - b
        // otherwise (the predecessor donates its slab tail; the antichain
        // keeps this non-negative). The ins_s_ == F1.end() test doubles as
        // last-pivot maintenance: it sees the new last only when erasures
        // emptied the front or erased the old last - an accepted append
        // with intact neighbors always took the fast path above.
        {
            ins_s_ = std::next(pos);
            if (ins_s_ == F1.end()) {   // p is the new last node
                sucA_ = refX_; sucB_ = refY_;   // no successor
                lastA_ = a; lastB_ = b;
            } else {
                sucA_ = ins_s_->first;
                sucB_ = ins_s_->second.first;
            }
            const double nextA = (ins_s_ == F1.end()) ? refX_ : sucA_;
            hv_ += (nextA - a) * (w == F1.end() ? (refY_ - b)
                                                : (w->second.first - b));
        }
        ins_pos_ = pos; posA_ = a; posB_ = b;

        cacheUpdate_(a, b);
        return true;
    }

    double hypervolume() const { return hv_; }

    // THE front, read-only: begin() -> end() is the ranking itself, ascending
    // S = from best-balanced toward the extremes.
    const Front& getFront() const { return F1; }

private:
    // Cache maintenance, acceptance-only. Soundness: a cached point can
    // leave the front only by being dominated by a later accepted point,
    // and that point passes the first replacement test below, so the cache
    // always names a live front member and the rejection test stays valid.
    // The corner sentinel is evicted by the first accepted point (it is
    // elementwise <= and has a smaller mean). Replacement policy: a point
    // dominating the current cache covers its rejection region and more;
    // otherwise keep the more central point (smaller a+b mean). The
    // heuristic affects rejection speed only, never membership or HV.
    void cacheUpdate_(double a, double b) {
        const double mean = a * 0.5 + b * 0.5; // half-sum: finite even if a+b overflows
        if ((a <= cacheA_ && b <= cacheB_) || mean < cacheMean_) {
            cacheA_ = a; cacheB_ = b; cacheMean_ = mean;
        }
    }

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

    Front F1;

    // Cached coordinates use the reference corner as the "no real state
    // yet" sentinel: every dependent test is constructed so it cannot fire
    // on the sentinel (the rejection cache fires only for a newcomer
    // exactly at the corner, which carries zero area).
    double cacheA_, cacheB_, cacheMean_;
    double lastA_, lastB_;  // last front node
    double posA_, posB_;    // ins_pos_ coordinates
    double sucA_, sucB_;    // ins_s_ coordinates
    double refX_, refY_;
    double hv_;
    // Last inserted node and its successor; refreshed at every acceptance
    // and read only before any mutation (sentinel: F1.end()).
    Front::iterator ins_pos_, ins_s_;
};

#endif // BITEHV2_INCLUDED
