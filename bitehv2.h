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

// The front F1 is a std::vector<std::tuple<double,double,double>>, sorted
// ascending by element 0 — the dominance-ranking statistic
//
//        S = s1 + s2,  s1 = a*(a+b),  s2 = b*(a+b)  ==>  S = a^2 + 2ab + b^2 .
//
// WHY KEYING BY S IS LEGAL (Pareto compliance, for a,b >= 0):
//   each s_i is non-decreasing in both objectives, so
//       x dominates y  ==>  S(x) <= S(y)   (with equality POSSIBLE).
//   Therefore, for a newcomer p with key S:
//     - only entries with key <  S can dominate p        (rejection zone),
//     - only entries with key >  S can be dominated by p (erasure zone),
//   and a single ordered scan decides acceptance and performs all erasures.
//
// THE TIE CAVEAT (the only place geometry is unavoidable):
//   dominance at EQUAL sums exists — (1,4) dominates (2,4), both S = 25 —
//   so inside the tie block (key == S) dominance may point either way and
//   is checked geometrically: a reject pass over the WHOLE block first
//   (a later tie may dominate p; no mutation before rejection is settled),
//   then an erase pass over the ties p dominates.
//
// HV: after each accepted insertion, HV is recomputed by one sweep of the
// front in its required geometric order (ascending a). F1 is a strict
// antichain by construction (accepted points were verified non-dominated;
// everything they dominate was erased), so along ascending a the b values
// are strictly decreasing and NO coverage test is needed: HV telescopes
// into   sum_{i=1}^{m-1} (a_{i+1} - a_i)*(refY - b_i) + (refX - a_m)*(refY - b_m),
// which the loop evaluates with preloaded (prevA, prevB) starting at i = 1.

#include <algorithm>
#include <cassert>
#include <tuple>
#include <vector>

class RankSumFrontHV {
public:
    // Front entry: (S, a, b). Sort order is by S; the tuple's own operator<
    // (lexicographic) is used for insertion to get a deterministic slot inside
    // a tie run.
    typedef std::tuple<double, double, double, std::vector<double> > Entry; // (S, a, b, x[])
    typedef std::vector<Entry> Front;                      // kept sorted by S

    // Key-only ordering: std::lower_bound / std::upper_bound with this
    // comparator behave exactly like their std::multimap counterparts
    // (first key >= S / first key > S), since the vector is ascending in S.
    struct ByS {
        bool operator()(const Entry& x, const Entry& y) const {
            return std::get<0>(x) < std::get<0>(y);
        }
    };

    // (refX, refY): upper-right reference corner; all points must satisfy
    // 0 <= a <= refX, 0 <= b <= refY.
    explicit RankSumFrontHV(double refX, double refY) : refX_(refX), refY_(refY), hv_(0.0) {}

    // Insert one point; returns true iff it entered the front.
    bool addPoint(double a, double b, const double *x = nullptr,
        const int xN = 0) {

        assert(a >= 0.0 && b >= 0.0 && "objectives must be non-negative");
        assert(a <= refX_ && b <= refY_ && "point outside the reference box");

        // ---- the statistic, expanded form to preserve precision ----
        const double S = a * a + 2.0 * a * b + b * b;

        // ---- zone boundaries by binary search on the key ----
        const Entry probe = std::make_tuple(S, 0.0, 0.0,
            std::vector<double>() );   // key = S, coords irrelevant

        // == Zone 1: keys < S — can ONLY dominate the newcomer ==
        // (If the newcomer dominated one of them, S would be <= their key < S:
        //  impossible.) Pure rejection checks; no erasure possible here.
        Front::iterator zone2Begin = std::lower_bound(F1.begin(), F1.end(), probe, ByS());
        for (Front::iterator it = F1.begin(); it != zone2Begin; ++it)
            if (std::get<1>(*it) <= a && std::get<2>(*it) <= b)
                return false;                           // dominated; F1 untouched

        // == Zone 2: tie block [zone2Begin, zone2End) — both directions possible ==
        Front::iterator zone2End = std::upper_bound(zone2Begin, F1.end(), probe, ByS());

        // 2a) Rejection pass over the WHOLE block before any mutation.
        for (Front::iterator t = zone2Begin; t != zone2End; ++t)
            if (std::get<1>(*t) <= a && std::get<2>(*t) <= b)
                return false;                           // a tie dominates p

        // 2b) Erasure pass: compact the ties p dominates over the survivors,
        //     then drop the tail in one erase. (We do not erase element by
        //     element here because the erase pass runs inside [zone2Begin,
        //     zone2End) where an iterator-invalidation footgun lives; a
        //     single erase of the tail range keeps it obvious.)
        Front::iterator newEnd = std::remove_if(zone2Begin, zone2End,
            [&](const Entry& e) {
                return std::get<1>(e) >= a && std::get<2>(e) >= b;   // p dominates e
            });
        F1.erase(newEnd, zone2End);

        // == Zone 3: keys > S — can ONLY be dominated by the newcomer ==
        // (If one dominated p, its key would be <= S.) Start recomputed:
        // the erase above shifted the tail, so old iterators are stale.
        for (Front::iterator it = std::upper_bound(F1.begin(), F1.end(), probe, ByS());
             it != F1.end(); )
            if (std::get<1>(*it) >= a && std::get<2>(*it) >= b)
                it = F1.erase(it);          // vector::erase returns next (C++11)
            else
                ++it;

        // == Commit the newcomer ==
        // Full-tuple lower_bound: any slot inside the S-tie run is valid.
        F1.insert(std::lower_bound(F1.begin(), F1.end(), probe),
                  std::make_tuple(S, a, b, (x==nullptr||xN<=0 ?
            std::vector<double>() : std::vector<double>(x, x + xN) )));

        // == HV: one sweep in the required order (ascending a) ==
        // F1 is an antichain, so no coverage test: preload point 0 and
        // telescope slabs from i = 1.
        std::vector<std::pair<double, double> > pts;
        pts.reserve(F1.size());
        for (Front::const_iterator e = F1.begin(); e != F1.end(); ++e)
            pts.push_back(std::make_pair(std::get<1>(*e), std::get<2>(*e)));
        std::sort(pts.begin(), pts.end());              // by a ascending

        hv_ = 0.0;
        if (!pts.empty()) {                             // degenerate guard
            double prevA = pts[0].first;                // preloaded sweep state
            double prevB = pts[0].second;
            for (std::size_t i = 1; i < pts.size(); ++i) {
                // slab between prevA and a_i, covered up to height prevB
                hv_ += (pts[i].first - prevA) * (refY_ - prevB);
                prevA = pts[i].first;
                prevB = pts[i].second;
            }
            hv_ += (refX_ - prevA) * (refY_ - prevB);   // trailing slab
        }
        return true;
    }

    double hypervolume() const { return hv_; }

    // THE front: vector sorted by S. begin() -> end() is the ranking itself,
    // ascending S = from best-balanced toward the extremes.
    Front F1;

private:
    double refX_, refY_;
    double hv_;
};

#endif // BITEHV2_INCLUDED
