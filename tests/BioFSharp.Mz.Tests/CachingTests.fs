module CachingTests

open Expecto
open BioFSharp.Mz

let mkCache () =
    let cache = Cache.createCache<int,string>
    Cache.addItem cache (20,"b")
    Cache.addItem cache (10,"a")
    Cache.addItem cache (30,"c")
    cache

[<Tests>]
let tests =
    testList "CachingTests" [
        testCase "addItem inserts keys in ascending sorted order" <| fun _ ->
            let cache = mkCache()
            Expect.equal cache.Count 3 "cache contains all three items"
            Expect.equal (List.ofSeq cache.Keys) [10; 20; 30] "keys are sorted in ascending order"
            // SortedList maintains ascending key order regardless of insertion order (BCL contract); addItem delegates to Add for new keys.

        testCase "addItem replaces the value of an existing key without changing the count" <| fun _ ->
            let cache = Cache.createCache<int,string>
            Cache.addItem cache (3,"c")
            Cache.addItem cache (3,"z")
            Expect.equal cache.Count 1 "re-adding a key keeps one entry"
            Expect.equal (Cache.getItemBy cache 3) (true, "z") "the newest value is retained"
            // documented add-or-update semantics - re-adding a key must keep one entry holding the newest value.

        testCase "getItemBy distinguishes present and missing keys" <| fun _ ->
            let cache = mkCache()
            Expect.equal (Cache.getItemBy cache 20) (true, "b") "present keys return their value"
            Expect.isFalse (fst (Cache.getItemBy cache 15)) "missing keys return a false success flag"
            // TryGetValue contract - success flag plus value for hits, false for misses.

        testCase "bulkInsertBy inserts every supplied pair" <| fun _ ->
            let cache = Cache.createCache<int,string>
            Cache.bulkInsertBy cache [(1,"x");(2,"y");(3,"z")] |> ignore
            Expect.equal cache.Count 3 "cache contains every supplied pair"
            Expect.equal (Cache.getItemBy cache 1) (true, "x") "key 1 is retrievable with its value"
            Expect.equal (Cache.getItemBy cache 2) (true, "y") "key 2 is retrievable with its value"
            Expect.equal (Cache.getItemBy cache 3) (true, "z") "key 3 is retrievable with its value"
            // documented bulk insert - every pair must land in the cache.

        testCase "binarySearch Border.Upper resolves each existing key to its sorted index" <| fun _ ->
            let cache = mkCache()
            Expect.equal (Cache.binarySearch Cache.Border.Upper cache 10) 0 "10 is at sorted index 0"
            Expect.equal (Cache.binarySearch Cache.Border.Upper cache 20) 1 "20 is at sorted index 1"
            Expect.equal (Cache.binarySearch Cache.Border.Upper cache 30) 2 "30 is at sorted index 2"
            // keys sorted ascending are [10;20;30], so their indices are 0,1,2 by definition of sorted
            // position. Note: only 20 is found by the search loop itself; 10 and 30 terminate as misses
            // and are resolved to the correct index by the Border.Upper miss fallback - the expectation is
            // contractually right for Border.Upper but does NOT extend to Border.Lower (see the pending
            // Border.Lower test).

        testCase "binarySearch Border.Upper clamps a probe above all keys to the last index" <| fun _ ->
            Expect.equal (Cache.binarySearch Cache.Border.Upper (mkCache()) 99) 2 "probe above all keys clamps to the last index"
            // this clamp is the contract relied upon by SearchEngineGeneric (which feeds the result directly into getValuesByIdx as an upper index when a query range extends past the cache); any sensible corrected implementation (last key <= probe, or first key >= probe clamped to Count-1) also returns 2 for probe 99, so the assertion is fix-compatible, not bug-locking.

        // PENDING: for keys [10;20;30;40;50], the existing key 20 sits at sorted index 1, and a
        // search that claims to locate keys must return it for either border. The early-terminating
        // loop misses the key and the Border.Lower fallback subtracts one, returning 0. On a 3-key
        // cache the loop happens to compare key 20 directly, so the bug does not manifest there -
        // the 5-key fixture is required.
        ptestCase "binarySearch Border.Lower resolves an existing key to its sorted index" <| fun _ ->
            let cache = Cache.createCache<int,string>
            Cache.addItem cache (10,"a")
            Cache.addItem cache (20,"b")
            Cache.addItem cache (30,"c")
            Cache.addItem cache (40,"d")
            Cache.addItem cache (50,"e")
            Expect.equal (Cache.binarySearch Cache.Border.Lower cache 20) 1 "20 is at sorted index 1"

        testCase "getItemsByIdx returns the inclusive index range" <| fun _ ->
            let cache = mkCache()
            Expect.equal (Cache.getItemsByIdx cache (0,1)) [(10,"a");(20,"b")] "indices 0 through 1 are included"
            // documented inclusive index range over the sorted list.

        testCase "getValuesByIdx returns the values of the inclusive index range" <| fun _ ->
            Expect.equal (Cache.getValuesByIdx (mkCache()) (1,2)) ["b";"c"] "values at indices 1 through 2 are included"
            // same inclusive range contract, values only.

        testCase "getItemsByRange with existing endpoint keys returns exactly the items in that key range" <| fun _ ->
            Expect.equal (Cache.getItemsByRange (mkCache()) (20,30)) [(20,"b");(30,"c")] "both existing endpoints and their range are included"
            // documented inclusive key range; both endpoints are existing keys so no gap-resolution semantics are involved.

        testCase "containsItemsBetween reports the index bounds of a range holding two keys" <| fun _ ->
            Expect.equal (Cache.containsItemsBetween (mkCache()) (15,35)) (Some (1,2)) "keys 20 and 30 occupy indices 1 and 2"
            // keys 20 and 30 lie inside [15,35] and sit at sorted indices 1 and 2 (hand-computed from the sorted key order [10;20;30]).

        testCase "containsItemsBetween index bounds are minimal and maximal for existing endpoints" <| fun _ ->
            let cache = Cache.createCache<int,string>
            Cache.addItem cache (10,"a")
            Cache.addItem cache (20,"b")
            Cache.addItem cache (30,"c")
            Cache.addItem cache (45,"d")
            Cache.addItem cache (60,"e")
            Cache.addItem cache (80,"f")
            Expect.equal (Cache.containsItemsBetween cache (20, 60)) (Some (1, 4)) "existing endpoints resolve to their exact sorted indices"
            Expect.equal (Cache.getValuesByIdx cache (1, 4)) ["b";"c";"d";"e"] "the inclusive bounds return every value in the endpoint range"
            // both endpoints are existing keys, so by inclusive-interval membership the bounds must be exactly their sorted indices - this pins lower-bound MINIMALITY, which the general invariant test (upper maximality only) cannot: an implementation skipping the first qualifying key would pass it. Consumers feed these indices straight into getValuesByIdx; a high lower bound silently drops candidates.

        testCase "a custom comparer governs ordering and search" <| fun _ ->
            let reverseComparer = { new System.Collections.Generic.IComparer<int> with member _.Compare(a, b) = compare b a }
            let cache = Cache.createCacheWith reverseComparer 1
            Cache.addItem cache (20,"b")
            Cache.addItem cache (10,"a")
            Cache.addItem cache (30,"c")
            Expect.equal (List.ofSeq cache.Keys) [30;20;10] "the supplied comparer orders keys in descending order"
            Expect.equal (Cache.binarySearch Cache.Border.Upper cache 30) 0 "30 is at sorted index 0 under the custom comparer"
            Expect.equal (Cache.binarySearch Cache.Border.Upper cache 20) 1 "20 is at sorted index 1 under the custom comparer"
            Expect.equal (Cache.binarySearch Cache.Border.Upper cache 10) 2 "10 is at sorted index 2 under the custom comparer"
            Expect.equal cache.Count 3 "all items are retained after automatic growth"
            // the supplied comparer dictates sorted order (SortedList contract) and binarySearch must consult the cache's own comparer; growing past the initial capacity of 1 also exercises documented automatic growth without pinning a capacity value.

        testCase "containsItemsBetween returns None for ranges holding no keys" <| fun _ ->
            let cache = mkCache()
            Expect.equal (Cache.containsItemsBetween cache (12,18)) None "a gap between keys contains no cached key"
            Expect.equal (Cache.containsItemsBetween cache (40,50)) None "a range above all keys contains no cached key"
            Expect.equal (Cache.containsItemsBetween cache (0,5)) None "a range below all keys contains no cached key"
            // these ranges contain no cached key, and production (SearchDB/SearchEngineGeneric) uses Some/None as the sole guard deciding whether cached data is served - a false Some silently yields wrong peptide lists.

        testCase "containsItemsBetween index bounds always point at keys inside the range" <| fun _ ->
            let cache = Cache.createCache<int,string>
            Cache.addItem cache (10,"a")
            Cache.addItem cache (20,"b")
            Cache.addItem cache (30,"c")
            Cache.addItem cache (45,"d")
            Cache.addItem cache (60,"e")
            Cache.addItem cache (80,"f")
            let ranges = [(15,50); (10,80); (25,70); (55,90)]
            for range in ranges do
                let lowerMass, upperMass = range
                match Cache.containsItemsBetween cache range with
                | Some (lo, hi) ->
                    for idx in lo .. hi do
                        let key = cache.Keys.[idx]
                        Expect.isTrue (key >= lowerMass && key <= upperMass) (sprintf "key %d at index %d lies inside range [%d,%d]" key idx lowerMass upperMass)
                    Expect.isTrue (cache.Keys.[hi] <= upperMass) (sprintf "key at upper index %d is at or below upper bound %d" hi upperMass)
                    Expect.isTrue (hi = cache.Count - 1 || cache.Keys.[hi + 1] > upperMass) (sprintf "upper index %d is maximal for upper bound %d" hi upperMass)
                | None ->
                    failtestf "range [%d,%d] must report at least two cached keys" lowerMass upperMass
            // SearchDB reads Keys.[upperMassIdx] and re-queries the database from that key upward, so a hi bound pointing at a key outside the range (or not the maximal in-range key) corrupts the re-query window - this invariant is what the caller stakes correctness on.

        testCase "containsItemsBetween returns None for an empty cache" <| fun _ ->
            let cache = Cache.createCache<int,string>
            Expect.equal (Cache.containsItemsBetween cache (0,100)) None "an empty cache has no contained items"
            // an empty cache contains no items in any range.

        // PENDING: the range [15,25] contains exactly the key 20, and the doc says border items are included, so the function should report it. The implementation demands strictly different lower/upper indices and returns None for single-key ranges.
        ptestCase "containsItemsBetween reports a range containing exactly one key" <| fun _ ->
            Expect.equal (Cache.containsItemsBetween (mkCache()) (15,25)) (Some (1,1)) "the single contained key 20 occupies both bounds"

        // PENDING: only key 20 lies inside [15,25], but the implementation returns ALL THREE items:
        // the lower probe resolves the missing bound to index 0 (leaking key 10 below the range) and
        // the upper probe clamps to the last index (leaking key 30 above it). Both ends leak.
        ptestCase "getItemsByRange excludes keys outside the range" <| fun _ ->
            Expect.equal (Cache.getItemsByRange (mkCache()) (15,25)) [(20,"b")] "keys below the lower bound are excluded"

        // PENDING: a range query over an empty collection must be empty (as containsItemsBetween
        // already answers None); the implementation dereferences index 0 and throws
        // ArgumentOutOfRangeException.
        ptestCase "getItemsByRange over an empty cache is empty" <| fun _ ->
            Expect.equal (Cache.getItemsByRange (Cache.createCache<int,string>) (0,5)) [] "an empty cache has no items in the range"

        // PENDING: doc says items with keys smaller than the cutoff are deleted, so [20;30] must survive. The implementation adds the survivors back into the ORIGINAL cache (duplicate-key ArgumentException on the first re-add) and returns the empty replacement cache.
        // The gap cutoff 15 must yield the same survivors - a fix that keys off binarySearch Upper 15 (= index 0) would wrongly retain key 10.
        ptestCase "bulkDeleteBy removes items below the cutoff and returns the surviving items" <| fun _ ->
            let cleaned = Cache.bulkDeleteBy (mkCache()) 20
            Expect.equal (List.ofSeq cleaned.Keys) [20;30] "items at and above the cutoff survive"
            let cleanedGap = Cache.bulkDeleteBy (mkCache()) 15
            Expect.equal (List.ofSeq cleanedGap.Keys) [20;30] "items at and above a gap cutoff survive"
    ]
