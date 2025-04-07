# Score Compounds

This module performs ranking of compound variants.

> [!WARNING]
> Ranking of compound variants is only done for the first family in the VCF.

During the ranking of these compounds the rank score might be modified in place.
See `genmod/score_variants/compound_scorer.py:L248`.

## Rankscore Capping
Since the rank scores are modified in place in this module, there's a risk
that the modified rank score might fall outside the valid range of normalization
bounds `(MIN, MAX)` that was established in the `score_variants` module.

This applies to variants belonging to the lower range of rank scores.

When this happens, the modified rank score is capped to  `(MIN, )` if it's of `RankScore` type
or `(0, 1)` if it's of `RankScoreNormalized` type.

In previous Genmod versions there were no such capping rule in effect.
Earlier ranked variants from `compounds` module might show lower rank
scores compared to this implementation.
