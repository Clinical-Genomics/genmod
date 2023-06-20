# Score Variant

## Rank Score Normalization

The rank score is MAXMIN normalized into range (0, 1) according to the following formula:

```
RankScoreNormalized = (RankScore - CategorySumMin) / (CategorySumMax - CategorySumMin)
```
where `RankScore` is the sum of rank score across categories (including rules such as min, max, sum etc)
`RankScore = SUM(Score_category_n) for 0...n categories`
and `CategorySumMin` is the sum of minimal score values for all categories,
i. e `CategorySumMin = SUM(CategoryMin_n) for 0...n categories`.
The same applies to `CategorySumMax = SUM(CategoryMax_n) for 0...n categories`.

Refer to `score_variants.py::score()` method for implementation details.

Additionally, also read in the `score-compounds.md` on compound scoring step that affects
final rank score values.