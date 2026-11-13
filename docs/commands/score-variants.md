# Score Variant

## Rank model rules

### Numeric ranges (`lower` / `upper`)

For numeric plugins, each rule covers a half-open interval `[lower, upper)`:

- `lower` is **inclusive** (`value >= lower`)
- `upper` is **exclusive** (`value < upper`)

When two rules meet at the same value, the rule that has it as `lower` is used. In the example below a value of exactly `0.02` scores `-12` (`common`), and `0.0199` scores `1` (`intermediate`):

```
  [[common]]
    score = -12
    lower = 0.02
    upper = 1.1

  [[intermediate]]
    score = 1
    lower = 0.005
    upper = 0.02
```

Intervals should not overlap. If they do, which rule is used is not defined.

### String rules and multiple values (`separators`)

A `data_type = string` rule is matched against whole values, not parts of them. The annotation is split into values using `separators`, applied in the order listed, and each rule `string` is compared to those values.

VEP joins multiple consequences of one transcript with `&`, for example `intron_variant&non_coding_transcript_variant`. Set `'&'` as the separator so that every consequence is considered. With `record_rule = max` the consequence with the highest `priority` / score is then used, not only the first one:

```
[most_severe_consequence]
  category = consequence
  data_type = string
  description = Most severe consequence for this variant
  field = INFO
  info_key = CSQ
  csq_key = Consequence
  record_rule = max
  separators = '&'
```

For `info_key = CSQ`, genmod itself splits the transcripts (`,`) and the VEP columns (`|`) and picks the column named by `csq_key`. `separators` only describes how to split the value *within* that column, so `,` and `|` should not be listed for a CSQ plugin.

If `'&'` is missing, `intron_variant&non_coding_transcript_variant` is one single value that matches none of the rules, and the variant gets no score for the plugin. Setting `exact_match = false` on the plugin restores the older substring matching, which finds the consequences without `&` as a separator, but then a rule such as `AD` also matches `AD_dn`.

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