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

VEP joins multiple consequences of one transcript with `&`, for example `intron_variant&non_coding_transcript_variant`. Set `'&'` as the separator so that every consequence is considered. With `record_rule = max` the consequence with the highest `priority` is then used, not only the first one:

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

### Several numeric values (`record_rule` / `record_rule_on`)

When a `data_type = integer` or `float` annotation has several values, for example `CLNSIG=4|5`, `record_rule` (`min` or `max`, default `max`) chooses which one to score. By default it compares the **values**: `max` takes the highest value and `min` the lowest. This suits annotations where a higher value means a different kind of variant, such as an allele frequency, where `max` scores the most common frequency.

`record_rule_on` changes what is compared:

- `record_rule_on = value` (default): the highest or lowest annotation value is chosen, then scored.
- `record_rule_on = score`: every value is scored first, and the highest or lowest **score** is used. This is what you want when the score does not follow the value, for example when value `4` scores `100` but value `5` scores `1`.

```
[CLNSIG]
  field = INFO
  data_type = integer
  category = Clinical_significance
  info_key = CLNSIG
  separators = '|',
  record_rule = max
  record_rule_on = score

  [[likely_pathogenic]]
    score = 100
    value = 4

  [[pathogenic]]
    score = 1
    value = 5
```

For `CLNSIG=4|5` this scores `100`. Without `record_rule_on = score` it scores `1`, since the highest value `5` is chosen.

`record_rule_on` only affects `integer` and `float` plugins. For `data_type = string`, `record_rule` always compares the `priority` of the matching rules.

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