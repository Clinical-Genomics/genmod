This PR adds a functionality or fixes a bug.
OR
This PR marks a new Genmod release. We apply semantic versioning. This is a major/minor/patch release for reasons.

Test on your local machine and describe how to do it.
Also consider if it would be useful to test the changes in a staging environment, on e.g. `hasta` or `hasta-login`.
If so, document the tests here, but think once more about if it is ok to paste details about case names, pseudonyms or variants.
This is fine for e.g. public test data, such as GIAB or NIST, or synthetic data, but likely not for real patient data.

```
uv run genmod annotate -o examples/test_vcf_annotated.vcf --annotate_regions examples/test_vcf.vcf
uv run genmod models -f examples/recessive_trio.ped -o examples/test_vcf_recessive_models_annotated.vcf examples/test_vcf_annotated.vcf
```

<details>
<summary>Testing docs locally with `uv`</summary>
1. Build and serve documents: `uv run --group docs mkdocs serve --strict`
1. Visit [http://127.0.0.1:4000/scout/](http://127.0.0.1:4000/scout/)
</details>

**How to test**:
1. how to test, possibly with real cases/data

**Expected outcome**:
The functionality should be working
Take a screenshot and attach or copy/paste the output.

**Review:**
- [ ] code approved by
- [ ] tests executed by
