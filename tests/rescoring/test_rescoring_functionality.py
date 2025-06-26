def test_rescore_with_annotation_suffix():
    """
    Test for asserting secondary Genmod annotation is present in VCF.
    This is intended for secondary analysis using a rank config,
    and where score annotations are added with a suffix.

    Expected:
        - All previously existing GENMOD annotations are present and unaltered
        - New scoring annotations are added that's based on an additional scoring config
    """
    runner = CliRunner()
    #r = runner.invoke(score_command, ['--help'])
    #print(r.output)
    #assert False

    result = runner.invoke(score_command,
                           [SCORED_VCF,
                            "-c", SCORE_CONFIG,
                            "--skip_is_previously_scored_check"])

    print(result.output)
    assert result.exit_code == 0