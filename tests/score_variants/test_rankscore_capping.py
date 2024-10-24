from genmod.score_variants.cap_rank_score_to_min_bound import cap_rank_score_to_min_bound, MIN_SCORE_NORMALIZED


MIN_SCORE: float = -5.0


def test_rankscore_normalized_capping():
    """
    Test the MIN normalization bounds capping of rankscore normalized.
    """
    # GIVEN a normalized rank score
    # WHEN running cap method
    # THEN expect rank score to be larger than min bound
    for rank_score_normalized in range(-10, 10):
        assert cap_rank_score_to_min_bound(rank_score_type='RankScoreNormalized',
                                           rank_score=float(rank_score_normalized),
                                           min_rank_score_value=MIN_SCORE_NORMALIZED) >= MIN_SCORE_NORMALIZED


def test_rankscore_capping():
    """
    Test the MIN normalization bounds capping of rankscore.
    """

    # GIVEN a rank score
    # WHEN running cap method
    # THEN expect rank score to be larger than min bound
    for rank_score in range(-10, 10):
        assert cap_rank_score_to_min_bound(rank_score_type='RankScore',
                                           rank_score=rank_score,
                                           min_rank_score_value=MIN_SCORE) >= MIN_SCORE
