from genmod.score_variants.score_variant import MIN_SCORE_NORMALIZED
from genmod.score_variants.rank_score_variant_definitions import RANK_SCORE_TYPE_NAMES


def cap_rank_score_to_min_bound(rank_score_type: str,
                                rank_score,
                                min_rank_score_value: float) -> float:
    """
    Caps rank_score to fall withing MIN bound of normalized rank score, if it's outside valid range.
    Args:
        rank_score_type: Type of rank score
        rank_score: The value to bounds check
        min_rank_score_value: Minimum allowed bound according to rank score normalization
    Returns:
        Bounds capped rank score, either to min_rank_score_value (if RankScore)
        or MIN_SCORE_NORMALIZED if RankScoreNormalized type.
    """

    if rank_score_type not in set(RANK_SCORE_TYPE_NAMES):
        raise ValueError(f'Unknown rank score type {rank_score_type}')

    if rank_score_type == 'RankScoreNormalized':
        min_rank_score_value = MIN_SCORE_NORMALIZED

    if rank_score < min_rank_score_value:
        return min_rank_score_value
    return rank_score
