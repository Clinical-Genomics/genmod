from typing import Dict, List

# Types of rank scores provided by scoring, raw and normalized rank scores
RANK_SCORE_TYPES: Dict[str, str]= {
    'RankScore': 'The rank score for this variant in this family. family_id:rank_score.',
    'RankScoreNormalized': 'The normalized rank score in range(0, 1) for this variant in this family. family_id:rank_score.'
}
RANK_SCORE_TYPE_NAMES: List[str] = list(RANK_SCORE_TYPES.keys())
