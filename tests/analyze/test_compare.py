import pytest
from sisana.analyze_networks.analyze import transform_edge_to_positive_val, calc_log2_fc

@pytest.mark.parametrize("edge_val, expected_output", [
    [800, 800],
    [1, 1.313],
    [0, 0.693],
    [-10, 0.0000454]
])

def test_transform_edge_to_positive_val(edge_val: float, expected_output):
    assert transform_edge_to_positive_val(edge_val) == pytest.approx(expected_output, 0.001)

@pytest.mark.parametrize("g1, g2, expected_output", [
    [[2, 3, 5], [9, 8 , 4], 1.07039]
])

def test_calc_log2_fc(g1: list, g2: list, expected_output):
    assert calc_log2_fc(g1, g2) == pytest.approx(expected_output, 0.001)

@pytest.mark.parametrize("g1, g2", [
    [[-2, 3, 4], [4, 5, 6]] 
])

def test_calc_log2_fc_neg_value(g1: list, g2: list):
    with pytest.raises(Exception):
        calc_log2_fc(g1, g2)
