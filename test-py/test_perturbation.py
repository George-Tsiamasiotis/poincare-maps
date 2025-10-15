def test_perturbation_fields(perturbation):
    harmonics = perturbation.get_harmonics()
    assert isinstance(harmonics, list)
    assert len(harmonics) == 2
