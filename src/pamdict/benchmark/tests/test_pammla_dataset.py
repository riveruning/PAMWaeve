from scripts.build_pammla_benchmark import components,splits


def test_hamming_components_keep_transitive_neighbors_together():
    names=['AAAAAA','CCAAAA','CCCCAA','TTTTTT']
    assert components(names)==[['AAAAAA','CCAAAA','CCCCAA'],['TTTTTT']]


def test_pilot_neighbors_stay_development_and_order_invariant():
    names=['AAAAAA','CCAAAA','CCCCAA','TTTTTT','SSSSSS']
    result=splits(names,{'AAAAAA'})
    assert result==splits(list(reversed(names)),{'AAAAAA'})
    assert all(result[n][1]=='development' for n in names[:3])


def test_duplicate_variant_codes_rejected():
    try:components(['AAAAAA','AAAAAA'])
    except ValueError:pass
    else:raise AssertionError('Duplicate accepted')


def test_quantitative_control_and_constant_predictions():
    from scripts.eval_pammla_benchmark import score_vectors
    target=list(range(64))
    identical=score_vectors(target,target,target)
    assert identical['delta_spearman']==0 and identical['delta_top5_overlap']==0
    constant=score_vectors([0]*64,target,[0]*64)
    assert constant['spearman'] is None and constant['delta_spearman'] is None
    assert abs(constant['top5_overlap']-5/64)<1e-12


def test_development_scope_excludes_reserve():
    from scripts.eval_pammla_benchmark import select_scope
    rows=[{'variant_id':'DSGERT','future_task_split':'development'},
          {'variant_id':'AAAAAA','future_task_split':'reserved_evaluation'}]
    assert select_scope(rows,'development')==rows[:1]
    assert select_scope(rows,'all')==rows
