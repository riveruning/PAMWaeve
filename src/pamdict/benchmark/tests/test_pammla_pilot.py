import numpy as np
from scripts.score_pammla_pilot import ranks, spearman, overlap, top_membership
from scripts.prepare_pammla_pilot import translate


def test_average_ranks_and_correlation():
    assert np.array_equal(ranks([3,1,1,2]), [4,1.5,1.5,3])
    assert abs(spearman([1,2,3],[10,20,30])-1)<1e-12
    assert abs(spearman([1,2,3],[30,20,10])+1)<1e-12


def test_constant_is_undefined():
    assert spearman([1,1,1],[1,2,3]) is None
    assert spearman([1,2,3],[0,0,0]) is None
    assert spearman(np.maximum([-7,-6,-5],-5),[1,2,3]) is None


def test_top5_fractional_ties():
    assert abs(overlap(list(range(10)),list(range(10)))-1)<1e-12
    assert overlap(list(range(10)),list(reversed(range(10))))==0
    assert abs(overlap([0]*64,[0]*64)-5/64)<1e-12
    membership = top_membership([2,2,1,1,1,1,0],5)
    assert np.allclose(membership,[1,1,.75,.75,.75,.75,0])


def test_invalid_vectors_fail_closed():
    for fn,args in [(ranks,([1,float('nan')],)),(spearman,([1],[1,2])),
                    (top_membership,([1,2],5)),(translate,('ATGN',))]:
        try:
            fn(*args)
        except ValueError:
            pass
        else:
            raise AssertionError('Expected rejection')


def test_standard_code():
    assert translate('ATGGATTCTGGAGAAAGAACTTAA')=='MDSGERT*'


def test_representation_change_tracks_scale_and_identity():
    from scripts.diagnose_pammla_precision import vector_change
    same=vector_change([3,4],[3,4])
    assert same['l2']==0 and same['relative_l2']==0 and same['cosine']==1
    scaled=vector_change([6,8],[3,4])
    assert scaled['l2']==5 and scaled['relative_l2']==1 and scaled['cosine']==1
    zero=vector_change([0,0],[0,0])
    assert zero['relative_l2'] is None and zero['cosine'] is None
