"""Tests for the shared binning/weighting layer (lidar_diff_icp.binstats).

The two behaviours worth protecting are the ones that were being got wrong by hand: that
the cluster-robust error is genuinely larger than the naive one when data are spatially
clustered, and that binning keeps the sparse extremes instead of quietly truncating them.
"""
import numpy as np
import pytest

from lidar_diff_icp import binstats as bs


def test_nmad_matches_sigma_for_a_normal_sample():
    rng = np.random.default_rng(0)
    assert bs.nmad(rng.normal(scale=3.0, size=200000)) == pytest.approx(3.0, rel=0.02)


def test_nmad_ignores_a_wild_tail():
    a = np.concatenate([np.zeros(999), [1e9]])
    assert bs.nmad(a) == pytest.approx(0.0, abs=1e-9)


def test_block_ids_group_by_spatial_block_not_by_cell():
    nx = 100
    cell = np.array([0, 1, 9, 10, nx, nx*10])          # ix = 0,1,9,10,0,0 ; iy = 0,0,0,0,1,10
    b = bs.block_ids(cell, nx=nx, res=5.0, block_m=50.0)   # 10 cells per block
    assert b[0] == b[1] == b[2] == b[4], "same 10x10 block must share an id"
    assert b[3] != b[0], "crossing the block edge in x must change the id"
    assert b[5] != b[0], "crossing the block edge in y must change the id"


def test_quantile_edges_span_every_observation():
    rng = np.random.default_rng(1)
    x = np.concatenate([np.zeros(6000), rng.uniform(0, 0.5, 3000), [0.95]])   # skewed + lone extreme
    e = bs.quantile_edges(x, 10, first_edge=0.02)
    assert e[0] <= x.min() and e[-1] > x.max(), "edges must cover the full range"
    st = bs.binned_stats(x, np.zeros_like(x), e, min_n=1)
    assert st.n.sum() == x.size, "no observation may be dropped"
    assert st.hi[-1] > 0.95, "the lone extreme must land inside a bin, not outside"


def test_quantile_edges_isolate_a_leading_spike():
    x = np.concatenate([np.zeros(8000), np.linspace(0.05, 0.6, 2000)])
    e = bs.quantile_edges(x, 8, first_edge=0.02)
    assert e[1] == pytest.approx(0.02), "the spike gets its own bin"
    st = bs.binned_stats(x, np.zeros_like(x), e, min_n=1)
    assert st.n[0] == 8000
    assert st.n[1:].min() > 0, "remaining bins are populated"


def test_cluster_robust_se_exceeds_the_naive_one_when_clustered():
    """Each block has its own offset: the naive SE sees n returns, the truth is n blocks."""
    rng = np.random.default_rng(2)
    nblk, per = 12, 500
    block = np.repeat(np.arange(nblk), per)
    y = np.repeat(rng.normal(scale=20.0, size=nblk), per) + rng.normal(scale=1.0, size=nblk*per)
    x = np.zeros_like(y)
    st = bs.binned_stats(x, y, np.array([-1.0, 1.0]), block=block, min_n=1)
    assert st.n[0] == nblk*per and st.n_block[0] == nblk
    assert st.se_block[0] > 5 * st.se_return[0], "clustering must inflate the error"
    assert st.se[0] == st.se_block[0], "the cluster-robust SE is the one reported"


def test_cluster_robust_se_matches_naive_when_unclustered():
    rng = np.random.default_rng(3)
    y = rng.normal(scale=10.0, size=12000)
    block = np.arange(y.size) % 400                    # blocks carry no offset of their own
    st = bs.binned_stats(np.zeros_like(y), y, np.array([-1.0, 1.0]), block=block, min_n=1)
    assert st.se_block[0] == pytest.approx(st.se_return[0], rel=0.5)


def test_weights_favour_bins_with_more_independent_data():
    rng = np.random.default_rng(4)
    x = np.concatenate([np.zeros(20000), np.ones(400)])
    y = np.concatenate([rng.normal(scale=10, size=20000), rng.normal(scale=10, size=400)])
    block = np.concatenate([np.arange(20000) % 500, np.arange(400) % 10])
    st = bs.binned_stats(x, y, np.array([-0.5, 0.5, 1.5]), block=block, min_n=1)
    assert st.weights[0] > st.weights[1], "the better-sampled bin must weigh more"
    assert np.all(np.isfinite(st.weights))


def test_sparse_bins_are_kept_not_dropped():
    x = np.concatenate([np.zeros(5000), [0.9, 0.91, 0.92]])
    y = np.concatenate([np.zeros(5000), [-100.0, -110.0, -120.0]])
    e = bs.quantile_edges(x, 4, first_edge=0.02)
    st = bs.binned_stats(x, y, e, min_n=1)
    assert st.y.min() < -50, "the sparse extreme bin must survive and carry its value"
    assert st.n.sum() == x.size


def test_falls_back_to_return_se_without_blocks():
    rng = np.random.default_rng(5)
    y = rng.normal(size=5000)
    st = bs.binned_stats(np.zeros_like(y), y, np.array([-1.0, 1.0]), min_n=1)
    assert st.n_block[0] == 0 and not np.isfinite(st.se_block[0])
    assert st.se[0] == st.se_return[0]
