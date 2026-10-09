import numpy as np

from liposome_pipeline import umbrella


def test_chunk_evenly_balances_and_drops_empty():
    chunks = umbrella.chunk_evenly(list(range(7)), 3)
    assert [len(c) for c in chunks] == [3, 2, 2]
    assert umbrella.chunk_evenly([1, 2], 5) == [[1], [2]]


def test_target_distances_cover_range_without_duplicates():
    targets = umbrella.target_distances(1.03, 5.0)
    assert targets == sorted(set(targets))
    assert targets[0] == 1.1
    assert 2.0 in targets and 4.0 in targets
    assert max(targets) <= 5.0


def test_select_frames_matches_nearest_and_limits_distance():
    times = np.arange(0, 1000.0)
    distances = np.linspace(1.0, 12.0, 1000)
    frames = umbrella.select_frames(times, distances)
    assert frames and all(d < umbrella.MAX_WINDOW_DISTANCE_NM for d, _ in frames)
    for target, time in frames:
        assert abs(distances[int(time)] - target) < 0.02
    assert len({t for _, t in frames}) == len(frames)


def test_read_xvg_skips_headers(tmp_path):
    path = tmp_path / "d.xvg"
    path.write_text("# c\n@ s\n0 1.5 9\n10 2.5 9\n")
    times, dists = umbrella.read_xvg(path)
    assert list(times) == [0, 10] and list(dists) == [1.5, 2.5]


def test_window_script_uses_minimisation_flags():
    em = umbrella.window_script([], "a.mdp", "t.top", "i.ndx", True)
    md = umbrella.window_script([], "a.mdp", "t.top", "i.ndx", False)
    assert "-bonded gpu" not in em and "-bonded gpu" in md
