"""Deterministic categorical colours shared by every cohort report view."""

import json
from functools import lru_cache
from pathlib import Path

import distinctipy

SEED = 42


@lru_cache(maxsize=16)
def _colors(count):
    if count < 0:
        raise ValueError("Colour count must be nonnegative")
    return tuple(distinctipy.get_colors(count, exclude_colors=[(1, 1, 1), (0, 0, 0)], rng=SEED))


def categorical_colors(count):
    return list(_colors(count))


def method_palette(labels):
    ordered = sorted(set(labels))
    return dict(zip(ordered, (distinctipy.get_hex(color) for color in _colors(len(ordered)))))


def colors_for_labels(labels, palette=None):
    palette = method_palette(labels) if palette is None else palette
    # Missing labels indicate an incomplete cohort map; never cycle or recolour.
    return [palette[label] for label in labels]


def save_palette(output_dir, palette):
    payload = dict(schema="amber-method-palette-v1", generator="distinctipy",
                   version=distinctipy.__version__, seed=SEED,
                   excluded_rgb=[[1, 1, 1], [0, 0, 0]], colors=dict(sorted(palette.items())))
    Path(output_dir, "method_palette.json").write_text(json.dumps(payload, indent=2) + "\n")
