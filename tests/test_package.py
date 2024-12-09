from __future__ import annotations

import importlib.metadata

import rnanotuple as m


def test_version():
    assert importlib.metadata.version("rnanotuple") == m.__version__
