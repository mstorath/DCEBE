"""Test configuration for the DCEBE Python port.

The repository ships its source under ``python/dcebe`` rather than at the
root, so ``[tool.pytest.ini_options].pythonpath`` in pyproject.toml adds
``python/`` to ``sys.path`` for test discovery without an editable
install.
"""
