# Pre-PyPI publication checklist

Hand-off list for taking `dcebe` from the current branch state to a published `pip install dcebe` package.

The agent (Claude) **does not** open PRs, push branches, edit `.github/`, or upload to PyPI — these all require the maintainer's hand on the wheel. The items below are framed as concrete maintainer actions.

## 1. Get the bug-fix and port branches reviewed and merged

The port branch `claude/dcebe-port-2026-05` is based on `claude/dcebe-prep-2026-05`. The latter is a prerequisite (it carries the search-interval bug fix that the Python's `default_interval` mirrors).

Recommended order:

1. `git push -u origin claude/dcebe-prep-2026-05` and open a PR against `master`. Merge.
2. `git push -u origin claude/dcebe-port-2026-05` and open a PR against `master`. Merge.

Both branches are local-only at the time of this checklist's writing.

## 2. Wire up CI (one-time)

Copy the two templates into `.github/workflows/`:

```bash
mkdir -p .github/workflows
cp docs/ci_templates/ci.yml .github/workflows/ci.yml
cp docs/ci_templates/release.yml .github/workflows/release.yml
git add .github/workflows
git commit -m "Add CI workflows"
```

Then push and verify the CI badge turns green on the next PR.

## 3. Configure PyPI OIDC trusted publishing (one-time)

The `release.yml` template uses OIDC instead of a PyPI API token. Setup:

1. Log in to https://pypi.org and visit *Account → Publishing*.
2. Add a new "pending publisher":
   - PyPI Project Name: `dcebe`
   - Owner: `mstorath`
   - Repository name: `DCEBE`
   - Workflow name: `release.yml`
   - Environment name: `pypi`
3. In GitHub: *Settings → Environments → New environment* → name `pypi` (no protection rules required, but optional reviewer approval recommended).

Once configured, the `Release` workflow will publish on any `v*.*.*` tag without an API token.

## 4. Reserve the `dcebe` name on PyPI (optional but recommended)

To prevent name-squatting, push a placeholder `0.1.0.dev0` build to TestPyPI before the real release:

```bash
# From a clean checkout, with a `release.yml.dryrun` variant pointing at
# https://test.pypi.org/legacy/, or via twine directly:
python -m build
twine upload --repository-url https://test.pypi.org/legacy/ dist/*
```

Then verify with `pip install -i https://test.pypi.org/simple/ dcebe` in a fresh venv.

## 5. Tag and push v0.1.0

Once steps 1–3 are done and CI is green:

```bash
git tag -a v0.1.0 -m "DCEBE Python v0.1.0 — initial release"
git push origin v0.1.0
```

The `Release` workflow will pick up the tag, build sdist + wheel, and publish to PyPI.

## 6. Smoke-test the public install

In a fresh venv:

```bash
pip install dcebe
python -c "from dcebe import estimate_bat; print(estimate_bat.__doc__[:200])"
python -c "import dcebe; print(dcebe.__version__)"
```

If these succeed, the publication path is healthy.

## 7. Update README.md to mention the Python package

Currently `README.md` (the MATLAB-facing one) doesn't mention `dcebe` on PyPI. Add a one-line pointer near the top:

> **Python users**: `pip install dcebe` — see [README_PYTHON.md](README_PYTHON.md).

## What the agent left in place

- `pyproject.toml` with hatchling backend, `dcebe` package name, version `0.1.0.dev0`, dependencies pinned to NumPy ≥ 1.23 and SciPy ≥ 1.10, classifiers for cp3.9–cp3.13.
- `python/dcebe/{__init__,_spline,_search,_barriers,_gcv,_factorize,_solve}.py` — full Python port.
- `tests_py/` — 81 unit tests (1 skipped placeholder), including 4 MATLAB parity tests against `matlab_fixtures/`.
- `matlab_fixtures/fixture_{small,explicit_si,common}.mat` — committed reference fixtures.
- `scripts/gen_fixtures.m` — regen tool for the fixtures.
- `demos_py/{demo,demo_common,demo_in_vivo}.py` — runnable demonstrations.
- `README_PYTHON.md` — package-level docs.
- `CITATION.cff` (added in `claude/dcebe-prep-2026-05`).
- `docs/ci_templates/{ci,release}.yml` — copy-and-paste CI.
- `docs/PRE_PYPI_CHECKLIST.md` — this file.

## What the agent deliberately left out

- Any file under `.github/` (CLAUDE.md forbids editing CI workflows directly).
- Any push to `origin` (CLAUDE.md requires the maintainer to push).
- Any PyPI upload, TestPyPI upload, or git tag (deliberate — irreversible, requires owner judgement).
- Performance optimisation (pure-Python is functionally correct but per-call overhead is higher than MATLAB; see the comment in `_solve.py` about `eps` and the demo runtime caveat. Future work: vectorise the per-(t, beta) inner loop in `gcv_score_qr`, or factor the matrix once per `t` and reuse across `beta`).
