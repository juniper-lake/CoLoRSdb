.DEFAULT_GOAL := help

.PHONY: help  ## Display this message
help:
	@grep -E \
		'^.PHONY: .*?## .*$$' $(MAKEFILE_LIST) | \
		sort | \
		awk 'BEGIN {FS = ".PHONY: |## "}; {printf "\033[36m%-19s\033[0m %s\n", $$2, $$3}'

.PHONY: install  ## Create a virtual environment and install dependencies
install:
	@echo "Installing project in dev mode, hope you're using python 3.8+ for formatting compatibility with walrus operator."
	@echo "Creating a virtualenv and installing dependencies ..."
	rm -rf .venv
	python3 -m venv .venv
	./.venv/bin/pip install -U pip
	./.venv/bin/pip install -r requirements-dev.txt
	./.venv/bin/pre-commit install --install-hooks

.PHONY: update  ## Update dependencies in virtual environment
update:
	./.venv/bin/pip install -r requirements-dev.txt

.PHONY: format  ## Auto-format staged python files with ruff
format:
	./.venv/bin/pre-commit run --hook-stage manual ruff-fix --all-files
	./.venv/bin/pre-commit run --hook-stage manual ruff-format-fix --all-files

.PHONY: lint  ## Lint staged python files with ruff
lint:
	./.venv/bin/pre-commit run ruff --all-files
	./.venv/bin/pre-commit run ruff-format --all-files

.PHONY: lint-wdl  ## Lint WDL files with miniwdl
lint-wdl:
	./.venv/bin/pre-commit run miniwdl-check --all-files

.PHONY: codespell  ## Spell-check staged files with codespell
codespell:
	./.venv/bin/pre-commit run codespell --all-files

.PHONY: mypy  ## Perform type-checking on staged files with mypy
mypy:
	./.venv/bin/pre-commit run mypy --all-files

.PHONY: test-wdl ## Run full WDL test
test-wdl:
	.venv/bin/miniwdl run --verbose --cfg miniwdl.cfg --dir miniwdl_execution wdl/workflows/colors_main.wdl cohort_id=test_small sample_sheet=../sample_sheets/test.small.tsv reference_bundle=../colorsdb.v1.0.1.resources.tgz backend=HPC preemptible=false

.PHONY: clean  ## Clean unused files
clean:
	@find ./ -name '*.pyc' -exec rm -f {} \;
	@find ./ -name '__pycache__' -exec rm -rf {} \;
	@find ./ -name 'Thumbs.db' -exec rm -f {} \;
	@find ./ -name '*~' -exec rm -f {} \;
	@rm -rf .cache
	@rm -rf .pytest_cache
	@rm -rf .ruff_cache
	@rm -rf .mypy_cache
	@rm -rf build
	@rm -rf dist
	@rm -rf *.egg-info
	@rm -rf htmlcov
	@rm -rf .coverage
	@rm -rf coverage.xml
	@rm -rf .tox/
	@rm -rf docs/_build
