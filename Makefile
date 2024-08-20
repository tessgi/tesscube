.PHONY: all clean pytest coverage flake8 black mypy isort

CMD:=poetry run
PYMODULE:=src
TESTS:=tests

# Run all the checks which do not change files
all: ruff pytest

# Run the unit tests using `pytest`
pytest:
	$(CMD) pytest $(PYMODULE) $(TESTS)

# run ruff linter
ruff:
	$(CMD) ruff format --check $(PYMODULE) $(TESTS)
