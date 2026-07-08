.PHONY: dev install-uv install-deps install-hooks test lint format build

# Development environment setup
dev: install-uv install-deps install-hooks
	@echo "Development environment setup complete!"

# Run pytest
test:
	uv run pytest

# Run linting
lint:
	uv run ruff check .

# Run formatting
format:
	uv run ruff format .

# Build the distribution
build:
	uv build

install-uv:
	@if command -v uv >/dev/null 2>&1; then \
		echo "uv is already installed"; \
	else \
		echo "Installing uv..."; \
		curl -LsSf https://astral.sh/uv/install.sh | sh; \
	fi

install-deps:
	@echo "Installing Python dependencies with uv..."
	uv sync --group dev

install-hooks:
	@echo "Installing pre-commit hooks..."
	uv run pre-commit install
	@echo "Pre-commit hooks installed and active!"
