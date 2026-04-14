"""Shared pytest fixtures."""
import pytest
from pathlib import Path

@pytest.fixture
def project_root() -> Path:
    return Path(__file__).parent.parent

@pytest.fixture
def test_data(project_root: Path) -> Path:
    d = project_root / "tests" / "data"
    d.mkdir(exist_ok=True)
    return d
