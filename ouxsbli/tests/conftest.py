"""Shared pytest fixtures and markers for OUxSBLI integration tests."""
import pathlib
import pytest


REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "integration: marks tests that build and run the CUDA Fortran solver "
        "(requires NVIDIA HPC SDK and a CUDA-capable GPU)",
    )


@pytest.fixture(scope="session")
def repo_root():
    return REPO_ROOT

