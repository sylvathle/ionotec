import importlib
import pkgutil

import pytest
import ionotec

MODULES = [m.name for m in pkgutil.walk_packages(ionotec.__path__, prefix="ionotec.")]


@pytest.mark.parametrize("name", MODULES)
def test_module_imports(name):
    importlib.import_module(name)
