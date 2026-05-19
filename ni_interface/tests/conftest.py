"""Shared fixtures for the interface.cpp pytest suite.

Brian2 owns compilation: every test asks for a fresh ``cpp_standalone``
device pointed at a per-test ``tmp_path``. There is no hand-rolled g++
fixture — see /memories/session/plan.md.
"""
import ctypes
import os
import platform

import pytest

NI_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
LIBNIDAQMX_PATH = os.path.join(NI_DIR, "libnidaqmx.so")


# ---------------------------------------------------------------------------
# Hardware probe
# ---------------------------------------------------------------------------

def _probe_ni_hardware():
    """Lightweight ctypes probe: try DAQmxResetDevice("Dev1") once.

    Returns True iff the call succeeds (status code 0). Any DAQmx error or
    missing library -> False. Cached at module scope via session fixture.
    """
    try:
        lib = ctypes.CDLL(LIBNIDAQMX_PATH)
    except OSError:
        return False
    try:
        lib.DAQmxResetDevice.argtypes = [ctypes.c_char_p]
        lib.DAQmxResetDevice.restype = ctypes.c_int
        status = lib.DAQmxResetDevice(b"Dev1")
    except Exception:
        return False
    return status == 0


@pytest.fixture(scope="session")
def has_ni_hardware():
    return _probe_ni_hardware()


def pytest_collection_modifyitems(config, items):
    """Skip @pytest.mark.hardware tests when no NI card is present."""
    if _probe_ni_hardware():
        return
    skip_hw = pytest.mark.skip(reason="NI Dev1 not reachable (no hardware)")
    for item in items:
        if "hardware" in item.keywords:
            item.add_marker(skip_hw)


# ---------------------------------------------------------------------------
# PREEMPT_RT detection (jitter test gating)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def on_rt_kernel():
    return "PREEMPT_RT" in platform.uname().version


# ---------------------------------------------------------------------------
# Brian2 device fixture
# ---------------------------------------------------------------------------

@pytest.fixture
def brian2_standalone(tmp_path):
    """Yield a fresh Brian2 cpp_standalone device rooted at ``tmp_path``.

    Tests should build their network, then call ``device.build()`` (or
    rely on ``run()`` which builds implicitly with ``build_on_run=True``).
    Between tests we ``reinit()`` + ``activate()`` so global Brian2 state
    is clean.
    """
    from brian2 import device, set_device, prefs

    build_dir = tmp_path / "b2_out"
    build_dir.mkdir()
    set_device("cpp_standalone", directory=str(build_dir), build_on_run=True)
    # Make sure the rig's libnidaqmx is discoverable at link/run time even
    # when Brian2 is invoked from a tmp cwd.
    prefs.codegen.cpp.include_dirs = [NI_DIR]
    prefs.codegen.cpp.library_dirs = [NI_DIR]
    prefs.codegen.cpp.headers = ['"interface.h"', '"NIDAQmx.h"']

    yield device

    try:
        device.reinit()
        device.activate()
    except Exception:
        pass


@pytest.fixture
def brian2_standalone_debug(tmp_path):
    """Same as ``brian2_standalone`` but compiles interface.cpp with -DDEBUG.

    Used by the read_times dump test.
    """
    from brian2 import device, set_device, prefs

    build_dir = tmp_path / "b2_out_debug"
    build_dir.mkdir()
    set_device("cpp_standalone", directory=str(build_dir), build_on_run=True)
    prefs.codegen.cpp.include_dirs = [NI_DIR]
    prefs.codegen.cpp.library_dirs = [NI_DIR]
    prefs.codegen.cpp.headers = ['"interface.h"', '"NIDAQmx.h"']
    # Inject -DDEBUG; copy any preexisting flags so we don't clobber them.
    existing = list(prefs.codegen.cpp.extra_compile_args_gcc or [])
    prefs.codegen.cpp.extra_compile_args_gcc = existing + ["-DDEBUG"]

    yield device

    try:
        device.reinit()
        device.activate()
    except Exception:
        pass
