"""Brian2-callable wrappers around the test-helper accessors exposed at the
bottom of ``interface.cpp`` (inside the ``extern "C"`` test-helper block).

These mirror the ``@implementation('cpp', ...)`` pattern used by
``ni_brian2.step_clamp`` so every test exercises the same codegen path
Brian2 uses in production. The Python bodies are never executed at
runtime — they exist only so Brian2 can register the function names.

Important name-mangling note
----------------------------
Brian2's cpp_standalone codegen emits the *exact* function name that
appears in the equation string. The wrapper names below therefore
**must not** collide with the underlying C symbols (``get_steps_taken``,
``reset_state``, ...) or Brian2 will bypass the wrapper and call the
zero-argument C function with a ``t`` argument, producing
``too many arguments`` compile errors. We use a ``ni_*`` prefix to keep
the wrappers in their own namespace.

Wrapper definitions are ``static inline`` so Brian2 may safely emit
them into multiple translation units without multiple-definition
linker errors. Forward declarations of the underlying C symbols are
injected via the ``@implementation`` body so we do not extend the
public ``interface.h`` ABI.
"""
import os

from brian2 import implementation, check_units, second

_NI_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
_INTERFACE_CPP = os.path.join(_NI_DIR, "interface.cpp")
_LIBNIDAQMX = os.path.join(_NI_DIR, "libnidaqmx.so")

_IMPL_KW = dict(
    sources=[_INTERFACE_CPP, _LIBNIDAQMX],
    headers=['"interface.h"', '"NIDAQmx.h"'],
    include_dirs=[_NI_DIR],
)

# Forward decls for the symbols defined inside interface.cpp's
# test-helper extern "C" block.
_FWD_DECLS = '''
extern "C" {
    long get_steps_taken(void);
    double get_total_debt(void);
    double get_last_read_t(void);
    double get_total_rate(void);
    double get_run_time(void);
    void reset_state(void);
}
'''


@implementation('cpp', _FWD_DECLS + '''
static inline double ni_steps_taken(double dummy) {
    (void)dummy; return (double)get_steps_taken();
}
''', **_IMPL_KW)
@check_units(dummy=second, result=1)
def ni_steps_taken(dummy):
    raise NotImplementedError("Linked in from interface.cpp")


@implementation('cpp', _FWD_DECLS + '''
static inline double ni_total_debt(double dummy) {
    (void)dummy; return get_total_debt();
}
''', **_IMPL_KW)
@check_units(dummy=second, result=second)
def ni_total_debt(dummy):
    raise NotImplementedError("Linked in from interface.cpp")


@implementation('cpp', _FWD_DECLS + '''
static inline double ni_last_read_t(double dummy) {
    (void)dummy; return get_last_read_t();
}
''', **_IMPL_KW)
@check_units(dummy=second, result=second)
def ni_last_read_t(dummy):
    raise NotImplementedError("Linked in from interface.cpp")


@implementation('cpp', _FWD_DECLS + '''
static inline double ni_total_rate(double dummy) {
    (void)dummy; return get_total_rate();
}
''', **_IMPL_KW)
@check_units(dummy=second, result=second)
def ni_total_rate(dummy):
    raise NotImplementedError("Linked in from interface.cpp")


@implementation('cpp', _FWD_DECLS + '''
static inline double ni_run_time(double dummy) {
    (void)dummy; return get_run_time();
}
''', **_IMPL_KW)
@check_units(dummy=second, result=second)
def ni_run_time(dummy):
    raise NotImplementedError("Linked in from interface.cpp")


@implementation('cpp', _FWD_DECLS + '''
static inline double ni_reset_state(double dummy) {
    (void)dummy; reset_state(); return 0.0;
}
''', **_IMPL_KW)
@check_units(dummy=second, result=1)
def ni_reset_state(dummy):
    raise NotImplementedError("Linked in from interface.cpp")
