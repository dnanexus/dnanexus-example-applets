"""Helper functions and classes
Frontmatter helper class.
SectionParser:src/code parser class
---
"""
from __future__ import print_function
import os
import sys
from shutil import move
import re
from subprocess import check_output, CalledProcessError
from contextlib import contextmanager
from importlib import import_module
from inspect import getsource


IMPORT_ERROR_MATCHER = re.compile(r'No module named \b(.+)\b')
DXPY_RUN_SEARCH = re.compile(r'.*dxpy\.run\(\).*')
DXPY_DECORATOR = re.compile(r'.*@dxpy\.entry_point\(\S+\).*')


class InvalidSection(Exception):
    pass


@contextmanager
def pushd_popd(abs_dir_context):
    """Just mimic bash pushd and popd"""
    curr_dir = os.getcwd()
    try:
        os.chdir(abs_dir_context)
        yield
    finally:
        os.chdir(curr_dir)


@contextmanager
def _temp_applet_src_alter(dirpath, module_path):
    """Temporarily comment out dxpy.run() and alter python path so applet src can be imported

    FIXME seek to dxpy.run() line the second time
    Inspiration: https://stackoverflow.com/questions/17211078/how-to-temporarily-modify-sys-path-in-python
    """
    def file_comment_lines(re_searches, source_path, comment=True):
        new_path = "{path}_extra".format(path=source_path)
        with open(new_path, 'w') as nf:
            with open(source_path, 'r') as f:
                for line in f:
                    for re_search in re_searches:
                        if re_search.match(line):
                            line = "{pre}{ori}".format(pre="#", ori=line) if comment else line.lstrip("#")
                    nf.write(line)
        os.remove(source_path)
        move(new_path, source_path)

    sys.path.insert(0, dirpath)
    file_comment_lines(
        re_searches=[DXPY_RUN_SEARCH, DXPY_DECORATOR], source_path=module_path)
    yield
    file_comment_lines(
        re_searches=[DXPY_RUN_SEARCH, DXPY_DECORATOR], source_path=module_path, comment=False)
    print("Removing temp PATH")  # How does this work with multiple processes?
    sys.path.remove(dirpath)


def resolve_module(module_name, depth=0):
    """Attempt to import a module, mock failed imports as they occur.

    FIXME: del imports once module resolved
    FIXME: High risk function, find some other way of parsing source code
    FIXME: Add proper Mock module and use that instead
    FIXME: After a duh moment, I realize I could just temporarily comment out r'^import.*' lines... duh
    Module mocking inspiration: https://stackoverflow.com/questions/8658043/how-to-mock-an-import
    """
    try:
        module_import = import_module(name=module_name)
    except ImportError as ie:
        module_to_mock = IMPORT_ERROR_MATCHER.match(ie.message).group(1)
        if module_to_mock == module_name:
            raise ImportError('Circular dependency.')
        sys.modules[module_to_mock] = import_module('string')  # placeholder module for now, issues?
        depth += 1
        return resolve_module(module_name, depth=depth)
    else:
        print('Import Depth {} reached before full module resolution: {}'.format(depth, module_name))
        return module_import


def get_python_function_by_path(path_to_src, func_name):
    """https://stackoverflow.com/questions/6677424/how-do-i-import-variable-packages-in-python-like-using-variable-variables-i

    Raises:
        AttributeError - function doesn't exist
        OSError
    """
    py_path, module_py = os.path.split(path_to_src)
    module_name = module_py.rstrip('.py')
    with _temp_applet_src_alter(dirpath=py_path, module_path=path_to_src):
        with pushd_popd(py_path):
            module_import = resolve_module(module_name=module_name)
            func_obj = getattr(module_import, func_name)
            return getsource(func_obj)


def get_bash_function_by_path(path_to_src, func_name):
    """https://stackoverflow.com/questions/6916856/can-bash-show-a-functions-definition"""
    cmd = "source {src_code} &>/dev/null; declare -f {func}".format(
        src_code=path_to_src, func=func_name)
    try:
        return check_output(cmd, shell=True)
    except CalledProcessError:
        return
