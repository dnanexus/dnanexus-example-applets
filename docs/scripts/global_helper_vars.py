from multiprocessing import cpu_count
import re

NUM_CORES = cpu_count()

TUTORIAL_TYPES_SEARCH = {
    "parallel": re.compile(r'[para]'),
    "distributed": re.compile(r'[distr|distributed]')  # distr is good enough
}

SUPPORTED_INTERPRETERS = {  # dxapp.json supported interpreters. switch case setup
    "bash": "bash",
    "python2.7": "python"
}
"""
LANG_COMMENT = {
    "bash": ["#"],
    "python2.7": ["#", "\"\"\""],
}
"""
LANG_COMMENT = {
    "bash": {"regular_comment": ["#"]},
    "python2.7": {"regular_comment": ["#"], "block_comment": ["\"\"\""]},
}

CODE_SECTION_SEARCH = {  # section (STARTre, ENDre), similar for now
    "bash": (re.compile(r'.*# CODE-SECTION:\s*\b(.+)$'), re.compile(r'.*#\s*CODE-SECTION-END')),
    "python2.7": (re.compile(r'.*# CODE-SECTION:\s*\b(.+)$'), re.compile(r'.*#\s*CODE-SECTION-END')),
}
