from multiprocessing import cpu_count
from collections import namedtuple
import re  # re-read docs and make re not greedy

NUM_CORES = cpu_count()

AppObj = namedtuple('AppObj', ['readme_md', 'src_code', 'app_name', 'title', 'interpreter', 'summary'])

TUTORIAL_TYPES_SEARCH = {
    "parallel": re.compile(r'.*para.*'),
    "distributed": re.compile(r'.*distr.*'),
    "web_app": re.compile(r'.*web.*')
}

SUPPORTED_INTERPRETERS = {  # dxapp.json supported interpreters. switch case setup
    "bash": "bash",
    "python2.7": "python"
}
LANG_COMMENT = {
    "bash": {"regular_comment": ["#"]},
    "python2.7": {"regular_comment": ["#"], "block_comment": ["\"\"\""]},
}

LANG_FUNC_SECTION_SEARCH = {  # multilple match groups
    "python2.7": re.compile(r'^(\s*)\bdef\s\b(.+)\(.+:'),
    "bash": re.compile(r'^\s*\b(.+)\(.+{|^\s*\bfunction\s(.+){')
}

CODE_SECTION_SEARCH = {  # section (STARTre, ENDre), similar for now
    "bash": (re.compile(r'.*#\s*SECTION:\s*\b(.+)$'), re.compile(r'.*#\s*SECTION-END')),
    "python2.7": (re.compile(r'.*#\s*SECTION:\s*\b(.+)$'), re.compile(r'.*#\s*SECTION-END')),
}
