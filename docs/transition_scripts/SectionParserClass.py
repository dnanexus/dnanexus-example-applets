from __future__ import print_function
from helpers import InvalidSection, get_bash_function_by_path, get_python_function_by_path
from global_helper_vars import LANG_COMMENT, CODE_SECTION_SEARCH, SUPPORTED_INTERPRETERS


class SectionParser(object):
    """Parse a src/code.* file for sections.
    TODO handle block comments. and block strings in python
    TODO Make matching less sensitive and catch correct block comment start and ends
    """
    func_finders = {
        'bash': get_bash_function_by_path,
        'python': get_python_function_by_path}

    def __init__(self, code_file_path, logger, keep_comments=True, interpreter=""):
        self.code = code_file_path
        self._active_section = None
        self._section_mapping = {}
        self.logger = logger
        self._func_finder = None
        self.applet_full_script = ""
        self.section_start, self.section_end = CODE_SECTION_SEARCH[interpreter]
        self.language = SUPPORTED_INTERPRETERS.get(interpreter, "")
        self.comment_chars = LANG_COMMENT.get(interpreter, [])
        self.keep_comments = keep_comments
        self._tempsection = []

    @classmethod
    def create_code_region(cls, code_str, language):
        code_str = code_str.strip('\n\r')
        content = code_str[:-len('dxpy.run()')].strip("\n\r") if code_str.endswith('dxpy.run()') else code_str  # final section edgecase
        code = "{start}{content}{end}".format(
            start="```{language}\n".format(language=language),
            content=content,
            end="\n```\n")
        return code

    @property
    def active_section(self):
        return self._active_section

    @active_section.setter
    def active_section(self, value):
        self.logger.info("Section change to: {0}".format(value))
        value = value.strip() if value is not None else None
        if value in self._section_mapping:
            raise InvalidSection("You Only Add Secions Once. Give sections unique names")
        if self._tempsection:
            self._section_mapping[self.active_section] = self.create_code_region(
                code_str="".join(self._tempsection), language=self.language)
        self._active_section = value
        self._tempsection = []

    @property
    def func_finder(self):
        if not self._func_finder:
            self._func_finder = self.func_finders.get(self.language)
        if self._func_finder is None:
            raise Exception('Language is not supported yet')
        return self._func_finder

    def get_parse_dict(self):
        """Maybe validate."""
        return self._section_mapping

    def get_func_code(self, func_name):
        """Import python source code and get function code through inspect"""
        if func_name == 'FULL SCRIPT':
            return self.applet_full_script
        return self.create_code_region(
            code_str=self.func_finder(self.code, func_name), language=self.language)

    def parse_file(self):
        """
        FIXME: Doc strings done on the same line not correctly parsed. Throws
               entire document parsing off.
        """
        def section_parse(line, flag):
            if self.active_section is None:
                match = self.section_start.match(line)
                if match is not None:
                    self.active_section = match.group(1)
            else:
                end = self.section_end.match(line)
                start_check = self.section_start.match(line)
                if start_check is not None:
                    self.active_section = start_check.group(1).strip()
                elif end is not None:
                    self.active_section = None
                else:
                    if not flag:
                        self._tempsection.append(line)

        self.logger.info("Code region search")
        comment_flag = False  # False,True, or comment flag. Block comment char if needed
        all_code = []
        with open(self.code, "r") as f:
            for line in f:
                comment_flag = self.is_comment(line, comment_flag)
                section_parse(line, comment_flag)
                if not comment_flag:  # move into c
                    all_code.append(line)
        if self.active_section is not None and self._tempsection:
            self.active_section = None
        self.logger.info("Creating code regions")
        self.applet_full_script = self.create_code_region(
            code_str="".join(all_code).strip('\n\r'), language=self.language)
        return self

    def is_comment(self, line, comment_flag=False):
        def startswith_prefixes(word, prefixes):
            pref_match = [word.startswith(prefix) for prefix in prefixes]
            if any(pref_match):
                return prefixes[pref_match.index(True)]
            return False
        if not self.keep_comments:
            return False

        line = line.strip()
        if startswith_prefixes(line, self.comment_chars.get("regular_comment", [])):
            self.logger.debug("regular comment, line: {line_is}".format(line_is=line))
            return comment_flag if comment_flag else True
        # block comment code block
        block_str = startswith_prefixes(line, self.comment_chars.get("block_comment", []))
        if type(block_str) is str:
            self.logger.debug("block comment, line: {0}".format(block_str))
            return True if block_str == comment_flag else block_str
        return comment_flag if type(comment_flag) is str else False
