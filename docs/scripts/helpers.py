"""Helper functions and classes
Frontmatter helper class.
SectionParser:src/code parser class
---
"""
from __future__ import print_function
from global_helper_vars import LANG_COMMENT, CODE_SECTION_SEARCH, SUPPORTED_INTERPRETERS
import logging


class InvalidSection(Exception):
    pass


class FrontMatter(object):
    """ Extend func as needed
    Front matter is the part of the markdown file that looks like:
    ---
    title: SAMtools distributed count by region
    tutorial_type: distributed
    source: samtools_count_dist_chr_region_sh
    language: bash
    ---
    """
    allowed_fields = ['title', 'tutorial_type', 'source', 'language']

    def __init__(self, **kwargs):
        logging.debug("init FM")
        self.mapping = {}
        for k, v in kwargs.iteritems():
            self.add_field(k, v)
        logging.debug("done FM")

    def add_field(self, field, value):
        logging.info("Front matter field: {field} value {val}".format(field=field, val=value))
        if field not in self.allowed_fields:
            raise Exception("{field} is not a defined frontmatter field. Do you need to add it?".format(field=field))
        if field in self.mapping:
            raise Exception("{field} already defined".format(field=field))
        self.mapping[field] = value

    def __str__(self):
        frontmatter = "---"
        for k, v in self.mapping.iteritems():
            frontmatter += "\n{key}: {value}".format(key=k, value=v)
        frontmatter += "\n---\n"
        return frontmatter


class SectionParser(object):
    """Parse a src/code.* file for sections.
    TODO handle block comments. and block strings in python
    """
    def __init__(self, code_file_path, ignore_comments=True, interpreter=""):
        self.code = code_file_path
        self._active_section = None
        self._section_mapping, self._func_mapping = {}, {}
        self.section_start, self.section_end = CODE_SECTION_SEARCH[interpreter]
        self.language = SUPPORTED_INTERPRETERS.get(interpreter, "")
        self.comment_chars = LANG_COMMENT.get(interpreter, [])
        self.ignore_comments = ignore_comments
        self._tempsection = []

    @property
    def active_section(self):
        return self._active_section

    @active_section.setter
    def active_section(self, value):
        logging.info("Section change to: {0}".format(value))
        value = value.strip() if value is not None else None
        if value in self._section_mapping:
            raise InvalidSection("You Only Add Secions Once. Give sections unique names")
        if self._tempsection:
            code = "{start}{content}{end}".format(
                start="```{language}\n".format(language=self.language),
                content="".join(self._tempsection).strip("\n\r"),
                end="\n```\n")
            self._section_mapping[self.active_section] = code
        self._active_section = value
        self._tempsection = []

    def get_parse_dict(self):
        """Maybe validate."""
        return self._section_mapping, self._func_mapping

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

        def func_parse(line, flag):
            pass

        logging.info("Code region search")
        comment_flag = False  # False,True, or comment flag. Block comment char if needed
        # func_flag = (0, 0)  # ( Int number of space chars, newline count )
        all_code = []
        with open(self.code, "r") as f:
            for line in f:
                comment_flag = self.is_comment(line, comment_flag)
                # func_flag = self.is_in_func(line, func_flag, comment_flag)
                section_parse(line, comment_flag)
                if not comment_flag:  # move into c
                    all_code.append(line)
        if self.active_section is not None and self._tempsection:
            self.active_section = None

        self._func_mapping["FULL SCRIPT"] = "{code_start}{code_blk}{code_end}".format(
            code_start="\n```{lang}\n".format(lang=self.language),
            code_blk="".join(all_code).strip(),
            code_end="\n```\n")
        logging.info("Code regions created")
        return self

    def is_comment(self, line, comment_flag=False):
        def startswith_prefixes(word, prefixes):
            pref_match = [word.startswith(prefix) for prefix in prefixes]
            if any(pref_match):
                return prefixes[pref_match.index(True)]
            return False
        if not self.ignore_comments:
            return False

        line = line.strip()
        if startswith_prefixes(line, self.comment_chars.get("regular_comment", [])):
            logging.debug("regular comment, line: {line_is}".format(line_is=line))
            return comment_flag if comment_flag else True
        # block comment code block
        block_str = startswith_prefixes(line, self.comment_chars.get("block_comment", []))
        if type(block_str) is str:
            logging.debug("block comment, line: {0}".format(block_str))
            return True if block_str == comment_flag else block_str
        return comment_flag if type(comment_flag) is str else False

        def is_in_func(self, line, func_flag, comment_flag):
            def _count_prefix_spaces(sentence):
                curr_spaces = 0
                for ch in sentence:
                    if ch == " ":
                        curr_spaces += 1
                    elif ch == "#":
                        return False
                    else:
                        break
                return curr_spaces
            func_space_count, newline_count = func_flag
            curr_spaces = _count_prefix_spaces(line)
        #  TODO logic to see when a file ends
