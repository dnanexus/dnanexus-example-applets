"""Helper functions and classes
Frontmatter helper class.
SectionParser:src/code parser class
---
"""
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
        self._section_mapping = {}
        self._tempsection = []
        self.section_start, self.section_end = CODE_SECTION_SEARCH[interpreter]
        self.language = SUPPORTED_INTERPRETERS.get(interpreter, "")
        self.comment_chars = LANG_COMMENT.get(interpreter, [])
        self.ignore_comments = ignore_comments

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
                content="".join(self._tempsection).strip(),
                end="\n```\n")
            self._section_mapping[self.active_section] = code
        self._active_section = value
        self._tempsection = []

    def get_section_dict(self):
        """Maybe validate."""
        return self._section_mapping

    def parse_file(self):
        logging.info("Code region search")
        comment_flag = False  # False,True, or comment flag. Block comment char if needed
        with open(self.code, "r") as f:
            for line in f:
                if self.active_section is None:
                    match = self.section_start.match(line)
                    if match is None:
                        continue
                    self.active_section = match.group(1)
                else:
                    end = self.section_end.match(line)
                    start_check = self.section_start.match(line)
                    if start_check is not None:
                        self.active_section = start_check.group(1).strip()
                        continue
                    elif end is not None:
                        self.active_section = None
                        continue
                    comment_flag = self.is_comment(line, comment_flag)
                    if end is None and not comment_flag:
                        self._tempsection.append(line)
            if self.active_section is not None and self._tempsection:
                self.active_section = None
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
        block_str = startswith_prefixes(line.strip(), self.comment_chars.get("block_comment", []))
        if type(block_str) is str:
            logging.debug("block comment, line: {0}".format(block_str))
            return True if block_str == comment_flag else comment_flag
        return False
