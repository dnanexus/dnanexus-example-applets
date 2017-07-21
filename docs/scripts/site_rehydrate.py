#!/usr/bin/env python
from __future__ import print_function
import argparse
import os
import json
import re
import logging
from multiprocessing import Pool
from itertools import izip
from helpers import FrontMatter, SectionParser
from global_helper_vars import TUTORIAL_TYPES_SEARCH, SUPPORTED_INTERPRETERS, NUM_CORES


def _create_code_block_dict(tutorial_dict, logger):
    code_block_dict, func_dict = {}, {}
    if SUPPORTED_INTERPRETERS.get(tutorial_dict["interpreter"]) is not None:
        code_block_dict, func_dict = SectionParser(
            code_file_path=tutorial_dict["src_code"],
            ignore_comments=tutorial_dict["ignore_comments"],
            interpreter=tutorial_dict["interpreter"], logger=logger).parse_file().get_parse_dict()
    return code_block_dict, func_dict


def _rehydrate_site(user_args):
    def update_dict(d1, d2):
        d1.update(d2)
        return d1
    user_input_dict = {
        "overwrite": user_args.overwrite,
        "ignore_comments": user_args.skip_comments,
        "site_pages_dir": user_args.site_pages_dir
    }
    dxapp_files = find_all_matches(user_args.tutorials_dir, "dxapp.json")
    tutorial_dicts = [update_dict(d, user_input_dict) for d in _resolve_applets(dxapp_files)]
    md_maker = Pool(processes=NUM_CORES)
    status = md_maker.map(create_jekyll_markdown_tutorial, tutorial_dicts)
    md_maker.close()
    md_maker.join()
    for i in xrange(len(status)):
        msg = "SUCCESS: {title}" if status[i][0] else "FAIL: {title}"
        print(msg.format(title=tutorial_dicts[i]["title"]))
        if not status[i][0]:
            print(status[i][1])
            #logging.info(status[i][1])


def _resolve_applet(dxapp_path):
    """Moved to top level due to Pool.map() inability to use nested func or output dict"""
    with open(dxapp_path, "r") as dxf:
        dxapp_obj = json.load(dxf)
    src_code = os.path.join(dxapp_path[:-11], dxapp_obj["runSpec"]["file"])  # pretty sure this is a required field
    readme_md = os.path.join(dxapp_path[:-11], "Readme.md")
    if os.path.exists(src_code) and os.path.exists(readme_md):
        app_name = dxapp_obj["name"]
        title = dxapp_obj["title"]
        interpreter = dxapp_obj["runSpec"]["interpreter"]
        return (readme_md, src_code, app_name, title, interpreter)
    return None


def _resolve_applets(dxapp_paths):
    """
    TODO: just using multiprocessing.processing since mapping pickling issue makes this complex
    Create app tutorial representations to be used downstream
    dictionary keys returned: "readme_md", "src_code", "name", "title", "interpreter"
    """
    workers = Pool(NUM_CORES)
    tutorial_tuples = [res for res in workers.map(_resolve_applet, dxapp_paths) if res is not None]
    workers.close()
    workers.join()
    mapping_list = [{k: v for k, v in izip(("readme_md", "src_code", "name", "title", "interpreter"), tupe)} for tupe in tutorial_tuples]
    return mapping_list


def _write_kmarkdown(fh_md, readme_md_path, logger, section_parser={}, func_parser={}):
    """Creates Kramdown webpage

    Notes: Add {match: func(line)} to dictionary in order to handle special matches
    """
    def _section_match(match):
        section = match.group(1).strip()
        logger.debug("Kramdown code region in {0} {1}".format(readme_md_path, section))
        return section_parser[section]

    def _force_line_match(match):
        return "\n{insert}\n".format(insert=match.group(1).strip())

    def _func_match(match):
        func = match.group(1).strip()
        return func_parser[func]

    def _opt_header2_match(match):
        fh_md.write('<hr>')

    special_match = {  # Only one of these matches can be used. TODO flag multiple matches
        re.compile(r'.*<!--\s*SECTION:\s*\b(.*)-->'): _section_match,
        re.compile(r'.*<!--\s*INCLUDE:\s*(\S.*)-->'): _force_line_match,
        re.compile(r'.*<!--\s*FUNCTION:\s*(\S.*)-->'): _func_match
    }

    optional_match = {  # As many of these matches can be present in one line
        re.compile(r'^\#\#\s.*|.*<!--\s*INCLUDE:\s*\#\#\s.*-->'): _opt_header2_match
    }

    with open(readme_md_path) as readme_md:
        for line in readme_md:
            for matcher, opt_func in optional_match.iteritems():
                match = matcher.match(line)
                if match is not None:
                    opt_func(match)
            for matcher, func in special_match.iteritems():
                match = matcher.match(line)
                if match is not None:
                    line = func(match)
                    break
            fh_md.write(line)


def find_all_matches(tutdir, filename, exclude_dirs=[]):
    """Return list of """
    matches = []
    for root, dirs, files in os.walk(tutdir):
        if filename in files:
            matches.append(os.path.join(root, filename))
        for dir_ex in exclude_dirs:
            if dir_ex in dirs:
                dirs.remove(dir_ex)
    return matches


def get_parser():
    """Return an argument parser for the script."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--tutorials", help="Directory where tutorials are located", default=None, metavar="Tutorials directory", dest="tutorials_dir")
    parser.add_argument("--site-pages", help="Directory where site pages are located", default=None, metavar="Site Pages directory", dest="site_pages_dir")
    parser.add_argument("--keep-comments", help="Ignore other comments in code when parsing sections.", action="store_false", dest="skip_comments")
    parser.add_argument("--overwrite-files", help="Overwrites old files with generated files", action="store_true", dest="overwrite")

    return parser


def create_jekyll_markdown_tutorial(tutorial_dict):
    """
    :param dxapp_obj: tuple describing an applet.
        dict keys: "readme_md", "src_code", "name", "title"
                    "interpreter" "overwrite" "ignore_comments" "site_pages_dir"
    :type dxapp_obj: dict
    :overwrite: if target markdown already exist determins if it is overriden.
    :type overwrite: boolean
    TODO switch back to mapping once multiprocessing is implementended using Processing
    """
    log_file_unq_name = "{tut_name}_{title}".format(
        tut_name=tutorial_dict["name"].strip(), title=tutorial_dict["name"].strip())
    setup_logger(
        logger_name=log_file_unq_name, log_dir='log_temp_dir',
        level=logging.DEBUG)
    proc_logger = logging.getLogger(log_file_unq_name)

    target_file = os.path.join(
        tutorial_dict["site_pages_dir"], tutorial_dict["name"].strip() + ".md")
    if os.path.exists(target_file):
        if tutorial_dict["overwrite"]:
            proc_logger.debug("removing: {0}".format(target_file))
            os.remove(target_file)
        else:
            proc_logger.info("Exist: {fn}".format(fn=target_file))
            return True

    code_block_dict, func_dict = _create_code_block_dict(tutorial_dict, logger=proc_logger)

    try:
        with open(target_file, "w") as tutorial_md:
            _write_front_matter(dxapp_obj=tutorial_dict, logger=proc_logger, file_handle=tutorial_md)
            _write_kmarkdown(
                fh_md=tutorial_md,
                readme_md_path=tutorial_dict["readme_md"],
                section_parser=code_block_dict,
                func_parser=func_dict,
                logger=proc_logger)
    except Exception as e:
        return False, "Failed with Error:\n{err}\n Review log {logname}".format(
            err=e.message, logname=log_file_unq_name)

    return True, ""


def setup_logger(logger_name, log_dir, level=logging.INFO):
    """Create logger instance
    Inspiration: https://stackoverflow.com/questions/17035077/python-logging-to-multiple-log-files-from-different-classes
    """
    try:
        os.mkdir(log_dir)
    except OSError:
        pass
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)
    log_file = os.path.join(log_dir, logger_name)
    log_instance = logging.getLogger(logger_name)
    formatter = logging.Formatter('%(asctime)s : %(message)s')
    fileHandler = logging.FileHandler(log_file, mode='w')
    fileHandler.setFormatter(formatter)
    # streamHandler = logging.StreamHandler()
    # streamHandler.setFormatter(formatter)

    log_instance.setLevel(level)
    log_instance.addHandler(fileHandler)
    # log_instance.addHandler(streamHandler)


def _write_front_matter(dxapp_obj, logger, file_handle=None):
    """Generate and write liquid front matter"""
    proc_logger = logger
    frontmatter = FrontMatter(
        logger=logger,
        title=dxapp_obj["title"],
        source=dxapp_obj["name"])
    tut_type = "basic"
    for k, v in TUTORIAL_TYPES_SEARCH.iteritems():
        if v.match(dxapp_obj["name"]):
            tut_type = k
            break
    frontmatter.add_field("tutorial_type", tut_type)
    lang = SUPPORTED_INTERPRETERS.get(dxapp_obj["interpreter"], "none")
    frontmatter.add_field("language", lang)
    if file_handle is None:
        return frontmatter.__str__()
    file_handle.write(frontmatter.__str__())
    proc_logger.info("front matter written")


def main():
    """Script Entry point."""
    parser = get_parser()
    args = parser.parse_args()

    if args.tutorials_dir is None:
        script_dir = os.path.dirname(os.path.realpath(__file__))
        assumed_tutorial_path = os.path.join(script_dir, "..", "..", "Tutorials")
        args.tutorials_dir = os.path.normpath(assumed_tutorial_path)
    if args.site_pages_dir is None:
        script_dir = os.path.dirname(os.path.realpath(__file__))
        assumed_site_page_path = os.path.join(script_dir, "..", "_tutorials")
        args.site_pages_dir = os.path.normpath(assumed_site_page_path)

    _rehydrate_site(args)


if __name__ == '__main__':
    main()
