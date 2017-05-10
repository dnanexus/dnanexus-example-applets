from __future__ import print_function
import argparse
import os
import json
import re
import logging
from multiprocessing import Pool
from itertools import izip
from helpers import FrontMatter, SectionParser, InvalidSection
from global_helper_vars import TUTORIAL_TYPES_SEARCH, SUPPORTED_INTERPRETERS, NUM_CORES


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
        msg = "SUCCESS: {title}" if status[i] else "FAIL: {title}"
        print(msg.format(title=tutorial_dicts[i]["title"]))


def _resolve_applet(dxapp_path):
    """Moved to top level due to Pool.map() inability to use nested func"""
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


def _write_kmarkdown(fh_md, readme_md_path, section_parser={}):
    """Creates Kramdown webpage"""
    def _section_match(match):
        section = match.group(1).strip()
        logging.debug("Kramdown code region {0}".format(section))
        return section_parser[section]

    def _force_line_match(match):
        return match.group(1).strip()
    special_match = {
        re.compile(r'.*<!--\s*SECTION:\s*\b(.*)-->'): _section_match,
        re.compile(r'.*<!--\s*INCLUDE:\s*(\S.*)-->'): _force_line_match,
    }
    with open(readme_md_path) as readme_md:
        for line in readme_md:
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
    parser.add_argument("-v", help="Verbosity", action="store_true", dest="verbose")

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
    def _write_front_matter(dxapp_obj, file_handle):
        """Generate and write liquid front matter"""
        frontmatter = FrontMatter(
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
        file_handle.write(frontmatter.__str__())
        logging.info("front matter written")

    target_file = os.path.join(
        tutorial_dict["site_pages_dir"], tutorial_dict["name"].strip() + ".md")
    if os.path.exists(target_file):
        if tutorial_dict["overwrite"]:
            logging.debug("removing: {0}".format(target_file))
            os.remove(target_file)
        else:
            logging.info("exist: {fn}".format(fn=target_file))
            return True

    # Handle JSON for assets at some point, should be easy.
    code_block_dict = {}
    if SUPPORTED_INTERPRETERS.get(tutorial_dict["interpreter"]) is not None:
        code_block_dict = SectionParser(
            code_file_path=tutorial_dict["src_code"],
            ignore_comments=tutorial_dict["ignore_comments"],
            interpreter=tutorial_dict["interpreter"]).parse_file().get_section_dict()
    try:
        with open(target_file, "w") as tutorial_md:
            _write_front_matter(tutorial_dict, tutorial_md)
            _write_kmarkdown(
                fh_md=tutorial_md,
                readme_md_path=tutorial_dict["readme_md"],
                section_parser=code_block_dict)
    except InvalidSection as IE:
        # Change to logging
        logging.info(IE.message)
        return False

    return True


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
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)

    _rehydrate_site(args)


if __name__ == '__main__':
    main()
