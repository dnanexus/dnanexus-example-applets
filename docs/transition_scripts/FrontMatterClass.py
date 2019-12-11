from datetime import date
import yaml


class FrontMatter(object):
    """ Extend functionality as needed
    Front matter is the part of the markdown file at the very top that looks like:
    ---
    date: 2017-08-02
    title: Pysam tools count
    categories:
      - python # language
      - basic  # Tutorial type
    description: SAMtools count on an input BAM using Pysam, a python wrapper for SAMtools
    type: Document
    ---

    Attributes:
        allowed_inputs: set of allowed class inputs, extend as needed
        frontmatter_tutorial_types: set of allowed tutorial types, extend as needed
        logger: Explicitly passed logger TODO: remove this and log in parallel correctly

    """
    allowed_inputs = {'title', 'description', 'categories', 'date', 'github_link', 'summary'}

    def __init__(self, logger, isdocument, **kwargs):
        self.logger = logger
        self.logger.debug('Entering Front Matter class')
        invalid_inputs = set(kwargs.keys()) - self.allowed_inputs
        if invalid_inputs:
            msg = "Following inputs are not valid: {inps}".format(inps="\n".join(invalid_inputs))
            self.logger.info(msg)
            raise Exception(msg)
        self._frontmatter_mapping = {
            'date': date.today().isoformat(),
            'type': 'Document' if isdocument else 'Video'  # TODO: Video not fully supported yet
        }

    def _add_direct_to_front(self, key, value):
        self._frontmatter_mapping[key] = value

    def _add_category(self, *kwargs):
        category = kwargs[1]
        try:
            self._frontmatter_mapping['categories'].add(category)
        except KeyError:
            self._frontmatter_mapping['categories'] = set([category])

    def add_field(self, field, value):
        """Add a field to the front matter
        TODO: Improve implementation, being more pythonic helps here. Try... except
        """
        self.logger.info("Front matter field: {field} value {val}".format(field=field, val=value))
        if field not in self.allowed_inputs:
            raise Exception("{field} is not a defined frontmatter field. Do you need to add it?".format(field=field))
        if field in self._frontmatter_mapping and not self.validate_iter(self._frontmatter_mapping[field]):
            raise Exception("{field} already defined".format(field=field))
        special_func_mapping = {
            'categories': self._add_category
        }
        handler_func = special_func_mapping.get(field, self._add_direct_to_front)
        handler_func(field, value)

    @property
    def front_matter(self):
        return self._frontmatter_mapping

    @classmethod
    def validate_iter(cls, obj):
        try:
            obj_iter = iter(obj)  # just don't want print to terminal
        except TypeError:
            return False
        return True

    def __str__(self):
        try:  # working around pyaml set collection issues.
            self._frontmatter_mapping['categories'] = list(self._frontmatter_mapping['categories'])
        except KeyError:
            pass

        frontmatter = "---\n{yaml_dumps}---\n".format(
            yaml_dumps=yaml.safe_dump(
                self._frontmatter_mapping, default_flow_style=False))
        return frontmatter
