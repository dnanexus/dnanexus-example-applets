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

    def __init__(self, logger, **kwargs):
        self.logger = logger
        self.mapping = {}
        for k, v in kwargs.iteritems():
            self.add_field(k, v)

    def add_field(self, field, value):
        self.logger.info("Front matter field: {field} value {val}".format(field=field, val=value))
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
