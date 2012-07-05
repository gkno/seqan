# encoding: utf-8
"""Trac text boxes macro.

This macro allows to create boxes for shell and IDE output.

Example:

{{{
{{{
#!ShellBox
# cd /tmp
# touch foo
}}}

{{{
#!IdeBox
Success!
}}}
}}}
"""

import StringIO

from trac.core import *
from trac.web.chrome import ITemplateProvider, add_stylesheet
from trac.wiki.macros import WikiMacroBase

import trac.wiki
import genshi.core
from genshi.builder import tag

CLASSES = {
    'ShellBox' : 'shell_box',
    'IdeBox' : 'ide_box',
    'WarningBox': 'system-message',  # Use trac class.
}

STYLES = {
    'ShellBox' : 'background-color: black; color: lightgray;',
    'IdeBox' : 'background-color: white; color: black; border:1px solid black;',
    'WarningBox' : None,
}

class TextBoxMacro(WikiMacroBase):
    implements(ITemplateProvider)
    
    def get_macros(self):
        yield 'ShellBox'
        yield 'IdeBox'
        yield 'WarningBox'

    def expand_macro(self, formatter, name, content, args):
        # TODO(holtgrew): Actually, we would like to add a style sheet but this does not work. Thus we use a style tag below.
        #add_stylesheet(formatter.req, 'text_boxes/css/text_boxes.css')
        className = CLASSES.get(name, 'ShellBox')
        if name == 'WarningBox':
            content_html = self.format_wiki(formatter, content)
            return tag.div(genshi.core.Markup(content_html), class_=className, style=STYLES.get(name))
        else:
            return tag.pre(content, class_='wiki ' + className, style=STYLES.get(name))

    def format_wiki(self, formatter, wiki_string):
        """Format the given string wiki_string to HTML."""
        out = StringIO.StringIO()
        trac.wiki.Formatter(self.env, formatter.context).format(wiki_string, out)
        return out.getvalue()

    ### ITemplateProvider methods

    def get_templates_dirs(self):
        return []

    def get_htdocs_dirs(self):
        from pkg_resources import resource_filename
        return [('text_boxes', resource_filename(__name__, 'htdocs'))]

