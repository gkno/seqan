#!/usr/bin/env python
"""Processed version of the documentation.

The documentation from the objects of raw_doc is further processed into
objects from the module proc_doc.  These objects can then be processed
into structured documents such as HTML more easily.
"""

# TODO(holtgrew): Location traceability for entries and text.

import os.path
import logging
import re
import sys
import xml.sax.saxutils

import inc_mgr
import sig_parser
import dox_tokens
import raw_doc


def escapeForXml(s):
    """Return escaped XML of s."""
    return xml.sax.saxutils.escape(s)


class DocumentationBuildException(Exception):
    """Thrown when there is a logical error on building the documentation."""


class LinkResolver(object):
    """Member of ProcDoc for resolving links."""
    
    def __init__(self, proc_doc):
        self.proc_doc = proc_doc

    def resolveEntry(self, link_str):
        """Resolve a name to an ProcEntry objects."""
        return self.proc_doc.entries[link_str]


def splitSecondLevelEntry(name):
    """Split second-level entry and return (first, second) pair.

    If name is not a second level entry then (None, name) is returned.
    """
    xs = None
    if name.count('::') > 1 and ' ' in name:
        xs = name.split(' ', 1)
    elif '#' in name:
        xs = name.split('#', 1)
    elif '::' in name:
        xs = name.rsplit('::', 1)
    if xs:
        return xs
    return (None, name)


class ProcDoc(object):
    """Collection of the top-level documentation entries."""

    def __init__(self):
        self.top_level_entries = {}
        self.second_level_entries = {}
        self.entries = {}
        self.link_resolver = LinkResolver(self)

    def addTopLevelEntry(self, x):
        """Add a top-level-entry."""
        self.registerEntry(x)
        self.top_level_entries[x.name] = x

    def addSecondLevelEntry(self, x):
        """Add a second-level entry."""
        self.registerEntry(x)
        self.second_level_entries[x.name] = x
        first, second = splitSecondLevelEntry(x.name)
        if first:
#            print '%s => %s as %s' % (x.name, second, x.kind)
            self.top_level_entries[first].registerSubentry(x)
        
    def addVariable(self, x):
        """Add a second-level entry."""
        self.registerEntry(x)
        #print 'x ==', x
        #print 'x.type ==', x.type
        if self.top_level_entries.get(x.type):
            self.second_level_entries[x.name] = x
            self.top_level_entries[x.type].registerSubentry(x)
        elif '::' in x.name:
            self.second_level_entries[x.name] = x
            first, second = splitSecondLevelEntry(x.name)
            self.top_level_entries[first].registerSubentry(x)
        else:
            self.top_level_entries[x.name] = x
        
    def registerEntry(self, x):
        """Register an entry."""
        if x.name in self.entries:
            raise DocumentationBuildException('%s already known! (%s)' % (x.name, x))
        self.entries[x.name] = x
        x.doc = self

    def runTextVisitor(self, v):
        """Run visitor v on all Text members of all entries and sub entries.
        """
        for e in self.entries.itervalues():
            e.runTextVisitor(v)


class TextNode(object):
    """A node represents a part of a processed text.

    Processed text is text generated from tokens lexed from the input file.
    For example, text in the paragraph of a entry's body can be representd by
    TextNode objects.

    TextNode objects are similar to DOM nodes, i.e. they can contain children
    and have attributes.  This means that we can have a link node that has a
    href/target attribute with a target URL and one or more child nodes that
    contain the link's label.

    Additionally, we store the source location (begin and end line/column) of
    the node in its source file.

    We represent plain links, i.e. where the label is the same as the target
    using the representation for "<a href="seqan:$target">$target</a>".

    We represent included code snippets as "<code type='.cpp'>$code</code>."

    @ivar type: The type of the node, as a string.  Reserved values are
                '<text>' for plain text nodes.
    @ivar attrs: A dict object mapping attribute names to string values.
    @ivar children: A list of TextNode objects.
    @ivar text: The text value of a node, a string.
    """

    def __init__(self, type='<text>', verbatim=False, text='', attrs={}):
        self.type = type
        self.attrs = dict(attrs)
        self.children = []
        if verbatim:
            self.text = text
        else:
            self.text = escapeForXml(text)

    def __str__(self):
        attrs = (repr(self.type), repr(self.text), repr(self.attrs), len(self.children))
        return 'TextNode(type=%s, text=%s, attrs=%s, len(children)=%d)' % attrs

    def __repr__(self):
        return str(self)
        
    def setAttr(self, key, value):
        self.attrs[escapeForXml(key)] = escapeForXml(value)

    def addChild(self, n):
        self.children.append(n)
        return self.children[-1]

    @property
    def X(self):
        """Returns first child, used to retrieve member of top-level <div>."""
        if self.type == '<text>':
            return self
        else:
            return self.children[0]
        
    def toHtmlLike(self, skip_top_tag=False, **kwargs):
        """Returns a string with a HTML-like representation for debuggin.

        @param skip_top_tag: Do ont output top-level tag.
        @param kwargs: Additional attributes to add.
        """
        if self.type == '<text>':
            if self.attrs:
                print >>sys.stderr, 'WARNING: Attributes on text node!'
            return self.text
        else:
            res = []
            if not skip_top_tag:
                res += ['<', self.type]
                for key, value in self.attrs.iteritems():
                    res += [' ', key, '=', '"', repr(value)[1:-1], '"']
                for key, value in kwargs.iteritems():
                    res += [' ', key, '=', '"', value, '"']
                res.append('>')
            res += [x.toHtmlLike() for x in self.children]
            if not skip_top_tag:
                res += ['</', self.type, '>']
            return ''.join(res)


class ProcEntry(object):
    """A processed representation of a documentation entry.

    A documentation entry has a kind (string), a name (string), a brief
    description (TextNode(type='<text>')), and a list of references/sees to
    other elements (list of TextNode(type='<link>')).  Also, it has a body
    which is a TextNode with children.

    @ivar kind: The kind of the entry, string.
    @ivar name: The name of the entry, string.
    @ivar brief: A brief description, a text-typed TextNode or None.
    @ivar body: A TextNode object with children for the documentation body.
    @ivar sees: A list of link-typed TextNode objects, can be empty.
    @ivar doc:  The owning, ProcDoc, set on ProcDoc.registerEntry
    @ivar subentries: Sub entries, dir, grouped by type.
    @ivar raw_entry: The RawEntry object that this ProcEntry was generated from.
    """

    def __init__(self, name, brief=None, body=None, sees=[]):
        self.name = name
        self.brief = brief
        self.body = body
        self.sees = list(sees)
        self.doc = None
        self.subentries = {}
        self.raw_entry = None

    def registerSubentry(self, proc_entry):
        self.subentries.setdefault(proc_entry.kind, []).append(proc_entry)

    def visitTextNodes(self, visitor):
        """Visit all text nodes using the given visitor."""
        visitor.visit(self.brief)
        visitor.visit(self.body)
        for see in self.sees:
            visitor.visit(see)

    @property
    def kind(self):
        return self.__class__.__name__.replace('Proc', '').lower()


class ProcCodeEntry(ProcEntry):
    """A processed code entry.

    @ivar signatures: A TextNode with the signatures of the entry.  They are
                      properly formatted to be displayed as verbatim text.
    @ivar signature_entries: A list of sig_parser.SigEntry objects.
    @ivar headerfiles: A list of str objects with the arguments to #include.
    @ivar deprecation_msgs: List of TextNode objects with deprecation messages.
    """

    def __init__(self, name, brief=None, body=None, sees=[]):
        ProcEntry.__init__(self, name, brief, body, sees)
        self.signatures = []
        self.signature_entries = []
        self.headerfiles = []
        self.deprecation_msgs = []

    def addSignature(self, s):
        self.signatures.append(s)

    def addSignatureEntry(self, e):
        self.signature_entries.append(e)

    def addHeaderfile(self, h):
        self.headerfiles.append(h)
        
    def addDeprecationMsg(self, m):
        self.deprecation_msgs.append(m)
        
    def subEntries(self, kind):
        return []


class ProcEnum(ProcCodeEntry):
    """A processed enum documentation.

    @ivar values: A list of ProcVariable entries that represent values
                  of this enum.
    """

    def __init__(self, name, brief=None, body=None, sees=[]):
        ProcCodeEntry.__init__(self, name, brief, body, sees)
        self.values = []


class ProcAdaption(ProcCodeEntry):
    """A processed adaption documentation.

    @ivar values: A list of ProcVariable entries that represent values
                  of this adaption.
    """

    def __init__(self, name, brief=None, body=None, sees=[]):
        ProcCodeEntry.__init__(self, name, brief, body, sees)
        self.values = []


class ProcTypedef(ProcCodeEntry):
    """A processed typedef documentation.

    @ivar values: A list of ProcVariable entries that represent values
                  of this typedef.
    """

    def __init__(self, name, brief=None, body=None, sees=[]):
        ProcCodeEntry.__init__(self, name, brief, body, sees)
        self.values = []

    @property
    def kind(self):
        if '#' in self.name:
            return 'grouped_typedef'
        elif '::' in self.name:
            return 'member_typedef'
        else:
            return 'global_typedef'


class ProcConcept(ProcCodeEntry):
    """A processed concept documentation.

    @ivar extends: A list of str values with the names of the extended
                   concepts.
    @ivar all_extended: A set of str values with the names of all extended
                        concepts, also transitively.
    @ivar all_extending: A set of str values with the names of all extending
                         concepts.
    @ivar all_implementing: A set of str values with the names of all
                            implementing classes.
    """

    def __init__(self, name, brief=None, body=None, sees=[]):
        ProcCodeEntry.__init__(self, name, brief, body, sees)
        self.extends = []
        self.all_extended = set()
        self.all_extending = set()
        self.all_implementing = set()

    def addExtends(self, s):
        self.extends.append(s)


class ProcClass(ProcCodeEntry):
    """A processed class documentation.

    @ivar extends: A list of str values with the names of the extended
                   classes.
    @ivar implements: A list of str values with the names of the implemented
                      concepts.
    @ivar all_implemented: Set of str values with the names of all implemented
                           concepts.
    @ivar all_extending: Set of str values with the names of all extending
                         classes.
    @ivar all_extended: Set of str values with the names of all extended classes.
    """

    def __init__(self, name, brief=None, body=None, sees=[]):
        ProcCodeEntry.__init__(self, name, brief, body, sees)
        self.extends = []
        self.implements = []
        self.all_implemented = set()
        self.all_extending = set()
        self.all_extended = set()
        self.tparams = []
        self.typedefs = []

    def addExtends(self, s):
        self.extends.append(s)

    def addImplements(self, s):
        self.implements.append(s)

    def addTParam(self, t):
        self.tparams.append(t)

    def addTypedef(self, t):
        self.typedefs.append(t)


class ProcTag(ProcCodeEntry):
    """A processed tag documentation.
    """

    @property
    def kind(self):
        if '#' in self.name:
            return 'grouped_tag'
        else:
            return 'tag'


class ProcParam(object):
    """Representation of a parameter.

    @ivar name: The name of the parameter. str.
    @ivar in_out: One of IN, OUT, IN_OUT, None.
    @ivar desc: Documentation of the parameter. TextNode.
    """

    def __init__(self):
        self.name = None
        self.in_out = None
        self.desc = TextNode()

    def visitTextNodes(self, visitor):
        """Visit all text nodes using the given visitor."""
        visitor.visit(self.desc)


ProcParam.IN = 'IN'
ProcParam.OUT = 'IN'
ProcParam.IN_OUT = 'IN_OUT'


class ProcTParam(object):
    """Documentation of a processed template parameter.

    @ivar type: The type of the parameter. str
    @ivar desc: Documentation of the parameter. TextNode.
    """

    def __init__(self):
        self.type = None
        self.desc = TextNode()

    def visitTextNode(self, visitor):
        """Visit all text nodes using the given visitor."""
        visitor.visit(self.desc)


class ProcReturn(object):
    """Documentation of a @return entry.

    @ivar type: The return type. str.
    @ivar desc: The documentation of the return value. TextNode.
    """

    def __init__(self):
        self.type = None
        self.desc = TextNode()

    def visitTextNode(self, visitor):
        """Visit all text nodes using the given visitor."""
        visitor.visit(self.desc)


class ProcFunction(ProcCodeEntry):
    """A processed function documentation.

    @ivar params: A list of str values with the names of the extended
                  concepts.
    @ivar tparams:
    @ivar returns:
    """

    def __init__(self, name, brief=None, body=None, sees=[]):
        ProcCodeEntry.__init__(self, name, brief, body, sees)
        self.params = []
        self.tparams = []
        self.returns = []

    @property
    def kind(self):
        if '#' in self.name:
            return 'interface_function'
        elif '::' in self.name:
            return 'member_function'
        else:
            return 'global_function'

    def visitTextNode(self, visitor):
        """Visit all text nodes using the given visitor."""
        ProcCodeEntry.visitTextNode(self, visitor)
        for p in self.params:
            p.visitTextNode(p)
        for p in self.tparams:
            p.visitTextNode(p)
        for p in self.returns:
            p.visitTextNode(p)
        
    def addParam(self, p):
        self.params.append(p)

    def addTParam(self, t):
        self.tparams.append(t)

    def addReturn(self, r):
        self.returns.append(r)


class ProcMacro(ProcCodeEntry):
    """A processed macro documentation.

    @ivar params: A list of str values with the names of the extended
                  concepts.
    @ivar returns:
    """

    def __init__(self, name, brief=None, body=None, sees=[]):
        ProcCodeEntry.__init__(self, name, brief, body, sees)
        self.params = []
        self.returns = []

    @property
    def kind(self):
        if '#' in self.name:
            return 'grouped_macro'
        else:
            return 'macro'

    def visitTextNode(self, visitor):
        """Visit all text nodes using the given visitor."""
        ProcCodeEntry.visitTextNode(self, visitor)
        for p in self.params:
            p.visitTextNode(p)
        for p in self.returns:
            p.visitTextNode(p)
        
    def addParam(self, p):
        self.params.append(p)

    def addReturn(self, r):
        self.returns.append(r)


class ProcMetafunction(ProcCodeEntry):
    """A processed function documentation.

    @ivar tparams: A list of str values with the names of the extended
                   concepts.
    @ivar returns: A list of ProcReturn values.
    """

    def __init__(self, name, brief=None, body=None, sees=[]):
        ProcCodeEntry.__init__(self, name, brief, body, sees)
        self.tparams = []
        self.returns = []

    @property
    def kind(self):
        if '#' in self.name:
            return 'interface_metafunction'
        else:
            return 'global_metafunction'

    def visitTextNode(self, visitor):
        """Visit all text nodes using the given visitor."""
        ProcCodeEntry.visitTextNode(self, visitor)
        for p in self.tparams:
            p.visitTextNode(p)
        for p in self.returns:
            p.visitTextNode(p)
        
    def addTParam(self, t):
        self.tparams.append(t)

    def addReturn(self, r):
        self.returns.append(r)


class ProcVariable(ProcCodeEntry):
    """A processed function documentation.

    @ivar type: A string with the name of a type.
    """

    def __init__(self, name, brief=None, body=None, sees=[]):
        ProcCodeEntry.__init__(self, name, brief, body, sees)
        self.type = None

    @property
    def kind(self):
        if '::' in self.name:
            return 'member_variable'
        else:
            return 'variable'


class ProcPage(ProcEntry):
    """A processed page."""

    def __init__(self, name, brief=None, body=None, sees=[]):
        ProcEntry.__init__(self, name, brief, body, sees)
        self.title = None

    def __str__(self):
        return 'Page(name=%s)' % repr(self.name)


class ProcGroup(ProcEntry):
    """A processed group."""

    def __init__(self, name, brief=None, body=None, sees=[]):
        ProcEntry.__init__(self, name, brief, body, sees)
        self.title = None
        self.tags = []
        self.typedefs = []

    def __str__(self):
        return 'Group(name=%s)' % repr(self.name)

    def addTypedef(self, t):
        self.typedefs.append(t)



class RawTextToTextNodeConverter(object):
    def __init__(self, strip_lt_line_space=False):
        self.tag_stack = []
        self.node_stack = []
        self.current = None
        self.strip_lt_line_space = strip_lt_line_space

    def handleTag(self, token):
        regex = r'<(/)?\s*([_a-zA-Z][_a-zA-Z0-9]*)(?:\s+([_a-zA-Z][_a-zA-Z0-9]*)="([^"]*)")*\s*(?:/)?>'
        m = re.match(regex, token.val)
        if m.group(1):  # </...
            #print >>sys.stderr, 'CLOSING %s' % token.val
            #print >>sys.stderr, '  B tag stack:  %s' % self.tag_stack
            #print >>sys.stderr, '  B node stack: %s' % self.node_stack
            if self.tag_stack and self.tag_stack[-1] == m.group(2):
                self.tag_stack.pop()
            elif self.tag_stack and self.tag_stack[-1] != m.group(2):
                args = (m.group(1), self.tag_stack[-1])
                print >>sys.stderr, 'WARNING: Closing wrong tag %s instead of %s' % args
                self.tag_stack.pop()
                return
            else:  # not self.tag_stack
                print >>sys.stderr, 'WARNING: Closing not opened tag %s' % m.group(1)
            if self.node_stack:
                self.current = self.node_stack[-1]
                self.node_stack.pop()
            else:
                print >>sys.stderr, 'WARNING: Having closed too many tags!'
            #print >>sys.stderr, '  A tag stack:  %s' % self.tag_stack
            #print >>sys.stderr, '  A node stack: %s' % self.node_stack
        else: # <$name
            #print >>sys.stderr, 'OPENING %s' % token.val
            #print >>sys.stderr, '  B tag stack:  %s' % self.tag_stack
            #print >>sys.stderr, '  B node stack: %s' % self.node_stack
            self.tag_stack.append(m.group(2))
            self.node_stack.append(self.current)
            tag = TextNode(type=m.group(2))
            if m.group(2) == 'a':
                tag.setAttr(m.group(3), m.group(4))
            self.current = self.current.addChild(tag)
            #print >>sys.stderr, '  A tag stack:  %s' % self.tag_stack
           # print >>sys.stderr, '  A node stack: %s' % self.node_stack

    def run(self, raw_text, verbatim=False):
        #print >>sys.stderr, '================== %s' % raw_text.text
        #print >>sys.stderr, [(t.type, t.val) for t in raw_text.tokens]
        self.current = TextNode(type='div')
        root = self.current
        at_line_start = True
        for i, t in enumerate(raw_text.tokens):
            if t.type in dox_tokens.WHITESPACE:
                if i == 0 or (i + 1) == len(raw_text.tokens):
                    continue  # Ignore leading and trailing whitespace.
                if t.type == 'SPACE' and at_line_start:
                    continue  # Ignore space at the beginning of a line.
                if t.type == 'BREAK':
                    self.current.addChild(TextNode(text='\n'))
                else:
                    self.current.addChild(TextNode(text=' '))
            elif not verbatim and t.type == 'HTML_TAG':
                at_line_start = False
                self.handleTag(t)
            else:
                at_line_start = False
                # TODO(holtgrew): Escape values.
                self.current.addChild(TextNode(text=t.val))
            at_line_start = t.type in ['EMPTY_LINE', 'BREAK']
        return root

    def process(self, raw_entry):
        raise Exception('Not implemented!')


class EntryConverter(object):
    """Base class for the conversion of raw entries processed entries.

    @ivar doc_proc: DocProcessor object.
    @ivar entry_class: The class of the ProcEntry type to create.
    """

    def __init__(self, doc_proc):
        self.doc_proc = doc_proc
        self.entry_class = None
    
    def rawTextToTextNode(self, raw_text, strip_lt_line_space=False, verbatim=False):
        """Convert RawText object into a TextNode object.

        The text node will have the type 'div'.

        @param strip_lt_breaks_lines: Whether or not to remove leading
                                      space for lines.
        @param verbatim: Whether or not to convert HTML tags.
        """
        converter = RawTextToTextNodeConverter(strip_lt_line_space)
        return converter.run(raw_text, verbatim)

    def bodyToTextNode(self, raw_body):
        """Convert a RawBody to a TextNode."""
        res = TextNode(type='div')
        for p in raw_body.paragraphs:
            if p.getType() == 'paragraph':
                if not p.text.text.strip():
                    continue  # Skip whitespace
                p = self.rawTextToTextNode(p.text)
                p.type = 'p'
                res.addChild(p)
            elif p.getType() == 'section':
                h = self.rawTextToTextNode(p.heading)
                h.type = 'h%d' % p.level
                res.addChild(h)
            elif p.getType() == 'include':
                ftype = os.path.splitext(p.path.text)[1]
                code_text = self.doc_proc.include_mgr.loadFile(p.path.text)
                proc_include = TextNode(type='code', attrs={'type': ftype})
                proc_include.addChild(TextNode(text=code_text, verbatim=True))
                res.addChild(proc_include)
            elif p.getType() == 'snippet':
                ftype = os.path.splitext(p.path.text)[1]
                code_text = self.doc_proc.include_mgr.loadSnippet(p.path.text, p.name.text)
                proc_snippet = TextNode(type='code', attrs={'type': ftype})
                proc_snippet.addChild(TextNode(text=code_text, verbatim=True))
                res.addChild(proc_snippet)
            elif p.getType() == 'code':
                code_text = p.text.text
                type = '.txt'
                m = re.match(r'^{[^}]+}', code_text)
                if m:
                    type = m.group(0)[1:-1]
                code_text = code_text[len(type) + 2:].strip()
                #print [repr(t.val) for t in p.text.tokens]
                x = TextNode(type='code', attrs={'type':type})
                x.addChild(TextNode(text=code_text, verbatim=True))
                res.addChild(x)
        return res

    def process(self, raw_entry):
        entry = self.entry_class(name=raw_entry.name.text)
        # Convert first brief member.  We already warned about duplicate ones
        # elsewhere.
        if raw_entry.briefs:
            entry.brief = self.rawTextToTextNode(raw_entry.briefs[0].text)
        # Convert the body
        if raw_entry.body:
            entry.body = self.bodyToTextNode(raw_entry.body)
        # Convert the sees entries.
        for see in raw_entry.sees:
            link = self.rawTextToTextNode(see.text)
            link.type = 'a'
            link.attrs['href'] = 'seqan:%s' % see.text.text
            entry.sees.append(link)
        # Store the raw entry in the processed ones.
        entry.raw_entry = raw_entry
        return entry

    
class CodeEntryConverter(EntryConverter):
    """Base for the processing RawCodeEntry objects into processed entries."""

    def __init__(self, doc_proc):
        EntryConverter.__init__(self, doc_proc)
        self.parse_signature = True
    
    def process(self, raw_entry):
        entry = EntryConverter.process(self, raw_entry)
        # Add headerfile paths as list of strings.
        for s in raw_entry.headerfiles:
            entry.addHeaderfile(s.text.text.strip())
        # Add deprecation messages as list of TextNodes.
        for s in raw_entry.deprecation_msgs:
            entry.addDeprecationMsg(self.rawTextToTextNode(s.text, strip_lt_line_space=True,
                                                           verbatim=True))
        # Add signatures as a text node with code.
        for s in raw_entry.signatures:
            entry.addSignature(self.rawTextToTextNode(s.text, strip_lt_line_space=True,
                                                      verbatim=True))
        # Use sig_parser to convert the signature texts to SigEntry objects.
        # They are used for the list of functions/metafunctions for a type.
        if self.parse_signature:
            for s in raw_entry.signatures:
                try:
                    sig_entry = sig_parser.SigParser(s.text.text).parse()
                    entry.addSignatureEntry(sig_entry)
                except sig_parser.SigParseException, e:
                    print >>sys.stderr, '\nWARNING: Could not parse signature: %s' % e
                    print >>sys.stderr, 'Signature is: %s' % s.text.text.strip()
                
        return entry


class EnumConverter(CodeEntryConverter):
    def __init__(self, doc_proc):
        CodeEntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcEnum
    
    def process(self, raw_entry):
        return CodeEntryConverter.process(self, raw_entry)


class AdaptionConverter(CodeEntryConverter):
    def __init__(self, doc_proc):
        CodeEntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcAdaption
        self.parse_signature = False
    
    def process(self, raw_entry):
        return CodeEntryConverter.process(self, raw_entry)


class TypedefConverter(CodeEntryConverter):
    def __init__(self, doc_proc):
        CodeEntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcTypedef
        self.parse_signature = False
    
    def process(self, raw_entry):
        return CodeEntryConverter.process(self, raw_entry)


class ConceptConverter(CodeEntryConverter):
    def __init__(self, doc_proc):
        CodeEntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcConcept
    
    def process(self, raw_entry):
        concept = CodeEntryConverter.process(self, raw_entry)
        for e in raw_entry.extends:
            concept.addExtends(e.text.text.strip())
        return concept


class ClassConverter(CodeEntryConverter):
    def __init__(self, doc_proc):
        CodeEntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcClass
    
    def process(self, raw_entry):
        klass = CodeEntryConverter.process(self, raw_entry)
        for e in raw_entry.extends:
            klass.addExtends(e.text.text.strip())
        for e in raw_entry.implements:
            klass.addImplements(e.text.text.strip())
        for t in raw_entry.tparams:
            proc_tparam = ProcTParam()
            proc_tparam.type = t.name.text
            proc_tparam.desc = self.rawTextToTextNode(t.text)
            klass.addTParam(proc_tparam)
        return klass


class TagConverter(CodeEntryConverter):
    def __init__(self, doc_proc):
        CodeEntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcTag
        self.parse_signature = False


class FunctionConverter(CodeEntryConverter):
    def __init__(self, doc_proc):
        CodeEntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcFunction
        self.in_out_map = {
            'in': ProcParam.IN,
            'out': ProcParam.OUT,
            'in,out': ProcParam.IN_OUT,
            }
    
    def process(self, raw_entry):
        function = CodeEntryConverter.process(self, raw_entry)
        for p in raw_entry.params:
            proc_param = ProcParam()
            proc_param.name = p.name.text
            if p.inout:
                proc_param.in_out = self.in_out_map.get(p.inout.val[1:-1])
            proc_param.desc = self.rawTextToTextNode(p.text)
            function.addParam(proc_param)
        for t in raw_entry.tparams:
            proc_tparam = ProcTParam()
            proc_tparam.type = t.name.text
            proc_tparam.desc = self.rawTextToTextNode(t.text)
            function.addTParam(proc_tparam)
        for r in raw_entry.returns:
            proc_return = ProcReturn()
            proc_return.type = r.name.text
            proc_return.desc = self.rawTextToTextNode(r.text)
            function.addReturn(proc_return)
        return function


class MacroConverter(CodeEntryConverter):
    def __init__(self, doc_proc):
        CodeEntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcMacro
        self.in_out_map = {
            'in': ProcParam.IN,
            'out': ProcParam.OUT,
            'in,out': ProcParam.IN_OUT,
            }
        self.parse_signature = False
        
    def process(self, raw_entry):
        macro = CodeEntryConverter.process(self, raw_entry)
        for p in raw_entry.params:
            proc_param = ProcParam()
            proc_param.name = p.name.text
            if p.inout:
                proc_param.in_out = self.in_out_map.get(p.inout.val[1:-1])
            proc_param.desc = self.rawTextToTextNode(p.text)
            macro.addParam(proc_param)
        for r in raw_entry.returns:
            proc_return = ProcReturn()
            proc_return.type = r.name.text
            proc_return.desc = self.rawTextToTextNode(r.text)
            macro.addReturn(proc_return)
        return macro


class MetafunctionConverter(CodeEntryConverter):
    def __init__(self, doc_proc):
        CodeEntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcMetafunction

    def process(self, raw_entry):
        metafunction = CodeEntryConverter.process(self, raw_entry)
        for t in raw_entry.tparams:
            proc_tparam = ProcTParam()
            proc_tparam.type = t.name.text
            proc_tparam.desc = self.rawTextToTextNode(t.text)
            metafunction.addTParam(proc_tparam)
        for r in raw_entry.returns:
            proc_return = ProcReturn()
            proc_return.type = r.name.text
            proc_return.desc = self.rawTextToTextNode(r.text)
            metafunction.addReturn(proc_return)
        return metafunction


class VariableConverter(CodeEntryConverter):
    def __init__(self, doc_proc):
        CodeEntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcVariable
    
    def process(self, raw_entry):
        variable = CodeEntryConverter.process(self, raw_entry)
        if raw_entry.type:
            variable.type = raw_entry.type.text
        return variable


class TagStack(object):
    """Helper class for processing nested HTML tags."""

    def __init__(self):
        self.stack = []
         
    def push(self, token):
        pass
             
    def pop(self, token):
        pass


class PageConverter(EntryConverter):
    """Process a RawPage into a Page object."""

    def __init__(self, doc_proc):
        EntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcPage
    
    def process(self, raw_entry):
        page = EntryConverter.process(self, raw_entry)
        # Convert the title
        if raw_entry.title:
            page.title = self.rawTextToTextNode(raw_entry.title)
        return page


class GroupConverter(EntryConverter):
    """Process a RawGroup into a Group object."""

    def __init__(self, doc_proc):
        EntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcGroup
    
    def process(self, raw_entry):
        page = EntryConverter.process(self, raw_entry)
        # Convert the title
        if raw_entry.title:
            page.title = self.rawTextToTextNode(raw_entry.title)
        return page


class TextNodeVisitor(object):
    """Interface/abstract base class for visiting text nodes of entries or such."""

    def visit(self, text_node):
        """Visit TextNode, possibly translating its content.

        @param text_node: TextNode object or None.
        """
        pass


class DocProcessor(object):
    """Convert a RawDoc object into a ProcDoc object.

    @ivar converters: Dict that maps RawEntry kinds to Converter objects.
    @ivar include_dir: The base path for including files  that can be used in
                       the @include and @snippet commands.
    @ivar include_mgr: inc_mgr.IncludeManager object for file/snippet
                       inclusion.
    """
    
    def __init__(self, logger=None, include_dir='.'):
        self.logger = logger
        self.include_dir = include_dir
        self.include_mgr = inc_mgr.IncludeManager(self.include_dir)
        self.converters = {
            'class': ClassConverter(self),
            'concept': ConceptConverter(self),
            'enum': EnumConverter(self),
            'adaption': AdaptionConverter(self),
            'global_typedef': TypedefConverter(self),
            'member_typedef': TypedefConverter(self),
            'grouped_typedef': TypedefConverter(self),
            'global_function': FunctionConverter(self),
            'global_metafunction': MetafunctionConverter(self),
            'defgroup': GroupConverter(self),
            'grouped_macro' : MacroConverter(self),
            'interface_function': FunctionConverter(self),
            'interface_metafunction': MetafunctionConverter(self),
            'macro' : MacroConverter(self),
            'member_function': FunctionConverter(self),
            'member_variable': VariableConverter(self),
            'page': PageConverter(self),
            'tag': TagConverter(self),
            'grouped_tag': TagConverter(self),
            'variable': VariableConverter(self),
            }

    def run(self, doc):
        res = ProcDoc()
        self.log('Processing Documentation...')
        self.convertTopLevelEntries(doc, res)
        self.convertSecondLevelEntries(doc, res)
        self.convertVariables(doc, res)
        self.checkLinks(doc, res)
        self.buildInheritanceLists(res)
        return res

    def convertTopLevelEntries(self, doc, res):
        """Convert top level entries.

        Variables are not converted yet.  They are converted in a separate
        step since they might encode enum values.
        """
        self.log('  1) Converting Top-Level Entries.')
        print 'doc.entries', [e.name.text for e in doc.entries]
        for raw_entry in doc.entries:
            # Get fitting converter or warn if there is none.
            kind = raw_entry.getType()
            if not kind in ['concept', 'class', 'global_function',
                            'global_metafunction', 'page', 'tag',
                            'defgroup', 'macro', 'adaption', 'global_typedef', 'enum']:
                continue  # Not a top-level entry.
            converter = self.converters.get(kind)
            if not converter:
                self.logWarning('Could not find converter for kind "%s".', kind)
                continue  # Skip if no converter could be found.
            # Perform conversion.
            proc_entry = converter.process(raw_entry)
            # Store object in ProcDoc.
            self.log('    * %s (%s)' % (proc_entry.name, proc_entry))
            res.addTopLevelEntry(proc_entry)

    def convertSecondLevelEntries(self, doc, res):
        self.log('  2) Converting Second-Level Entries.')
        for raw_entry in doc.entries:
            # Get fitting converter or warn if there is none.
            kind = raw_entry.getType()
            if not kind in ['member_function', 'interface_function',
                            'interface_metafunction',
                            'grouped_tag', 'grouped_macro', 'member_typedef',
                            'grouped_typedef']:
                continue  # Not a top-level entry.
            converter = self.converters.get(kind)
            if not converter:
                self.logWarning('Could not find converter for kind "%s".', kind)
                continue  # Skip if no converter could be found.
            # Perform conversion.
            proc_entry = converter.process(raw_entry)
            # Store object in ProcDoc.
            self.log('    * %s' % proc_entry.name)
            res.addSecondLevelEntry(proc_entry)

    def convertVariables(self, doc, res):
        self.log('  3) Converting Variable entries.')
        var_types = ['member_variable', 'grouped_variable', 'variable']
        for raw_entry in [e for e in doc.entries if e.getType() in var_types]:
            kind = raw_entry.getType()
            converter = self.converters.get(kind)
            if not converter:
                self.logWarning('Could not find converter for kind "%s".', kind)
                continue  # Skip if no converter could be found.
            # Perform conversion.
            proc_entry = converter.process(raw_entry)
            # Store object in ProcDoc.
            self.log('    * %s %s' % (proc_entry.type, proc_entry.name))
            res.addVariable(proc_entry)

    def checkLinks(self, doc, res):
        """Check <link> items of text nodes and references.

        References are given either explicitely in items like @extends and
        @implements.
        """
        self.log('  3) Checking References.')
        self.logWarning('    WARNING: Not implemented yet!')

    def buildInheritanceLists(self, doc):
        """Build lists regarding the inheritance in the classes and concepts in doc.

        We will build the equivalent to what Javadoc builds.

        For concepts, this is the list of (a) all extended concepts, (b) all
        known extending concepts, (c) all known implementing classes.

        For classes, this is the list of (a) all implemented concepts, (b) all
        direct known subclasses, (c) all extended classes.

        @param doc: The ProcDoc object with the classes and concept.

        """
        self.log('  4) Building Inheritance Lists.')
        # Process concepts: All extended and all extending.
        concepts = [x for x in doc.top_level_entries.values()
                    if x.kind == 'concept']
        # Get all concepts that c extends into c.all_extended.
        for c in concepts:
            q = list(c.extends)  # Queue for recursion
            while q:
                name = q[0]
                q.pop(0)
                if name in c.all_extended:
                    continue  # Skip to break loops.
                c.all_extended.add(name)
                q += doc.top_level_entries[name].extends
        # Now, build list of all extending concepts into c.all_extending.
        for c in concepts:
            for name in c.all_extended:
                doc.top_level_entries[name].all_extending.add(name)
        # Process classes: All extended and all extending classes.
        classes = [x for x in doc.top_level_entries.values()
                   if x.kind == 'class']
        # Get all classes that c extends into c.all_extended.
        for c in classes:
            q = list(c.extends)  # Queue for recursion
            while q:
                name = q[0]
                q.pop(0)
                if name in c.all_extended:
                    continue  # Skip to break loops.
                c.all_extended.add(name)
                q += doc.top_level_entries[name].extends
        # Now, build list of all extending clsses into c.all_extending.
        for c in classes:
            for name in c.all_extended:
                doc.top_level_entries[name].all_extending.add(c.name)
        # Build list of all implementing classes for all concepts.
        for cl in classes:
            for name in cl.implements:
                if '\u0001' in name:
                    continue  # Skip transitive inheritance.
                co = doc.top_level_entries[name]
                co.all_implementing.add(cl.name)
                co.all_implementing.update(cl.all_extending)
        # Build list of all implemented concepts for all classes.
        for co in concepts:
            for name in co.all_implementing:
                cl = doc.top_level_entries[name]
                cl.all_implemented.add(co.name)
            
        
    def log(self, msg, *args, **kwargs):
        """Print the given message to the configured logger if any.
        """
        if not self.logger:
            return
        self.logger.info(msg, *args, **kwargs)
        
    def logWarning(self, msg, *args, **kwargs):
        """Print the given message to the configured logger if any.
        """
        if not self.logger:
            return
        self.logger.warning('WARNING: ' + msg, *args, **kwargs)
