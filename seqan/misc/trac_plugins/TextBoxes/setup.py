from setuptools import setup

setup(
    name='TextBoxesMacros',
    version='0.1.3',
    packages=['text_boxes'],
    package_data={'text_boxes' : ['htdocs/css/*.css',
                                  'htdocs/dialog-information.png',
                                  'htdocs/dialog-warning.png',
                                  'htdocs/emblem-important.png',
                                  ]},
    entry_points = {'trac.plugins': ['textboxesmacros = text_boxes.macro']},
    )
