from setuptools import setup

setup(
    name='TextBoxesMacros',
    version='0.1',
    packages=['text_boxes'],
    package_data={'text_boxes' : ['htdocs/css/*.css']},
    entry_points = {'trac.plugins': ['textboxesmacros = text_boxes.macro']},
    )
