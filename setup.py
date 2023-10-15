#!/usr/bin/env python

import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
	long_description = fh.read()

setuptools.setup(
	name = 'phylosophos',
	version = '1.1.1', 
	author = 'Min Hyung Cho', 
	author_email = 'mhcho@bmdrc.org', 
	description = 'PhyloSophos scientific name mapper', 
	long_description = long_description,
	long_description_content_type = 'text/markdown', 
	url = "https://github.com/mhcho4096/phylosophos", 
	python_requires = ">=3.8",
	install_requires= ["numpy"],
	py_modules=['phylosophos'],
	license='MIT',
)
