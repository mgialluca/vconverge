# -*- coding: utf-8 -*-
import os
from setuptools import setup

# Setup!
setup(
    name="vconverge",
    description="VPLANET parameter sweep and convergence helper",
    long_description=open("README.md", "r").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/mgialluca/vconverge",
    author="Megan Gialluca",
    author_email="gialluca@uw.edu",
    license="MIT",
    packages=["vconverge"],
    include_package_data=True,
    use_scm_version={
        "write_to": os.path.join("vconverge", "vconverge_version.py"),
        "write_to_template": '__version__ = "{version}"\n',
    },
    install_requires=[
        "numpy",
        "matplotlib",
        "argparse",
        "astropy"
    ],
    entry_points={
        "console_scripts": [
            "vconverge=vconverge.vconverge:main",
        ],
    },
    setup_requires=["setuptools_scm"],
    zip_safe=False,
)