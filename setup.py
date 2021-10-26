# SPDX-License-Identifier: LGPL-2.1-or-later
import setuptools

setuptools.setup(
     name='channelunit',
     version='1.0',
     author="Joanna Jędrzejewska-Szmek",
     author_email="j.jedrzejewska-szmek@nencki.edu.pl",
     description="Compare nmodl channel activation and inactivation with experimental data ",
     url="https://github.com/asiaszmek/channelunit",
     classifiers=[
         "Programming Language :: Python :: 3",
         "Operating System :: OS Independent",
     ],
    project_urls={
        "Bug Tracker": "https://github.com/asiaszmek/channelunit/issues",
    },
    package_dir={"": "src"},
    python_requires=">=3.6",
    include_package_data=False,
    license='LGPL-2.1-or-later',
 )
