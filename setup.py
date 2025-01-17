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
    packages=["channelunit", "channelunit.capabilities", "channelunit.tests",
              "channelunit.tests", "channelunit.scores", "channelunit_tests",
              "demo_CA1"],
    include_package_data=True,
    package_data={
        'channelunit': [
            "mechanisms/*mod",
        ],
        "demo_CA1":  ["ion_channels/*mod",
                      "data/*csv",]

    },
    
    python_requires=">=3.11",
    license='LGPL-2.1-or-later',
    install_requires=['sciunit>=0.2.1']
 )
